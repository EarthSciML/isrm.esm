#!/usr/bin/env python3
"""Offline oracle for the plume-rise addition to ``isrm.esm``: how is each EGU
emission record's mass split across the SR matrix's three emission layers?

Not "which layer": InMAP's ``sr.Reader.layerFracs`` splits a record across TWO
layers whenever its plume's model layer falls strictly between two entries of
``sr.layers`` = ``[0, 3, 6]`` — model layers 1-2 between SR 0 and 1, model
layers 4-5 between SR 1 and 2 — interpolating on cell-CENTRE heights. So the
answer is three per-record weights summing to 1, plus the integer index of the
lower layer.

Why this exists
---------------
``isrm.esm`` states plume rise itself, in its own observeds, the way the
InMAP source-receptor tutorial (https://inmap.run/blog/2019/04/20/sr/) does —
it moves most of the emitted mass up into SR layers 1 and 2. Plume rise
changes exactly one intermediate quantity — the per-record
``(source cell, SR-layer weights)`` assignment — and everything downstream
follows from it. This script is the independent second opinion on that one
quantity.

Checking it does not need the 330 GB SR matrix: it needs 13 one-
dimensional meteorology/geometry arrays off the ISRM zarr, ~16 MB compressed
total. So this is the *cheap* check on the new physics, and it is the one to
run first when the document's deathsK moves. If the digest here matches what a
runner reports, the plume-rise physics is right and any remaining difference is
in the contraction; if it does not match, stop before touching the SR matrix.

The algorithm is InMAP's, not ours:

  * ``github.com/ctessum/atmos/plumerise`` — ``ASMEPrecomputed``,
    ``calcDeltaHPrecomputed``, ``findLayer``
  * ``github.com/spatialmodel/inmap/plumerise.go`` — ``(*Cell).IsPlumeIn``
  * ``github.com/spatialmodel/inmap/sr/srreader.go`` — ``(*Reader).Concentrations``,
    ``(*Reader).layerFracs``

One deliberate deviation from InMAP — see ABOVE_L7_NOTE below and the
``above_model_layer_7`` block in the emitted record.

Usage
-----
    python3 contract/plume_oracle.py [--out PATH] [--zip PATH] [--cache DIR]
                                     [--firstn N]

``--firstn N`` produces the REDUCED target: the first N delivered records, the
same truncation ``ISRM_FIRSTN=N`` applies to the three shims, so a reduced run
has something to be checked against. The full-scale expected values do not
apply to a truncation and are not asserted there.

Needs only the system python3 (3.9) plus numpy and numcodecs. The zarr chunks
are cached under ``$ISRM_SCRATCH`` (default ``/scratch.local/$USER/isrm-esm``,
the same resolution ``run-py/paths.py`` does) — NOT under /tmp, which is a
tmpfs on this cluster and comes out of the memory cgroup.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import subprocess
import sys
import zipfile

import numpy as np
import numcodecs

CONTRACT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(CONTRACT_DIR)
sys.path.insert(0, CONTRACT_DIR)
sys.path.insert(0, os.path.join(REPO, "run-py"))

from compare_results import Report  # noqa: E402  (contract/compare_results.py)
from results import field_summary, int_seq_sha256, ppl_sha256  # noqa: E402

# ---------------------------------------------------------------------------
# The ISRM grid, as published
# ---------------------------------------------------------------------------
S3 = "https://inmap-model.s3.amazonaws.com/isrm_v1.2.1.zarr"

NS = 52411          # cells per model layer, layers 0..7 (the ground grid)
N_HIGH = 9324       # cells per model layer, layers 8..26
N_GROUND_LAYERS = 8
N_LAYERS = 27
NALL = 596444       # 8*52411 + 19*9324

# Every one of these is a `[596444] float64` zarr array stored as a single
# blosc chunk keyed "0", laid out layer-major over the 27 model layers.
MET = ("Temperature", "WindSpeed", "S1", "SClass", "WindSpeedInverse",
       "WindSpeedMinusThird", "WindSpeedMinusOnePointFour", "LayerHeight", "Dz",
       "W", "S", "E", "N")
# `layers` is the SR matrix's own list of the MODEL layers the SR calculation was
# performed for. It is read only to be REPORTED — see SR_MODEL_LAYERS.
INT_ARRAYS = ("layers",)

# The model layers the three SR emission layers were computed for. DECLARED
# HERE, and deliberately NOT taken from the store, exactly as isrm.esm declares
# them in its SR_MODEL_L1 / SR_MODEL_L2 metaparameters.
#
# isrm_v1.2.1.ncf on Zenodo — the authoritative deposition — has
# layers = [0, 3, 6]. The isrm_v1.2.1.zarr this repo reads serves [0, 1, 2]
# instead: a machine-generated arange that displaced the real array during the
# zarr conversion. Corroboration independent of the NetCDF: sr/sr.go attaches
# `description: "Layer indices for which the SR calculation was performed"` to
# `layers`, and the zarr's copy has NO attributes at all while neighbouring
# InMAP variables (LayerHeight, Dz) keep theirs; the same conversion pass lost
# the per-cell `Layer` variable and stripped .zattrs off the SR arrays, which is
# why the shims need seed_empty_zattrs. Zenodo calls the three heights
# "ground-level, low-stack, and high-stack", which [0,3,6] (0-58 m / 253-391 m /
# 786-1049 m) is and [0,1,2] (three adjacent layers, all under 245 m) is not.
# The matrix DATA is fine: Zenodo's zip is md5-identical across depositions and
# its rows are byte-identical to the zarr's. One index array was replaced,
# nothing was rescaled.
SR_MODEL_LAYERS = [0, 3, 6]
SR_LAYERS_CORRUPT = [0, 1, 2]      # what the zarr serves, and must never be used

G = 9.80665                       # m/s2, plumerise.g
FT_TO_M = 0.3048

PATHWAYS = ("VOC", "NOx", "NH3", "SOx", "PM25")
# What the document calls the same five pathways (the SR matrix's pollutants).
PATHWAY_SR_NAME = {"VOC": "SOA", "NOx": "pNO3", "NH3": "pNH4",
                   "SOx": "pSO4", "PM25": "PrimaryPM25"}

ABOVE_L7_NOTE = (
    "Records whose plume rises above the top of model layer 7 — i.e. out of the "
    "52411-cell ground grid and into the 9324-cell high-altitude grid. InMAP has "
    "a latent defect here: sr.Reader.layerFracs clamps the emission into SR layer "
    "2 (correct — every model layer above 6 clamps), but the horizontal index it pairs with "
    "that layer is sr.indices[c], which srreader.go resets to 0 at every layer "
    "boundary. For a cell in layer >= 8 that index is a position in the COARSE "
    "9324-cell grid, and Reader.source then reinterprets it as a position in the "
    "52411-cell ground grid — so these emissions land on the wrong source cell. "
    "THIS ORACLE IMPLEMENTS THE CORRECT BEHAVIOUR: clamp to SR layer 2 at the "
    "source cell the emission actually came from. Any residual disagreement with "
    "the blog's published deathsK/deathsL is attributable to this group."
)

# ---------------------------------------------------------------------------
# Expected outputs. Verified during planning; the script fails loudly if the
# world drifts out from under them.
# ---------------------------------------------------------------------------
EXPECT = {
    "n_rec": 43650,
    "n_zero_height": 1100,
    "n_emitting_cells": 1520,
    "max_stack_height_m": 316.3824,
    "stack_layer_hist": [29313, 8460, 3899, 878],
    "plume_layer_hist": [12108, 9155, 11017, 6392, 2156, 1177, 553, 438, 654],
    "branch": {"momentum": 2073, "stable_buoyant": 868,
               "unstable_buoyant": 29601, "no_buoyancy": 11108},
    "max_plume_height_m": 9437.1,        # to 0.1 m
    "n_above_layer7": 654,
    # Pathway totals, short tons/yr, and the emitting-cell support set below:
    # both are the document's own values, from
    # contract/records/ground-level-only/python.json. Agreeing with them proves
    # the oracle ingests the FF10 zip and contains points in cells exactly the
    # way isrm.esm does — so a disagreement in the SR-layer assignment is about
    # plume rise and nothing else.
    #
    # They survived the plume-rise re-baselining that retired those records:
    # plume rise moves mass between LAYERS, never into or out of a pathway, and
    # changes which layer a record emits into, never which cell. A live
    # plume-era record must still carry these same five totals and this same
    # `ppl` digest, and compare_results.py checks that.
    "emis_total": {"VOC": 33452.80359612757, "NOx": 1314462.8823038577,
                   "NH3": 25012.479214860767, "SOx": 1571216.8812541936,
                   "PM25": 140822.70985515337},
    "ppl_sha256": "6f784d7e66f63872901126dabb2dd7354a96cdcd3d4585b2f52d20b6105a875b",
    # Record counts by LOWER SR layer. Not a partition of mass: a record whose
    # model layer is 1, 2, 4 or 5 is split across sr_lower and sr_lower + 1.
    "sr_lower_hist": [32280, 9725, 1645],
    # Percent of each pathway's mass landing in SR layer 0/1/2, to 0.1 pt.
    # These moved a long way when the document switched from the clamp
    # `sr_layer = min(plume_layer, 2)` — which followed the zarr's corrupt
    # `layers = [0, 1, 2]` and shoved ~90% of everything into SR layer 2 — to
    # sr.Reader.layerFracs over the authoritative [0, 3, 6]. SOx used to read
    # 0.5/4.3/95.2 here and now reads 12.7/51.2/36.1; the mass is the same, the
    # placement is not.
    "sr_split_pct": {
        "VOC":  [26.1, 54.4, 19.5],
        "NOx":  [22.4, 53.5, 24.1],
        "NH3":  [36.9, 54.7, 8.4],
        "SOx":  [12.7, 51.2, 36.1],
        "PM25": [21.9, 51.9, 26.2],
    },
}


# ---------------------------------------------------------------------------
# Chunk cache
# ---------------------------------------------------------------------------
def resolve_cache(explicit: str = "") -> str:
    """<scratch>/plume-oracle-zarr, with scratch resolved the run-py way."""
    if explicit:
        os.makedirs(explicit, exist_ok=True)
        return explicit
    try:
        import paths  # run-py/paths.py — the existing scratch convention
        scratch = paths.SCRATCH
    except Exception:                                        # pragma: no cover
        user = os.environ.get("USER", os.environ.get("LOGNAME", "user"))
        scratch = os.environ.get("ISRM_SCRATCH",
                                 os.path.join(f"/scratch.local/{user}", "isrm-esm"))
    d = os.path.join(scratch, "plume-oracle-zarr")
    os.makedirs(d, exist_ok=True)
    return d


def fetch(cache: str, name: str) -> str:
    """Download one single-chunk zarr array, once. ~16 MB for all of them."""
    p = os.path.join(cache, name + ".blosc")
    if not os.path.exists(p):
        tmp = p + ".part"
        subprocess.run(["curl", "-sSf", "-o", tmp, f"{S3}/{name}/0"], check=True)
        os.replace(tmp, p)
    return p


def rd(cache: str, name: str, dtype: str = "<f8") -> np.ndarray:
    with open(fetch(cache, name), "rb") as fh:
        return np.frombuffer(numcodecs.Blosc().decode(fh.read()), dtype=dtype)


# ---------------------------------------------------------------------------
# FF10 ingest — the blog's row handling, mirrored
# ---------------------------------------------------------------------------
VOCSET = {'VOC', 'VOC_INV', 'XYL', 'TOL', 'TERP', 'PAR', 'OLE', 'NVOL', 'MEOH', 'ISOP',
          'IOLE', 'FORM', 'ETOH', 'ETHA', 'ETH', 'ALD2', 'ALDX', 'CB05_ALD2', 'CB05_ALDX',
          'CB05_BENZENE', 'CB05_ETH', 'CB05_ETHA', 'CB05_ETOH', 'CB05_FORM', 'CB05_IOLE',
          'CB05_ISOP', 'CB05_MEOH', 'CB05_OLE', 'CB05_PAR', 'CB05_TERP', 'CB05_TOL',
          'CB05_XYL', 'ETHANOL', 'NHTOG', 'NMOG'}
PMSET = {'PM25-PRI', 'PM2_5', 'DIESEL-PM25', 'PAL', 'PCA', 'PCL', 'PEC', 'PFE', 'PK',
         'PMG', 'PMN', 'PMOTHR', 'PNH4', 'PNO3', 'POC', 'PSI', 'PSO4', 'PTI'}
NOXSET = {'NOX', 'HONO', 'NO', 'NO2'}


def read_ff10(zip_path: str):
    """Every recognised EGU record: (pathway, tons/yr, height m, diam m,
    temp K, velocity m/s, lon, lat).

    Blank stack fields become 0.0; a blank annual total drops the row; an
    unrecognised pollutant drops the row. Matches the document's ingest.
    """
    rows, members = [], []
    z = zipfile.ZipFile(zip_path)
    for nm in sorted(n for n in z.namelist() if "egu" in n and n.endswith(".csv")):
        members.append(nm)
        with z.open(nm) as f:
            for row in csv.reader(io.TextIOWrapper(f, newline='')):
                if not row or not row[0] or row[0][0] == '#':
                    continue
                pol, em = row[12], row[13]
                if em == '':
                    continue
                if pol in VOCSET:
                    k = 0
                elif pol in NOXSET:
                    k = 1
                elif pol == 'NH3':
                    k = 2
                elif pol == 'SO2':
                    k = 3
                elif pol in PMSET:
                    k = 4
                else:
                    continue
                rows.append((
                    k, float(em),
                    float(row[17]) * FT_TO_M if row[17] != '' else 0.0,   # stack height
                    float(row[18]) * FT_TO_M if row[18] != '' else 0.0,   # stack diameter
                    (float(row[19]) - 32) * 5.0 / 9.0 + 273.15 if row[19] != '' else 0.0,
                    float(row[21]) * FT_TO_M if row[21] != '' else 0.0,   # exit velocity
                    float(row[23]), float(row[24])))                      # lon, lat
    return rows, members


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def lcc_forward(lon, lat):
    """The InMAP/ISRM Lambert conformal conic: R=6370997, std 33/45, ref 40/-97."""
    R0, p1, p2 = 6370997.0, math.radians(33.0), math.radians(45.0)
    p0, l0 = math.radians(40.0), math.radians(-97.0)
    n = (math.log(math.cos(p1) / math.cos(p2))
         / math.log(math.tan(math.pi / 4 + p2 / 2) / math.tan(math.pi / 4 + p1 / 2)))
    F = R0 * math.cos(p1) * math.tan(math.pi / 4 + p1 / 2) ** n / n
    rho0 = F / math.tan(math.pi / 4 + p0 / 2) ** n
    rho = F / np.tan(np.pi / 4 + np.radians(lat) / 2) ** n
    theta = n * (np.radians(lon) - l0)
    return rho * np.sin(theta), rho0 - rho * np.cos(theta)


def containing_cell(W, S, E, N, x, y, rep: Report):
    """0-based ground-grid cell containing each point, half-open [W,E) x [S,N).

    Every cell edge in the ISRM ground grid is a multiple of 1 km (cell widths
    run 1, 2, 4, 12, 24, 48 km), so a dense 1 km raster is an exact lookup and
    costs ~86 MB. Asserted, not assumed.
    """
    fine = 1000.0
    for nm, a in (("W", W), ("S", S), ("E", E), ("N", N)):
        rep.check(bool(np.all(np.mod(a, fine) == 0)), f"grid/{nm}-on-1km",
                  "ground-grid edges are not all multiples of 1000 m")
    x0, y0 = float(W.min()), float(S.min())
    ncx = int(round((float(E.max()) - x0) / fine))
    ncy = int(round((float(N.max()) - y0) / fine))

    grid = np.full(ncx * ncy, -1, dtype=np.int32)
    ci = np.rint((W - x0) / fine).astype(np.int64)
    cj = np.rint((S - y0) / fine).astype(np.int64)
    cw = np.rint((E - W) / fine).astype(np.int64)
    ch = np.rint((N - S) / fine).astype(np.int64)
    for s in range(len(W)):
        base = ci[s] * ncy + cj[s]
        for a in range(cw[s]):
            grid[base + a * ncy: base + a * ncy + ch[s]] = s
    rep.check(not np.any(grid < 0), "grid/raster-complete",
              "the 1 km raster has holes — the ground grid is not a full cover")

    gi = np.floor((x - x0) / fine).astype(np.int64)
    gj = np.floor((y - y0) / fine).astype(np.int64)
    ok = (gi >= 0) & (gi < ncx) & (gj >= 0) & (gj < ncy)
    cell = np.full(len(x), -1, dtype=np.int64)
    cell[ok] = grid[gi[ok] * ncy + gj[ok]]
    return cell


# ---------------------------------------------------------------------------
# ASME plume rise
# ---------------------------------------------------------------------------
def plume_layers(cache, cell, hs, ds, Ts, vs, rep: Report):
    """ctessum/atmos ASMEPrecomputed + inmap IsPlumeIn + sr.Reader.layerFracs.

    Returns (stack_layer, plume_height, plume_layer, weights, sr_lower, branch)
    where `plume_layer == 8` is the sentinel "at or above the top of model layer
    7" (InMAP's plumerise.ErrAboveModelTop, on a column truncated at layer 7 —
    harmless for the SR answer because every model layer above SR_MODEL_LAYERS
    [-1] gets the same weights), `weights` is a (3, R) array of the per-record
    SR-layer shares and `sr_lower` is the integer index of the lower layer each
    record's mass goes to.
    """
    def col(a, l):
        return a[l * NS:(l + 1) * NS]

    LH, Dz = rd(cache, "LayerHeight"), rd(cache, "Dz")
    R = len(hs)
    s = np.where(cell >= 0, cell, 0)

    # IsPlumeIn builds layerHeights as a staggered [0, top of L0, top of L1, ...]
    # by accumulating Dz up the column. LayerHeight is that same staggered
    # coordinate, so lh[l] is the bottom of model layer l and lh[0] == 0.
    rep.check(bool(np.all(col(LH, 0) == 0.0)), "met/LayerHeight-L0-zero",
              "LayerHeight is not the staggered bottom-of-layer coordinate")
    lh = np.zeros((9, R))
    for l in range(1, N_GROUND_LAYERS):
        lh[l] = col(LH, l)[s]
    lh[8] = col(LH, 7)[s] + col(Dz, 7)[s]          # top of model layer 7

    # findLayer(layerHeights, h) = sort.SearchFloat64s(layerHeights, h), then
    # -1 unless the result is 0. Because lh[0] == 0, for h > 0 that is exactly
    # the count of strictly-lower boundaries above index 0; for h == 0 it is 0.
    stack_layer = np.zeros(R, dtype=np.int64)
    for l in range(1, N_GROUND_LAYERS):
        stack_layer += (lh[l] < hs).astype(np.int64)
    rep.check(bool(np.all(hs < lh[8])), "plume/stack-inside-ground-grid",
              "a stack is taller than the top of model layer 7")

    def met(name):
        """Value of `name` in the stack's own cell (column s, layer stack_layer)."""
        a = rd(cache, name)
        o = np.zeros(R)
        for l in range(N_GROUND_LAYERS):
            m = stack_layer == l
            o[m] = col(a, l)[s[m]]
        return o

    Ta, u = met("Temperature"), met("WindSpeed")
    sc, s1 = met("SClass"), met("S1")
    uinv = met("WindSpeedInverse")
    u13 = met("WindSpeedMinusThird")
    u14 = met("WindSpeedMinusOnePointFour")

    # calcDeltaHPrecomputed, branch for branch.
    mom = ((Ts - Ta) < 50.) & (vs > u) & (vs > 10.)
    td = np.where(Ts - Ta == 0, 0.0,
                  2 * (Ts - Ta) / np.where(Ts + Ta == 0, 1.0, Ts + Ta))
    F = G * td * vs * (ds / 2) ** 2                # buoyancy flux, m4/s3

    buoy = ~mom
    stable = buoy & (sc > 0.5) & (s1 != 0) & (F > 0)
    unstable = buoy & ~stable & (F > 0)
    none = buoy & ~stable & ~unstable              # F <= 0: deltaH = 0

    # Each branch is evaluated only on its own records — Go evaluates one
    # branch per record, and the discarded branches here would take cube roots
    # of negative numbers on records the branch does not own.
    dH = np.zeros(R)
    dH[mom] = ds[mom] * np.power(vs[mom], 1.4) * u14[mom]
    dH[stable] = 29. * np.power(F[stable] / s1[stable], 0.333333333) * u13[stable]
    dH[unstable] = (7.4 * np.power(F[unstable] * np.power(hs[unstable], 2.),
                                   0.333333333) * uinv[unstable])
    rep.check(not np.any(np.isnan(dH)), "plume/deltaH-finite", "deltaH has NaNs")

    plume_height = hs + dH
    plume_layer = np.zeros(R, dtype=np.int64)
    for l in range(1, 9):
        plume_layer += (lh[l] < plume_height).astype(np.int64)

    # Concentrations(): `if e.Height != 0 { IsPlumeIn(...) } else { ground }`.
    # A zero stack height skips plume rise entirely and stays in layer 0.
    plume_layer[hs == 0] = 0

    # ---- sr.Reader.layerFracs, with sr.layers == SR_MODEL_LAYERS ----------
    #
    #   for i := range sr.layers            { if c.Layer == sr.layers[i]  -> {i}, {1} }
    #   for i := range sr.layers[:len-1]    { if sr.layers[i] < c.Layer < sr.layers[i+1] {
    #                                             below := layerHeights[sr.layers[i]]
    #                                             above := layerHeights[sr.layers[i+1]]
    #                                             frac  := (plumeHeight - below)/(above - below)
    #                                             -> {i, i+1}, {frac, 1 - frac} } }
    #   if c.Layer > sr.layers[len-1]       { -> {len-1}, {1}, AboveTopErr }
    #
    # `layerHeights` comes from sr.d.VerticalProfile, which fills
    # `height[i] = c.LayerHeight + c.Dz/2.` — CELL-CENTRE heights, not layer
    # bottoms. Getting that wrong is worth ~0.4% on the national total.
    #
    # THE INTERPOLATION IS INVERTED: the LOWER layer gets `frac`, so a HIGHER
    # plume puts MORE weight on the LOWER layer. That is an InMAP bug. It
    # conserves mass, and this oracle reproduces it rather than fixing it,
    # because the question it answers is what InMAP computes.
    #
    # Every model layer above the top entry gets the same weights as an exact
    # match on it, so capping the model layer there loses nothing — which is
    # exactly the saturation isrm.esm's `plume_model_layer` performs.
    top = SR_MODEL_LAYERS[-1]
    M = np.minimum(plume_layer, top)

    def cent(l):
        """Cell-CENTRE height of model layer l over each record's source cell."""
        return col(LH, l)[s] + 0.5 * col(Dz, l)[s]

    weights = np.zeros((3, R))
    sr_lower = np.zeros(R, dtype=np.int64)
    for i in range(len(SR_MODEL_LAYERS)):
        L = SR_MODEL_LAYERS[i]
        exact = M == L if L != top else M >= L
        weights[i][exact] = 1.0
        sr_lower[exact] = i
    for i in range(len(SR_MODEL_LAYERS) - 1):
        lo, hi = SR_MODEL_LAYERS[i], SR_MODEL_LAYERS[i + 1]
        between = (M > lo) & (M < hi)
        below, above = cent(lo), cent(hi)
        span = above - below
        f = (plume_height - below) / np.where(span == 0, 1.0, span)
        weights[i][between] = f[between]
        weights[i + 1][between] = 1.0 - f[between]
        sr_lower[between] = i
    dev = float(np.max(np.abs(weights.sum(axis=0) - 1.0))) if R else 0.0
    rep.check(dev <= 1e-12, "plume/weights-sum-to-one",
              f"max |w0+w1+w2 - 1| = {dev} — layerFracs conserves mass, so a "
              "record's weights must sum to 1")

    branch = {"momentum": int(mom.sum()), "stable_buoyant": int(stable.sum()),
              "unstable_buoyant": int(unstable.sum()), "no_buoyancy": int(none.sum())}
    return stack_layer, plume_height, plume_layer, weights, sr_lower, branch


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(CONTRACT_DIR, "records",
                                                  "plume_oracle.json"))
    ap.add_argument("--zip", default="", help="FF10 point-source zip (default: "
                    "$EGU_ZIP, else <repo>/data/2016fd_inputs_point.zip)")
    ap.add_argument("--cache", default="", help="zarr chunk cache dir "
                    "(default: $ISRM_SCRATCH/plume-oracle-zarr)")
    ap.add_argument("--no-assert", action="store_true",
                    help="emit the record even if the expected values drift")
    ap.add_argument("--firstn", type=int, default=0, metavar="N",
                    help="REDUCED target: keep only the first N delivered "
                         "records, the same truncation ISRM_FIRSTN applies to "
                         "the shims. The full-scale expected values do not "
                         "apply to a truncation and are not asserted; the "
                         "structural facts about the grid and the algorithm "
                         "still are.")
    args = ap.parse_args(argv)
    firstn = args.firstn if args.firstn and args.firstn > 0 else 0
    full = firstn == 0
    if not full and args.out == ap.get_default("out"):
        args.out = os.path.join(CONTRACT_DIR, "records",
                                "plume_oracle_first%d.json" % firstn)

    zip_path = args.zip or os.environ.get("EGU_ZIP", "") or os.path.join(
        REPO, "data", "2016fd_inputs_point.zip")
    if not os.path.isfile(zip_path):
        print(f"error: EGU zip not found at {zip_path}; set EGU_ZIP or pass --zip",
              file=sys.stderr)
        return 2
    cache = resolve_cache(args.cache)
    rep = Report()
    print(f"zip:   {zip_path}\ncache: {cache}")

    # ---- grid ------------------------------------------------------------
    W, S, E, N = (rd(cache, v) for v in ("W", "S", "E", "N"))
    rep.check(len(W) == NALL, "grid/n-all", f"{len(W)} != {NALL}")
    rep.check(N_GROUND_LAYERS * NS + (N_LAYERS - N_GROUND_LAYERS) * N_HIGH == NALL,
              "grid/layer-arithmetic", "8*52411 + 19*9324 != 596444")
    for l in range(1, N_GROUND_LAYERS):
        same = all(np.array_equal(a[:NS], a[l * NS:(l + 1) * NS]) for a in (W, S, E, N))
        rep.check(same, f"grid/L{l}-matches-L0",
                  "model layers 0..7 must share a byte-identical horizontal grid")
    # `layers` is READ ONLY TO BE REPORTED. This oracle interpolates on the
    # DECLARED SR_MODEL_LAYERS, because the store's copy is the corrupt arange.
    # The check that used to live here asserted `layers == [0, 1, 2]` — i.e. it
    # asserted the corruption, and would have gone green forever while the
    # physics was wrong. Inverted, it is two checks:
    #   * the layers this oracle actually USES are never the corrupt arange —
    #     a tripwire on the constant, so a future edit that "simplifies"
    #     SR_MODEL_LAYERS back to the store's value fails here rather than
    #     silently returning the old wrong answer;
    #   * the store's own value is one of the two states we have evidence for.
    #     Anything else means something changed upstream and is a hard failure,
    #     not a note.
    # The store serving [0, 1, 2] is therefore reported loudly and recorded in
    # the emitted record (grid.sr_layers_in_store), never assumed.
    srl = rd(cache, "layers", "<i4")  # zarr dtype for `layers` is <i4
    store_layers = srl.tolist()
    rep.check(SR_MODEL_LAYERS != SR_LAYERS_CORRUPT, "sr/layers-not-the-arange",
              "the SR model layers this oracle interpolates on are the zarr's "
              "corrupt arange [0,1,2] — that is the bug, not the baseline")
    rep.check(store_layers in (SR_LAYERS_CORRUPT, SR_MODEL_LAYERS),
              "sr/layers-store-state",
              f"the store serves layers = {store_layers}, which is neither the "
              f"known-corrupt arange {SR_LAYERS_CORRUPT} nor the authoritative "
              f"{SR_MODEL_LAYERS} — something changed upstream; find out what "
              "before trusting anything downstream of it")
    if store_layers == SR_LAYERS_CORRUPT:
        print(f"NOTE: the store still serves the corrupt `layers` = "
              f"{store_layers}; interpolating on the declared "
              f"{SR_MODEL_LAYERS} instead (see SR_MODEL_LAYERS).")
    elif store_layers == SR_MODEL_LAYERS:
        print(f"NOTE: the store now serves `layers` = {store_layers} — the "
              "zarr conversion bug has been FIXED upstream and the declared "
              "value could be read from the store again.")

    W, S, E, N = W[:NS], S[:NS], E[:NS], N[:NS]

    # ---- emissions -------------------------------------------------------
    rows, members = read_ff10(zip_path)
    if not full:
        # The shims truncate the DELIVERED record list — the rows that survive
        # the loader's record_filter — so truncating here, after read_ff10 has
        # applied the same filter, picks exactly the same records.
        print(f"REDUCED: first {firstn} of {len(rows)} delivered records")
        rows = rows[:firstn]
    n_rec = len(rows)
    pol = np.array([r[0] for r in rows], dtype=np.int64)
    em = np.array([r[1] for r in rows])
    hs = np.array([r[2] for r in rows])
    ds = np.array([r[3] for r in rows])
    Ts = np.array([r[4] for r in rows])
    vs = np.array([r[5] for r in rows])
    lon = np.array([r[6] for r in rows])
    lat = np.array([r[7] for r in rows])
    n_zero = int((hs == 0).sum())
    if full:
        rep.check(n_rec == EXPECT["n_rec"], "emis/n_rec",
                  f"{n_rec} != {EXPECT['n_rec']}")
        rep.check(n_zero == EXPECT["n_zero_height"], "emis/n_zero_height",
                  f"{n_zero} != {EXPECT['n_zero_height']}")
        rep.check(abs(hs.max() - EXPECT["max_stack_height_m"]) < 1e-6,
                  "emis/max_stack_height",
                  f"{hs.max()} != {EXPECT['max_stack_height_m']}")
    else:
        rep.check(n_rec == firstn, "emis/n_rec_reduced",
                  f"{n_rec} != the requested {firstn}")

    x, y = lcc_forward(lon, lat)
    cell = containing_cell(W, S, E, N, x, y, rep)
    n_outside = int((cell < 0).sum())
    rep.check(n_outside == 0, "emis/all-in-grid",
              f"{n_outside} records fall outside the ISRM ground grid")
    emitting = sorted(set(int(c) for c in cell if c >= 0))
    # 1-based, matching the document's `ppl`.
    ppl_ids = [c + 1 for c in emitting]
    ppl_hash = ppl_sha256(ppl_ids)
    if full:
        rep.check(len(emitting) == EXPECT["n_emitting_cells"], "emis/n_emitting_cells",
                  f"{len(emitting)} != {EXPECT['n_emitting_cells']}")
        rep.check(ppl_hash == EXPECT["ppl_sha256"], "emis/ppl_sha256",
                  f"{ppl_hash} != the document's ppl digest")

    # ---- plume rise ------------------------------------------------------
    stack_layer, ph, plume_layer, weights, sr_lower, branch = plume_layers(
        cache, cell, hs, ds, Ts, vs, rep)

    pos = hs > 0
    sl_hist = np.bincount(stack_layer[pos], minlength=4).tolist()
    pl_hist = np.bincount(plume_layer, minlength=9).tolist()
    if full:
        rep.check(sl_hist == EXPECT["stack_layer_hist"], "plume/stack_layer_hist",
                  f"{sl_hist} != {EXPECT['stack_layer_hist']}")
        rep.check(pl_hist == EXPECT["plume_layer_hist"], "plume/model_layer_hist",
                  f"{pl_hist} != {EXPECT['plume_layer_hist']}")
        rep.check(branch == EXPECT["branch"], "plume/branch_usage",
                  f"{branch} != {EXPECT['branch']}")
        rep.check(abs(round(float(ph.max()), 1) - EXPECT["max_plume_height_m"]) < 5e-2,
                  "plume/max_height", f"{ph.max()} != ~{EXPECT['max_plume_height_m']}")
        sr_hist = np.bincount(sr_lower, minlength=3).tolist()
        rep.check(sr_hist == EXPECT["sr_lower_hist"], "plume/sr_lower_hist",
                  f"{sr_hist} != {EXPECT['sr_lower_hist']}")

    # ---- the above-layer-7 group (InMAP's latent defect) -----------------
    above = plume_layer >= N_GROUND_LAYERS
    n_above = int(above.sum())
    if full:
        rep.check(n_above == EXPECT["n_above_layer7"], "plume/n_above_layer7",
                  f"{n_above} != {EXPECT['n_above_layer7']}")
    em_total_all = float(em.sum())
    above_mass = float(em[above].sum())

    # ---- per-pathway split by SR layer -----------------------------------
    pathways = {}
    for k, name in enumerate(PATHWAYS):
        m = pol == k
        tot = float(em[m].sum())
        # A record split across two SR layers contributes a SHARE to each, so
        # this is a weighted sum, not a partition of the records.
        by = [float((em[m] * weights[L][m]).sum()) for L in (0, 1, 2)]
        pathways[name] = {
            "sr_name": PATHWAY_SR_NAME[name],
            "emis_sum": tot,
            "n_rec": int(m.sum()),
            "by_sr_layer": by,
            "frac_by_sr_layer": [b / tot if tot else 0.0 for b in by],
            "above_model_layer_7": {"n_rec": int((m & above).sum()),
                                    "emis_sum": float(em[m & above].sum())},
        }
        if full:
            rep.check(abs(tot - EXPECT["emis_total"][name]) <= 1e-6 * abs(tot),
                      f"emis/total/{name}",
                      f"{tot} != {EXPECT['emis_total'][name]}")
            pct = [round(100.0 * b / tot, 1) for b in by]
            rep.check(pct == EXPECT["sr_split_pct"][name], f"emis/sr_split/{name}",
                      f"{pct} != {EXPECT['sr_split_pct'][name]}")

    # ---- record ----------------------------------------------------------
    phs = field_summary(ph)
    phs["mean"] = phs["sum"] / len(ph)
    rec = {
        "kind": "plume_oracle",
        "firstn": firstn,
        "binding": "python",
        "purpose": ("record -> SR layer assignment for the plume-rise addition to "
                    "isrm.esm; checkable without fetching the SR matrix"),
        "source": {
            "emissions_zip": os.path.basename(zip_path),
            "members": members,
            "isrm_zarr": S3,
            "arrays": list(MET) + list(INT_ARRAYS),
        },
        "grid": {
            "n_cells_all": NALL,
            "n_cells_ground": NS,
            "n_cells_high": N_HIGH,
            "n_model_layers": N_LAYERS,
            "n_ground_grid_layers": N_GROUND_LAYERS,
            # Two different things, deliberately reported side by side.
            "sr_layers_used": list(SR_MODEL_LAYERS),
            "sr_layers_in_store": store_layers,
            "sr_layers_note": (
                "sr_layers_used is DECLARED (isrm_v1.2.1.ncf on Zenodo, and "
                "isrm.esm's SR_MODEL_L1/SR_MODEL_L2 metaparameters). "
                "sr_layers_in_store is what isrm_v1.2.1.zarr serves, which is a "
                "machine-generated arange that displaced the real array during "
                "the zarr conversion; it is reported, never used."),
        },
        "n_rec": n_rec,
        "n_emitting_cells": len(emitting),
        "n_records_outside_grid": n_outside,
        "ppl": {"count": len(ppl_ids), "sha256": ppl_hash},
        "sr_lower": {
            "sha256": int_seq_sha256(sr_lower),
            "encoding": ("sha256 over the per-record LOWER SR layer as ASCII "
                         "decimals joined by ',' in record order "
                         "(contract/results.py int_seq_sha256)"),
            "histogram": np.bincount(sr_lower, minlength=3).tolist(),
            "note": ("the lower of the at most two SR layers layerFracs charges "
                     "the record to; the upper, when there is one, is always "
                     "this + 1. It replaces the old per-record `sr_layer`, "
                     "which stopped describing the model once a record could be "
                     "split across two layers."),
        },
        "weights": {
            "count": n_rec,
            "max_sum_error": float(np.max(np.abs(weights.sum(axis=0) - 1.0)))
            if n_rec else 0.0,
            **{"w_sr%d" % L: field_summary(weights[L]) for L in (0, 1, 2)},
            "note": ("sr.Reader.layerFracs's per-record shares, interpolated on "
                     "CELL-CENTRE heights (sr.d.VerticalProfile fills "
                     "height[i] = c.LayerHeight + c.Dz/2), with the LOWER layer "
                     "taking `frac` — InMAP's inverted interpolation, "
                     "reproduced rather than corrected. FLOATS, so these are "
                     "compared within RTOL_FIELD and their sha256 is "
                     "informative, not asserted: they descend from "
                     "plume_height, whose cube roots differ by an ulp between "
                     "languages."),
        },
        "stack_layer": {
            "n_zero_height": n_zero,
            "max_stack_height_m": float(hs.max()),
            # Two histograms, because they answer two questions. The
            # height>0 one is the physics: how many STACKS top out in each
            # model layer, with the zero-height records (which never enter
            # plume rise) excluded. The all-records one is what a runner can
            # emit from the document's `stack_layer` observed alone — that
            # observed is 0 for a zero-height record, so bin 0 carries them —
            # and is therefore the one the contract compares.
            "histogram_height_gt_0": sl_hist,
            "histogram": np.bincount(stack_layer, minlength=4).tolist(),
            "sha256": int_seq_sha256(stack_layer),
        },
        "plume_model_layer": {
            "histogram": pl_hist,
            "note": ("index 8 is the sentinel 'at or above the top of model layer "
                     "7' (plumerise.ErrAboveModelTop on a column truncated at "
                     "layer 7); it is not a model layer index"),
            "sha256": int_seq_sha256(plume_layer),
        },
        "plume_height_m": phs,
        "branch_usage": branch,
        "pathways": pathways,
        "emis_total": em_total_all,
        "above_model_layer_7": {
            "n_rec": n_above,
            "emis_sum": above_mass,
            "emis_frac": above_mass / em_total_all,
            "sr_layer_assigned": 2,
            "oracle_behaviour": "clamp to SR layer 2 at the correct source cell",
            "inmap_behaviour": ("clamp to SR layer 2 but keep the 9324-cell "
                                "high-layer index, reinterpreted against the "
                                "52411-cell ground grid"),
            "note": ABOVE_L7_NOTE,
        },
        "units": {"emissions": "short tons/yr", "heights": "m"},
    }

    # ---- report ----------------------------------------------------------
    print(f"\nrecords: {n_rec}, zero-height records: {n_zero}, "
          f"max stack height: {hs.max():.4f} m")
    print(f"distinct emitting cells: {len(emitting)}")
    print(f"stack-layer histogram (records with height>0): {sl_hist}")
    print(f"plume model-layer histogram: {pl_hist}")
    print("   (index 8 = \"at or above model layer 8\")")
    print("branch usage: momentum {momentum}, stable-buoyant {stable_buoyant}, "
          "unstable-buoyant {unstable_buoyant}, F<=0 {no_buoyancy}".format(**branch))
    print(f"max plume height: {ph.max():.1f} m")
    print("\nSR-layer emission split (short tons/yr, % of pathway total):")
    for name in PATHWAYS:
        p = pathways[name]
        f = p["frac_by_sr_layer"]
        print(f"  {name:<6} total {p['emis_sum']:>12.1f}   "
              f"L0 {100*f[0]:>4.1f}%  L1 {100*f[1]:>4.1f}%  L2 {100*f[2]:>4.1f}%")
    print(f"\nabove model layer 7: {n_above} records, {above_mass:.1f} tons/yr "
          f"({100 * above_mass / em_total_all if em_total_all else 0.0:.2f}% "
          f"of emitted mass) — "
          f"clamped to SR layer 2 at the correct source cell")
    print(f"lower-SR-layer histogram: {rec['sr_lower']['histogram']}")
    print(f"sr_lower sha256: {rec['sr_lower']['sha256']}")
    print("weight sums (w_sr0/w_sr1/w_sr2): " + " / ".join(
        "%.9f" % rec["weights"]["w_sr%d" % L]["sum"] for L in (0, 1, 2)) +
        "   max|Σw - 1| = %g" % rec["weights"]["max_sum_error"])

    ok = rep.close()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    if ok or args.no_assert:
        with open(args.out, "w") as fh:
            json.dump(rec, fh, indent=2, sort_keys=True)
        print(f"wrote {args.out}")
    else:
        print("NOT writing the record: expected values drifted "
              "(pass --no-assert to write anyway)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
