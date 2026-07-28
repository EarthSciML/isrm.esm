"""HARNESS INPUT PREPARATION for the Python ISRM runner.

Everything in this file is *outside* the model: it produces exactly the arrays
`isrm_pushdown.esm` declares as `parameter`s, and nothing more. Specifically:

    X, Y            emission-point coordinates in the ISRM LCC plane
    emis_annual     annual emission value per record (tons/yr, unconverted)
    pollutant       pathway CODE per record (the .esm's is_* masks select by band)
    W, S, E, N      full-grid cell rectangles
    TotalPop        population per cell
    MortalityRate   baseline mortality rate per cell

The strict bar for this project is that every *model* number comes from the
binding's own evaluation of the .esm graph. Reading bytes and projecting
lon/lat into the grid's plane is the sanctioned impure boundary on the input
side — the .esm declares X/Y as parameters, so the projection is upstream of the
spec, and the Julia runner does the identical thing with the identical
constants (``run-model-jl-pushdown/l3_common.jl``). Everything DOWNSTREAM of
these arrays — the spatial join, the binning, the contraction, the deaths — is
evaluated from the graph.

The two readers are EarthSciIO's: the FF10 point reader over the EGU zip (the
`member=` kwarg extracts a zip member without disturbing the content-addressed
cache key) and the zarr reader over the live ISRM store.
"""

from __future__ import annotations

import math
import os
import zipfile
from dataclasses import dataclass

import numpy as np

import earthsciio as eio

from paths import EGU_ZIP, N_SRC, SCRATCH, ZARR_URL

# --------------------------------------------------------------------------- #
# Snyder 1987 spherical Lambert Conformal Conic forward.
# Matches lambert_conformal.esm and run-model-jl's validated implementation
# constant for constant; the two runners MUST agree here or their X/Y differ
# and the derived support set diverges for a reason that is not the model.
# --------------------------------------------------------------------------- #
LAT_1, LAT_2, LAT_0, LON_0, LCC_R = 33.0, 45.0, 40.0, -97.0, 6370997.0
D2R = 0.017453292519943295


def _lcc_t(lat: float) -> float:
    return math.tan(0.7853981633974483 + lat * 0.008726646259971648)


LCC_N = math.log(math.cos(LAT_1 * D2R) / math.cos(LAT_2 * D2R)) / math.log(
    _lcc_t(LAT_2) / _lcc_t(LAT_1)
)
LCC_F = math.cos(LAT_1 * D2R) * _lcc_t(LAT_1) ** LCC_N / LCC_N


def _lcc_rho(lat):
    return LCC_R * LCC_F / _lcc_t_arr(lat) ** LCC_N


def _lcc_t_arr(lat):
    return np.tan(0.7853981633974483 + np.asarray(lat) * 0.008726646259971648)


LCC_RHO0 = LCC_R * LCC_F / _lcc_t(LAT_0) ** LCC_N


def lcc_forward(lon, lat):
    """Vectorised forward projection -> (X, Y) in metres."""
    theta = LCC_N * (np.asarray(lon) - LON_0) * D2R
    rho = _lcc_rho(lat)
    return rho * np.sin(theta), LCC_RHO0 - rho * np.cos(theta)


# --------------------------------------------------------------------------- #
# POLID -> pathway (STEP 0.6). Verbatim from run-model-jl/run-model.jl so both
# runners recognise exactly the same 43,650 records.
# --------------------------------------------------------------------------- #
VOC_SET = {
    "VOC", "VOC_INV", "XYL", "TOL", "TERP", "PAR", "OLE", "NVOL", "MEOH", "ISOP", "IOLE",
    "FORM", "ETOH", "ETHA", "ETH", "ALD2", "ALDX", "CB05_ALD2", "CB05_ALDX", "CB05_BENZENE",
    "CB05_ETH", "CB05_ETHA", "CB05_ETOH", "CB05_FORM", "CB05_IOLE", "CB05_ISOP", "CB05_MEOH",
    "CB05_OLE", "CB05_PAR", "CB05_TERP", "CB05_TOL", "CB05_XYL", "ETHANOL", "NHTOG", "NMOG",
}
PM25_SET = {
    "PM25-PRI", "PM2_5", "DIESEL-PM25", "PAL", "PCA", "PCL", "PEC", "PFE", "PK", "PMG",
    "PMN", "PMOTHR", "PNH4", "PNO3", "POC", "PSI", "PSO4", "PTI",
}
NOX_SET = {"NOX", "HONO", "NO", "NO2"}
NH3_SET = {"NH3"}
SOX_SET = {"SO2"}

#: Pathway CODE per record. The .esm's is_VOC/is_NOx/... observeds select by code
#: BAND (SOA 1-35, pNO3 36-39, pNH4 40, pSO4 41, PrimaryPM25 42-59), so the code
#: — not the pathway name — is the model's input. 0.0 = unrecognised (dropped).
PATH_CODE = {"VOC": 1.0, "NOx": 36.0, "NH3": 40.0, "SOx": 41.0, "PM25": 42.0}


def pathway_code(polid: str) -> float:
    p = str(polid).strip().upper()
    if p in VOC_SET:
        return PATH_CODE["VOC"]
    if p in PM25_SET:
        return PATH_CODE["PM25"]
    if p in NOX_SET:
        return PATH_CODE["NOx"]
    if p in NH3_SET:
        return PATH_CODE["NH3"]
    if p in SOX_SET:
        return PATH_CODE["SOx"]
    return 0.0


# --------------------------------------------------------------------------- #


@dataclass
class Inputs:
    """The .esm's declared parameter arrays, ready to bind as const arrays."""

    X: np.ndarray
    Y: np.ndarray
    emis_annual: np.ndarray
    pollutant: np.ndarray
    W: np.ndarray
    S: np.ndarray
    E: np.ndarray
    N: np.ndarray
    TotalPop: np.ndarray
    MortalityRate: np.ndarray
    N_REC: int
    N_SRC: int


def _egu_members(zippath: str) -> list[str]:
    """The EGU point CSV members of the 2016fd inputs zip (same selection rule as
    run-model.jl: name contains 'egu' and ends in '.csv')."""
    with zipfile.ZipFile(zippath) as z:
        return sorted(
            n for n in z.namelist() if "egu" in n.lower() and n.lower().endswith(".csv")
        )


def _read_ff10_member(zippath: str, member: str, variables, workdir: str):
    """One EGU zip member -> a NativeDataset, with the column-name header row
    stripped first.

    These 2016fd members carry a ``country_cd,region_cd,...`` header row that is
    NOT a ``#`` comment, and neither EarthSciIO's Python nor its Julia FF10
    reader skips one (the fixed 77-column schema comes from the spec, not from
    the file). run-model.jl works around it the same way — filter the header
    line, then hand the reader a clean CSV — so both runners see identical
    records. Kept in the runner rather than patched into the reader so the two
    bindings' input paths stay equivalent.
    """
    os.makedirs(workdir, exist_ok=True)
    with zipfile.ZipFile(zippath) as z:
        text = z.read(member).decode("utf-8")
    lines = [ln for ln in text.splitlines() if not ln.strip().lower().startswith("country_cd")]
    tmp = os.path.join(workdir, os.path.basename(member))
    with open(tmp, "w") as fh:
        fh.write("\n".join(lines))
    try:
        return eio.FF10Reader().read_native(tmp, variables=variables, kind="point")
    finally:
        os.remove(tmp)


def read_emissions(zippath: str = EGU_ZIP, cache_root: str | None = None):
    """Read every EGU FF10 point record through EarthSciIO's FF10 reader.

    Returns ``(lon, lat, ann, code)`` for the RECOGNISED-pathway, finite records
    only — the same filter run-model.jl applies, so N_REC matches (43,650 at full
    scale).
    """
    workdir = os.path.join(cache_root or SCRATCH, "py_ff10_tmp")
    lon_l, lat_l, ann_l, code_l = [], [], [], []
    members = _egu_members(zippath)
    if not members:
        raise RuntimeError(f"no EGU csv members found in {zippath}")
    for m in members:
        ds = _read_ff10_member(
            zippath, m, ["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"], workdir
        )
        polid = np.asarray(ds["POLID"].data)
        ann = np.asarray(ds["ANN_VALUE"].data, dtype=float)
        lon = np.asarray(ds["LONGITUDE"].data, dtype=float)
        lat = np.asarray(ds["LATITUDE"].data, dtype=float)
        code = np.array([pathway_code(p) for p in polid], dtype=float)
        lon_l.append(lon)
        lat_l.append(lat)
        ann_l.append(ann)
        code_l.append(code)
    lon = np.concatenate(lon_l)
    lat = np.concatenate(lat_l)
    ann = np.concatenate(ann_l)
    code = np.concatenate(code_l)
    # Drop unrecognised pathways and any non-finite coordinate/value. This also
    # removes a stray header row, whose numeric columns parse to NaN.
    keep = (code > 0.0) & np.isfinite(lon) & np.isfinite(lat) & np.isfinite(ann)
    return lon[keep], lat[keep], ann[keep], code[keep]


def read_grid(cache_root: str | None = None):
    """W/S/E/N + TotalPop/MortalityRate for the first N_SRC cells of the ISRM
    zarr — the container grid, exactly the prefix run-model.jl bins against
    (pop_cells == SR source cells)."""
    root = cache_root or os.path.join(SCRATCH, "py_cache_meta")
    cache = eio.Cache(root=root)
    out = {}
    for group in (["W", "S", "E", "N"], ["TotalPop", "MortalityRate"]):
        loader = eio.DataLoader(
            name="ISRM_grid", format="zarr", url=ZARR_URL, variables=group
        )
        ds = eio.Provider(loader, cache).materialize()
        for v in group:
            out[v] = np.asarray(ds[v].data, dtype=float)[:N_SRC].copy()
    return out


def build_inputs(firstn: int | None = None) -> Inputs:
    """Assemble every declared parameter array. ``firstn`` truncates the emission
    record list for a reduced, fast run (the Python mirror of L3_FIRSTN)."""
    lon, lat, ann, code = read_emissions()
    if firstn is not None:
        lon, lat, ann, code = lon[:firstn], lat[:firstn], ann[:firstn], code[:firstn]
    X, Y = lcc_forward(lon, lat)
    g = read_grid()
    return Inputs(
        X=np.ascontiguousarray(X),
        Y=np.ascontiguousarray(Y),
        emis_annual=np.ascontiguousarray(ann),
        pollutant=np.ascontiguousarray(code),
        W=g["W"], S=g["S"], E=g["E"], N=g["N"],
        TotalPop=g["TotalPop"], MortalityRate=g["MortalityRate"],
        N_REC=int(len(X)), N_SRC=N_SRC,
    )
