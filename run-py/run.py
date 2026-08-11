#!/usr/bin/env python3
# =============================================================================
# run.py — Phase 3 (clean consolidation): the PYTHON binding drives the single
# clean `isrm.esm` end to end through the PUBLIC earthsci_ast surface.
#
#   * `prepare(doc, providers=…, const_arrays=…, pushdown_rewrite=True)` — the
#     automatic projection-pushdown rewrite runs inside the engine; the SR
#     provider gates are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * the SR / grid / pop / mortality providers come FROM THE DOCUMENT
#     (`providers_from_document`: format = `metadata.esio_format`, URL =
#     `source.url_template`);
#   * the EGU FF10 records are read through EarthSciIO's OWN ff10 reader with
#     the document-declared container options (`metadata.x_esd`: the zip
#     `member_filter` glob and `skip_header_row`) — no hand unzip, no hand
#     header-strip; the POLID→code map is the document's `pollutant_codes`;
#   * every reported number is the binding's evaluation of the document's
#     observed graph (`observed_field`) — NO hand-written STEP-0 math here,
#     and NO hand LCC projection (raw emis_lon/emis_lat are the parameters;
#     X/Y are in-model observeds the engine evaluates).
#
#   FULL run  (default)          → assert sum(deathsK/L) ≈ 7524.92 / 16979.63
#   REDUCED   (ISRM_FIRSTN=n)    → first n emission records, totals reported
#
# Emits the cross-language contract record (contract/results_schema.json) with
# model="isrm.esm", mode="runtime_observed_graph", binding="python".
# =============================================================================

from __future__ import annotations

import json
import os
import sys
import time

import paths

paths.ensure_importable()

import numpy as np  # noqa: E402

import earthsciio as eio  # noqa: E402

sys.path.insert(0, os.path.join(paths.REPO, "contract"))
import results as contract  # noqa: E402

T0 = time.time()
N_SRC = 52411  # metaparameter default; SR source/receptor cells
ORACLE_K = 7524.918845602511
ORACLE_L = 16979.632171487083

#: zarr array name -> (emissions observed, concentration observed)
PW_OBS = [
    ("SOA", "E_VOC", "conc_SOA"),
    ("pNO3", "E_NOx", "conc_pNO3"),
    ("pNH4", "E_NH3", "conc_pNH4"),
    ("pSO4", "E_SOx", "conc_pSO4"),
    ("PrimaryPM25", "E_PM25", "conc_PrimaryPM25"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def peak_rss_bytes() -> int:
    """VmHWM (peak resident set) from /proc, in bytes."""
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024
    return 0


# --------------------------------------------------------------------------- #
# EGU FF10 ingest — through EarthSciIO's ff10 reader, options FROM THE DOCUMENT.
#
# The document's EGU_Emis loader declares the zip container (`x_esd.container`),
# the member glob (`x_esd.member_filter`, "*egu*"), the EPA column-header quirk
# (`x_esd.skip_header_row`), and the POLID→pathway-code integer map
# (`x_esd.pollutant_codes`) the is_* observeds compare against. This function
# contributes input PLUMBING only: read the table, map the codes, drop the
# unrecognised/non-finite records (per the document's note), truncate for a
# reduced run. The reader concatenates members in sorted order — identical to
# the zip order for this archive, so ISRM_FIRSTN selects the same records as
# every previous runner.
# --------------------------------------------------------------------------- #
def read_egu(egu_meta: dict, firstn: int | None):
    xe = egu_meta["metadata"]["x_esd"]
    eio.register_format_readers()
    if os.path.isfile(paths.EGU_ZIP):
        url = "file://" + os.path.abspath(paths.EGU_ZIP)
    else:
        # Fall back to the document's declared source, fetched once through the
        # EarthSciIO cache (content-addressed under the cache root).
        url = str(egu_meta["source"]["url_template"])
        log(f"  EGU zip not found at {paths.EGU_ZIP} — fetching {url} via the cache ...")
    loader = eio.DataLoader(
        name="EGU_Emis",
        format=str(egu_meta["metadata"]["esio_format"]),
        url=url,
        variables=["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"],
        reader_kwargs={
            "member_glob": str(xe["member_filter"]),
            "skip_header_row": bool(xe.get("skip_header_row", False)),
        },
    )
    cache = eio.Cache(root=os.path.join(paths.ESIO_CACHE, "EGU_Emis"))
    nds = eio.Provider(loader, cache).materialize()

    codes_map = {str(k).upper(): float(v) for k, v in xe["pollutant_codes"].items()}
    polid = [str(p).strip().upper() for p in nds["POLID"].data]
    code = np.array([codes_map.get(p, 0.0) for p in polid])
    ann = np.asarray(nds["ANN_VALUE"].data, dtype=float)
    lon = np.asarray(nds["LONGITUDE"].data, dtype=float)
    lat = np.asarray(nds["LATITUDE"].data, dtype=float)
    n_all = len(code)
    keep = (code > 0.0) & np.isfinite(lon) & np.isfinite(lat) & np.isfinite(ann)
    lon, lat, ann, code = lon[keep], lat[keep], ann[keep], code[keep]
    log(f"  EGU FF10 records: {n_all} read, {len(lon)} recognised-pathway kept")
    if firstn is not None:
        lon, lat, ann, code = lon[:firstn], lat[:firstn], ann[:firstn], code[:firstn]
    return lon, lat, ann, code


# zarr workaround (unchanged from the validated runners): the SR arrays carry NO
# `.zattrs` object (live GET → 404), so seed an empty `{}` as a cache HIT
# before any SR fetch. A no-op when the blobs are already cached.
def seed_empty_zattrs(cache_root: str, base: str, arrays) -> None:
    cache = eio.Cache(root=cache_root)
    base = base.rstrip("/")
    for arr in arrays:
        key = eio.cache_key(f"{base}/{arr}/.zattrs")
        if cache.store.get_blob(key) is None:
            staged = cache.store.staging_path("json")
            staged.write_text("{}")
            cache.store.put_blob(key, staged)


def main() -> int:
    firstn_env = os.environ.get("ISRM_FIRSTN")
    firstn = int(firstn_env) if firstn_env else None
    reduced = firstn is not None
    os.environ.setdefault("ESS_OBSERVED_PROGRESS", "1")  # per-observed hoist logs

    from earthsci_ast.data_loaders.esio_provider import providers_from_document
    from earthsci_ast.prepare import observed_field, prepare
    from earthsci_ast.simulation_array import BuildInspection
    from earthsci_ast.simulation_loaders import _provider_sample_field

    log(
        f"{'REDUCED' if reduced else 'FULL'} run"
        + (
            f" — first {firstn} emission records"
            if reduced
            else f" — whole domain (target deathsK≈{ORACLE_K:.2f}, deathsL≈{ORACLE_L:.2f})"
        )
    )
    log(f"model:   {paths.MODEL}")
    log(f"scratch: {paths.SCRATCH}")
    log(f"cache:   {paths.ESIO_CACHE}")

    with open(paths.MODEL) as fh:
        doc = json.load(fh)

    # ---- inputs: EGU records via the document-declared FF10 reader options --
    log("reading EGU emissions (EarthSciIO ff10: member_glob + skip_header_row) ...")
    t = time.time()
    lon, lat, ann, code = read_egu(doc["data_loaders"]["EGU_Emis"], firstn)
    n_rec = len(lon)
    t_egu = time.time() - t
    log(f"  N_REC = {n_rec}  in {t_egu:.1f} s")

    # ---- providers FROM THE DOCUMENT ----------------------------------------
    sr_meta = doc["data_loaders"]["ISRM_SR"]["metadata"]["x_esd"]
    sr_arrays = [
        a
        for a in sr_meta["arrays"]
        if a not in ("TotalPop", "MortalityRate", "W", "S", "E", "N")
    ]
    seed_empty_zattrs(
        os.path.join(paths.ESIO_CACHE, "ISRM_SR"),
        doc["data_loaders"]["ISRM_SR"]["source"]["url_template"],
        sr_arrays,
    )
    log("building providers from the document (ISRM_SR) ...")
    providers = providers_from_document(
        doc, cache_root=paths.ESIO_CACHE, loaders=["ISRM_SR"]
    )
    log(f"  providers: {sorted(providers)}")

    # ---- src-cell rectangles: the [0:N_SRC] prefix of the full grid ---------
    # The document declares src_W/src_S/src_E/src_N ([src_cells]) alongside the
    # full-grid W/S/E/N ([pop_cells]); the loader note pins the prefix
    # relationship. Slicing is input plumbing, not model math.
    log("materializing grid prefix (src_W/S/E/N) ...")
    t = time.time()
    ca: dict[str, np.ndarray] = {
        f"src_{k}": np.asarray(
            _provider_sample_field(providers[f"ISRM_SR.{k}"], 0.0), dtype=float
        )[:N_SRC]
        for k in ("W", "S", "E", "N")
    }
    t_grid = time.time() - t
    log(f"  grid prefix in {t_grid:.1f} s")

    # ---- const arrays: the EGU loader's variables, src rects ----------------
    ca["EGU_Emis.lon"] = lon
    ca["EGU_Emis.lat"] = lat
    ca["EGU_Emis.annual"] = ann
    ca["EGU_Emis.pollutant"] = code
    for v in ("stkhgt", "stkdiam", "stktemp", "stkvel"):
        ca[f"EGU_Emis.{v}"] = np.zeros(n_rec)  # declared, unused downstream

    # ---- PREPARE (rewrite → VI → gated SR fetch → observed-graph hoist) -----
    log(f"prepare(pushdown_rewrite=True) — N_REC={n_rec} ...")
    insp = BuildInspection()
    t = time.time()
    prep = prepare(
        doc,
        providers=providers,
        const_arrays=ca,
        inspect=insp,
        metaparameters={"N_REC": n_rec},
        pushdown_rewrite=True,
    )
    t_prep = time.time() - t
    log(
        f"PREPARE done in {t_prep:.1f} s  (peak RSS so far: "
        f"{peak_rss_bytes() / 2**30:.2f} GiB)"
    )

    # ---- the engine-derived support set (for the contract record) -----------
    mf_keys = [k for k in insp.const_arrays if str(k).startswith("pd_member_factor__")]
    if not mf_keys:
        raise RuntimeError("no pd_member_factor__* const array — did the rewrite fire?")
    members = sorted(int(m) for m in np.asarray(insp.const_arrays[mf_keys[0]]))
    n_ppl = len(members)
    log(f"engine-derived support set: |members| = {n_ppl} of {N_SRC} source cells")
    if not reduced and n_ppl != 1520:
        log("  WARNING: expected 1520 emission-bearing cells at full scale")

    # ---- results through the prepared document's own graph ------------------
    # (already evaluated by the const-geometry hoist inside prepare; the
    # per-observed timings are in the [ess-observed] log lines above)
    t = time.time()
    dK = np.asarray(observed_field(prep, "deathsK"), dtype=float)
    dL = np.asarray(observed_field(prep, "deathsL"), dtype=float)
    tp = np.asarray(observed_field(prep, "TotalPM25"), dtype=float)
    pathways = {}
    for arr, evar, cvar in PW_OBS:
        ep = np.asarray(observed_field(prep, evar), dtype=float)
        cp = np.asarray(observed_field(prep, cvar), dtype=float)
        pathways[arr] = {
            "emis_sum": float(ep.sum()),
            "conc_sum": float(cp.sum()),
            "conc_max": float(cp.max()),
        }
    t_eval = time.time() - t
    log(f"field readback in {t_eval:.1f} s")

    sK, sL = float(dK.sum()), float(dL.sum())
    log("\n" + "=" * 70)
    log(f"  sum(deathsK) = {sK!r}")
    log(f"  sum(deathsL) = {sL!r}")
    log(f"  Σ TotalPM25  = {float(tp.sum())!r}")
    ok_k = abs(sK - ORACLE_K) <= 1e-4 * ORACLE_K
    ok_l = abs(sL - ORACLE_L) <= 1e-4 * ORACLE_L
    if not reduced:
        log(f"  target deathsK={ORACLE_K}  rel.err {100 * (sK - ORACLE_K) / ORACLE_K:.6f}%")
        log(f"  target deathsL={ORACLE_L} rel.err {100 * (sL - ORACLE_L) / ORACLE_L:.6f}%")
        log(f"PHASE 3 FULL: {'PASS' if ok_k and ok_l else 'FAIL'}")
    log("=" * 70)

    # ---- contract record ----------------------------------------------------
    import earthsci_ast

    out = os.path.join(
        paths.RUNPY_DIR, "results_reduced.json" if reduced else "results.json"
    )
    contract.write_results(
        out,
        binding_version=(
            f"python {sys.version.split()[0]} / earthsci_ast "
            f"{getattr(earthsci_ast, '__version__', '?')}"
        ),
        model=paths.MODEL,
        mode="runtime_observed_graph",
        n_src=N_SRC,
        n_rcv=N_SRC,
        n_rec=n_rec,
        ppl=members,
        pathways=pathways,
        total_pm25=tp.tolist(),
        deathsK=dK.tolist(),
        deathsL=dL.tolist(),
        timing={
            "wall_seconds": time.time() - T0,
            "egu_seconds": t_egu,
            "grid_seconds": t_grid,
            "prepare_seconds": t_prep,
            "peak_rss_bytes": peak_rss_bytes(),
        },
    )
    if not reduced and not (ok_k and ok_l):
        raise SystemExit("full-scale totals off oracle (see above)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
