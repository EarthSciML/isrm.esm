#!/usr/bin/env python3
# =============================================================================
# run.py — the PYTHON binding drives the single clean `isrm.esm` end to end
# through the PUBLIC earthsci_ast surface. NOTHING MODEL-SHAPED LIVES HERE:
# this file names no pollutant, no column, no grid extent and no record count.
#
#   * `prepare(doc, providers=…, pushdown_rewrite=True)` — the automatic
#     projection-pushdown rewrite runs inside the engine; the SR provider gates
#     are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * EVERY provider comes FROM THE DOCUMENT (`providers_from_document`:
#     format = `metadata.esio_format`, URL = `source.url_template`) — the SR
#     slabs, the grid, the population, AND the EGU FF10 table, whose ingest the
#     loader now declares in full (esm-spec §8.9): `reader_options` (the zip
#     member glob + header row), `codes` (POLID text -> the pathway enum, an
#     unrecognised code dropping the record), `record_filter` (no coordinate /
#     no annual total is not a record) and `extent` (the surviving count binds
#     N_REC). The src-cell rectangles are the `select` range `W[0:N_SRC]` on
#     their own loader variables;
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
    """Peak resident set size, in bytes.

    ``VmHWM`` from ``/proc/self/status`` is the number the cross-language table
    reports (and the one this repo's notes insist on over the cgroup's, which
    the page cache saturates). ``getrusage`` is the portable fallback for a
    machine with no procfs — ``ru_maxrss`` is kibibytes on Linux and BYTES on
    macOS/BSD, hence the platform branch.
    """
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak) if sys.platform == "darwin" else int(peak) * 1024


def metaparam(doc: dict, name: str) -> int:
    """A metaparameter's declared default, read from the document (so no grid
    extent is written down here)."""
    mp = doc.get("metaparameters", {}).get(name, {})
    return int(mp.get("default", 0))


def record_loaders(doc: dict) -> list[str]:
    """The loaders that DISCOVER their own extent (``extent.metaparameter``) —
    the record-bearing tables of the document, whatever they happen to be
    called. The two knobs below are scale/locality concerns of a RUN, not of
    the model, and both are expressed in the document's own vocabulary."""
    return [
        name
        for name, ld in (doc.get("data_loaders") or {}).items()
        if isinstance(ld, dict) and isinstance(ld.get("extent"), dict)
        and ld["extent"].get("metaparameter")
    ]


def truncate_records(doc: dict, n: int) -> None:
    """REDUCED runs: truncate every record-discovering loader to its first ``n``
    DELIVERED records with a loader-level ``select`` range (esm-spec §8.9.2).
    Because the selection follows the loader's own ``record_filter``, this picks
    the same records the previous runners' post-filter ``[:n]`` did — and
    ``extent`` then re-discovers the smaller N_REC by itself."""
    for name in record_loaders(doc):
        doc["data_loaders"][name]["select"] = {"axes": [{"range": {"start": 0, "stop": n}}]}


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
    if reduced:
        truncate_records(doc, firstn)

    # ---- providers FROM THE DOCUMENT — ALL of them ---------------------------
    # Including the FF10 table: the loader declares its own reader options, code
    # map, record filter and extent, so there is nothing left here to read, map,
    # filter or count.
    log("building providers from the document ...")
    t = time.time()
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
    # A local copy of a record loader's source is a LOCALITY choice of this run
    # (gaftp.epa.gov is slow and flaky), so it is a url_override rather than an
    # edit to the document.
    url_overrides: dict[str, str] = {}
    if os.path.isfile(paths.EGU_ZIP):
        mirror = "file://" + os.path.abspath(paths.EGU_ZIP)
        url_overrides = {name: mirror for name in record_loaders(doc)}
        log(f"  record source mirrored from {paths.EGU_ZIP}")
    providers = providers_from_document(
        doc, cache_root=paths.ESIO_CACHE, url_overrides=url_overrides
    )
    log(f"  providers: {sorted(providers)}")
    t_providers = time.time() - t

    # ---- PREPARE (extent → rewrite → VI → gated fetch → observed-graph) -----
    log("prepare(pushdown_rewrite=True) — N_REC discovered by the loader ...")
    insp = BuildInspection()
    t = time.time()
    prep = prepare(
        doc,
        providers=providers,
        inspect=insp,
        pushdown_rewrite=True,
    )
    n_rec = int(np.asarray(observed_field(prep, "X")).size)
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
    n_src, n_rcv = metaparam(doc, "N_SRC"), metaparam(doc, "N_RCV")
    log(f"engine-derived support set: |members| = {n_ppl} of {n_src} source cells")
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
        n_src=n_src,
        n_rcv=n_rcv,
        n_rec=n_rec,
        ppl=members,
        pathways=pathways,
        total_pm25=tp.tolist(),
        deathsK=dK.tolist(),
        deathsL=dL.tolist(),
        timing={
            "wall_seconds": time.time() - T0,
            "providers_seconds": t_providers,
            "prepare_seconds": t_prep,
            "peak_rss_bytes": peak_rss_bytes(),
        },
    )
    if not reduced and not (ok_k and ok_l):
        raise SystemExit("full-scale totals off oracle (see above)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
