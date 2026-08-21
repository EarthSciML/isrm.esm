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
#   FULL run  (default)          → report sum(deathsK/L) against the tutorial
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
# The InMAP source-receptor tutorial's published national totals
# (https://inmap.run/blog/2019/04/20/sr/), which account for plume rise.
# A TARGET, not an assertion: the document deliberately does not reproduce
# InMAP's high-plume source-index defect (a plume above model layer 7 keeps an
# index built in the coarse 9324-cell grid, then read against the 52411-cell
# ground grid), which misplaces 654 of 43650 records - 0.43% of emitted mass -
# onto the wrong source cell. So a run lands NEAR rather than ON these, and the
# deviation is printed rather than failed.
ORACLE_K = 6928.959583
ORACLE_L = 15623.924632
#: Beyond this the deviation is more than the clean-physics choice can explain.
#: MEASURED at full scale on 2026-08-19: deathsK 6983.9385617781645 (+0.79%)
#: and deathsL 15752.315804140908 (+0.82%) against the published totals. The
#: misplaced group is 0.43% of emitted mass but buys about twice that in
#: deaths, because putting it back on the cells the emissions came from puts it
#: back over people. This threshold sits just above what was measured.
#: STALE, and left in place only as a loose upper bound: it was fitted to a
#: full-scale run that predates sr.Reader.layerFracs and the [0, 3, 6] fix, and
#: no full-scale run of the current document has been made. At reduced scale the
#: document now reproduces the live inmap cloud service to 9e-9, so the real
#: threshold is expected to be far tighter — see SERVICE_DEATHS in
#: contract/compare_results.py.
ORACLE_NOTABLE_REL = 8.3e-3

#: zarr array name -> (the per-SR-layer emissions observeds, concentration
#: observed). The emissions side grew a LAYER dimension when the document
#: started stating plume rise: a record is charged to the SR layer its plume
#: reaches, so a pathway's emissions are three arrays, not one. The
#: concentration side did not - ``conc_<p>`` is the document's own sum over the
#: three layered contractions.
PW_OBS = [
    ("SOA", ["E_VOC_L0", "E_VOC_L1", "E_VOC_L2"], "conc_SOA"),
    ("pNO3", ["E_NOx_L0", "E_NOx_L1", "E_NOx_L2"], "conc_pNO3"),
    ("pNH4", ["E_NH3_L0", "E_NH3_L1", "E_NH3_L2"], "conc_pNH4"),
    ("pSO4", ["E_SOx_L0", "E_SOx_L1", "E_SOx_L2"], "conc_pSO4"),
    ("PrimaryPM25", ["E_PM25_L0", "E_PM25_L1", "E_PM25_L2"], "conc_PrimaryPM25"),
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
    """The data sources that DISCOVER their own extent (``extent.metaparameter``) —
    the record-bearing tables of the document, whatever they happen to be
    called. The two knobs below are scale/locality concerns of a RUN, not of
    the model, and both are expressed in the document's own vocabulary."""
    return [
        name
        for name, ld in (doc.get("data_sources") or {}).items()
        if isinstance(ld, dict) and isinstance(ld.get("extent"), dict)
        and ld["extent"].get("metaparameter")
    ]


def gated_loader_arrays(doc: dict) -> dict:
    """The sources that declare a ``gated_select``, and the arrays each gates -
    the document's own statement of how many model arrays the pushdown rewrite
    must end up gating. Derived rather than written down, so splitting or
    merging a gated source keeps the check honest."""
    out: dict[str, list[str]] = {}
    for name, ld in (doc.get("data_sources") or {}).items():
        gs = ((ld or {}).get("metadata") or {}).get("x_esd", {}).get("gated_select")
        if isinstance(gs, dict):
            out[name] = [str(a) for a in gs.get("applies_to", [])]
    return out


def truncate_records(doc: dict, n: int) -> None:
    """REDUCED runs: truncate every record-discovering loader to its first ``n``
    DELIVERED records with a loader-level ``select`` range (esm-spec §8.9.2).
    Because the selection follows the loader's own ``record_filter``, this picks
    the same records the previous runners' post-filter ``[:n]`` did — and
    ``extent`` then re-discovers the smaller N_REC by itself."""
    for name in record_loaders(doc):
        doc["data_sources"][name]["select"] = {"axes": [{"range": {"start": 0, "stop": n}}]}


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

    from earthsci_ast.data_sources.esio_provider import providers_from_document
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
    # Only the GATED loaders' arrays: those are the SR slabs, and they are the
    # ones the store has no `.zattrs` for. The 1-D grid arrays are read whole
    # through the ordinary path and must keep whatever attrs the store has.
    for lname, arrs in gated_loader_arrays(doc).items():
        seed_empty_zattrs(
            os.path.join(paths.ESIO_CACHE, lname),
            doc["data_sources"][lname]["source"]["url_template"],
            arrs,
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
    t_prep = time.time() - t

    # ---- the gate covers EVERY declared SR array ---------------------------
    # A malformed E_* or conc_* body does not fail: the pathway simply drops
    # out of the rewrite's `applies_to` list, the rest of the rewrite reports
    # success, and the un-gated array is then fetched WHOLE - 330 GB, which
    # surfaces hours later as a memory failure rather than an error. The
    # document says how many arrays it declared for gating; anything less here
    # is that silent drop, so stop on it.
    expect_gated = sum(len(v) for v in gated_loader_arrays(doc).values())
    gated = [
        str(a)
        for a in ((prep.doc or {}).get("metadata", {}).get("x_esd", {})
                  .get("pushdown", {}).get("gated_select", {})
                  .get("applies_to", []))
    ]
    log(f"gated arrays: {len(gated)} of {expect_gated} declared")
    if len(gated) != expect_gated:
        raise SystemExit(
            f"the pushdown rewrite gated {len(gated)} arrays but the document "
            f"declares {expect_gated} ({gated}) - a pathway dropped out of the "
            "gate silently, and its SR array would be fetched UNGATED. Check "
            "that each conc_*_L* body is a plain two-factor SR*E product and "
            "that the containment ifelse is the FIRST ifelse in every E_* body."
        )

    n_rec = int(np.asarray(observed_field(prep, "X")).size)
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
    emis_by_layer = {}
    for arr, evars, cvar in PW_OBS:
        # A pathway's emissions are now spread over the three SR layers; the
        # record's `emis_sum` is the pathway TOTAL, so it stays comparable to
        # the ground-level-only baselines (plume rise moves mass between
        # layers, never into or out of a pathway).
        es = [np.asarray(observed_field(prep, v), dtype=float).ravel() for v in evars]
        ep = np.concatenate(es)
        cp = np.asarray(observed_field(prep, cvar), dtype=float)
        pathways[arr] = {
            "emis_sum": float(ep.sum()),
            "conc_sum": float(cp.sum()),
            "conc_max": float(cp.max()),
        }
        # How much mass plume rise put in each SR layer — the physics made
        # visible as tons, per pathway.
        emis_by_layer[arr] = [float(e.sum()) for e in es]
    # The layer assignment itself — now a SPLIT, not a single layer: InMAP's
    # sr.Reader.layerFracs charges a record to two SR layers whenever its model
    # layer falls between two entries of `layers`. `sr_lower` is the lower index
    # (integer, compared exactly) and w_sr0/1/2 are the three shares. These are
    # the document's OWN observeds, read through the same `observed_field` path
    # as everything else — this runner does not know what ASME is, and must not:
    # the point of the contract's `plume` block is that the ENGINE produced the
    # assignment from the spec. contract/plume_oracle.py computes the same
    # quantity independently, from the meteorology arrays and without the SR
    # matrix, and compare_results.py checks the two against each other.
    def per_record(name):
        return np.asarray(observed_field(prep, name), dtype=float).ravel()

    plume = contract.plume_block(
        sr_lower=per_record("sr_lower"),
        stack_layer=per_record("stack_layer"),
        weights={w: per_record(w) for w in ("w_sr0", "w_sr1", "w_sr2")},
        emis_by_sr_layer=emis_by_layer,
    )
    t_eval = time.time() - t
    log(f"field readback in {t_eval:.1f} s")

    sK, sL = float(dK.sum()), float(dL.sum())
    log("\n" + "=" * 70)
    log(f"  sum(deathsK) = {sK!r}")
    log(f"  sum(deathsL) = {sL!r}")
    log(f"  Σ TotalPM25  = {float(tp.sum())!r}")
    if not reduced:
        rk = (sK - ORACLE_K) / ORACLE_K
        rl = (sL - ORACLE_L) / ORACLE_L
        log(f"  tutorial deathsK={ORACLE_K}  deviation {100 * rk:.6f}%")
        log(f"  tutorial deathsL={ORACLE_L} deviation {100 * rl:.6f}%")
        if abs(rk) > ORACLE_NOTABLE_REL or abs(rl) > ORACLE_NOTABLE_REL:
            log(
                f"  WARNING: deviation exceeds {100 * ORACLE_NOTABLE_REL:.2f}% - "
                "more than the above-layer-7 group has been measured to be worth "
                "(0.43% of emitted mass, +0.79%/+0.82% of deaths at full scale), "
                "so something else differs."
            )
    log("  lower-SR-layer histogram (records per layer 0/1/2) = "
        f"{plume['sr_lower']['histogram']}")
    log(f"  sr_lower sha256 = {plume['sr_lower']['sha256']}")
    log("  Σ w_sr0/w_sr1/w_sr2 = "
        + " / ".join(repr(plume["weights"][w]["sum"])
                     for w in ("w_sr0", "w_sr1", "w_sr2"))
        + f"   max|Σw - 1| = {plume['weights']['max_sum_error']!r}")
    log(
        "    (check it against `python3 contract/plume_oracle.py"
        + (f" --firstn {firstn}`" if reduced else "`")
        + " — no SR matrix needed)"
    )
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
        plume=plume,
        timing={
            "wall_seconds": time.time() - T0,
            "providers_seconds": t_providers,
            "prepare_seconds": t_prep,
            "peak_rss_bytes": peak_rss_bytes(),
        },
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
