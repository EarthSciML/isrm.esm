#!/usr/bin/env python3
"""ISRM through the PYTHON binding's own evaluation of isrm_pushdown.esm.

The strict bar for this project: every number in the result record comes from
the binding's own evaluation of the .esm graph — no hand-written STEP-0
arithmetic anywhere in this file. Concretely:

  * the emission-bearing source-cell support set (`ppl`) is DERIVED by
    ``earthsci_ast.value_invention`` running the model's own producer aggregate,
    gated by the spatial `join.overlap` broad phase (CONFORMANCE_SPEC §5.5.6);
  * every field — the per-source binned emissions E_*, the per-receptor
    concentrations conc_*, TotalPM25, deathsK/deathsL — is produced by
    ``earthsci_ast.numpy_interpreter.eval_expr`` on the .esm's own observed
    expressions;
  * this file contributes ORCHESTRATION only: which observed to evaluate when,
    and where the bytes come from.

The one impure boundary is I/O, which is what EarthSciIO is for: the emission
records and the grid come from ``inputs.py``, and the SR slab is fetched with a
projection pushdown built from the members value-invention just derived
(``SR[layer 0, ppl, :]`` — 1,520 rows, not 52,411).

Usage
-----
    python run_model.py                 # full run -> results.json
    ISRM_FIRSTN=200 python run_model.py # reduced run -> results_reduced.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

import earthsciio as eio
from earthsci_ast.numpy_interpreter import EvalContext, eval_expr
from earthsci_ast.parse import _parse_expression
from earthsci_ast.value_invention import materialize_value_invention

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "contract"))
import results as contract  # noqa: E402

from inputs import build_inputs  # noqa: E402
from paths import MODEL, N_SRC, SCRATCH, ZARR_URL  # noqa: E402

#: zarr array name -> (model SR variable, emissions observed, concentration observed)
PATHWAYS = [
    ("SOA", "SR_SOA", "E_VOC", "conc_SOA"),
    ("pNO3", "SR_pNO3", "E_NOx", "conc_pNO3"),
    ("pNH4", "SR_pNH4", "E_NH3", "conc_pNH4"),
    ("pSO4", "SR_pSO4", "E_SOx", "conc_pSO4"),
    ("PrimaryPM25", "SR_PrimaryPM25", "E_PM25", "conc_PrimaryPM25"),
]

ORACLE_K = 7524.918845602511
ORACLE_L = 16979.632171487083


def log(msg: str) -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------- #
# Document preparation
# --------------------------------------------------------------------------- #


def resolve_sizes(doc: dict, mp: dict[str, int]) -> dict:
    """Replace string index-set sizes ("N_SRC", ...) with integers — the Python
    mirror of l3_common.jl's ``resolve_sizes!``. A `derived` set has no size:
    value invention gives it one."""
    for iset in (doc.get("index_sets") or {}).values():
        if isinstance(iset, dict) and isinstance(iset.get("size"), str):
            name = iset["size"]
            if name in mp:
                iset["size"] = int(mp[name])
    return doc


def model_view(doc: dict, name: str = "ISRM") -> dict:
    """The model dict carrying the document-scoped index-set registry, which is
    the shape ``materialize_value_invention`` consumes."""
    model = dict(doc["models"][name])
    model.setdefault("index_sets", doc.get("index_sets", {}))
    return model


def scalar_params(model: dict) -> dict[str, float]:
    """Every declared 0-D parameter's default — `fact`, the risk ratios, the LCC
    constants. These are load-time constants the graph reads directly."""
    out: dict[str, float] = {}
    for n, v in (model.get("variables") or {}).items():
        if v.get("type") == "parameter" and not v.get("shape") and v.get("default") is not None:
            out[str(n)] = float(v["default"])
    return out


def observed_defs(model: dict) -> dict[str, dict]:
    """name -> raw expression dict, for every `observed` variable."""
    return {
        str(n): v["expression"]
        for n, v in (model.get("variables") or {}).items()
        if v.get("type") == "observed" and v.get("expression") is not None
    }


def observed_order(defs: dict[str, dict]) -> list[str]:
    """Dependency order over the observeds.

    Each aggregate declares the factors it reads in its own ``args`` list, so the
    edges come from the MODEL, not from a re-derivation here: an observed follows
    every observed named in its args. A cycle is a malformed model and is a hard
    error rather than an arbitrary order.
    """
    pending = dict(defs)
    ordered: list[str] = []
    done: set[str] = set()
    while pending:
        ready = [
            n
            for n, e in pending.items()
            if all(a in done for a in (e.get("args") or []) if a in defs and a != n)
        ]
        if not ready:
            raise RuntimeError(f"cyclic observed dependency among {sorted(pending)}")
        for n in sorted(ready):
            ordered.append(n)
            done.add(n)
            del pending[n]
    return ordered


# --------------------------------------------------------------------------- #
# The gated SR fetch (projection pushdown, §5.5.6 "gated_select")
# --------------------------------------------------------------------------- #


def fetch_sr(members1: np.ndarray, evict: bool = True) -> dict[str, np.ndarray]:
    """Fetch ``SR[layer 0, members, :]`` for each pathway.

    ``members1`` are the 1-based full-grid source-cell ids value invention just
    derived; EarthSciIO indexes 0-based, and the selection goes down to the zarr
    reader so only the intersecting chunks are fetched. Pathways are fetched and
    EVICTED one at a time so peak disk is one pathway (~2.7 GiB), not five.
    """
    root = os.path.join(SCRATCH, os.environ.get("ISRM_SR_DIR", "py_l3_cache_sr"))
    idx0 = [int(m) - 1 for m in members1]
    select = {"axes": [{"indices": [0]}, {"indices": idx0}, "all"]}
    out: dict[str, np.ndarray] = {}
    for zname, mname, _, _ in PATHWAYS:
        if evict:
            blobs = os.path.join(root, "v1", "blobs")
            if os.path.isdir(blobs):
                import shutil

                shutil.rmtree(blobs, ignore_errors=True)
        cache = eio.Cache(root=root)
        loader = eio.DataLoader(
            name=f"ISRM_SR[{zname}]", format="zarr", url=ZARR_URL, variables=[zname]
        )
        t = time.time()
        ds = eio.Provider(loader, cache).materialize(select=select)
        # (1, |ppl|, rcv) -> (|ppl|, rcv): the model's SR is [emis_src_cells, rcv_cells];
        # emisLayer 0 is already sliced out by the selection.
        arr = np.asarray(ds[zname].data, dtype=float)[0]
        out[mname] = np.ascontiguousarray(arr)
        log(f"    [gated] fetched {zname} -> {arr.shape}  in {time.time() - t:.1f} s")
    return out


# --------------------------------------------------------------------------- #


def main() -> int:
    t0 = time.time()
    firstn_env = os.environ.get("ISRM_FIRSTN")
    firstn = int(firstn_env) if firstn_env else None
    reduced = firstn is not None
    log(
        f"{'REDUCED' if reduced else 'FULL'} run"
        + (f" — first {firstn} emission records" if reduced else
           f" — whole domain (target deathsK≈{ORACLE_K:.2f}, deathsL≈{ORACLE_L:.2f})")
    )

    log("building inputs ...")
    inp = build_inputs(firstn=firstn)
    log(f"  N_REC={inp.N_REC}  N_SRC={inp.N_SRC}")

    with open(MODEL) as fh:
        doc = json.load(fh)
    doc = resolve_sizes(
        doc,
        {"N_SRC": inp.N_SRC, "N_RCV": inp.N_SRC, "N_POP": inp.N_SRC,
         "N_LAYER": 3, "N_REC": inp.N_REC},
    )
    model = model_view(doc)
    params = scalar_params(model)

    # The .esm's declared parameter arrays. Nothing derived lives here.
    const_arrays: dict[str, np.ndarray] = {
        "X": inp.X, "Y": inp.Y,
        "emis_annual": inp.emis_annual, "pollutant": inp.pollutant,
        "W": inp.W, "S": inp.S, "E": inp.E, "N": inp.N,
        "TotalPop": inp.TotalPop, "MortalityRate": inp.MortalityRate,
    }

    # ---- VALUE INVENTION: the graph derives its own support set -------------
    log("value invention (spatial overlap gate -> emis_src_cells) ...")
    t = time.time()
    vi = materialize_value_invention(model, const_arrays=const_arrays, params=params)
    faq = "emis_src_cells_faq"
    members = np.asarray(vi.members[faq], dtype=float)  # 1-based full-grid cell ids
    n_ppl = len(members)
    log(f"  |emis_src_cells| = {n_ppl}   in {time.time() - t:.1f} s")
    if not reduced and n_ppl != 1520:
        log(f"  WARNING: expected 1520 emission-bearing cells at full scale, got {n_ppl}")

    # Hook 1: the derived set's `member_factor` is fed back as a const factor so
    # the downstream cell_W/cell_S/... gathers read the compact derived axis.
    mfactor = (doc["index_sets"]["emis_src_cells"] or {}).get("member_factor")
    if mfactor:
        const_arrays[str(mfactor)] = members

    # ---- Hook 2: the gated fetch, driven by the members just derived --------
    log(f"gated SR fetch: layer 0, {n_ppl} of {inp.N_SRC} source cells, all receptors")
    const_arrays.update(fetch_sr(members))

    # ---- Evaluate the observed graph ---------------------------------------
    defs = observed_defs(model)
    order = observed_order(defs)
    ctx = EvalContext(
        state_layout={}, state_shapes={}, param_values=dict(params),
        observed_values={}, y=np.empty((0,), dtype=float), t=0.0,
        index_sets=dict(doc.get("index_sets") or {}),
        derived_extents={k: int(v) for k, v in vi.extents.items()},
        input_arrays=dict(const_arrays),
    )

    fields: dict[str, np.ndarray] = {}
    log(f"evaluating {len(order)} observeds through the .esm graph ...")
    for name in order:
        t = time.time()
        val = eval_expr(_parse_expression(defs[name]), ctx)
        arr = np.asarray(val, dtype=float)
        fields[name] = arr
        # Feed each result back so downstream observeds resolve it by name — the
        # same dependency-ordered substitution the array/PDE build performs.
        ctx.input_arrays[name] = arr
        log(f"  {name:18s} shape={arr.shape}  {time.time() - t:6.1f} s")

    dK, dL, tp = fields["deathsK"], fields["deathsL"], fields["TotalPM25"]
    sK, sL = float(dK.sum()), float(dL.sum())
    log("\n" + "=" * 70)
    log(f"  sum(deathsK) = {sK!r}")
    log(f"  sum(deathsL) = {sL!r}")
    log(f"  Σ TotalPM25  = {float(tp.sum())!r}")

    ok = True
    if not reduced:
        ok = (abs(sK - ORACLE_K) <= 1e-4 * ORACLE_K) and (abs(sL - ORACLE_L) <= 1e-4 * ORACLE_L)
        log(f"  target deathsK={ORACLE_K}  rel.err {100 * (sK - ORACLE_K) / ORACLE_K:.4f}%")
        log(f"  target deathsL={ORACLE_L} rel.err {100 * (sL - ORACLE_L) / ORACLE_L:.4f}%")
        log(f"PYTHON FULL: {'PASS' if ok else 'FAIL'}")
    log("=" * 70)

    pathways = {
        zname: {
            "emis_sum": float(fields[evar].sum()),
            "conc_sum": float(fields[cvar].sum()),
            "conc_max": float(fields[cvar].max()),
        }
        for zname, _, evar, cvar in PATHWAYS
    }
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "results_reduced.json" if reduced else "results.json")
    import earthsci_ast

    contract.write_results(
        out,
        binding_version=f"python {sys.version.split()[0]} / earthsci_ast "
                        f"{getattr(earthsci_ast, '__version__', '?')}",
        model=MODEL, mode="runtime_observed_graph",
        n_src=inp.N_SRC, n_rcv=inp.N_SRC, n_rec=inp.N_REC,
        ppl=[int(m) for m in members],
        pathways=pathways,
        total_pm25=tp.tolist(), deathsK=dK.tolist(), deathsL=dL.tolist(),
        timing={"wall_seconds": time.time() - t0},
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
