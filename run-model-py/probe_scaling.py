#!/usr/bin/env python3
"""How does evaluating conc_* / E_* scale with |ppl|?

The reduced run spent 3.0 s on a 9 x 52,411 contraction that einsum should do in
milliseconds, so before committing an hour to a full run, measure the growth
directly on the REAL .esm expressions with synthetic factors. Isolated and cheap
(a few hundred MB), so it is safe to run alongside another job.

    python probe_scaling.py 9 50 200
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np

from earthsci_ast.numpy_interpreter import EvalContext, eval_expr
from earthsci_ast.parse import _parse_expression

from paths import MODEL, N_SRC
from run_model import model_view, observed_defs, observed_order, resolve_sizes, scalar_params

N_RCV = N_SRC


def build_ctx(n_ppl: int, n_rec: int):
    with open(MODEL) as fh:
        doc = json.load(fh)
    doc = resolve_sizes(doc, {"N_SRC": N_SRC, "N_RCV": N_RCV, "N_POP": N_SRC,
                              "N_LAYER": 3, "N_REC": n_rec})
    model = model_view(doc)
    params = scalar_params(model)
    rng = np.random.default_rng(12345)
    ca = {
        "X": rng.uniform(-2e6, 2e6, n_rec), "Y": rng.uniform(-2e6, 1e6, n_rec),
        "emis_annual": rng.uniform(0, 100, n_rec),
        "pollutant": rng.choice([1.0, 36.0, 40.0, 41.0, 42.0], n_rec),
        "W": rng.uniform(-2e6, 2e6, N_SRC), "S": rng.uniform(-2e6, 1e6, N_SRC),
        "E": rng.uniform(-2e6, 2e6, N_SRC), "N": rng.uniform(-2e6, 1e6, N_SRC),
        "TotalPop": rng.uniform(0, 1e4, N_SRC),
        "MortalityRate": rng.uniform(0, 1e3, N_SRC),
        "src_cell_of_ppl": np.arange(1, n_ppl + 1, dtype=float),
    }
    for v in ("SR_SOA", "SR_pNO3", "SR_pNH4", "SR_pSO4", "SR_PrimaryPM25"):
        ca[v] = rng.uniform(0, 1e-9, (n_ppl, N_RCV))
    ctx = EvalContext(
        state_layout={}, state_shapes={}, param_values=dict(params),
        observed_values={}, y=np.empty((0,), dtype=float), t=0.0,
        index_sets=dict(doc.get("index_sets") or {}),
        derived_extents={"emis_src_cells_faq": n_ppl},
        input_arrays=dict(ca),
    )
    return ctx, observed_defs(model)


def main() -> int:
    sizes = [int(a) for a in sys.argv[1:]] or [9, 50, 200]
    n_rec = int(sys.argv[0] and 2000)
    print(f"n_rcv={N_RCV}  n_rec={n_rec}  (synthetic factors; real .esm expressions)")
    watch = ["E_VOC", "conc_SOA", "deathsK"]
    prev: dict[str, float] = {}
    for n_ppl in sizes:
        ctx, defs = build_ctx(n_ppl, n_rec)
        order = observed_order(defs)
        times: dict[str, float] = {}
        for name in order:
            t = time.time()
            ctx.input_arrays[name] = np.asarray(
                eval_expr(_parse_expression(defs[name]), ctx), dtype=float
            )
            times[name] = time.time() - t
        total = sum(times.values())
        parts = []
        for w in watch:
            g = f" ({times[w] / prev[w]:.1f}x)" if w in prev and prev[w] > 0 else ""
            parts.append(f"{w}={times[w]:.2f}s{g}")
        print(f"  |ppl|={n_ppl:5d}  total={total:6.2f}s   " + "  ".join(parts))
        prev = times
    return 0


if __name__ == "__main__":
    sys.exit(main())
