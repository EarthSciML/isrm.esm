#!/usr/bin/env python3
"""Full-scale value invention ONLY — no SR fetch, no observed evaluation.

De-risks the long full run: the overlap broad phase is 43,650 query points x
52,411 cell rectangles, and the derived support set must come out at exactly
1,520 to match Julia. Cheap (seconds, a few hundred MB), so it is safe to run
next to another job.
"""
import json, time
import numpy as np
from earthsci_ast.value_invention import materialize_value_invention, _VI_ENUM_VISITS
from earthsci_ast import broad_phase
from inputs import build_inputs
from paths import MODEL
from run_model import model_view, resolve_sizes, scalar_params

t0 = time.time()
inp = build_inputs()
print(f"inputs: N_REC={inp.N_REC} N_SRC={inp.N_SRC}  ({time.time()-t0:.1f} s)", flush=True)

with open(MODEL) as fh:
    doc = json.load(fh)
doc = resolve_sizes(doc, {"N_SRC": inp.N_SRC, "N_RCV": inp.N_SRC, "N_POP": inp.N_SRC,
                          "N_LAYER": 3, "N_REC": inp.N_REC})
model = model_view(doc)
ca = {"X": inp.X, "Y": inp.Y, "emis_annual": inp.emis_annual, "pollutant": inp.pollutant,
      "W": inp.W, "S": inp.S, "E": inp.E, "N": inp.N,
      "TotalPop": inp.TotalPop, "MortalityRate": inp.MortalityRate}

t = time.time()
cands = broad_phase.broad_phase_candidates(
    broad_phase.envelope_vectors(["X", "Y"], ca),
    broad_phase.envelope_vectors(["W", "S", "E", "N"], ca))
print(f"broad phase: {len(cands)} candidate pairs from "
      f"{inp.N_REC}x{inp.N_SRC}={inp.N_REC*inp.N_SRC:,} "
      f"({time.time()-t:.1f} s)", flush=True)

_VI_ENUM_VISITS[0] = 0
t = time.time()
vi = materialize_value_invention(model, const_arrays=ca, params=scalar_params(model))
mem = vi.members["emis_src_cells_faq"]
print(f"value invention: |ppl|={len(mem)}  visits={_VI_ENUM_VISITS[0]}  "
      f"({time.time()-t:.1f} s)", flush=True)
print(f"  expected 1520 -> {'MATCH' if len(mem)==1520 else 'MISMATCH'}")
print(f"  first/last members: {mem[:5]} ... {mem[-5:]}")
np.save("ppl_members.npy", np.asarray(mem, dtype=np.int64))
print(f"total {time.time()-t0:.1f} s")
