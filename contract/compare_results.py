#!/usr/bin/env python3
"""Compare ISRM result records emitted by the Julia / Rust / Python runners.

    python3 contract/compare_results.py run-model-jl/results.json \
                                        run-model-py/results.json \
                                        run-model-rs/results.json

Validates each file against ``contract/results_schema.json`` (when ``jsonschema``
is importable), then cross-compares every pair against the tolerances below and
against the tutorial oracle. Exits non-zero if any check fails.

Why these tolerances
--------------------
``ppl`` is INTEGER-valued: the emission-bearing support set is a distinct member
set, which CONFORMANCE_SPEC §5.5 requires to be byte-identical across bindings
regardless of candidate-generation backend. So it is compared EXACTLY — any
difference is a real disagreement, never float noise.

Float fields are NOT compared exactly. Three engines contract a 1520x52411 sum
in different orders, and reassociating a float sum changes the last bits; the
Julia wall2 work measured 6.2e-15 max relative difference between its own bit-
exact and BLAS-reassociated paths on this very contraction. A cross-language
tolerance below that would be asserting on reduction order, not on the model.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
RTOL_FIELD = 1e-12       # cross-binding, per-field scalars (sum/min/max/sample)
RTOL_ORACLE = 1e-3       # vs the published tutorial totals (they are rounded)

# The tutorial / run-model.jl oracle. It is a FULL-SCALE number: it is only
# asserted on records whose grid matches the full problem (see check_oracle).
ORACLE_DEATHS_K = 7524.918845602511
ORACLE_DEATHS_L = 16979.632171487083
FULL_N_SRC = 52411
FULL_N_REC = 43650

SAMPLE_N = 25


def sample_indices(n_rcv: int) -> list[int]:
    """The fixed 1-based sample indices every runner must report.

    Deliberately PURE INTEGER arithmetic (round-half-up via +d//2, then floor
    division) so the index set cannot drift between languages' float rounding
    modes. Mirrored in each runner; this is the definition of record.
    """
    d = SAMPLE_N - 1
    return [1 + (k * (n_rcv - 1) + d // 2) // d for k in range(SAMPLE_N)]


def rel_diff(a: float, b: float) -> float:
    """Relative difference, falling back to absolute near zero."""
    if a == b:
        return 0.0
    scale = max(abs(a), abs(b))
    if scale < 1e-300:
        return abs(a - b)
    return abs(a - b) / scale


class Report:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.checks = 0

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        self.checks += 1
        if not ok:
            self.failures.append(f"{label}: {detail}")
        return ok

    def close(self, kept_quiet_about: int = 0) -> bool:
        print(f"\n{self.checks} checks, {len(self.failures)} failed")
        for f in self.failures:
            print("  FAIL", f)
        return not self.failures


def validate_schema(path: str, doc: dict, rep: Report) -> None:
    schema_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "results_schema.json")
    try:
        import jsonschema
    except ImportError:
        print(f"  (jsonschema not importable — skipping schema check for {path})")
        return
    with open(schema_path) as fh:
        schema = json.load(fh)
    errs = sorted(jsonschema.Draft202012Validator(schema).iter_errors(doc),
                  key=lambda e: list(e.path))
    rep.check(not errs, f"schema/{os.path.basename(path)}",
              "; ".join(f"/{'/'.join(str(p) for p in e.path)}: {e.message}" for e in errs[:4]))


def field_summary_pairs(doc: dict):
    """Yield (label, summary-dict) for every FieldSummary in a record."""
    yield "total_pm25", doc["total_pm25"]
    yield "deaths.krewski", doc["deaths"]["krewski"]
    yield "deaths.lepeule", doc["deaths"]["lepeule"]


def compare_pair(name_a: str, a: dict, name_b: str, b: dict, rep: Report) -> None:
    tag = f"{name_a} vs {name_b}"
    print(f"\n--- {tag} ---")

    # Dimensions must match or nothing else is meaningful.
    if not rep.check(a["grid"] == b["grid"], f"{tag} grid",
                     f"{a['grid']} != {b['grid']}"):
        print("  grid mismatch — skipping the rest of this pair")
        return

    # --- ppl: exact ---------------------------------------------------------
    rep.check(a["ppl"]["count"] == b["ppl"]["count"], f"{tag} ppl.count",
              f"{a['ppl']['count']} != {b['ppl']['count']}")
    rep.check(a["ppl"]["sha256"] == b["ppl"]["sha256"], f"{tag} ppl.sha256",
              "support sets differ — this is a real disagreement, not float noise")
    print(f"  ppl: count={a['ppl']['count']} sha256 "
          f"{'MATCH' if a['ppl']['sha256'] == b['ppl']['sha256'] else 'DIFFER'}")

    # --- per-pathway --------------------------------------------------------
    keys = sorted(set(a["pathways"]) | set(b["pathways"]))
    for k in keys:
        if k not in a["pathways"] or k not in b["pathways"]:
            rep.check(False, f"{tag} pathways.{k}", "present in only one record")
            continue
        for field in ("emis_sum", "conc_sum"):
            d = rel_diff(a["pathways"][k][field], b["pathways"][k][field])
            rep.check(d <= RTOL_FIELD, f"{tag} pathways.{k}.{field}",
                      f"rel diff {d:.3e} > {RTOL_FIELD:.0e}")
        print(f"  pathway {k:<12} emis rel {rel_diff(a['pathways'][k]['emis_sum'], b['pathways'][k]['emis_sum']):.2e}"
              f"   conc rel {rel_diff(a['pathways'][k]['conc_sum'], b['pathways'][k]['conc_sum']):.2e}")

    # --- fields -------------------------------------------------------------
    for (label, sa), (_, sb) in zip(field_summary_pairs(a), field_summary_pairs(b)):
        worst = 0.0
        for scalar in ("sum", "min", "max"):
            if scalar in sa and scalar in sb:
                d = rel_diff(sa[scalar], sb[scalar])
                worst = max(worst, d)
                rep.check(d <= RTOL_FIELD, f"{tag} {label}.{scalar}",
                          f"rel diff {d:.3e} > {RTOL_FIELD:.0e}")
        if len(sa["sample"]) == len(sb["sample"]):
            for i, (x, y) in enumerate(zip(sa["sample"], sb["sample"])):
                d = rel_diff(x, y)
                worst = max(worst, d)
                rep.check(d <= RTOL_FIELD, f"{tag} {label}.sample[{i}]",
                          f"rel diff {d:.3e} > {RTOL_FIELD:.0e}")
        else:
            rep.check(False, f"{tag} {label}.sample", "differing sample lengths")
        bitid = sa.get("sha256") and sa.get("sha256") == sb.get("sha256")
        print(f"  {label:<16} max rel diff {worst:.3e}"
              f"{'   (bit-identical)' if bitid else ''}")


def check_oracle(name: str, doc: dict, rep: Report) -> None:
    # The tutorial totals are a FULL-SCALE result. A reduced record (a runner
    # driven with a truncated emission-record list, to exercise the pipeline
    # cheaply) is a different problem, not a failing one — asserting the oracle
    # against it reports a bogus ~99% error and buries the cross-binding
    # comparison that IS meaningful. Cross-binding checks still run on it.
    grid = doc.get("grid") or {}
    if (grid.get("n_src"), grid.get("n_rec")) != (FULL_N_SRC, FULL_N_REC):
        print(f"\n--- {name} vs tutorial oracle: SKIPPED ---")
        print(f"  reduced record (n_src={grid.get('n_src')}, n_rec={grid.get('n_rec')}); "
              f"the tutorial totals only apply at n_src={FULL_N_SRC}, n_rec={FULL_N_REC}")
        return
    print(f"\n--- {name} vs tutorial oracle ---")
    for label, key, target in (("krewski", "krewski", ORACLE_DEATHS_K),
                               ("lepeule", "lepeule", ORACLE_DEATHS_L)):
        got = doc["deaths"][key]["sum"]
        d = rel_diff(got, target)
        rep.check(d <= RTOL_ORACLE, f"{name} deaths.{label} vs oracle",
                  f"{got!r} vs {target!r}, rel diff {d:.3e}")
        print(f"  sum(deaths {label:<8}) = {got:.6f}   target {target:.6f}   rel {d:.2e}")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    docs = {}
    rep = Report()
    loaded = []
    for path in argv[1:]:
        with open(path) as fh:
            doc = json.load(fh)
        # `contract/records/` also holds records that are not results records —
        # e.g. plume_oracle.json, which answers a different question and has a
        # different shape. They are tagged with `kind`; results records are not.
        # Skip them so `contract/records/*.json` stays a usable glob.
        if "kind" in doc:
            print(f"skipping {path}  (kind={doc['kind']!r} — not a results record)")
            continue
        validate_schema(path, doc, rep)
        loaded.append((path, doc))

    # Label by binding, but DISAMBIGUATE when a binding appears more than once.
    # Two records can share a binding and still be the most interesting pair in
    # the set: the Julia oracle_step0 reference vs the Julia runtime_observed_graph
    # run is exactly what shows the runtime reproduces the hand-written oracle.
    # Keying on `binding` alone would silently drop one of them and compare nothing.
    bindings = [d.get("binding", os.path.basename(p)) for p, d in loaded]
    for (path, doc), b in zip(loaded, bindings):
        name = b if bindings.count(b) == 1 else f"{b}[{doc.get('mode', '?')}]"
        while name in docs:                       # same binding AND same mode
            name += "'"
        docs[name] = doc
        print(f"loaded {path}  as {name}  model={doc.get('model')} "
              f"mode={doc.get('mode')}")
        if doc.get("mode") != "runtime_observed_graph":
            print(f"  NOTE: mode={doc.get('mode')!r} — this record is a reference oracle, "
                  f"not a demonstration that the binding executes the .esm.")

    for name, doc in docs.items():
        check_oracle(name, doc, rep)

    for (na, a), (nb, b) in itertools.combinations(docs.items(), 2):
        compare_pair(na, a, nb, b, rep)

    ok = rep.close()
    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
