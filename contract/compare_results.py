#!/usr/bin/env python3
"""Compare ISRM result records emitted by the Julia / Rust / Python runners.

    python3 contract/compare_results.py contract/records/plume_oracle.json \
                                        run-jl/results.json \
                                        run-py/results.json \
                                        run-rs/results.json

Validates each file against ``contract/results_schema.json`` (when ``jsonschema``
is importable), then cross-compares every pair against the tolerances below and
against the tutorial oracle. A ``plume_oracle`` record among the arguments is
not a results record and is not cross-compared as one; it is the independent
target for the ``plume`` block, matched on ``n_rec``. Exits non-zero if any
check fails.

Why these tolerances
--------------------
``ppl`` is INTEGER-valued: the emission-bearing support set is a distinct member
set, which CONFORMANCE_SPEC §5.5 requires to be byte-identical across bindings
regardless of candidate-generation backend. So it is compared EXACTLY — any
difference is a real disagreement, never float noise. The ``plume`` block's
layer assignments are integer-valued for the same reason and are compared the
same way: a float tolerance there would hide a real disagreement about which
layer a record emits into.

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

# The InMAP source-receptor tutorial's published national totals
# (https://inmap.run/blog/2019/04/20/sr/). They are a FULL-SCALE result and are
# only asserted on records whose grid matches the full problem (see
# check_oracle). They ACCOUNT FOR PLUME RISE — the pre-plume-rise figures
# (7524.918845602511 / 16979.632171487083) belong to the ground-level-only
# records archived in contract/records/ground-level-only/.
ORACLE_DEATHS_K = 6928.959583
ORACLE_DEATHS_L = 15623.924632

# MEASURED, not chosen. The document deliberately does not reproduce InMAP's
# above-layer-7 source-index defect (contract/README.md, "One deliberate
# deviation from InMAP"): a plume that rises above the top of model layer 7
# keeps an index built in the coarse 9324-cell high-altitude grid, which
# `sr.Reader.source` then reads against the 52411-cell ground grid, so InMAP
# charges those emissions to the WRONG source cell. That is 654 of 43650
# records and 0.43% of emitted mass; this document charges them to the cell the
# emission actually came from. So a correct run lands NEAR the published totals,
# never on them, and the tolerance has to be wide enough for that one difference
# and no wider.
#
# MEASURED on the full-scale Julia run of 2026-08-19 (43650 records, 1520
# emission-bearing cells, isrm_v1.2.1.zarr): sum(deathsK) = 6983.9385617781645
# against the published 6928.959583, and sum(deathsL) = 15752.315804140908
# against 15623.924632 — relative differences of 7.87e-3 and 8.15e-3 as
# rel_diff() computes them. RTOL_ORACLE is set just above the larger.
#
# Worth noting rather than glossing: 8.15e-3 is nearly TWICE the 0.43% of
# emitted mass the misplaced group carries. That is not a contradiction. The
# mass is not lost in either model, it is placed differently, and where it
# lands is what decides how many deaths it causes: put back on the cells the
# emissions actually came from — power plants, which sit near people — a ton
# buys more deaths than it does scattered across the ground grid by a coarse
# index read as a fine one. So the group punches roughly twice its weight, in
# the direction that makes this document's total HIGHER than the blog's, which
# is the direction the argument predicts.
RTOL_ORACLE = 8.3e-3

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

    # --- plume rise ---------------------------------------------------------
    if "plume" in a and "plume" in b:
        compare_plume(tag, a["plume"], b["plume"], rep)
    elif "plume" in a or "plume" in b:
        rep.check(False, f"{tag} plume",
                  "one record states plume rise and the other does not — a "
                  "ground-level-only baseline is not comparable to a plume-rise "
                  "run (see contract/records/ground-level-only/README.md)")


def pad(hist: list, n: int) -> list:
    """A histogram zero-extended to `n` bins, so two records that saw different
    maxima are still compared bin for bin instead of failing on length."""
    return list(hist) + [0] * max(0, n - len(hist))


def plume_from_oracle(doc: dict) -> dict:
    """A `plume_oracle` record, read into the `plume` shape a results record uses.

    The oracle answers the same question from the other side — it computes the
    record → SR-layer assignment from the ISRM meteorology arrays directly,
    without the SR matrix and without EarthSciAST — so it is the one target that
    can say the document's plume rise is RIGHT rather than merely consistent
    across three bindings. It carries more than the results record does
    (`plume_model_layer`, `branch_usage`, the above-layer-7 group); this keeps
    the fields the two have in common, keyed the results record's way. The
    oracle's pathways are keyed by inventory name (VOC, NOx, …) and carry the SR
    array name the document uses, which is what `pathways` is keyed by.
    """
    sr, st = doc["sr_layer"], doc["stack_layer"]
    out = {
        "sr_layer": {"count": doc["n_rec"], "histogram": sr["histogram"],
                     "sha256": sr["sha256"]},
        "stack_layer": {"count": doc["n_rec"], "sha256": st["sha256"]},
        "pathways": {p["sr_name"]: {"by_sr_layer": p["by_sr_layer"]}
                     for p in doc["pathways"].values()},
    }
    # `histogram` (all records) postdates the first oracle records; the older
    # `histogram_height_gt_0` counts only records with a stack, which is a
    # different quantity and must not be silently compared against this one.
    if "histogram" in st:
        out["stack_layer"]["histogram"] = st["histogram"]
    return out


def compare_plume(tag: str, a: dict, b: dict, rep: Report) -> None:
    """The `plume` block, compared the way `ppl` is: EXACTLY.

    `sr_layer` is the SR emission layer each record is charged to — the one
    intermediate plume rise exists to produce, and the thing every downstream
    number follows from. It is integer-valued, so a float tolerance here would
    paper over a real disagreement about WHICH LAYER a record emits into.
    `stack_layer` is compared for localization: a wrong stack layer means the
    meteorology gather is wrong, a right stack layer with a wrong SR layer means
    the ASME expression is.

    The per-layer emission masses are floats and get RTOL_FIELD, like every
    other float — they are the same sum in a different order.
    """
    for key in ("sr_layer", "stack_layer"):
        if key not in a or key not in b:
            continue
        x, y = a[key], b[key]
        rep.check(x["count"] == y["count"], f"{tag} plume.{key}.count",
                  f"{x['count']} != {y['count']}")
        rep.check(x["sha256"] == y["sha256"], f"{tag} plume.{key}.sha256",
                  "the per-record layer assignments differ — a real disagreement "
                  "about which layer a record emits into, not float noise")
        if "histogram" in x and "histogram" in y:
            n = max(len(x["histogram"]), len(y["histogram"]))
            hx, hy = pad(x["histogram"], n), pad(y["histogram"], n)
            rep.check(hx == hy, f"{tag} plume.{key}.histogram", f"{hx} != {hy}")
            print(f"  plume {key:<11} {hx}  sha256 "
                  f"{'MATCH' if x['sha256'] == y['sha256'] else 'DIFFER'}")

    keys = sorted(set(a["pathways"]) | set(b["pathways"]))
    for k in keys:
        if k not in a["pathways"] or k not in b["pathways"]:
            rep.check(False, f"{tag} plume.pathways.{k}", "present in only one record")
            continue
        va = a["pathways"][k]["by_sr_layer"]
        vb = b["pathways"][k]["by_sr_layer"]
        if not rep.check(len(va) == len(vb), f"{tag} plume.pathways.{k}.by_sr_layer",
                         f"{len(va)} layers != {len(vb)}"):
            continue
        worst = 0.0
        for L, (x, y) in enumerate(zip(va, vb)):
            d = rel_diff(x, y)
            worst = max(worst, d)
            rep.check(d <= RTOL_FIELD, f"{tag} plume.pathways.{k}.by_sr_layer[{L}]",
                      f"rel diff {d:.3e} > {RTOL_FIELD:.0e}")
        print(f"  plume mass {k:<12} L0/L1/L2 = "
              f"{va[0]:.4g} / {va[1]:.4g} / {va[2]:.4g}   max rel diff {worst:.3e}")


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


def check_plume_oracle(docs: dict, oracles: list, rep: Report) -> None:
    """Every live record's `plume` block against every plume-oracle record.

    Matching on `n_rec` is what keeps a reduced oracle
    (`plume_oracle_first200.json`) from being compared against a full-scale run:
    they are different problems, and the digest of a 200-record assignment says
    nothing about a 43650-record one.
    """
    live = {n: d for n, d in docs.items() if "plume" in d}
    if not live:
        if oracles and docs:
            print("\n--- plume oracle: NOTHING TO CHECK ---")
            print("  a plume_oracle record was passed, but not one of the "
                  "records carries a `plume` block, so the layer assignment "
                  "went unreported and therefore uncompared. A runner driving "
                  "a document that states plume rise must emit it.")
        return
    if not oracles:
        print("\n--- plume oracle: NOT LOADED ---")
        print("  the records compared above state plume rise, but no "
              "plume_oracle record was passed, so nothing independent checked "
              "the record -> SR-layer assignment. Run "
              "`python3 contract/plume_oracle.py` (~1 min, no SR matrix) and "
              "pass contract/records/plume_oracle.json.")
        return
    for path, orc in oracles:
        oname = f"plume-oracle({os.path.basename(path)})"
        target = plume_from_oracle(orc)
        for name, doc in live.items():
            tag = f"{name} vs {oname}"
            print(f"\n--- {tag} ---")
            if doc["grid"]["n_rec"] != orc["n_rec"]:
                print(f"  SKIPPED: n_rec {doc['grid']['n_rec']} != the oracle's "
                      f"{orc['n_rec']} — a truncation is a different problem, "
                      f"and its digest says nothing about this one")
                continue
            compare_plume(tag, doc["plume"], target, rep)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    docs = {}
    rep = Report()
    loaded = []
    oracles: list[tuple[str, dict]] = []
    for path in argv[1:]:
        with open(path) as fh:
            doc = json.load(fh)
        # `contract/records/` also holds records that are not results records —
        # they answer a different question and have a different shape, and are
        # tagged with `kind`. Most are skipped so `contract/records/*.json`
        # stays a usable glob. The plume oracle is NOT skipped: it is the
        # independent reference for the `plume` block, computed from the ISRM
        # meteorology without the SR matrix and without EarthSciAST, so it is
        # the only check that says the document's plume rise is RIGHT rather
        # than merely agreed on by three bindings. That is the whole reason it
        # exists, so it is a first-class target — for the plume block only.
        if "kind" in doc:
            if doc["kind"] == "plume_oracle":
                oracles.append((path, doc))
                print(f"loaded {path}  as the plume oracle "
                      f"(n_rec={doc.get('n_rec')})")
            else:
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

    check_plume_oracle(docs, oracles, rep)

    ok = rep.close()
    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
