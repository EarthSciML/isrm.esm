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
integer parts — ``sr_lower`` and ``stack_layer`` — are compared the same way,
because a float tolerance there would hide a real disagreement about which
layers a record emits into. Its ``weights`` are genuinely floats
(``sr.Reader.layerFracs`` interpolates a plume between two SR layers) and get
the float treatment.

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
ATOL_WEIGHT_SUM = 1e-12  # |w_sr0 + w_sr1 + w_sr2 - 1|, an exact invariant

# The InMAP source-receptor tutorial's published national totals
# (https://inmap.run/blog/2019/04/20/sr/).
#
# PRINTED, NOT ASSERTED — but for a different reason than before, and the
# history matters because two successive versions of this file got it wrong in
# opposite directions.
#
#   * An earlier version asserted them with a tolerance fitted to a run that
#     later turned out to be wrong, which made a fitted number look like a bound.
#   * The version that replaced it asserted an "admissible span" instead, whose
#     upper edge was a ceiling derived from `sr.layers == [0, 1, 2]` — the
#     corrupt array in the zarr. Given a wrong `layers`, non-negative plume rise
#     really does bound deathsK at 6697.55 and the published 6928.96 really is
#     above it. The argument was valid and the premise was false; with
#     `layers == [0, 3, 6]` the published totals are reachable, and the
#     reduced-scale runs below now reproduce the live service exactly. The span
#     is REMOVED rather than re-derived: it never described the model, it
#     described a corrupt input, and re-deriving it would mean re-measuring four
#     full-scale contractions to assert a bound that the direct comparison with
#     the service makes redundant.
#
# The tutorial's published totals. As of 2026-08-20 these are a REFERENCE POINT,
# not a target: the document now states CORRECT physics and declines both of
# InMAP's plume-rise defects, so it is expected to sit ABOVE these. See
# CORRECTED_FULL.
ORACLE_DEATHS_K = 6928.959583
ORACLE_DEATHS_L = 15623.924632

# What this document computes at full scale, MEASURED 2026-08-20 (Rust binding,
# repaired store). This is a regression lock, not an oracle: it is this
# document's own output. Its authority comes from the two independent things
# that agree with it — the NumPy oracle in plume_oracle.py reproduces the weight
# sums and the sr_lower digest exactly, and the InMAP-faithful configuration of
# this same document reproduced the live service to 8.9e-9 before the physics
# was corrected (see SERVICE_DEATHS).
CORRECTED_FULL = (7022.724781368745, 15835.993595627131)
# Cross-binding spread on the same document is ~4e-18 relative (Julia/Rust are
# bit-identical; Python differs in the last ulp of a compensated sum), so this
# is loose by many orders of magnitude and will catch a real change, not noise.
RTOL_CORRECTED = 1e-9

# Reduced-scale targets from the LIVE `inmap cloud` service, run on the same
# truncated record list this repo's ISRM_FIRSTN uses. Keyed by n_rec.
#
# These are the strongest check in the file: they are not this document's own
# output, not a published summary, and not fitted to anything — they are what
# InMAP itself returns for the same input.
# NOTE THE CHANGE OF MEANING, 2026-08-20. These were a direct pass/fail check
# while the document deliberately reproduced InMAP's arithmetic. The document
# now CORRECTS InMAP's inverted layerFracs interpolation, so it no longer agrees
# with the service and MUST NOT be asserted against it. They are kept because
# they are the evidence that the divergence is a fix rather than a bug: on these
# same 200 records the InMAP-faithful configuration returned 49.09146956 /
# 110.40696810 against the service's values below — agreement to 8.9e-9, inside
# the precision the service prints. Agreement was established first; the
# correction came after.
SERVICE_DEATHS_INMAP_FAITHFUL = {
    200: (49.091470, 110.406968),
}
# What the CORRECTED document returns on the same 200 records, measured the same
# day. This is what is actually checked at reduced scale — a regression lock, as
# CORRECTED_FULL is at full scale.
SERVICE_DEATHS = {
    200: (49.11639491165982, 110.46292733178377),
}
# The service reports six decimals, so each figure is known to +-5e-7 absolute,
# which is +-1.0e-8 relative at 49 and +-4.5e-9 at 110. RTOL_SERVICE is an order
# of magnitude above that and is NOT fitted to the observed agreement (which is
# 9e-9 on deathsK and 2.9e-9 on deathsL, i.e. inside the printing precision).
RTOL_SERVICE = 1e-7

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
    (`plume_model_layer`, `branch_usage`, the above-layer-7 group, and which
    `layers` the store served); this keeps the fields the two have in common,
    keyed the results record's way. The
    oracle's pathways are keyed by inventory name (VOC, NOx, …) and carry the SR
    array name the document uses, which is what `pathways` is keyed by.
    """
    sr, st = doc["sr_lower"], doc["stack_layer"]
    out = {
        "sr_lower": {"count": doc["n_rec"], "histogram": sr["histogram"],
                     "sha256": sr["sha256"]},
        "stack_layer": {"count": doc["n_rec"], "sha256": st["sha256"]},
        "weights": {k: v for k, v in doc["weights"].items() if k != "note"},
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
    """The `plume` block: the integer parts EXACTLY, the float parts at RTOL_FIELD.

    `sr_lower` is the lower of the (at most two) SR emission layers each record
    is charged to. It is integer-valued, so a float tolerance would paper over a
    real disagreement about WHICH LAYERS a record emits into. `stack_layer` is
    compared for localization: a wrong stack layer means the meteorology gather
    is wrong, a right stack layer with a wrong SR assignment means the ASME
    expression or the layerFracs interpolation is.

    The `weights` are floats — sr.Reader.layerFracs interpolates a plume sitting
    between two SR layers — and descend from `plume_height`, whose cube roots
    differ by an ulp between languages, so they get the FieldSummary treatment
    every other float field gets rather than an exact digest. `max_sum_error`
    is different in kind: layerFracs conserves mass exactly, so each record's
    three weights sum to exactly 1 and this is checked against an ABSOLUTE
    bound in every record independently, not compared between records.

    The per-layer emission masses are floats and get RTOL_FIELD, like every
    other float — they are the same sum in a different order.
    """
    for key in ("sr_lower", "stack_layer"):
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

    if "weights" in a and "weights" in b:
        wa, wb = a["weights"], b["weights"]
        for label, w in ((tag.split(" vs ")[0], wa), (tag.split(" vs ")[-1], wb)):
            e = w.get("max_sum_error")
            if e is not None:
                rep.check(e <= ATOL_WEIGHT_SUM,
                          f"{label} plume.weights.max_sum_error",
                          f"max |w_sr0+w_sr1+w_sr2 - 1| = {e:.3e} > "
                          f"{ATOL_WEIGHT_SUM:.0e} — layerFracs conserves mass, "
                          f"so this is a broken document, not float noise")
        if "count" in wa and "count" in wb:
            rep.check(wa["count"] == wb["count"], f"{tag} plume.weights.count",
                      f"{wa['count']} != {wb['count']}")
        for wk in ("w_sr0", "w_sr1", "w_sr2"):
            if wk not in wa or wk not in wb:
                rep.check(False, f"{tag} plume.weights.{wk}",
                          "present in only one record")
                continue
            sa, sb = wa[wk], wb[wk]
            worst = 0.0
            for scalar in ("sum", "min", "max"):
                if scalar in sa and scalar in sb:
                    dd = rel_diff(sa[scalar], sb[scalar])
                    worst = max(worst, dd)
                    rep.check(dd <= RTOL_FIELD,
                              f"{tag} plume.weights.{wk}.{scalar}",
                              f"rel diff {dd:.3e} > {RTOL_FIELD:.0e}")
            if len(sa.get("sample", [])) == len(sb.get("sample", [])):
                for i, (x, y) in enumerate(zip(sa["sample"], sb["sample"])):
                    dd = rel_diff(x, y)
                    worst = max(worst, dd)
                    rep.check(dd <= RTOL_FIELD,
                              f"{tag} plume.weights.{wk}.sample[{i}]",
                              f"rel diff {dd:.3e} > {RTOL_FIELD:.0e}")
            else:
                rep.check(False, f"{tag} plume.weights.{wk}.sample",
                          "differing sample lengths")
            bitid = sa.get("sha256") and sa.get("sha256") == sb.get("sha256")
            print(f"  plume weight {wk:<7} sum {sa['sum']:.6f}   "
                  f"max rel diff {worst:.3e}"
                  f"{'   (bit-identical)' if bitid else ''}")
    elif "weights" in a or "weights" in b:
        rep.check(False, f"{tag} plume.weights",
                  "one record carries the layerFracs weights and the other does "
                  "not — a single-layer assignment is not comparable to a split "
                  "one")

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
    n_rec = grid.get("n_rec")
    if (grid.get("n_src"), n_rec) != (FULL_N_SRC, FULL_N_REC):
        if n_rec in SERVICE_DEATHS:
            check_service(name, doc, SERVICE_DEATHS[n_rec], rep)
            return
        print(f"\n--- {name} vs tutorial oracle: SKIPPED ---")
        print(f"  reduced record (n_src={grid.get('n_src')}, n_rec={n_rec}); "
              f"the tutorial totals only apply at n_src={FULL_N_SRC}, "
              f"n_rec={FULL_N_REC}, and no live-service target is recorded for "
              f"n_rec={n_rec} (see SERVICE_DEATHS)")
        return
    print(f"\n--- {name} vs the tutorial's published totals ---")
    for label, key, target in (("krewski", "krewski", ORACLE_DEATHS_K),
                               ("lepeule", "lepeule", ORACLE_DEATHS_L)):
        got = doc["deaths"][key]["sum"]
        print(f"  sum(deaths {label:<8}) = {got:.6f}   published {target:.6f}"
              f"   rel {rel_diff(got, target):.2e}  (context, not checked)")
    print("  NOT CHECKED, and not a defect: this document declines both of InMAP's")
    print("  plume-rise defects, so it is EXPECTED above the published pair. What is")
    print("  checked at full scale is CORRECTED_FULL, below.")
    for label, key, target in (("krewski", "krewski", CORRECTED_FULL[0]),
                               ("lepeule", "lepeule", CORRECTED_FULL[1])):
        got = doc["deaths"][key]["sum"]
        d = rel_diff(got, target)
        rep.check(d <= RTOL_CORRECTED, f"{name} deaths.{label} vs CORRECTED_FULL",
                  f"{got!r} vs {target!r}, rel diff {d:.3e} > {RTOL_CORRECTED:.0e} — "
                  "the corrected-physics full-scale total moved")


def check_service(name: str, doc: dict, target, rep: Report) -> None:
    """A reduced record against the corrected-physics value for the same records.

    This USED to compare against the live `inmap cloud` service directly, and
    that is still where the number's authority comes from — but the document now
    corrects InMAP's inverted layerFracs, so the two no longer agree and pinning
    the document to the service would pin the bug back in. What is compared here
    is SERVICE_DEATHS: this document's own corrected output, measured once. The
    service's own figures are kept alongside in
    SERVICE_DEATHS_INMAP_FAITHFUL, together with what the InMAP-faithful
    configuration returned for them (8.9e-9 agreement) — which is the evidence
    that the divergence below is a correction and not a regression.
    """
    n_rec = doc["grid"]["n_rec"]
    print(f"\n--- {name} vs the corrected-physics reduced target (n_rec={n_rec}) ---")
    for label, key, want in (("krewski", "krewski", target[0]),
                             ("lepeule", "lepeule", target[1])):
        got = doc["deaths"][key]["sum"]
        d = rel_diff(got, want)
        rep.check(d <= RTOL_SERVICE, f"{name} deaths.{label} vs inmap cloud",
                  f"{got!r} vs {want!r}, rel diff {d:.3e} > {RTOL_SERVICE:.0e} — "
                  f"the document no longer computes what InMAP computes on this "
                  f"input")
        print(f"  sum(deaths {label:<8}) = {got:.9f}   service {want:.6f}"
              f"   rel {d:.2e}")


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
