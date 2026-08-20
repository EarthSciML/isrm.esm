"""Emit a ``contract/results_schema.json`` record from a Python runner.

The Python mirror of ``contract/results.jl``. The hashing and sampling rules
here MUST match ``contract/compare_results.py`` (which is the definition of
record) and ``contract/results.jl`` — the three implementations of
``sample_indices`` are deliberately pure integer arithmetic so the index set
cannot drift between languages' float rounding.
"""

from __future__ import annotations

import hashlib
import json
import os
import struct
from typing import Any, Iterable, Mapping, Sequence

SAMPLE_N = 25


def sample_indices(n_rcv: int) -> list[int]:
    """The fixed 1-based sample indices every runner reports."""
    d = SAMPLE_N - 1
    return [1 + (k * (n_rcv - 1) + d // 2) // d for k in range(SAMPLE_N)]


def int_seq_sha256(values: Iterable[int]) -> str:
    """sha256 over an ORDERED integer sequence as ASCII decimals joined by ','.

    The same wire format as :func:`ppl_sha256` — ASCII decimals, ``,`` separator,
    no spaces — but order-preserving. ``ppl`` is a member *set*, so it is sorted
    before hashing; a per-record integer field (e.g. the plume-rise SR-layer
    assignment, one value per emission record) is a *sequence* whose order is
    part of the value, so it must not be sorted. Same convention, two uses.
    """
    s = ",".join(str(int(v)) for v in values)
    return hashlib.sha256(s.encode("ascii")).hexdigest()


def ppl_sha256(ids: Iterable[int]) -> str:
    """sha256 over sorted 1-based ids as ASCII decimals joined by ',' (no spaces)."""
    return int_seq_sha256(sorted(int(i) for i in ids))


def _as_int_seq(values: Iterable[float], label: str) -> list[int]:
    """An integer-valued observed, read back off the graph as float64, as ints.

    The document's `sr_lower` / `stack_layer` are sums of 0.0/1.0 indicators,
    so every value is an exact integer in float64 and this rounding is lossless.
    A value that is NOT integral means the observed is no longer the indicator
    sum it is supposed to be — that is a real disagreement about the physics, so
    it raises rather than rounding it away.
    """
    out: list[int] = []
    for i, x in enumerate(values):
        f = float(x)
        n = round(f)
        if abs(f - n) > 0:
            raise ValueError(
                f"{label}[{i}] = {f!r} is not integral; an integer-valued observed "
                "came back fractional, so the layer assignment is not what the "
                "document states"
            )
        out.append(int(n))
    return out


def histogram(values: Sequence[int], min_bins: int) -> list[int]:
    """Counts of 0, 1, 2, ... over a non-negative integer sequence.

    At least ``min_bins`` bins so a reduced run, where a layer may simply be
    empty, still emits the same shape a full run does; more if the data needs
    them, because a value the schema does not expect must be VISIBLE rather
    than dropped off the end of a fixed-width histogram.
    """
    if any(v < 0 for v in values):
        raise ValueError("negative value in a layer assignment")
    bins = [0] * max(min_bins, (max(values) + 1) if values else 0)
    for v in values:
        bins[v] += 1
    return bins


def weight_sum_error(w0: Sequence[float], w1: Sequence[float],
                     w2: Sequence[float]) -> float:
    """max over records of |w0 + w1 + w2 - 1|.

    InMAP's ``layerFracs`` conserves mass: whether a record lands wholly in one
    SR layer or is split across two, its weights sum to exactly 1. So this is
    float noise or a broken document, and nothing in between.
    """
    n = len(w0)
    if not (len(w1) == len(w2) == n):
        raise ValueError("the three SR-layer weight fields differ in length")
    return max((abs(float(a) + float(b) + float(c) - 1.0)
                for a, b, c in zip(w0, w1, w2)), default=0.0)


def plume_block(
    *,
    sr_lower: Iterable[float],
    stack_layer: Iterable[float],
    weights: Mapping[str, Sequence[float]],
    emis_by_sr_layer: Mapping[str, Sequence[float]],
) -> dict[str, Any]:
    """The schema's `plume` block, from the document's OWN observeds.

    `sr_lower`, `stack_layer` and the three `weights` fields are per-record
    observeds read straight off the graph — nothing here recomputes plume rise,
    and nothing here knows what ASME is. `weights` maps ``"w_sr0"`` / ``"w_sr1"``
    / ``"w_sr2"`` to those observeds; `emis_by_sr_layer` maps each SR array name
    to the three ``sum(E_<pathway>_L<layer>)`` totals, in layer order.

    Two digests are integer sequences in record order, hashed the way `ppl` is,
    so they can be compared EXACTLY — against the other bindings and against
    `contract/records/plume_oracle.json`, which computes the same assignment
    from the meteorology arrays without touching the SR matrix. The weights
    themselves are floats (``sr.Reader.layerFracs`` interpolates a plume between
    two SR layers) and get the FieldSummary treatment every other float here
    gets.
    """
    sr = _as_int_seq(sr_lower, "sr_lower")
    sl = _as_int_seq(stack_layer, "stack_layer")
    w = {k: [float(x) for x in weights[k]] for k in ("w_sr0", "w_sr1", "w_sr2")}
    return {
        "sr_lower": {
            "count": len(sr),
            "histogram": histogram(sr, 3),
            "sha256": int_seq_sha256(sr),
        },
        "stack_layer": {
            "count": len(sl),
            "histogram": histogram(sl, 4),
            "sha256": int_seq_sha256(sl),
        },
        "weights": {
            "count": len(w["w_sr0"]),
            "max_sum_error": weight_sum_error(w["w_sr0"], w["w_sr1"], w["w_sr2"]),
            **{k: field_summary(v) for k, v in w.items()},
        },
        "pathways": {
            str(k): {"by_sr_layer": [float(x) for x in v]}
            for k, v in emis_by_sr_layer.items()
        },
    }


def field_sha256(v: Sequence[float]) -> str:
    """sha256 over a float field as little-endian IEEE-754 float64 bytes."""
    h = hashlib.sha256()
    for x in v:
        h.update(struct.pack("<d", float(x)))
    return h.hexdigest()


def field_summary(v: Sequence[float]) -> dict[str, Any]:
    """Summarize one length-n_rcv field into the schema's FieldSummary shape."""
    vv = [float(x) for x in v]
    idx = sample_indices(len(vv))
    return {
        "sum": sum(vv),
        "min": min(vv),
        "max": max(vv),
        "sample": [vv[i - 1] for i in idx],
        "sha256": field_sha256(vv),
    }


def write_results(
    path: str,
    *,
    model: str,
    mode: str,
    n_src: int,
    n_rcv: int,
    n_rec: int,
    ppl: Iterable[int],
    pathways: Mapping[str, Mapping[str, float]],
    total_pm25: Sequence[float],
    deathsK: Sequence[float],
    deathsL: Sequence[float],
    binding_version: str = "",
    include_ppl_ids: bool = True,
    plume: Mapping[str, Any] | None = None,
    timing: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if mode not in ("runtime_observed_graph", "oracle_step0"):
        raise ValueError(
            f'mode must be "runtime_observed_graph" or "oracle_step0", got {mode!r}'
        )
    ids = sorted(int(i) for i in ppl)

    pw: dict[str, Any] = {}
    for k, v in pathways.items():
        entry = {"emis_sum": float(v["emis_sum"]), "conc_sum": float(v["conc_sum"])}
        if "conc_max" in v:
            entry["conc_max"] = float(v["conc_max"])
        pw[str(k)] = entry

    ppl_rec: dict[str, Any] = {"count": len(ids), "sha256": ppl_sha256(ids)}
    if include_ppl_ids:
        ppl_rec["ids"] = ids

    rec: dict[str, Any] = {
        "binding": "python",
        "binding_version": binding_version,
        "model": os.path.basename(model),
        "mode": mode,
        "grid": {"n_src": int(n_src), "n_rcv": int(n_rcv), "n_rec": int(n_rec)},
        "ppl": ppl_rec,
        "pathways": pw,
        "total_pm25": field_summary(total_pm25),
        "deaths": {
            "krewski": field_summary(deathsK),
            "lepeule": field_summary(deathsL),
        },
    }
    if plume is not None:
        rec["plume"] = dict(plume)
    if timing is not None:
        rec["timing"] = {str(k): v for k, v in timing.items()}

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    print(
        f"wrote {path}  (mode={mode}, |ppl|={len(ids)}, "
        f"sum deathsK={rec['deaths']['krewski']['sum']})"
    )
    return rec
