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


def ppl_sha256(ids: Iterable[int]) -> str:
    """sha256 over sorted 1-based ids as ASCII decimals joined by ',' (no spaces)."""
    s = ",".join(str(int(i)) for i in sorted(int(i) for i in ids))
    return hashlib.sha256(s.encode("ascii")).hexdigest()


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
