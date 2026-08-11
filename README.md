# isrm.esm — the InMAP ISRM as a language-agnostic model

> **TRANSITION (2026-08-10).** `isrm.esm` is now the former `isrm_clean.esm` — the
> fully-automatic variant (in-model Lambert projection, engine-derived pushdown,
> no runner pre-pass). The runners below still target the retired
> `isrm_pushdown.esm` and are being replaced by thin shims phase by phase; until
> that lands they do not run from this tree. The last fully-working pushdown
> state is tagged **`pushdown-era`**, and its four full-scale result records are
> frozen in `contract/records/` as the comparison baseline. The numbers quoted in
> this README are from that baseline.

`isrm.esm` states the InMAP source–receptor tutorial as an
[EarthSciAST](https://github.com/EarthSciML/EarthSciAST) `.esm` document, and
computes it three times — once through each of EarthSciAST's Julia, Python, and
Rust bindings — from that one document.

The point is not that three programs produce the same number. It is that **one
spec drives three independent engines**, and the engines agree.

## The result

Full scale: 52,411 receptor cells × 1,520 emission-bearing source cells × 43,650
EGU FF10 emission records, against the live `s3://inmap-model/isrm_v1.2.1.zarr`.

| binding | `sum(deathsK)` | vs tutorial |
|---|---|---|
| tutorial oracle | 7524.918845602511 | — |
| Julia | 7524.918845602511 | 0.0 |
| Python | 7524.9188456024895 | 2.9e-15 |
| Rust | 7524.918845602512 | 1.3e-16 |

Every pair of records, compared field by field
([`contract/compare_results.py`](contract/compare_results.py)):

```
ppl: count=1520 sha256 MATCH          (all three bindings)
total_pm25   max rel diff 0.000e+00   (BIT-IDENTICAL, every pair)
deaths       max rel diff <= 3.0e-15  (julia vs rust: 1.2e-16 — one ulp)
pathway emis/conc rel <= 2.1e-16
300 checks, 0 failed
```

## Why `ppl` is the number that matters

`ppl` is the set of grid cells that contain at least one emission source — 1,520
of 52,411. It is the output of a spatial join, and everything downstream is
shaped by it.

**No runner computes it.** Each binding *derives* it by evaluating the model's
own producer aggregate, gated by the spatial `join.overlap` broad phase
(CONFORMANCE_SPEC §5.5.6). All three arrive at a byte-identical set:

```
sha256 = 6f784d7e66f63872901126dabb2dd7354a96cdcd3d4585b2f52d20b6105a875b
```

That is required, not incidental: `ppl` is integer-valued, so §5.5 demands byte
identity regardless of which candidate-generation backend an engine uses (Julia's
STRtree, Rust's rstar R*-tree, Python's brute-force oracle).

This is also why the runners drive **`isrm_pushdown.esm`** and not `isrm.esm`:

| model | `emis_src_cells` | `ppl` derived by |
|---|---|---|
| `isrm_clean.esm` | `src_cells` (52,411) | n/a — SR would be 52411² ≈ 330 GB |
| `isrm.esm` | `interval` of size `N_PPL` | **the runner**, passed in as `src_cell_of_ppl` |
| `isrm_pushdown.esm` | `derived` via `join.overlap` | **the graph** |

Driving `isrm.esm` would have produced exactly the same death totals while
proving nothing: the spatial join would have happened in three hand-written
runners rather than in the spec.

## The projection pushdown

The full SR matrix is 52,411 × 52,411 per pathway — about 330 GB across five
pathways. None of it is downloaded.

Value invention runs *first* and derives the 1,520 members; the SR fetch is then
built from them and pushed down to the zarr reader, which fetches only the
intersecting chunks:

```
gated SR fetch: layer 0, 1520 of 52411 source cells, all receptors
```

The gate also drives *enumeration*, not just filtering. At full scale the
producer visits **43,668 tuples out of a 2,287,740,150-tuple product** — it
iterates the broad-phase candidate pairs directly rather than testing every
`(record, cell)` combination.

## Layout

```
isrm.esm                  the model (isrm.esm, isrm_clean.esm, isrm_pushdown.esm)
contract/                 the shared result record: schema, emitters, comparator
run-model-jl/             Julia STEP-0 reference oracle    (mode=oracle_step0)
run-model-jl-pushdown/    Julia runner through the graph   (mode=runtime_observed_graph)
run-model-py/             Python runner through the graph  (mode=runtime_observed_graph)
run-model-rs/             Rust runner through the graph    (mode=runtime_observed_graph)
```

Each runner emits a record conforming to
[`contract/results_schema.json`](contract/results_schema.json). The `mode` field
is load-bearing: `oracle_step0` marks numbers produced by hand-written STEP-0
arithmetic over the loader outputs — a reference, *not* evidence that a binding
executes the spec. Only `runtime_observed_graph` records make that claim.

Comparing the two Julia records is the sharpest check available, since it holds
the engine fixed and varies only whether the math came from the spec:

```
--- julia[oracle_step0] vs julia[runtime_observed_graph] ---
  ppl: count=1520 sha256 MATCH
  pathway emis rel 0.00e+00   conc rel 0.00e+00   (all five, bit-identical)
  deaths.krewski   max rel diff 3.183e-16
```

## Running them

```bash
# Julia
cd run-model-jl-pushdown && julia -t 2 --heap-size-hint=12G --project=. L3_full.jl

# Python (needs >= 3.11; zarr 3.x requires it)
cd run-model-py && python run_model.py

# Rust
cd run-model-rs && cargo build --release
LD_LIBRARY_PATH=$(dirname $(find target -name libs2bindings_shim.so | head -1)) \
    ./target/release/run-model-rs

# compare whatever records exist
python3 contract/compare_results.py run-model-jl/results.json \
    run-model-jl-pushdown/results.json run-model-py/results.json run-model-rs/results.json
```

Every runner takes `ISRM_FIRSTN` / `L3_FIRSTN` to truncate the emission-record
list for a fast reduced run, and honours `ISRM_SCRATCH`. **Point scratch at a
disk-backed path** — on a cluster whose root filesystem is `tmpfs`, SR blobs
written under `/tmp` consume the memory the model needs.

## Timing

Same machine, full scale, one core each:

| binding | wall | dominated by |
|---|---|---|
| Rust | 569 s | 382 s gated SR fetch |
| Julia | 2,579 s | 586 s build + ~300 s/observed |
| Python | 15,938 s | 15,223 s in five `E_*` joins |

The spread is one known gap, not three engines of differing quality. The `E_*`
emission-binning join is **3.3 s in Rust and ~3,000 s in Python** — the same
aggregate over the same 66.3M `(cell, record)` pairs, ~900× apart. Python's
interpreter vectorizes pure maps and scaled-product contractions (einsum), but
this body is a contraction with an `ifelse` containment predicate, which matches
neither, so it falls to a per-pair Python loop. The compiled evaluators handle it
directly.

## Tolerances

`ppl` is compared **exactly**. Float fields use `RTOL_FIELD = 1e-12`.

That is not arbitrary: three engines contract a 1,520-term sum in different
orders, and reassociating a float sum changes the last bits. The measured
cross-binding spread is ≤ 3.0e-15. Asserting tighter would be asserting on
summation order rather than on the model.

One subtlety worth recording, because it nearly produced a false failure: the
record's `sum` must be a property of the *data*, not of the summing language.
Julia's `sum` is pairwise and CPython 3.12's is Neumaier-compensated, but Rust's
`Iterator::sum` is a naive fold — with **bit-identical** `total_pm25` fields
(same sha256, same samples) the Rust total differed by 2.9e-13 purely from
accumulation error. At full scale that could have exceeded the tolerance and
reported a disagreement between provably identical fields. All emitters now use
compensated summation.
