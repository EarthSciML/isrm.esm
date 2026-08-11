# isrm.esm — the InMAP ISRM as a language-agnostic model

`isrm.esm` states the InMAP source–receptor tutorial as **one**
[EarthSciAST](https://github.com/EarthSciML/EarthSciAST) `.esm` document, and
computes it three times — once through each of EarthSciAST's Julia, Python, and
Rust bindings — from that one document.

The point is not that three programs produce the same number. It is that **one
spec drives three independent engines**, and the engines agree.

The document is fully automatic — nothing model-shaped lives in the runners:

* **In-model projection.** The raw EGU `lon`/`lat` are the parameters; the
  Lambert conformal conic `X`/`Y` are observeds the engine evaluates at build
  time. No runner projects a coordinate.
* **Engine-derived pushdown.** `prepare(…, pushdown_rewrite=true)` runs the
  projection-pushdown rewrite inside the engine: it derives the
  emission-bearing support set (`pd_support__src_cells`) from the model's own
  spatial `join.overlap`, and records the SR gate in its own rewrite record
  (`metadata.x_esd.pushdown.gated_select`). No runner hand-authors a gate.
* **Engine-side gated SR fetch.** The zarr reads are built from the rewrite's
  gate and pushed down to the store, so only the intersecting chunks of the
  ~330 GB SR matrix are fetched.
* **Document-declared data sources.** Every provider — the SR zarr, grid,
  population, mortality, and the EGU FF10 inventory — comes from the
  document's `data_loaders` via `providers_from_document` (format =
  `metadata.esio_format`, URL = `source.url_template`). The FF10-in-zip
  quirks are declared too: the `*egu*` member glob (`x_esd.member_filter`)
  and the EPA column-header row (`x_esd.skip_header_row`) are reader options
  EarthSciIO implements in all three languages, and the POLID→pathway map is
  the document's `x_esd.pollutant_codes`.

Each shim contributes input plumbing and orchestration only; every reported
number is the binding's evaluation of the document's observed graph.

## The result

Full scale: 52,411 receptor cells × 1,520 emission-bearing source cells ×
43,650 EGU FF10 emission records, against the live
`s3://inmap-model/isrm_v1.2.1.zarr`. Same machine, same cache discipline:

| binding | wall | `sum(deathsK)` | notes |
|---|---|---|---|
| Julia | 2,538 s (688 s prepare/build + 1,777 s observed eval) | 7524.918845602511 (exact vs oracle) | peak RSS 10.6 GiB; `julia -t 2 --heap-size-hint=12G` |
| Python | 14,455 s (~4.0 h; five `E_*` joins 2,752–2,891 s each — known interpreter gap) | 7524.9188456024895 | peak RSS 5.74 GiB |
| Rust | 476 s (474 s prepare, dominated by the cold gated zarr fetch; ~2 s/contraction) | 7524.918845602512 | peak RSS 6.22 GiB |

The tutorial oracle is `sum(deathsK) = 7524.918845602511`,
`sum(deathsL) = 16979.632171487083` — Julia and Rust hit both exactly.

Every record, compared field by field against every other — the three live
records **and** the four frozen pushdown-era baselines in `contract/records/`
([`contract/compare_results.py`](contract/compare_results.py)):

```
ppl: count=1520 sha256 MATCH          (every record)
total_pm25   max rel diff 0.000e+00   (BIT-IDENTICAL, every pair)
deaths       max rel diff <= 3.1e-15
2058 checks, 0 failed
```

## Why `ppl` is the number that matters

`ppl` is the set of grid cells that contain at least one emission source —
1,520 of 52,411. It is the output of a spatial join, and everything downstream
is shaped by it.

**No runner computes it.** The engine's pushdown rewrite derives it as
`pd_support__src_cells` by evaluating the model's own producer aggregate,
gated by the spatial `join.overlap` broad phase (CONFORMANCE_SPEC §5.5.6).
All three engines — and the frozen pushdown-era records, where the same set
was derived from hand-authored pushdown constructs — arrive at a
byte-identical set:

```
sha256 = 6f784d7e66f63872901126dabb2dd7354a96cdcd3d4585b2f52d20b6105a875b
```

That is required, not incidental: `ppl` is integer-valued, so §5.5 demands
byte identity regardless of which candidate-generation backend an engine uses
(Julia's STRtree, Rust's rstar R*-tree, Python's brute-force oracle). A float
tolerance here would paper over a real disagreement about *which cells emit*.

## The projection pushdown

The full SR matrix is 52,411 × 52,411 per pathway — about 330 GB across five
pathways. None of it is downloaded.

The rewrite runs *first* and derives the 1,520 members; the SR fetch is then
built from them and pushed down to the zarr reader, which fetches only the
intersecting chunks:

```
gated SR fetch: 1520 of 52411 source cells, all receptors
```

The gate also drives *enumeration*, not just filtering: the producer visits
the broad-phase candidate pairs directly rather than testing every
`(record, cell)` combination of the 2.3-billion-tuple product.

## Layout

```
isrm.esm          the model — one document, no variants
contract/         the shared result record: schema, emitters (results.jl/.py),
                  comparator, and records/ — four frozen pushdown-era
                  full-scale baselines the live runs are compared against
run-jl/           Julia shim   (run.jl; setup.jl instantiates the project)
run-py/           Python shim  (run.py; requirements.txt is the venv recipe)
run-rs/           Rust shim    (cargo project)
data/             untracked — the EGU FF10 zip lives here
```

Each shim emits a record conforming to
[`contract/results_schema.json`](contract/results_schema.json) with
`mode="runtime_observed_graph"` — the mode that claims the numbers came from
the binding's evaluation of the spec, not from hand-written arithmetic.

## Prerequisites

* **Sibling checkouts** (none of these packages are published) next to this
  repo, i.e. `../EarthSciAST` and `../EarthSciIO`:
  * EarthSciAST at `ea64f510` or later — all three bindings need `prepare` +
    `pushdown_rewrite`;
  * EarthSciIO at `68d544e` or later — the FF10 reader needs
    `members`/`member_glob` + `skip_header_row`.

  The shims resolve them relative to the repo; override with `EA_PATH` /
  `IO_PATH`.
* **The EGU zip**: `data/2016fd_inputs_point.zip` (69 MB), from
  `https://gaftp.epa.gov/air/emismod/2016/alpha/2016fd/emissions/2016fd_inputs_point.zip`
  (the document's `EGU_Emis.source.url_template`). Override with `EGU_ZIP`.
  When the file is absent each shim falls back to fetching that URL through
  the EarthSciIO cache — but gaftp.epa.gov is slow and flaky, so a manual
  download into `data/` is the reliable path.
* **Julia** (developed on 1.12): `julia --project=run-jl run-jl/setup.jl`
  dev-tracks the two checkouts and instantiates.
* **Python ≥ 3.11** (zarr 3.x requires it). Build the venv on disk-backed
  scratch per the recipe in [`run-py/requirements.txt`](run-py/requirements.txt):

  ```bash
  python3.12 -m venv /scratch.local/$USER/isrm-py-venv
  /scratch.local/$USER/isrm-py-venv/bin/pip install -r run-py/requirements.txt
  ```
* **Rust** (stable toolchain). The s2geometry shim is a cdylib that is not on
  the runtime linker path, so the built binary needs `LD_LIBRARY_PATH`
  pointing at the build output (see below).

## Running the shims

```bash
cd isrm.esm   # this repo

# Julia
julia -t 2 --heap-size-hint=12G --project=run-jl run-jl/run.jl

# Python
/scratch.local/$USER/isrm-py-venv/bin/python run-py/run.py

# Rust
( cd run-rs && cargo build --release )
LD_LIBRARY_PATH=$(dirname $(find run-rs/target -name libs2bindings_shim.so | head -1)) \
    ./run-rs/target/release/run-rs

# compare whatever records exist against the frozen baselines
python3 contract/compare_results.py contract/records/*.json \
    run-jl/results.json run-py/results.json run-rs/results.json
```

Each shim writes `results.json` (or `results_reduced.json`) next to itself.

* `ISRM_FIRSTN=200` truncates the emission-record list for a fast reduced
  run (~2 min in Julia after precompile).
* `ISRM_SCRATCH` **must be disk-backed**. On this cluster `/tmp` is a tmpfs:
  SR chunk blobs written there consume the same memory cgroup the model
  needs. The shims default to `/scratch.local/$USER/isrm-esm` for exactly
  that reason.
* Size `--heap-size-hint` to what is actually free in the memory cgroup —
  here a 40 GB cgroup shared with other jobs, hence 12G.
* `ISRM_ESIO_CACHE` can point one shim at another's cache (the format is
  cross-language) to skip the ~6-minute cold SR fetch.

## Timing

The wall-clock spread in the result table is **one known gap, not three
engines of differing quality**. The `E_*` emission-binning join is ~3.3 s in
Rust and ~2,800 s in Python — the same aggregate over the same 66.3M
`(cell, record)` pairs, ~850× apart. Python's interpreter vectorizes pure maps
and scaled-product contractions (einsum), but this body is a contraction with
an `ifelse` containment predicate, which matches neither, so it falls to a
per-pair Python loop. The compiled evaluators handle it directly. Everything
else — the gated fetch, the SR contractions, the deaths math — is within
small factors across the three.

## Tolerances

`ppl` is compared **exactly**. Float fields use `RTOL_FIELD = 1e-12`.

That is not arbitrary: three engines contract a 1,520-term sum in different
orders, and reassociating a float sum changes the last bits. The measured
cross-binding spread is ≤ 3.1e-15. Asserting tighter would be asserting on
summation order rather than on the model.

One subtlety worth recording, because it nearly produced a false failure: the
record's `sum` must be a property of the *data*, not of the summing language.
Julia's `sum` is pairwise and CPython 3.12's is Neumaier-compensated, but
Rust's `Iterator::sum` is a naive fold — with **bit-identical** `total_pm25`
fields (same sha256, same samples) the Rust total differed by 2.9e-13 purely
from accumulation error. At full scale that could have exceeded the tolerance
and reported a disagreement between provably identical fields. All emitters
use compensated summation.

## History

Tag **`pushdown-era`** holds the previous state of this repo: three `.esm`
variants (`isrm.esm` / `isrm_clean.esm` / `isrm_pushdown.esm`), four
runner-mediated `run-model-*` directories whose pushdown constructs were
hand-authored in the document rather than derived by the engine, and the
hand-written STEP-0 Julia oracle (`mode=oracle_step0`) that first reproduced
the tutorial totals. Its four full-scale result records are frozen in
`contract/records/` and remain part of every comparator run — the current
tree's records are bit-identical to them where the spec demands it.
