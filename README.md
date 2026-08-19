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
* **In-model plume rise.** The ASME plume-rise algorithm — stack layer,
  buoyancy flux, the four-branch ΔH, the SR-layer clamp — is stated in the
  document's own observeds, and each emission record is charged to the SR
  emission layer its plume reaches. No runner computes a plume height; the
  shims do not contain the word ASME.
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
  `metadata.esio_format`, URL = `source.url_template`).
* **Document-declared ingest** (esm-spec §8.9; the Rust binding today, see the
  binding note below). The EGU loader carries its own `reader_options` (the
  `*egu*` member glob and the EPA column-header row), a `codes` map turning the
  POLID text column into the pathway enum, a `record_filter` that drops a
  record with no coordinate or no annual total, and an `extent` that binds the
  surviving count to `N_REC`. The source-cell rectangles are a `select` range —
  `W[0:N_SRC]` on their own loader variables — rather than a prefix the caller
  slices. The FF10 table is read, mapped, filtered, counted and (for a reduced
  run) truncated **by the engine**, from the declaration: `run-rs` names no
  pollutant, no column, no grid extent and no record count.

Each shim contributes orchestration only; every reported number is the
binding's evaluation of the document's observed graph.

> **Binding note.** The §8.9 ingest fields are consumed by **all three**
> bindings' `providers_from_document` + `prepare`, and all three shims dropped
> their `read_egu` in this change. Each was checked the same way: the record a
> binding emits with the imperative ingest deleted, against the one it emitted
> with it, and against the other bindings'. `ppl` matches exactly — count and
> sha256 — the totals are bit-identical within a binding, and the bindings
> agree with each other to 5.1e-15.

## The result

Full scale: 52,411 receptor cells × 1,520 emission-bearing source cells ×
43,650 EGU FF10 emission records, against the live
`s3://inmap-model/isrm_v1.2.1.zarr` — and now with each record charged to the
SR emission layer its plume reaches, which makes the contraction fifteen SR
slabs (five pathways × three layers) rather than five.

The tutorial's published totals are `sum(deathsK) = 6928.959583` and
`sum(deathsL) = 15623.924632`. A correct run lands **near** them rather than
on them, by construction — see [Plume rise](#plume-rise) for the one
deliberate difference and what it costs.

| binding | wall | `sum(deathsK)` | vs the blog |
|---|---|---|---|
| Julia | 6,316 s (1,816 s prepare/build + 4,437 s observed eval) | 6983.9385617781645 | +0.79% |
| Python | — | not re-run at full scale since plume rise | |
| Rust | — | does not build today | |

`sum(deathsL) = 15752.315804140908`, +0.82%. Peak RSS 16.6 GiB under
`julia -t 2 --heap-size-hint=12G`. The two missing rows are not modesty about
the numbers; both bindings are currently unable to run this document at all,
for reasons unrelated to plume rise — see [What runs today](#what-runs-today).

`RTOL_ORACLE` in
[`contract/compare_results.py`](contract/compare_results.py) is set from these
two measurements and nothing wider.

**The residual is attributed but not proven.** The deliberate deviation below
is the expected cause, its sign is right, and the argument for why 0.43% of
misplaced mass buys ~0.8% of deaths is written out in that file — but nothing
here *measures* the split between that cause and any other. Two candidates are
untested: the hosted ISRM matrix may have changed since the 2019 tutorial, and
the tutorial's own path writes its emissions through a shapefile, whose DBF
numeric fields could truncate the stack parameters that plume rise reads. The
decisive test is to implement the bug-compatible variant and see whether it
lands on `6928.959583` exactly; that was considered and deliberately declined,
because reproducing a defect is a worse artifact than documenting one. Until
someone runs it, read the +0.79% as *consistent with* the explanation below
rather than as evidence for it.

Two numbers came back **unchanged** from the ground-level-only era, which is
the check that says plume rise moved what it was supposed to and nothing else:
`ppl` is still `6f784d7e…` over the same 1,520 cells, and all five per-pathway
`emis_sum` totals are identical to the last digit. Plume rise changes which
*layer* a record emits into, never which *cell*, and moves mass between layers,
never into or out of a pathway.

What **is** verified today is the physics the full-scale numbers rest on. The
record → SR-layer assignment is checkable without the 330 GB SR matrix — it
needs thirteen 1-D meteorology arrays, ~16 MB — and
[`contract/plume_oracle.py`](contract/plume_oracle.py) computes it
independently, from InMAP's Go source rather than from this document. Every
runner now reports the same quantity in its record's `plume` block, straight
off the document's own `plume_layer` observed, and the comparator checks the
two **exactly**:

| truncation | what it exercises | `plume.sr_layer` sha256 | vs the oracle |
|---|---|---|---|
| `ISRM_FIRSTN=200` | 2 of ASME's 4 branches; no zero-height stack, nothing above model layer 7 | `cd9796b8…` | exact match |
| `ISRM_FIRSTN=2000` | all 4 branches, 145 zero-height stacks, 286 above-layer-7 records | `574b3c3f…` | exact match |

Julia, both truncations: the per-record SR-layer digest **and** the per-record
stack-layer digest match the oracle exactly, and the fifteen per-layer emission
masses agree to ≤ 2.7e-16. Seeding a single record into the wrong layer makes
the digest, the histogram and every pairwise comparison fail — the check has
teeth.

The three emitters agree on the wire format too: `results.jl`, `results.py` and
`run-rs/src/contract.rs` produce byte-identical `sr_layer` / `stack_layer`
digests from the same input, and all three refuse a non-integral value rather
than rounding it away.

What is **not** yet closed: the full-scale run in the table above was launched
before the `plume` block existed, so it never emitted one, and the full-scale
`sr_layer` digest has therefore not been checked against
`contract/records/plume_oracle.json`'s
`808e0971a2eda1de1ffc53e242f7ea3fd9bbbda85c3b61702a56f71dd12b434b`. The next
full-scale run closes that, at no extra cost — the block comes off observeds
the run already evaluates.

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

Plume rise does not touch it. It changes which *layer* a record emits into,
never which *cell*, so that digest is as much a check on the current document
as it was on the pushdown-era one — a change there would be a bug, not the new
physics.

## Plume rise

A power plant does not emit at ground level. Its plume leaves a stack tens to
hundreds of metres up, and buoyancy and momentum carry it further — often
into the next model layer or the one above that. The ISRM matrix has three
emission layers for exactly this reason, and the InMAP tutorial uses all
three.

`isrm.esm` states the rise itself. Per emission record, in the document's own
observeds: `stack_layer` (which model layer the stack top sits in, hence which
cell's meteorology applies), `buoy_flux`, the four-branch `delta_h` —
momentum-dominated, stable-buoyant, unstable-buoyant, or none when `F <= 0` —
then `plume_height`, and finally `plume_layer`, the SR emission layer the
plume reaches. The `E_<pathway>_L<layer>` aggregates bin each record's mass
into that layer. The algorithm is `ctessum/atmos/plumerise`
(`ASMEPrecomputed`, `calcDeltaHPrecomputed`, `findLayer`) plus InMAP's
`IsPlumeIn` and `sr.Reader.layerFracs`, transcribed as expressions, not as
code in a runner.

It moves most of the mass upward. At full scale, per pathway, across SR layers
0/1/2:

| pathway | L0 | L1 | L2 |
|---|---|---|---|
| VOC | 6.3% | 13.1% | 80.5% |
| NOx | 3.1% | 5.3% | 91.5% |
| NH3 | 12.3% | 6.6% | 81.2% |
| SOx | 0.5% | 4.3% | 95.2% |
| PM2.5 | 3.5% | 5.8% | 90.7% |

Mass emitted higher up disperses before it reaches anyone, so deaths fall:
the ground-level-only document reported `sum(deathsK) = 7524.9188…`, and the
tutorial reports `6928.959583`.

### The one place this document does not follow InMAP

654 of the 43,650 records — 0.43% of emitted mass — have plumes that rise
above the top of model layer 7, out of the 52,411-cell ground grid and into
the 9,324-cell high-altitude grid.

InMAP has a latent defect there. `sr.Reader.layerFracs` clamps such an
emission into SR layer 2, which is right — model layers 3–26 all clamp. But
the horizontal index it pairs with that layer is `sr.indices[c]`, and
`srreader.go` builds that map by resetting the counter to 0 at every layer
boundary. For a cell in layer ≥ 8 the index is therefore a position in the
**coarse** 9,324-cell grid, which `Reader.source` then reinterprets as a
position in the ground grid — so those emissions are charged to the wrong
source cell. Not a rounding difference: a plant in Ohio can end up emitting
somewhere else entirely.

**This document implements the correct behaviour** — clamp to SR layer 2 at
the cell the emission actually came from — and is therefore *not*
bug-compatible with InMAP. That is a deliberate choice, and it is why the run
lands near `6928.959583` rather than on it: **+0.79% on `deathsK`, +0.82% on
`deathsL`**, measured at full scale.

That is nearly twice the group's 0.43% share of emitted mass, which is worth
saying plainly rather than hiding behind the mass number. The mass is not lost
in either model — it is *placed* differently, and placement is what decides how
many deaths a ton causes. Put back on the cells the emissions actually came
from, it sits over power plants, which sit near people; scattered across the
ground grid by a coarse index read as a fine one, it does not. So the group
punches about twice its weight, in the direction the argument predicts: this
document's totals come out **higher** than the blog's, because it puts that
mass back where the people are.

The `above_model_layer_7` block in `contract/records/plume_oracle.json` counts
and quantifies the group, per pathway, so the residual is attributable rather
than mysterious. `contract/compare_results.py`'s `RTOL_ORACLE` is set from
those two measurements and nothing wider; if a future run drifts past it, the
cause is something else and worth finding rather than accommodating.

## The projection pushdown

The full SR matrix is 52,411 × 52,411 per pathway — about 330 GB across five
pathways. None of it is downloaded.

The rewrite runs *first* and derives the 1,520 members; the SR fetch is then
built from them and pushed down to the zarr reader, which fetches only the
intersecting chunks:

```
gated SR fetch: 1520 of 52411 source cells, all receptors
```

Plume rise triples the gate's reach: fifteen SR arrays (five pathways × three
emission layers) instead of five. The 1,520 members fall in 416 of the 525
source chunks, so the fetch grows from ~43 GB to ~128 GB uncompressed and the
materialised SR from 1.6 GB to ~4.8 GB (measured while planning the change).
That is the dominant new cost and it is not avoidable — the emissions
genuinely land in all three layers.

The gate also drives *enumeration*, not just filtering, and — since the
EarthSciAST change that made plume rise expressible — it does so for **any**
aggregate carrying an `overlap` join, not only for the `distinct` producer
that derives the support set. That matters more than it sounds: after the
rewrite each `E_*` binning aggregate ranges over
`(1,520 support cells × 43,650 records)` and used to visit all 66.3M pairs to
find the roughly 43,650 that a point-in-rectangle containment actually admits
— about one per record. The engine now walks the broad-phase candidates
instead. See [Timing](#timing) for what that was worth.

## Layout

```
isrm.esm          the model — one document, no variants
contract/         the shared result record: schema, the three emitters
                  (results.jl, results.py, run-rs/src/contract.rs), the
                  comparator, and plume_oracle.py — the offline record ->
                  SR-layer check on plume rise (see contract/README.md)
contract/records/ plume_oracle.json, the SR-free reference the plume block is
                  checked against, and ground-level-only/ — the four
                  pre-plume-rise full-scale baselines, retired with a README
                  saying why they are kept
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
  * EarthSciAST at `53cc91f6` or later — all three bindings need `prepare` +
    `pushdown_rewrite`, and plume rise additionally needs the overlap gate to
    drive enumeration on **ordinary** aggregates, not only on `distinct`
    producers (that is what lets the per-record `plume_layer` observed carry
    the same spatial join the `E_*` aggregates do);
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

# the SR-free plume-rise reference (~1 min, ~16 MB fetched, no venv)
python3 contract/plume_oracle.py

# compare the live records against each other and against that reference
python3 contract/compare_results.py contract/records/plume_oracle.json \
    run-jl/results.json run-py/results.json run-rs/results.json
```

Each shim writes `results.json` (or `results_reduced.json`) next to itself.

The four pre-plume-rise baselines now live in
`contract/records/ground-level-only/` and are deliberately **not** in that
command: every emission in them is charged to layer 0, so they are answers to
a different model and the comparator will say so. They are kept, and still
asserted on, for the two numbers plume rise cannot change — see the README
there.

* `ISRM_FIRSTN=200` truncates the emission-record list for a fast reduced
  run (~2 min in Julia after precompile). Check it against
  `python3 contract/plume_oracle.py --firstn 200`, which truncates the same
  way. Worth knowing: the first 200 records exercise only two of ASME's four
  branches and contain no zero-height stack and nothing above model layer 7 —
  `ISRM_FIRSTN=2000` reaches all four branches, 145 zero-height records and
  286 above-layer-7 records, and still runs in minutes.
* `ISRM_SCRATCH` **must be disk-backed**. On this cluster `/tmp` is a tmpfs:
  SR chunk blobs written there consume the same memory cgroup the model
  needs. The shims default to `/scratch.local/$USER/isrm-esm` for exactly
  that reason.
* Size `--heap-size-hint` to what is actually free in the memory cgroup —
  here a 40 GB cgroup shared with other jobs, hence 12G. It is a hint, not a
  cap: plume rise triples the materialised SR, and the full-scale Julia
  `prepare` peaked at 16.6 GiB RSS under that same 12G hint (it was 10.6 GiB
  before).
* `ISRM_ESIO_CACHE` can point one shim at another's cache (the format is
  cross-language) to skip the ~6-minute cold SR fetch.

### What runs today

Written 2026-08-19, and worth checking before trusting it: of the three
shims, only Julia has been driven end to end since plume rise landed.

* **Rust does not build.** `earthsci-ast-rs` calls
  `earthsciio::DataLoader::reader_options`, which exists on no EarthSciIO
  branch. The break predates this work and is being fixed elsewhere. The Rust
  emitter's `plume` block is written, and it type-checks and produces digests
  byte-identical to the Python emitter's on the same input — but no Rust
  record has been produced from it.
* **Python stops inside `prepare`.** The build-time hoist drops `stack_layer`
  with *"join 'overlap' envelope factor 'src_W' is not bound as build-time
  const-array data"*, and every downstream observed goes unresolved with it.
  A reduced Python run of this same document succeeded a few hours earlier, so
  this looks like a binding in motion rather than a document error — but that
  is a guess, not a diagnosis.
* **The document's unit strings are outside the ESM unit table.** With `pint`
  installed, the Python binding's structural validation rejects `ton/yr`,
  `ft`, `ft/s` and `(m/s)^-1/3` (esm-spec §4.8.1). Julia does not enforce it,
  so this surfaces in one binding only — but the document, not the binding, is
  what is out of line, and `ton/yr` predates plume rise.

## Timing

The old wall-clock spread — Python four hours against Rust eight minutes —
was reported here as a "known interpreter gap": the `E_*` emission-binning
join took ~2,800 s in Python and ~3.3 s in Rust, and the explanation offered
was that a contraction with an `ifelse` containment predicate matches neither
of Python's vectorized shapes and falls to a per-pair loop.

That explanation was mostly wrong, and the fix that came with plume rise
proves it. The join was not slow because each pair was expensive; it was slow
because the engine visited **every** pair — all 1,520 × 43,650 = 66.3M of them
— to admit the ~43,650 that a point-in-rectangle containment actually
contains. The `overlap` join clause was already in the document, and the
broad-phase machinery was already there, but the gate only drove enumeration
for `distinct` producers; on an ordinary aggregate Rust dropped it as inert,
Julia used it as a filter, and Python rejected it outright. Making it drive
enumeration everywhere (an EarthSciAST change, needed anyway so the per-record
plume observeds could carry the same join) turns the join from `|support| ×
|records|` into roughly one candidate per record.

Measured on the Python binding at `ISRM_FIRSTN=8000`, before and after:

| | five `E_*` joins | wall |
|---|---|---|
| before | 628 s | 16:30 |
| after | 2.2 s | 2:08 |

`sum(deathsK)` was identical to the last digit — the gate changes which pairs
are *visited*, never which are *admitted*.

A compiled-vs-interpreted gap surely remains, but it is now small enough that
it has not been separately measured at full scale, and the honest statement is
that the four hours were mostly this. Full-scale timings for all three
bindings under plume rise are still to be collected.

## Tolerances

`ppl` is compared **exactly**, and so are the `plume` block's two layer
assignments — the per-record SR layer and the per-record stack layer. All
three are integer-valued, and a float tolerance on an integer would hide a
real disagreement about which cell emits, or which layer it emits into. Float
fields use `RTOL_FIELD = 1e-12`.

`RTOL_ORACLE`, the tolerance against the tutorial's published totals, is a
different kind of number: it is not float noise but a physics difference, the
above-layer-7 group this document deliberately places correctly and InMAP does
not. It is set from the measured deviation and no wider, with the measurement
recorded in the comment beside it.

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
`contract/records/ground-level-only/`. They were part of every comparator run
until plume rise; they are not any more, because every emission in them is
charged to layer 0 and a live run no longer answers that model. What survives
the re-baselining is asserted on still: the `ppl` digest and the five
per-pathway emission totals, which plume rise cannot change and which
`contract/plume_oracle.py` checks itself against.
