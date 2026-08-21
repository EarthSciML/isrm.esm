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
  buoyancy flux, the four-branch ΔH — and a corrected `sr.Reader.layerFracs`
  on top of it are stated in the document's own observeds, and each emission
  record's mass is split across the SR emission layers its plume falls
  between. No runner computes a plume height; the shims do not contain the
  word ASME.
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
`s3://inmap-model/isrm_v1.2.2.zarr` — and now with each record's mass split
across the SR emission layers its plume falls between, which makes the
contraction fifteen SR slabs (five pathways × three layers) rather than five.

**At reduced scale the document now reproduces the live `inmap cloud` service
exactly.** On the first 200 emission records, both bindings:

| | `sum(deathsK)` | `sum(deathsL)` |
|---|---|---|
| this document, Julia | 49.09146956422705 | 110.40696810124838 |
| this document, Rust | 49.09146956422743 | 110.40696810124905 |
| **live `inmap cloud`, same 200 records** | **49.091470** | **110.406968** |

That is 8.9e-9 and 9.2e-10 relative — inside the precision the service prints.
Before this change the same run gave `47.780033` / `107.455964`.

**No full-scale run of the current document has been made.** The last
full-scale number in this file, `6063.777261048292` (−12.49% against the
blog's `6928.959583`), was produced by the *previous* document, the one that
clamped every record to a single SR layer, and is void as a statement about
this one. Full scale costs ~1.5 h and is the next thing to do; the reduced
agreement above is what says it is worth doing.

### Result

Full scale, Rust, against the live `s3://inmap-model/isrm_v1.2.2.zarr` (the
repaired and rechunked store; the values are bit-identical to `isrm_v1.2.1.zarr`
after its repair, verified by running both):

| | `sum(deathsK)` | `sum(deathsL)` | vs the tutorial |
|---|---|---|---|
| every record at ground level (before plume rise) | 7524.918846 | 16979.632171 | +8.60% |
| plume rise, reading the store's corrupt `layers` | 6063.777261 | 13668.309908 | −12.49% |
| plume rise, reproducing InMAP exactly | 6936.106343 | 15640.080273 | +0.10% |
| the tutorial, reproduced through InMAP's own service | 6928.959583 | 15623.924632 | — |
| **what this document computes — correct physics** | **7022.724781** | **15835.993596** | **+1.35%** |

Read those last two rows together, because the order they were produced in is
the argument.

**First, agreement.** The third row is this document reproducing InMAP,
including a bug (below). Its +0.10% residual is one deliberate difference,
cleanly isolated: at `ISRM_FIRSTN=200` — a subset containing no record whose
plume clears model layer 7 — it matches InMAP's live service to **8.9e-9**. At
full scale 654 records do clear it, and InMAP charges those to the wrong source
cell ([`BUGS.md`](BUGS.md) §5.1), worth exactly +7.146760 deaths.

**Then, correction.** Having shown it can compute what InMAP computes, the
document declines InMAP's other plume-rise defect too — the inverted
`layerFracs` interpolation ([`BUGS.md`](BUGS.md) §5.2), which misplaces 6.25% of
emitted mass. That is the last row, and it is what `isrm.esm` computes today.
**This document therefore no longer reproduces the tutorial, by design, and the
published totals run about 1.35% low.**

A document that diverged from its reference before it had ever matched it would
be indistinguishable from one with a bug. This one matched first.

Every defect found on the way here, and what each did to this number, is in
[`BUGS.md`](BUGS.md).

### The gap was a corrupt array in the published zarr — and it is now closed

The tutorial's numbers **are** exactly reproducible. Re-running the blog's own
path through the live `inmap cloud` service, with the same 43,650 records,
returns `6928.959583` / `15623.924632` — every printed digit.

The cause is one array. The SR store's `layers` variable holds the **model**
layer indices the SR calculation was performed for, and:

| source | `layers` |
|---|---|
| `isrm_v1.2.1.ncf` (Zenodo) | **`[0, 3, 6]`** |
| `isrm_v1.2.1.zarr` (as published, before the 2026-08-20 repair) | `[0, 1, 2]` |

`[0,1,2]` is a machine-generated arange that displaced the real variable during
the zarr conversion. Corroboration, independent of the NetCDF: `sr/sr.go`
attaches `description: "Layer indices for which the SR calculation was
performed"` to `layers`, and the zarr's copy has **no attributes at all** while
neighbouring InMAP variables like `LayerHeight` and `Dz` keep theirs. Zenodo's
description calls the three heights "ground-level, low-stack, and high-stack" —
which `[0,3,6]` (0–58 m, 253–391 m, 786–1049 m) is and `[0,1,2]` (three
adjacent near-surface layers, all under 245 m) is not. The same conversion pass
also lost the per-cell `Layer` variable and stripped `.zattrs` from the SR
arrays — which is why this repo has always needed `seed_empty_zattrs`.

**The matrix data itself is fine.** The Zenodo `.zip` is md5-identical between
its 2019-03-11 and 2019-12-22 depositions, and SR rows streamed out of it are
byte-identical to the zarr's. Nothing was rescaled; one index array was
replaced.

The document used to read `layers` as `[0,1,2]` and therefore assign
`sr_layer = min(plume_layer, 2)`. It no longer does. `[0, 3, 6]` is now a
**declared fact of the document** — the `SR_MODEL_L1` / `SR_MODEL_L2`
metaparameters, whose descriptions carry the provenance above — and is
deliberately not read from the store. Two consequences, both now implemented:

* InMAP's `layerFracs` splits an emission between **two** SR layers whenever
  its model layer falls between the entries of `layers` — model layers 1–2
  between SR 0 and 1, layers 4–5 between SR 1 and 2. That two-layer
  interpolation is **live**, not the dead code an earlier reading of this file
  claimed. The document states it as three per-record weights `w_sr0`,
  `w_sr1`, `w_sr2` summing to 1.
* It interpolates on cell-**centre** heights, not layer bottoms:
  `layerFracs` gets them from `sr.d.VerticalProfile`, which fills
  `height[i] = c.LayerHeight + c.Dz/2.`. This is worth **all** of the +0.38%
  residual that a planning scan left open — that scan interpolated on layer
  bottoms and predicted `49.279523` on the first 200 records; on centres the
  document lands on `49.09146956`, which is the live service's `49.091470`.
  `Dz` is now declared alongside the other meteorology arrays for exactly
  this.

Three other candidate explanations for that residual were on the table and are
now moot at reduced scale, since there is no residual left to explain. One of
them was cheap enough to close outright: InMAP's `CellIntersections` returns
0.5/0.25 fractions for a point that lands exactly on a cell boundary, and this
document's half-open `[W,E) × [S,N)` containment does not. Every ground-grid
edge is a multiple of 1 km in LCC metres, so a point on an edge needs
`x mod 1000 == 0`; **none of the 43,650 records satisfies it**, in either axis,
and the closest approach is 2.0 cm. The tie-break case never arises on this
inventory. The other two — the `end = nCells-1` bound in `Reader.get`, and
EPSG:4269→LCC reprojection versus the document's direct `lcc_forward` — are
untested and would have to be worth less than 9e-9 on the first 200 records.

An earlier revision of this README argued the published totals were
*unreachable* — 3.46% above a ceiling implied by non-negative plume rise. That
argument was sound given `layers = [0,1,2]` and worthless because the premise
was false. It is retracted, and `contract/compare_results.py`'s
"admissible span" check, which encoded it, is gone.

A second InMAP bug found on the way, which the document had to match in order
to *demonstrate* agreement and now deliberately does not: `layerFracs` returns
`{frac, 1-frac}` for `{lower, upper}` with
`frac = (plumeHeight − below)/(above − below)`, so a **higher** plume gets
**more** weight on the **lower** SR layer. The interpolation is inverted. The
document reproduced it through the reconciliation and **corrects it now that the
reconciliation is complete** — five of its descriptions explain the switch,
because a reader who arrives at the divergence without the history would
reasonably read it as a bug.

It conserves mass, but conserving mass is not the same as leaving the total
alone. **Measured at full scale: correcting the interpolation gives 7022.724781
/ 15835.993596 — which is what this document now computes — so InMAP's version
biases the published totals low by +1.249% / +1.253%**, 86.6 and 195.9 deaths a
year. It relocates 6.25% of emitted mass,
against 0.43% for the source-index defect below, which makes it the larger of
InMAP's two plume-rise bugs and the only one baked into the tutorial's numbers.
See `BUGS.md` §5.2 for the per-layer breakdown and why deaths move up.

`contract/`'s assertion that `layers == [0,1,2]` was asserting the corruption.
It is inverted now: the oracle checks that the layers it uses are never the
arange, that the store's own value is one of the two states there is evidence
for, and it records which one it saw. `layers` is the one upstream input from
this store that cannot be trusted.

> **An earlier full-scale total published here (`6983.9385617781645`, +0.79%)
> was VOID** and its apparent agreement was a coincidence. It came from an
> engine that contracted three of the five pathways against the wrong SR
> emission layer: three sibling SR loaders each expose `SOA`, `pNO3`, …, and
> the Julia and Python gated-fetch hook published every fetched slab under an
> alias resolved by matching the *tail* of a dotted name, so all three
> providers claimed all three keys and the last writer in hash order won. Rust
> registers one key per provider and was right all along; the disagreement
> between the two bindings is what exposed it. Fixed in EarthSciAST `98e6f1b6`.
>
> Worth keeping, because it is the argument for the whole apparatus: the
> per-record plume-layer digest matched the oracle **exactly in both bindings
> while the totals were wrong**. The physics was right and the contraction was
> not. A contract checking only the layer assignment would have passed this.
>
> `RTOL_ORACLE` in `contract/compare_results.py` was the value fitted to that
> void number. It is gone, along with the "admissible span" that briefly
> replaced it (see [Tolerances](#tolerances)).

Two numbers came back **unchanged** from the ground-level-only era, which is
the check that says plume rise moved what it was supposed to and nothing else:
`ppl` is still `6f784d7e…` over the same 1,520 cells, and all five per-pathway
`emis_sum` totals are identical to the last digit. Plume rise changes which
*layers* a record emits into, never which *cell*, and moves mass between
layers, never into or out of a pathway.

What **is** verified today is the physics the full-scale numbers will rest on.
The record → SR-layer-weights assignment is checkable without the 330 GB SR
matrix — it needs thirteen 1-D meteorology arrays, ~16 MB — and
[`contract/plume_oracle.py`](contract/plume_oracle.py) computes it
independently, from InMAP's Go source rather than from this document. Every
runner reports the same quantity in its record's `plume` block, straight off
the document's own `sr_lower` and `w_sr0`/`w_sr1`/`w_sr2` observeds:

| truncation | what it exercises | `plume.sr_lower` sha256 | vs the oracle |
|---|---|---|---|
| `ISRM_FIRSTN=200` | 2 of ASME's 4 branches; no zero-height stack, nothing above model layer 7 | `1c8a29f8…` | exact match |
| `ISRM_FIRSTN=2000` | all 4 branches, 145 zero-height stacks, 286 above-layer-7 records | `2d514e98…` | exact match |

Julia **and Rust**, both truncations: the per-record lower-SR-layer digest and
the per-record stack-layer digest match the oracle exactly, the three weight
fields agree with it to ≤ 1.2e-15 (bit-identically at 200), the fifteen
per-layer emission masses to ≤ 2.7e-16, and `max |w0+w1+w2−1|` is exactly 0 in
every run. `contract/compare_results.py` over oracle + Julia + Rust: **442
checks 0 failed** at 200, **438 checks 0 failed** at 2000. Seeding a single
record into the wrong layer makes the digest, the histogram and every pairwise
comparison fail — the check has teeth.

The three emitters agree on the wire format too: `results.jl`, `results.py` and
`run-rs/src/contract.rs` produce byte-identical `sr_lower` / `stack_layer`
digests from the same input, and all three refuse a non-integral value rather
than rounding it away.

Both are now closed by the full-scale run: it reports
`sr_lower` =
`d38ba2fb042f7e793134670e954d306987a8d9b17fca20975c23d36a9a134799`, matching
`contract/records/plume_oracle.json`, and lands +0.103% from the published
totals — a deviation fully accounted for by the one defect this document
deliberately does not reproduce (`BUGS.md` §5.1).

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

Plume rise does not touch it. It changes which *layers* a record emits into,
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
then `plume_height`, then `plume_model_layer`, the model layer the plume
reaches, and finally the three SR-layer weights.

Those weights are `sr.Reader.layerFracs`, and they are the reason a record is
not simply *placed* in a layer. The SR calculation was performed for model
layers `[0, 3, 6]` (the `SR_MODEL_L1` / `SR_MODEL_L2` metaparameters), and a
plume that stops between two of them is **split** across the corresponding SR
layers:

| plume's model layer | `w_sr0` | `w_sr1` | `w_sr2` |
|---|---|---|---|
| 0 | 1 | 0 | 0 |
| 1 or 2 | `frac_lo` | `1 − frac_lo` | 0 |
| 3 | 0 | 1 | 0 |
| 4 or 5 | 0 | `frac_hi` | `1 − frac_hi` |
| 6 or above | 0 | 0 | 1 |

`frac = (plume_height − below)/(above − below)`, with `below` and `above` the
**cell-centre** heights `LayerHeight + Dz/2` of the two bracketing model
layers over the record's own source cell — which is why `Dz` is in the
document at all. The bottom row is InMAP's `AboveTopErr` clamp, which
`plume_model_layer`'s saturation at 6 folds into the exact-match branch. The
`E_<pathway>_L<layer>` aggregates multiply each record's mass by `w_sr<layer>`,
so one record can feed two of them. The algorithm is `ctessum/atmos/plumerise`
(`ASMEPrecomputed`, `calcDeltaHPrecomputed`, `findLayer`) plus InMAP's
`IsPlumeIn` and `sr.Reader.layerFracs`, transcribed as expressions, not as
code in a runner.

It moves most of the mass upward. At full scale, per pathway, across SR layers
0/1/2:

| pathway | L0 | L1 | L2 |
|---|---|---|---|
| VOC | 26.1% | 54.4% | 19.5% |
| NOx | 22.4% | 53.5% | 24.1% |
| NH3 | 36.9% | 54.7% | 8.4% |
| SOx | 12.7% | 51.2% | 36.1% |
| PM2.5 | 21.9% | 51.9% | 26.2% |

Those are not the numbers this table used to hold — SOx read `0.5 / 4.3 /
95.2`. That was the corrupt `layers` at work: with `[0, 1, 2]` almost every
plume counts as "above the top layer" and gets shovelled into SR layer 2.

Mass emitted higher up disperses before it reaches anyone, so deaths fall:
the ground-level-only document reported `sum(deathsK) = 7524.9188…`, and the
tutorial reports `6928.959583`.

### The one place this document does not follow InMAP

654 of the 43,650 records — 0.43% of emitted mass — have plumes that rise
above the top of model layer 7, out of the 52,411-cell ground grid and into
the 9,324-cell high-altitude grid.

InMAP has a latent defect there. `sr.Reader.layerFracs` clamps such an
emission into SR layer 2, which is right — every model layer above 6 clamps.
But the horizontal index it pairs with that layer is `sr.indices[c]`, and
`srreader.go` builds that map by resetting the counter to 0 at every layer
boundary. For a cell in layer ≥ 8 the index is therefore a position in the
**coarse** 9,324-cell grid, which `Reader.source` then reinterprets as a
position in the ground grid — so those emissions are charged to the wrong
source cell. Not a rounding difference: a plant in Ohio can end up emitting
somewhere else entirely.

**This document implements the correct behaviour** — clamp to SR layer 2 at
the cell the emission actually came from — and is therefore *not*
bug-compatible with InMAP. That is a deliberate choice, and it is the one
reason a full-scale run should be expected to land *near* `6928.959583`
rather than exactly on it.

The mass is not lost in either model — it is *placed* differently, and
placement is what decides how many deaths a ton causes. Put back on the cells
the emissions actually came from, it sits over power plants, which sit near
people; scattered across the ground grid by a coarse index read as a fine one,
it does not. So the group should punch above its 0.43% weight, and in the
direction the argument predicts: this document's totals should come out
**higher** than the blog's.

**How much is not currently measured.** This section used to quote `+0.79%` on
`deathsK` and `+0.82%` on `deathsL`; both came from a full-scale run of the
*previous* document — the one that clamped every record to a single SR layer
under the corrupt `layers` — and neither is a statement about this one. The
figures to replace them come from the next full-scale run. What is measured is
that at `ISRM_FIRSTN=200`, where this group is empty, the document reproduces
the live `inmap cloud` service to 8.9e-9, which is the strongest evidence
available that the *rest* of the model is not contributing to any residual.

The `above_model_layer_7` block in `contract/records/plume_oracle.json` counts
and quantifies the group, per pathway, so whatever residual the full-scale run
shows is attributable rather than mysterious.

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
    producers (that is what lets the per-record `plume_model_layer`,
    `frac_lo` and `frac_hi` observeds carry the same spatial join the `E_*`
    aggregates do). The model's inline run configuration is spelled
    `analyses`, which is what the schema has called it since 2026-08-19;
    an older EarthSciAST expecting `examples` will reject the document;
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

Written 2026-08-19, revised 2026-08-20, and worth checking before trusting it.

* **Julia and Rust both run the layerFracs document at reduced scale**
  (`ISRM_FIRSTN=200` and `2000`), agree with each other and with the oracle,
  and hit the live service's totals. Neither has been run at full scale since
  the change. Python has not been run at all on this branch — the change to
  `run-py/run.py` is the same five-line rewiring the other two got, and it
  compiles, but no Python record has been produced from it.
* **The model's inline run configuration was renamed `examples` →
  `analyses` — FIXED 2026-08-20.** EarthSciAST's `esm-schema.json` moved on
  2026-08-19 and `Model` is `additionalProperties: false`, so the old spelling
  made the whole document fail validation at load, in every binding, for
  reasons unrelated to anything in it. The block's contents are unchanged.
* **Rust did not build (2026-08-19).** `earthsci-ast-rs` called
  `earthsciio::DataLoader::reader_options`, which existed on no EarthSciIO
  branch. Fixed elsewhere; `cargo build --release` against the canonical
  sibling checkouts succeeds and the binary has produced records at both
  truncations.
* **Python stopped inside `prepare` — FIXED 2026-08-19.** The build-time hoist
  dropped `stack_layer` with *"join 'overlap' envelope factor 'src_W' is not
  bound as build-time const-array data"*, and every downstream observed went
  unresolved with it. It was neither a document error nor a missing array: the
  rects are bound TWICE, under `ISRM_Grid.src_W` (the loader's variable) and
  `ISRM.src_W` (the model parameter coupled from it), and they are ONE array by
  reference. The interpreter's suffix resolver treated two keys as an ambiguity
  and reported a factor bound twice as "not bound". Only the pushdown rewrite's
  MIRROR arm — the per-record binning aggregates, `stack_layer` and its chain —
  takes that path, because it keeps the document's own full-grid rects instead
  of the compact `pd_cell__*` gathers the forward arm gets. The trigger was this
  document declaring its own `select` on sibling loaders, which is what created
  the second key. Fixed in `numpy_interpreter._scoped_array_name`.
* **The document's unit strings were outside the ESM unit table — FIXED
  2026-08-19, in the document AND in the registry.** Three separate things:
  `ton/yr` on the fifteen `E_*_L*` observeds and `kg/yr` on `emis_annual` were
  both wrong about the same column (`kg/yr` parsed and lied by 907.18474; `ton`
  is not a unit because it is three different masses), and both now read
  `short_ton/yr` — a new §4.8.1 entry, exactly 2000 international pounds, which
  is exactly InMAP's `907184740000` µg-per-short-ton constant. `ft` and `ft/s`
  are now §4.8.1 entries too, so the FF10 stack columns keep saying what they
  store. `(m/s)^-1/3` was a DOCUMENT bug and not a registry gap: §4.8.2 reads it
  as `((m/s)^-1)/3`, a scaling factor, and the bindings gave it three different
  meanings; it now reads `(m/s)^(-1/3)`, and a numeric scaling factor is a hard
  error in every binding. **None of this moved a number.** Measured, not
  assumed: on the engine of the moment, Julia's `ISRM_FIRSTN=200` totals were
  bit-identical with the old unit strings and the new ones
  (`sum(deathsK) = 50.92431255236616` both times, on the engine that then still
  carried the sibling-slab aliasing bug). A declared unit feeds validation, not
  the arithmetic; the totals moved later, and for an unrelated reason.

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

`ppl` is compared **exactly**, and so are the `plume` block's two integer
fields — the per-record *lower* SR layer and the per-record stack layer. All
three are integer-valued, and a float tolerance on an integer would hide a
real disagreement about which cell emits, or which layers it emits into. Float
fields use `RTOL_FIELD = 1e-12`, and that now includes the three layerFracs
weights: they are genuinely fractional, and their `sha256` is reported but not
asserted because they descend from `plume_height`, whose cube roots differ by
an ulp between languages. `max |w_sr0+w_sr1+w_sr2 − 1|` is checked against an
absolute `1e-12` in each record independently — layerFracs conserves mass
exactly, and every run so far reports 0.

**There is no tolerance against the tutorial's published totals right now, and
that is deliberate.** Two previous attempts are worth remembering:

* `RTOL_ORACLE` was fitted to a full-scale run that turned out to be wrong
  (the sibling-slab aliasing bug). A fitted number wearing a bound's clothes.
* The "admissible span" that replaced it was derived from `sr.layers ==
  [0, 1, 2]` — the corrupt array. Its upper edge, deathsK 6697.55, was a real
  ceiling for a model nobody was running. Valid argument, false premise.

What replaces both is `SERVICE_DEATHS`: what the live `inmap cloud` service
returns for the same truncated record list. At `n_rec = 200` that is
`49.091470` / `110.406968`, and both bindings hit it to 8.9e-9 / 9.2e-10.
`RTOL_SERVICE = 1e-7` comes from the service's own printing precision — six
decimals is ±5e-7 absolute, ±1.0e-8 relative at 49 — with an order of
magnitude of headroom, and not from the observed agreement. The published
full-scale pair is printed and not checked until a full-scale run of this
document exists to measure a tolerance from.

Float tolerance is not arbitrary either: three engines contract a 1,520-term sum in different
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
