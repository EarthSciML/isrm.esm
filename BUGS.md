# Defects found while adding plume rise, and what each did to this analysis

Adding plume rise to `isrm.esm` meant a model consumed inputs it had only ever
declared — stack heights, exit temperatures, the SR matrix's layer axis — and
that flushed out a long list of defects across four codebases and one published
dataset. This file records each one and, specifically, **what it did to the
national mortality estimate this repository computes.**

The estimate moved like this:

| | `sum(deathsK)` | vs the tutorial |
|---|---|---|
| every record at ground level (before this work) | 7524.918846 | +8.60% |
| plume rise, reading the store's `layers` | 6063.777261 | −12.49% |
| plume rise, `layerFracs` over the true `[0,3,6]` | 6936.106343 | +0.10% |
| the InMAP tutorial, reproduced through its own service | 6928.959583 | — |
| **what this document computes today — correct physics** | **7022.724781** | **+1.35%** |

The third row is the reproduction of InMAP, and its +0.10% residual is one
deliberate difference (§5.1). The last row is what the document computes now:
having demonstrated it could reproduce InMAP, it declines InMAP's other
plume-rise defect too (§5.2). **The document is no longer a reproduction of the
tutorial, and the tutorial's totals run about 1.35% low.**

Two claims made during this work were **wrong and are retracted**: that the
~8.6% ground-level gap "is the plume rise" (it is not — most of it was the
corrupt layer axis), and that the tutorial's totals were unreachable by any
layer assignment (a proof that was valid given the corrupt axis and worthless
because of it).

---

## 1. The one that dominated everything: a corrupt array in the published store

**`s3://inmap-model/isrm_v1.2.1.zarr`, variable `layers`.** FIXED IN THE PUBLISHED STORE on 2026-08-20 (see §6) — was, at the time of writing,
`REPAIR.md` once the patch lands; the fix needs write access to the bucket.

`layers` holds the **model** layer indices the SR calculation was performed for.
The authoritative NetCDF says `[0, 3, 6]`. The zarr says `[0, 1, 2]` — a
machine-generated arange that displaced the real variable during conversion.

**Effect here: −12.6 percentage points.** Reading `[0,1,2]` makes SR layer *i*
look like model layer *i*, so this document assigned `sr_layer = min(plume_layer, 2)`
and pushed almost all emitted mass into the top SR layer. With `[0,3,6]`,
InMAP's `layerFracs` splits an emission across **two** SR layers whenever its
model layer falls between entries — model layers 1–2 between SR 0 and 1, layers
4–5 between SR 1 and 2 — interpolating on cell-**centre** heights. Correcting
only this array moved `sum(deathsK)` from `6063.777261` to `6936.106343`. The
SOx mass split went from 0.5 / 4.3 / 95.2 % to 12.7 / 51.2 / 36.1 %.

A fourth defect in the same store, found while drafting the repair:
**every array declares `fill_value: 0`**, which xarray promotes to a CF
`_FillValue` — so every genuine zero in the matrix reads back as NaN. It is
latent only because the missing `.zattrs` currently stop xarray opening the
store at all; measured once that is fixed,
`ds['SOA'][0,0:100,:].values.sum()` is `nan`. Repairing the attributes without
also clearing `fill_value` would turn a store nothing can open into one that
opens and silently returns NaN, which is worse.

The same conversion pass also **lost the per-cell `Layer` variable** (a
case-insensitive collision with the 3-element `layer` dimension coordinate) and
**stripped `.zattrs` from the SR arrays**, which is why this repo has always
needed `seed_empty_zattrs` — a workaround that was evidence of a lossy
conversion sitting in plain sight for years.

The matrix *data* appears fine: Zenodo's zip is md5-identical across both 2019
depositions, and sampled SR rows streamed from it are byte-identical to the
zarr's. That sampling is weaker than first reported — the extraction tool used
had an off-by-one in its buffer accounting (`buf_start += cut` advancing further
than `del buf[:cut]` removed whenever targets were far apart), so three of the
eight rows originally compared were garbage produced by the tool rather than by
the store. Five rows stand. One index array was replaced; nothing was rescaled.

**Why it went unnoticed:** the conversion had no round-trip verification. And
`[0,1,2]` is exactly what a plausible arange looks like, so nothing about it
reads as wrong until you notice that Zenodo describes the three heights as
"ground-level, low-stack, and high-stack" — which `[0,3,6]` (0–58 m, 253–391 m,
786–1049 m) is and `[0,1,2]` (three adjacent layers, all under 245 m) is not.

---

## 2. Defects that produced wrong numbers in this analysis

### 2.1 Gated providers published each other's slabs
**EarthSciAST, Julia and Python.** Fixed.

`_fetch_gated_providers` published each fetched SR slab under
`_const_factor_aliases`, which resolves a bare name against every model variable
sharing the same dotted **tail**. This document uses three sibling SR loaders
(`ISRM_SR_L0/L1/L2`), each exposing `SOA`, `pNO3`, …, so all three providers
claimed all three keys and the last writer in `Dict` hash order won.

**Effect here: three of the five pathways were contracted against the wrong SR
emission layer.** It produced `6983.938562`, which was *coincidentally* within
0.8% of the tutorial and was briefly reported as a success. Rust was immune —
it registers one key per provider, looked up by that provider's own key — and
the Julia/Rust disagreement is what exposed it.

Not pre-existing: the pre-plume document has a single SR loader, so no two
providers could collide. This document's design made a latent defect reachable.

### 2.2 `unit_conversion` was declared everywhere and applied nowhere
**EarthSciAST + EarthSciIO, all three bindings, on the ESIO provider path.** Fixed.

`esm-spec.md:2220` is normative: the runtime applies it when producing values in
the declared units. It existed only on a typed variable-mapping path that
nothing on the `providers_from_document` route reaches.

**Effect here: stack heights arrived in feet and temperatures in Fahrenheit**
(max 1038 ft, max 3103 °F) into a plume-rise formula expecting metres and
kelvin. At reduced scale the layer histogram was `[0, 0, 200]` — every record
slammed into the top SR layer — instead of `[0, 32, 168]`. Invisible for as long
as those columns went unused.

### 2.3 Julia's unit parser silently dropped scale factors
**EarthSciAST, Julia.** Fixed.

`parse_units("(m/s)^-1/3")` returned `s·m⁻¹`, **dropping the `1/3`**, and
`1000/s` returned `s⁻¹`, dropping the 1000. One string, three meanings across
bindings, two differing only in magnitude — invisible to a dimensional checker.

**Effect here:** the document's `(m/s)^-1/3` (for `WindSpeedMinusThird`) was
accepted by Julia with the wrong meaning and rejected by Python. The document now
spells it `(m/s)^(-1/3)`, and a numeric atom other than `1` is a hard error in
all five bindings.

### 2.4 `fact` was a truncated constant
**This document.** Fixed.

`28766.639` where InMAP's `emisConversionFactor` gives
`907184740000/(3600*8760)` = `28766.639396245562`.

**Effect here: 1.4e-8 relative** — no effect on any conclusion, but it capped
agreement with any external reference at ~8 significant figures, which matters
when the whole point is bit-level cross-binding comparison.

---

## 3. Defects that were silent and would have produced wrong numbers

None of these changed a published number here, but each could have, without
raising anything.

### 3.1 Expression templates silently suppressed the pushdown rewrite
**EarthSciAST, all three bindings.** Fixed.

The pushdown detector treated a surviving `apply_expression_template` node as a
leaf, so a join body factored through a template was unrecognisable. The rewrite
then **silently did not fire** — no derived support set, no gate, and the full
330 GB SR matrix in the fetch plan. Demonstrated: factoring three of the fifteen
`E_*` bodies through one template gated 12 of 15 arrays with no error and no
warning. The consequence is a memory failure hours later, not a diagnostic.

### 3.2 Rust returned `0.0` for an out-of-range const-array gather
**EarthSciAST, Rust.** Fixed.

The homogeneous-Dirichlet zero-ghost convention was applied to *every* gather,
including const-array gathers, contrary to `CONFORMANCE_SPEC.md:844-856`.
Measured: an off-the-end gather raised `BoundsError` in Julia, `IndexError` in
Python, and returned `[0.0, 0.0, 0.0, 0.0]` in Rust. This document's met gather
computes a flat offset into a 596,444-element array; an off-by-one would have
been a hard error in two bindings and silently zeroed emissions in the third.

Worth noting for anyone fixing something similar: the silent zero on the
measured path came from the **whole-array overlay's own ghost-0 fill**, not from
the obvious site. Fixing only the obvious one left the reproducer passing.

### 3.3 Python silently truncated complex results to their real part
**EarthSciAST, Python.** Fixed.

`^` evaluates as `a ** b`, so a negative base with a fractional exponent yields a
complex value, and the cast to float discarded the imaginary part with only a
`ComplexWarning`. Measured on `(-2.5)^0.333333333`: Julia raises, Rust gives
`NaN`, Python returned `0.6786044051723126` — a plausible wrong number, the worst
of the three outcomes. The plume-rise expression takes cube roots of a buoyancy
flux that is negative for cold plumes, so this was directly in the path; the
document clamps every fractional power with `max(·, 0)` regardless.

### 3.4 Reader options were swallowed rather than refused
**EarthSciIO, Python and Julia.** Fixed.

Every reader had a `**_` catch-all so a shared kwarg set could be handed to any
of them, so a mistyped option was silently ignored. A misspelled `member_filter`
meant the FF10 reader **read the entire zip** instead of the EGU members.
EarthSciAST already had tests asserting an `UnknownReaderOption` that existed
nowhere.

### 3.5 Python's `1/0` aborted the build instead of returning IEEE
**EarthSciAST, Python.** Fixed.

`CONFORMANCE_SPEC.md:689` is explicit — "division is IEEE, never an exception, so
every binding agrees" — and the interpreter raised `ZeroDivisionError`. Because
`ifelse` is eager, a guarded division `ifelse(d != 0, n/d, fallback)` aborted the
whole build on a branch the model never took. The binding also disagreed with
itself: `value_invention._vi_eval` had implemented IEEE division all along,
citing the same rule.

### 3.6 Three more cross-binding unit divergences
**Fixed.** Surfaced only once the registry pinned *scales* rather than dimensions:
Rust's `year` was a 365-day year (0.0685% short); Go lacked the `degree` alias
and falsely rejected a shared *valid* fixture; Julia never called
`Unitful.register`, so its own custom units threw `KeyError` on any conversion.

---

### 3.7 Splitting the document across files re-opened §3.1, three more ways
**EarthSciAST, all three bindings.** Fixed 2026-08-21.

§3.1 taught the pushdown detector to expand a surviving
`apply_expression_template`. That fix is complete only for a template declared
**inside the model**. Moving eleven shared bodies out into a template-library
file (`isrm_base.esm`, esm-spec §9.7.1) and importing them broke the build
again, in three places that had never had to agree before — and two of the three
failed exactly the way §3.1 did, silently.

* **Python and Rust handed the recogniser an unresolved document.** Both
  `prepare`s read the raw `.esm` JSON — `json.load` / `&serde_json::Value` — and
  pass it straight to `desugar_pushdown`. An import edge resolves at LOAD, which
  has not happened yet on that path, so the imported templates are not in scope,
  the binning body cannot expand, and the rewrite declines with §3.1's
  consequence precisely: no derived set, no gate, the whole SR matrix in the
  fetch plan, no diagnostic. Julia never had this because its `prepare` loads
  first and desugars the serialized form. Fixed by resolving the §9.7 machinery
  ahead of the rewrite in both bindings, gated on the document actually carrying
  an edge so a document without one still reaches the recogniser in the bytes it
  does today.

* **Python discovered loader extents on the wrong side of the rewrite.** Julia
  discovers first, then loads with the closed metaparameters; Python discovered
  *after* the rewrite. Harmless for as long as nothing in between folded a
  metaparameter — and the §9.7 resolution above does, which would have folded
  `N_REC` to its default of 0 and sized `emis_records` to nothing before the FF10
  loader had counted a row. Reordered to match Julia, with the gate/extent
  contradiction check moved to after the rewrite, where the derived gate exists.

* **Julia's build-time pre-passes had no expansion arm.** Four passes run at the
  `build_evaluator` front door, *before* the boundary where `_build_evaluator_impl`
  expands references: the binning-coordinate derivation, value invention, the
  `member_factor` feedback and the overlap-env derivation. They evaluate bodies
  through `_eval_cellwise`, which has no `apply_expression_template` arm. So the
  LCC projection — projected at build time by that first pass — died with
  `E_TREEWALK_UNSUPPORTED_OP: apply_expression_template` the moment it was
  factored into a template, while the same document built in Rust and Python.
  This one at least failed loudly. Fixed by expanding once for those passes
  (§9.6.4 rule 2: a reference denotes its expansion, and every consumer may
  expand); the impl keeps expanding its own copy with its sites.

The shape is worth naming: **each binding resolves the document at a different
depth before handing it to the same pass**, and nothing tests that they agree,
because until now every document was one file.

### 3.8 A library cannot own the axes, because two consumers read them raw
**EarthSciAST.** Not fixed — worked around, and the reason is recorded here.

esm-spec §9.7.5 lets a template-library file carry `index_sets` and
`metaparameters`, which is exactly where the ISRM grid's axes belong. They could
not move. Two consumers read those blocks out of the RAW document, before §9.7
resolution: the esio bridge's `providers_from_document`, which folds a loader
select's `range.stop: "N_SRC"` while building providers *ahead of* `prepare`, and
the pushdown prepass. The first fails outright —

```
ISRM.src_E.update.from.select: range.stop names "N_SRC", which is not a
metaparameter with an integer default
```

— and it cannot simply resolve first, because at that point the FF10 loader has
not discovered `N_REC` and a resolve would default it to 0. So `isrm_base.esm`
is templates only, and `src_cells` / `rcv_cells` / `all_cells` / `emis_layer` /
`pathways` and their six sizing metaparameters are declared in each consuming
document — about 2 KB duplicated per geometry sibling.

## 4. Defects in the engines, not in the physics

Two cost only performance (§4.0, §4.1); two stopped a binding from running at
all (§4.2, §4.3). None of them changed a number that came out.

### 4.0 Every warm cache hit pays a network round-trip — FIXED

**`EarthSciIO/earthsciio/validate.py`, `decide()`.** The freshness ladder is
first-wins, and it is ordered:

```python
if expected_checksum: ...          # 1. content hash
if manifest.etag or manifest.last_modified:
    return REVALIDATE              # 2. conditional GET
if temporal is None or temporal.immutable:
    return HIT                     # 3. never reached for an S3 store
```

**S3 always returns an ETag**, so rule 2 always fires and rule 3 is dead code
for every S3-backed store. A cached chunk that is already on local disk issues a
conditional GET to us-east-2 before it can be used — including for a source
explicitly declared immutable, which is precisely the declaration that should
make the round-trip unnecessary.

Measured on the ISRM store, 60 chunks, cache already warm:

| | per chunk |
|---|---|
| as shipped | **85.9 ms** |
| `immutable` checked before validators | **0.078 ms** |
| (raw `open()` + `read()` of the same file) | 0.078 ms |

**1,100x**, and the fixed cost falls to exactly the file-read floor. At ~6,240
SR chunks for a full run against the 100-row store that is ~536 s of the 704 s
PREPARE — most of the run is waiting on revalidation of data it already has.

It also silently distorts every chunk-size decision, which is how it was found.
The 5-row store (§6) read 11.4x fewer bytes and still ran 28% slower, because
3x more chunks meant 3x more round-trips; bytes were never the bottleneck. With
this fixed the smaller chunking wins on both axes, measured: 233.1 s against
262.2 s, and the document reads it again.

**It needed a second fix to be safe**: `temporal` was declared in `.esm` sources
but never passed to `DataSource` by `providers_from_document`, so no document
could say `immutable` at all — and worse, once an absent `temporal` means
immutable, a dropped cadence stops being ignored and starts pinning stale bytes
permanently. The two had to ship together, and did.

**Both fixed 2026-08-20**: EarthSciIO `fec6875` reorders the ladder (spec §4,
Python and Rust in lockstep, regression guards in both languages); EarthSciAST
`ca10f1214` passes `temporal` through, and teaches `_to_esio_temporal` to read
the raw-JSON mapping shape `providers_from_document` actually holds — reading
only attributes had made a declared cadence degrade to CONST indistinguishably
from a legitimate "no anchor".

**Rust caught up 2026-08-21** (EarthSciAST `b41ceefba`): it had no ISO-8601
converter at all, so its bridge built every `DataSource` with no cadence. It now
parses the duration and the instant itself — deliberately by hand, against the
same mean-Gregorian constants Python's `approximate_seconds` uses (365.2425 /
30.436875 days), because the two bindings must resolve one document's `P1M` to
the same second. Cross-checked spec by spec against the Python parser, including
which malformed spellings each rejects. Two spellings that are NOT errors: a
block with no `start` stays CONST (nothing anchors the schedule), and either
`frequency` or `file_period` alone fills in for the other. A block that anchors
and then names neither IS an error, because reading its first file forever would
be a wrong answer rather than a slow one — the one place where refusing to run
beats guessing.



### 4.1 The overlap gate did not drive enumeration on ordinary aggregates
**EarthSciAST, all three bindings.** Fixed.

The broad-phase driver ran only on `distinct` producers. Every rewritten `E_*`
aggregate still visited all 1,520 × 43,650 = 66,348,000 pairs, where a
point-in-rectangle broad phase yields roughly one candidate per record.

**Effect here:** measured on the Python binding at `ISRM_FIRSTN=8000`, the five
`E_*` joins went **628 s → 2.2 s** and wall clock **16:30 → 2:08**, with
`sum(deathsK)` identical to the last digit. Visit counts fell 20,000× on a
record-oriented aggregate and 200× on a forward one. The README's "known
interpreter gap" was largely this.

### 4.2 The Rust binding had not built for months
**EarthSciAST + EarthSciIO.** Fixed.

`earthsci-ast-rs` called `earthsciio::DataLoader::reader_options`, which existed
nowhere in EarthSciIO on any branch, from EarthSciAST commit `8b473947` onward.
It broke a second time the same way: `aaa2fd0c0` renamed the import to
`earthsciio::DataSource`, a type the published crate did not have either.

**Fixed 2026-08-20.** The name is settled as `DataSource` (and `SourceTemporal`)
in EarthSciIO `7a2788f`, shipped as 0.1.2 with the 0.1.1 spellings kept as
deprecated aliases so a consumer is not stranded mid-upgrade. `run-rs/Cargo.toml`
carries a `[patch.crates-io]` onto the sibling checkout, because `earthsci-ast`
resolves EarthSciIO from crates.io and the shim reads the same types the bridge
does — they must be one crate, not a registry copy plus a path copy.

**Effect here:** the repo's headline three-binding claim could not be reproduced
for one of the three, and Rust had never executed the plume-rise document — which
is also why §2.1 went undetected until it did, and why §4.3 below survived the
whole 1.0.0 migration.

### 4.3 The esio bridge read 0.x source variables, so it built no providers
**EarthSciAST (Rust binding).** Fixed.

`providers_from_document` looked up `data_sources[l].variables` — the map esm
1.0.0 deleted, replacing it with `update: {kind: "data", source, from}` on the
CONSUMING parameter. Against the migrated `isrm.esm` it matched nothing and
returned an empty provider list.

The failure surfaced nowhere near the cause: `prepare` went on to evaluate the
graph with no data, an aggregate ranged over a zero-length axis, and the run
died inside `col_major_to_arrayd` with `ShapeError/IncompatibleShape`. It named
neither the source, nor the observed, nor the shape.

It survived the migration because §4.2 meant the file had not compiled since the
rename, and EarthSciAST's CI never builds `--features esio` (it cannot: the
feature resolves EarthSciIO from a sibling checkout that does not exist on a
runner). The 1.0.0 fixture and the Python twin test had both moved; the Rust
test still poked `data_sources[X].variables[...]`, so it could not have caught
this even if it had run.

**Fixed 2026-08-20** in EarthSciAST `d1ceb5bdb`, which also ports the 16 ingest
tests to 1.0.0, resolves `record_filter.require_finite` in the reader's file-
variable vocabulary, and makes the `col_major_to_arrayd` panic name the shape,
its product and the element count.

---

## 5. Defects in InMAP itself

These are upstream and are **not** fixed in InMAP by this repo. Neither is
reproduced here: `isrm.esm` states correct physics. §5.1 was never reproduced;
§5.2 was, deliberately, until agreement with InMAP had been demonstrated, and
was corrected on 2026-08-20.

### 5.1 High plumes are charged to the wrong source cell — not reproduced
`sr/srreader.go` builds `sr.indices` by resetting a counter at every layer
boundary, so a cell in model layer ≥ 8 gets a position in the **coarse**
9,324-cell grid, which `Reader.source` then reads against the 52,411-cell ground
grid. Those emissions are charged to the wrong source cell — not a rounding
difference; a plant in Ohio can emit somewhere else entirely.

**Effect here: +7.146760 deaths (+0.103%), and it is the entire remaining gap.**
This document charges them to the cell the emission came from. The attribution is
clean rather than argued: at `ISRM_FIRSTN=200` — a subset containing **no**
records whose plume clears model layer 7, so the defect cannot fire — this
document matches InMAP's live service to **8.9e-9**. At full scale, 654 records
trigger it and the totals differ by exactly this.

### 5.2 `layerFracs` interpolates backwards — reproduced, then corrected
```go
frac := (plumeHeight - below) / (above - below)
return []int{i, i + 1}, []float64{frac, 1 - frac}, nil
```
`frac` rises with plume height and is applied to `i`, the **lower** SR layer. So
a plume at the lower reference puts all its weight on the upper layer, and vice
versa. Correct interpolation gives the lower layer `1 - frac`; the two weights
are simply swapped.

It conserves mass, so no emission is created or destroyed, and the reference
heights themselves are right (cell centres, matching how `sr/sr.go` places
sources at `LayerHeight + Dz/2`). But it systematically anti-correlates height
with placement, and it is not a rare path: **~61% of emitted mass** in this
inventory goes through it (29.6% in model layers 1–2, 31.7% in layers 4–5).

**Effect here: +86.618439 deaths (+1.249%) Krewski, +195.913322 (+1.253%)
Lepeule — the tutorial's published totals are biased low by that much.**
MEASURED at full scale 2026-08-20, by running this document twice through the
Rust binding against the same SR matrix, with only the four AST nodes that place
`frac_lo`/`frac_hi` swapped:

| | `sum(deathsK)` | `sum(deathsL)` |
|---|---|---|
| tutorial, published | 6928.959583 | 15623.924632 |
| this document, reproducing the defect | 6936.106343 | 15640.080273 |
| same document, interpolation corrected | 7022.724781 | 15835.993596 |

Mass conservation is not the same as an unbiased total, which is what an earlier
revision of this section assumed. The defect relocates **192,918 t/yr — 6.25% of
emitted mass**, 14× more than §5.1 and, unlike §5.1, fully baked into the
published numbers. It is the larger of InMAP's two plume-rise defects on both
counts.

The direction is not the one "inverted" suggests. Correcting it moves mass
*inward* toward the middle SR layer from both sides — L1 +192 kt, L0 −54 kt,
L2 −138 kt — because giving each split record more weight on the layer *farther*
from its plume pushes mass outward from L1 in both directions. Deaths rise
because the L2→L1 flow (138 kt, mostly SOx and NOx) outweighs the L0→L1 flow
(54 kt): more mass ends up nearer the surface, where it reaches people. Neither
`ppl` (1,520 emitting cells) nor the `sr_lower` assignment changes — only how a
split record's mass is divided between two layers.

Cross-checked before any SR contraction: the corrected run's own weight sums and
`sr_lower` digest (23532.168795934616 / 17289.645044880024 / 2828.1861591853617,
`d38ba2fb…`) match an independently patched NumPy oracle exactly, so the engine
and the reference implementation agree on what "corrected" means.

**As of 2026-08-20 this document CORRECTS the defect and no longer reproduces
it.** The order matters and is the whole argument: agreement with InMAP was
established *first* — 8.9e-9 against the live service on 200 records, +0.103% at
full scale, the exact worth of the one defect already declined — and only then
was the physics corrected. A document that diverged before it had ever agreed
would be indistinguishable from one with a bug.

So `isrm.esm` now declines **both** of InMAP's plume-rise defects, and computes
**7022.724781 / 15835.993596**. It is no longer a reproduction of the tutorial
and does not claim to be; `contract/compare_results.py` demotes the published
totals to context and checks `CORRECTED_FULL` instead. Five variable
descriptions now say why, where four used to say the opposite.

---

## 6. Still open

* ~~**`temporal` is ignored** by `providers_from_document`.~~ **DONE** — Python
  2026-08-20 (`ca10f1214`), Rust 2026-08-21 (`b41ceefba`). See §4.0. Nothing in
  this document moves (its five sources are all static, and the reduced record
  re-runs bit-for-bit), which is exactly why it went unnoticed for so long.
* **EarthSciAST's CI never builds `--features esio`** — the feature resolves
  EarthSciIO from a sibling checkout that is not on a runner, so the whole
  provider bridge is compiled by nothing but a human running `isrm.esm`. That is
  how §4.2 and §4.3 each survived for months. Now that EarthSciIO publishes the
  names the bridge imports (0.1.2), a CI job could build the feature against the
  registry instead.
* ~~**The pushdown rewrite re-points only the gate ENVELOPE factors.**~~
  **FIXED 2026-08-21**, all three bindings. When the rewrite compacted a binning
  aggregate's cell axis onto the derived support set, its rect map was built
  from `tgt_env` alone. Any OTHER array the body indexed by the cell symbol was
  left pointing at the full grid while the axis was now compact — full-grid
  values read at support positions, **wrong numbers with no diagnostic**;
  `_pd_assert_rects_rebound` checked the envelope factors only. Nothing in
  `isrm_point.esm` hit it, because every binning body there reads exactly the
  four rect factors; `isrm_polygon.esm` hits it on its first line, since its
  allocation weight is `polygon_intersection_area(cell_ring[c], rec_ring[k])`
  and `cell_ring` is a rank-3 `[cells, verts, 2]` array that is not an envelope
  factor. The gather family is now **every** array whose declared `shape[0]` is
  the cell set and that the body subscripts with the cell symbol, and each
  gather is **rank-preserving** — `pd_cell__<C>__<F>` is
  `[pd_support__<C>, …F's trailing axes]`, defined by a map whose `output_idx`
  names one generated symbol (`pd_t0`, `pd_t1`, …) per trailing axis. So the
  sliced polygon-operand spelling `index(F, c)` survives the substitution of the
  name unchanged, and `_pd_assert_rects_rebound` now covers the whole family.
  A cell-axis array read at a COMPUTED position (`index(F, c + 1)`) is refused
  with a hard error naming the array and the subscript — that one cannot be
  re-pointed at all, and neither declining silently (an ungated fetch, §3.1) nor
  emitting anyway (wrong numbers) is acceptable. Normative in CONFORMANCE_SPEC
  §5.5.7 "Cell-axis arrays"; pinned by the `pushdown/polygon_area` golden and by
  `pushdown_cell_geometry` tests in all three bindings, each comparing a gated
  compact run against a dense one.
* **Three Julia engine gaps sat behind it**, all found by making the fix
  evaluate rather than merely emit, all fixed 2026-08-21. (1) `and`/`or`/`not`
  lacked the `geo` registry flag, so a containment predicate — the canonical
  `ifelse(and(cmp, cmp, …), …)`, which every pushdown fixture uses — could not
  appear in a setup-time geometry body at all, though the runtime ladder these
  lower to already had the arms. (2) Value-invention extents were merged into
  `derived_extents` AFTER `_materialize_geometry_setup` ran, so a geometry body
  ranging over a derived support set could not resolve its own extent. (3)
  `_geo_index_extent` did not follow a `kind:"derived"` set through its
  `from_faq`, which is how `derived_extents` is keyed. Rust needed none of these
  — its dense evaluator computes the geometry inline.
* **A nested `aggregate` is not a portable geometry operand.** Building the cell
  ring inline inside the binning body was one of the three escapes probed while
  the item above was open; the rank-preserving gather makes it unnecessary, but
  the divergence is real and still there. Python evaluates it correctly (the
  inner aggregate does see the enclosing cell symbol). Julia refuses:
  `E_TREEWALK_GEOMETRY_SETUP: operand must be a build-time array (a const/setup
  array name, an index slice of one, or an intersect_polygon clip)`. So the
  spelling that works in one binding is not the spelling a document may use.
* **An indirect subscript on a geometry operand is wrong in two bindings.**
  `polygon_intersection_area(index(cell_ring, index(members, p)), …)` — the
  cheap alternative to a new gather — returns silently wrong values in Python
  and throws `BoundsError … at index [0, …]` in Julia, which reads the member id
  but never converts 0-based to 1-based. Two separate defects, one of them
  silent.
* **Index symbols are not usable as scalar VALUES.** `ifelse(d < 0.5, …)` inside
  an `aggregate` whose `output_idx` carries `d` evaluates to garbage rather than
  erroring — it neither reads the loop index nor refuses. A ring built from
  scalars must route every index symbol through an `index(...)` position
  instead (const selector arrays: `isx[d]`, `fx[v]`), which does work.
* **`source.mirrors` are dropped** — declared failover URLs never reach the loader.
* **`determinism` is read by nothing**, in either repo, despite `esm-spec.md`
  §8.9 being a normative MUST.
* **`load()` enforcement diverges**: Python enforces §4.8.4 unit errors at load;
  Julia, Rust, Go and TypeScript do not, as `esm-libraries-spec` §2.1a requires.
  So an out-of-spec document runs in four bindings and fails in one. Python's
  non-conformance is what caught §2.3. This needs a spec decision about which
  entry point owns §4.8.4, not an engine patch.
* ~~**The zarr repair** itself, which needs bucket write access.~~ **DONE
  2026-08-20.** `s3://inmap-model/isrm_v1.2.1.zarr` was repaired in place: it
  now serves `layers = [0, 3, 6]`, carries the restored per-cell `Layer`
  variable, has `.zattrs` on all three SR arrays that lacked them, and declares
  `fill_value: null` on all 70 arrays (§1's fourth defect — with the `.zattrs`
  restored but `fill_value: 0` left in place, xarray would have masked every
  genuine zero to NaN, so those two had to ship together). Root attribute
  `store_version` is now `1.2.1+repair1`; every overwritten byte is backed up.
  Verified against the live store: 39 checks, 0 failed, 0 skipped, plus an
  independent anonymous xarray read.
* A **rechunked and recompressed** copy was published alongside it as
  `s3://inmap-model/isrm_v1.2.2.zarr` — same values, `[1, 5, 52411]` chunks and
  zstd instead of `[1, 100, 52411]` and lz4, plus consolidated `.zmetadata` so
  the store can be opened over plain HTTPS at all. The document reads it.
  It first looked like a trade rather than a win — 3.5 GB fetched instead of
  40 GB, but PREPARE 28% slower (902 s vs 704 s, both warm) — and that slowdown
  turned out to be §4.0, not the chunk size: 3x the chunks meant 3x the cache
  revalidation round-trips. With §4.0 fixed, re-measured full scale on the same
  warm caches: **233.1 s against 262.2 s, and 3.5 GB against 40 GB**. See
  `RECHUNK.md`, including the part where the recommendation optimised bytes and
  asserted a speedup it had not measured — right in the end, but not for the
  reason given at the time.

---

## 7. The pattern

Thirteen of these were **silent**. Not one announced itself; every one produced a
plausible number, an ignored declaration, or a quiet performance cliff. The three
that changed this analysis's answer — the corrupt `layers`, the slab aliasing,
the unapplied `unit_conversion` — all had the same shape: **something was
declared, and something else did not honour the declaration.**

§3.7 adds a second shape, and it took a two-file document to expose it: **the
same pass, given the document at a different stage of resolution in each
binding.** As long as a model was one file with no import edge, the three
bindings' differing resolution depths were indistinguishable.

What actually found them:

* **Cross-binding disagreement.** §2.1 was found because Julia and Rust
  disagreed. Nothing else would have caught it; both bindings' internal checks
  passed.
* **Consuming a declared input for the first time.** §2.2 and §1 had been wrong
  for years and were invisible while nothing read those values.
* **An independent reimplementation.** The NumPy oracle, written from InMAP's Go
  source rather than from this document, is what let a layer assignment be
  checked without the 330 GB matrix — and what proved the physics was right while
  the totals were wrong.
* **Testing that a declared behaviour changes an answer**, not merely that a
  field parses. Every one of these fields parsed, validated and round-tripped
  perfectly.

The last point is the generalisable one. A conformance fixture for a
declared-behaviour field should assert the value actually *moved*, and include a
control where it should not — otherwise an implementation that applies nothing
passes.
