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
| plume rise, `layerFracs` over the true `[0,3,6]` | **6936.106343** | **+0.10%** |
| the InMAP tutorial, reproduced through its own service | 6928.959583 | — |

The remaining +0.10% is one deliberate difference, §5.1.

Two claims made during this work were **wrong and are retracted**: that the
~8.6% ground-level gap "is the plume rise" (it is not — most of it was the
corrupt layer axis), and that the tutorial's totals were unreachable by any
layer assignment (a proof that was valid given the corrupt axis and worthless
because of it).

---

## 1. The one that dominated everything: a corrupt array in the published store

**`s3://inmap-model/isrm_v1.2.1.zarr`, variable `layers`.** Not fixed yet — see
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

The same conversion pass also **lost the per-cell `Layer` variable** (a
case-insensitive collision with the 3-element `layer` dimension coordinate) and
**stripped `.zattrs` from the SR arrays**, which is why this repo has always
needed `seed_empty_zattrs` — a workaround that was evidence of a lossy
conversion sitting in plain sight for years.

The matrix *data* is fine: Zenodo's zip is md5-identical across both 2019
depositions and its rows are byte-identical to the zarr's. One index array was
replaced; nothing was rescaled.

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

## 4. Defects that cost only performance

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
**Effect here:** the repo's headline three-binding claim could not be reproduced
for one of the three, and Rust had never executed the plume-rise document — which
is also why §2.1 went undetected until it did.

---

## 5. Defects in InMAP itself

These are upstream and are **not** fixed here. The first is deliberately not
reproduced; the second deliberately is.

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

### 5.2 `layerFracs` interpolates backwards — reproduced deliberately
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

**Effect here: unmeasured, and this document reproduces it**, because matching
InMAP's published numbers requires matching its arithmetic. Four variable
descriptions say so, to stop anyone silently "correcting" it. If InMAP fixes
this, the tutorial's published totals change.

---

## 6. Still open

* **`temporal` is ignored** by `providers_from_document`: a loader declaring an
  hourly cadence is served CONST — reads the first file once, forever, no warning.
* **`source.mirrors` are dropped** — declared failover URLs never reach the loader.
* **`determinism` is read by nothing**, in either repo, despite `esm-spec.md`
  §8.9 being a normative MUST.
* **`load()` enforcement diverges**: Python enforces §4.8.4 unit errors at load;
  Julia, Rust, Go and TypeScript do not, as `esm-libraries-spec` §2.1a requires.
  So an out-of-spec document runs in four bindings and fails in one. Python's
  non-conformance is what caught §2.3. This needs a spec decision about which
  entry point owns §4.8.4, not an engine patch.
* **The zarr repair** itself, which needs bucket write access.

---

## 7. The pattern

Eleven of these were **silent**. Not one announced itself; every one produced a
plausible number, an ignored declaration, or a quiet performance cliff. The three
that changed this analysis's answer — the corrupt `layers`, the slab aliasing,
the unapplied `unit_conversion` — all had the same shape: **something was
declared, and something else did not honour the declaration.**

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
