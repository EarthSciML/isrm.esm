# `contract/` — the shared record

Everything here exists so that three independent runtimes can be shown to have
evaluated the *same* document and got the *same* answer.

| file | what it is |
| --- | --- |
| `results_schema.json` | the schema every runner's `results_<model>.json` conforms to |
| `results.py`, `results.jl` | the emitters — Python and Julia mirrors of the same sampling and hashing rules (`run-rs/src/contract.rs` is the third) |
| `compare_results.py` | the comparator, and the definition of record where the emitters disagree |
| `plume_oracle.py` | the offline plume-rise oracle (below) |
| `records/plume_oracle.json` | its full-scale output — the independent target for the `plume` block |
| `records/plume_oracle_first{200,2000}.json` | the same at the two documented truncations, frozen so a reduced run can be checked without network access |
| `records/ground-level-only/` | the four pre-plume-rise full-scale baselines, retired; see the README there |

Records are named after the DOCUMENT (`results_isrm_point.json`,
`results_isrm_polygon.json`, `results_isrm_line.json`), because the same three
shims drive all three models and they must not share one record file. The
schema is the same for all of them: a document that states no plume rise —
neither an area source nor a line source has a stack — simply carries no
`plume` block, and one that claims no published total is reported rather than
graded. Which observeds a runner reads comes from the document's own
`metadata.x_esd.report`; a record from a sibling document is comparable across
bindings, and not comparable to a point record, which the comparator's `grid`
check already enforces.

`compare_results.py` skips any record carrying a `kind` field, so
`contract/records/*.json` stays a usable glob as non-results records land
there — with one exception, `kind: "plume_oracle"`, which is loaded as a
first-class target for the `plume` block. Skipping it would have thrown away
the only check in the set that says the physics is *right* rather than merely
*agreed on*.

## The `plume` block

The schema's `plume` block is what makes plume rise a reported, compared
result rather than an unexplained shift in `deathsK`. Each runner fills it
from `observed_field` on the document's own `sr_lower`, `stack_layer` and
`w_sr0`/`w_sr1`/`w_sr2` observeds and its `E_<pathway>_L<layer>` aggregates —
no runner recomputes plume rise, and none of them contains the word ASME. That
is the claim being made: the *engine* derived the assignment from the spec.

A record is not charged to *a* layer. InMAP's `sr.Reader.layerFracs` splits it
across **two** whenever its plume's model layer falls strictly between two
entries of `sr.layers` = `[0, 3, 6]` — model layers 1–2 between SR 0 and 1,
model layers 4–5 between SR 1 and 2. So the block reports a split, not an
assignment:

* `sr_lower` — digest and histogram of the *lower* of the (at most two) SR
  layers each record's mass goes to. Integer-valued, and the upper one, when
  there is one, is always `sr_lower + 1`.
* `weights` — `w_sr0`/`w_sr1`/`w_sr2` as FieldSummaries, plus `max_sum_error`,
  the largest `|w0 + w1 + w2 − 1|` over records. layerFracs conserves mass
  exactly, so that number is 0 or the document is broken; it is checked
  against an absolute bound in each record independently rather than compared
  between records.
* `stack_layer` — the model layer the stack top sits in, which is what selects
  the meteorology the rise is computed from. It is there to localize a
  disagreement: a wrong stack layer means the met gather is wrong, a right
  stack layer with a wrong split means the ASME expression or the layerFracs
  interpolation is.
* `pathways.<SR name>.by_sr_layer` — emitted mass in each layer, per pathway,
  the weights integrated against mass. The physics, in tons.

`sr_lower` and `stack_layer` are **integer-valued**, so `compare_results.py`
compares them EXACTLY, the way it compares `ppl` — across every pair of live
records, and against `records/plume_oracle.json`. A float tolerance there
would hide a real disagreement about which layers a record emits into. The
weights and the per-layer masses are genuinely floats and get the ordinary
`RTOL_FIELD = 1e-12`; a weight field's `sha256` is reported but not asserted,
because the fractions descend from `plume_height`, whose cube roots differ by
an ulp between languages. (Measured: bit-identical across Julia, Rust and the
numpy oracle on the first 200 records, and within 1.2e-15 on the first 2000,
where all four ASME branches are in play.)

The oracle comparison is matched on `n_rec`, so a reduced oracle
(`plume_oracle_first200.json`) is never held against a full-scale run: a
truncation is a different problem, and the digest of a 200-record assignment
says nothing about a 43,650-record one.

## Hashing

Two conventions, one wire format. Integer sequences hash as ASCII decimals
joined by `,` with no spaces (`int_seq_sha256`); `ppl` is a member *set*, so
`ppl_sha256` sorts first, while a per-record field is a *sequence* and keeps
its record order. Float fields hash as little-endian IEEE-754 float64 bytes
(`field_sha256`). `sample_indices` is pure integer arithmetic in all three
languages so the sampled index set cannot drift on float rounding.

## `plume_oracle.py` — the cheap check on plume rise

`isrm_point.esm` reproduces the InMAP source-receptor tutorial
(<https://inmap.run/blog/2019/04/20/sr/>), **plume rise included**: the
document states the ASME rise itself and charges each emission record to the
SR emission layer its plume reaches. Emitting everything at ground level
instead gives `sum(deathsK) = 7524.918845602511`; the blog gets `6928.959583`.
Plume rise changes exactly one intermediate quantity — the
`(source cell, SR-layer weights)` assignment each record gets — and every
downstream number follows from it.

That assignment is checkable **without the 330 GB SR matrix**. It needs thirteen
1-D meteorology/geometry arrays off the ISRM zarr, ~16 MB compressed, each a
single blosc chunk fetched over plain HTTPS. So this is the fast test, and the
one to run first when the document's `deathsK` moves:

* digest matches → the new physics is right, look at the contraction;
* digest differs → stop before fetching a single SR chunk.

```sh
python3 contract/plume_oracle.py            # ~1 min, writes records/plume_oracle.json
python3 contract/plume_oracle.py --firstn 200   # the REDUCED target
python3 contract/plume_oracle.py --help     # --out / --zip / --cache / --no-assert
```

`--firstn N` truncates to the first N **delivered** records — the same
truncation `ISRM_FIRSTN=N` applies to the shims, taken after the loader's
`record_filter`, so it selects the same records — and writes
`records/plume_oracle_first<N>.json`. That is what a reduced run is checked
against while the full-scale one is too slow to iterate on. The expected
values above are full-scale facts and are not asserted on a truncation; the
structural checks on the grid and the algorithm still are.

System `python3` (3.9) plus `numpy` and `numcodecs`; no venv. Zarr chunks are
cached under `$ISRM_SCRATCH/plume-oracle-zarr`, resolved exactly the way
`run-py/paths.py` resolves scratch (default `/scratch.local/$USER/isrm-esm`) —
never `/tmp`, which is a tmpfs here and eats the memory cgroup. Emissions come
from `$EGU_ZIP`, default `data/2016fd_inputs_point.zip`.

The script re-implements nothing of its own: the algorithm is
`ctessum/atmos/plumerise` (`ASMEPrecomputed`, `calcDeltaHPrecomputed`,
`findLayer`), `inmap/plumerise.go` (`IsPlumeIn`) and `inmap/sr/srreader.go`
(`Concentrations`, `layerFracs`). It asserts its own expected outputs and every
structural fact it relies on — that the ISRM grid is 596444 = 8·52411 +
19·9324 cells laid out layer-major, that model layers 0–7 share a byte-identical
horizontal grid — and refuses to write the record if any of them drifts.

It used to assert that the SR matrix's `layers` array is `[0, 1, 2]`. That
check was **asserting the corruption**: the authoritative `isrm_v1.2.1.ncf` on
Zenodo holds `[0, 3, 6]`, and the zarr's `[0, 1, 2]` is a machine-generated
arange that displaced it during the conversion. It is inverted now. The oracle
interpolates on a *declared* `SR_MODEL_LAYERS = [0, 3, 6]`, mirroring the
document's `SR_MODEL_L1` / `SR_MODEL_L2` metaparameters, and checks instead
that the layers it uses are never the corrupt arange, and that the store's own
value is one of the two states there is evidence for — anything else is a hard
failure. Which state the store served is printed and recorded in the record as
`grid.sr_layers_in_store`, beside `grid.sr_layers_used`.

Two of its checks are free cross-validation against the document itself: the
five per-pathway emission totals and the `ppl` digest must equal those in
`records/ground-level-only/python.json`. Those two survived the plume-rise
re-baselining that retired that record — plume rise moves mass between layers,
never into or out of a pathway, and changes which layer a record emits into,
never which cell — so a live record must still carry them, and the comparator
checks that. Agreeing there proves the oracle ingests the FF10 zip
the way `isrm_point.esm` does, so a disagreement in the SR-layer assignment is about
plume rise and nothing else.

### The record

`records/plume_oracle.json`. `kind: "plume_oracle"`, so the comparator does not
read it as a results record — but it does not skip it either: it is the target
the `plume` block is checked against.

The headline fields are `sr_lower.sha256` — sha256 over the per-record *lower*
SR layer as ASCII decimals joined by `,` in record order, the same wire format
as `ppl_sha256` — and `weights`, the three per-record shares. Alongside them:
`stack_layer` and `plume_model_layer` histograms and digests, `plume_height_m`
(a `field_summary` plus `mean`), `branch_usage` (which of ASME's four branches
each record took), and the per-pathway emission mass split across SR layers
0/1/2 in short tons/yr.

`stack_layer` carries two histograms because they answer two questions.
`histogram_height_gt_0` is the physics — how many *stacks* top out in each
model layer, with the 1,100 zero-height records that never enter plume rise
excluded. `histogram` counts all records, which is what a runner can emit from
the document's `stack_layer` observed alone (that observed is 0 for a
zero-height record, so bin 0 carries them), and is therefore the one the
contract compares.

Current values:

```
records: 43650, zero-height records: 1100, max stack height: 316.3824 m
distinct emitting cells: 1520
stack-layer histogram (records with height>0): [29313, 8460, 3899, 878]
plume model-layer histogram: [12108, 9155, 11017, 6392, 2156, 1177, 553, 438, 654]
   (index 8 = "at or above model layer 8")
branch usage: momentum 2073, stable-buoyant 868, unstable-buoyant 29601, F<=0 11108
max plume height: 9437.1 m

SR-layer emission split (short tons/yr, % of pathway total):
  VOC    total      33452.8   L0 26.1%  L1 54.4%  L2 19.5%
  NOx    total    1314462.9   L0 22.4%  L1 53.5%  L2 24.1%
  NH3    total      25012.5   L0 36.9%  L1 54.7%  L2  8.4%
  SOx    total    1571216.9   L0 12.7%  L1 51.2%  L2 36.1%
  PM25   total     140822.7   L0 21.9%  L1 51.9%  L2 26.2%

lower-SR-layer histogram: [32280, 9725, 1645]
sr_lower sha256: d38ba2fb042f7e793134670e954d306987a8d9b17fca20975c23d36a9a134799
weight sums (w_sr0/w_sr1/w_sr2): 20855.831204066 / 18999.354955120 / 3794.813840815   max|Σw - 1| = 0
```

Everything above the SR split is **unchanged to the last digit** from the era
when the document clamped each record to one layer — the stack-layer
histogram, the plume model-layer histogram, the four branch counts, the five
emission totals, the above-layer-7 group and the `ppl` digest. Only the split
moved, and it moved a long way: SOx read `0.5 / 4.3 / 95.2` under the clamp
and reads `12.7 / 51.2 / 36.1` under `layerFracs`. That is what the corrupt
`layers` was doing — with `[0, 1, 2]` almost every plume counts as "above the
top" and gets shovelled into SR layer 2.

### One deliberate deviation from InMAP

654 records (13269.9 short tons/yr, **0.43% of emitted mass**) have plumes
above the top of model layer 7 — out of the 52411-cell ground grid and into
the 9324-cell high-altitude grid.

InMAP has a latent defect there. `sr.Reader.layerFracs` clamps such an emission
into SR layer 2, which is right (every model layer above 6 clamps). But the
horizontal index it pairs with that layer is `sr.indices[c]`, and `srreader.go`
builds that map by resetting the counter to 0 at every layer boundary. For a
cell in layer ≥ 8 the index is therefore a position in the **coarse** 9324-cell
grid, which `Reader.source` then reinterprets as a position in the 52411-cell
ground grid — so those emissions are charged to the wrong source cell.

**This oracle implements the correct behaviour**: clamp to SR layer 2 at the
source cell the emission actually came from. It is not bug-compatible with
InMAP. The `above_model_layer_7` block in the record counts and quantifies the
group, per pathway, so that any residual gap between the document's number and
the blog's published `6928.959583` is attributable rather than mysterious.

How much it is worth is **not currently measured**. The `+0.79% / +0.82%`
figures this file used to quote came from a full-scale run that predates both
the sibling-slab aliasing fix and `layerFracs`, and no full-scale run of the
current document has been made. What *is* measured is that at
`ISRM_FIRSTN=200` — where the group is empty — the document reproduces the
live `inmap cloud` service to 8.9e-9, and at `ISRM_FIRSTN=2000`, where 286 of
2000 records are in the group, no service target has been collected to compare
against.
