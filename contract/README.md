# `contract/` — the shared record

Everything here exists so that three independent runtimes can be shown to have
evaluated the *same* document and got the *same* answer.

| file | what it is |
| --- | --- |
| `results_schema.json` | the schema every runner's `results.json` conforms to |
| `results.py`, `results.jl` | the emitters — Python and Julia mirrors of the same sampling and hashing rules (`run-rs/src/contract.rs` is the third) |
| `compare_results.py` | the comparator, and the definition of record where the emitters disagree |
| `plume_oracle.py` | the offline plume-rise oracle (below) |
| `records/plume_oracle.json` | its full-scale output — the independent target for the `plume` block |
| `records/plume_oracle_first{200,2000}.json` | the same at the two documented truncations, frozen so a reduced run can be checked without network access |
| `records/ground-level-only/` | the four pre-plume-rise full-scale baselines, retired; see the README there |

`compare_results.py` skips any record carrying a `kind` field, so
`contract/records/*.json` stays a usable glob as non-results records land
there — with one exception, `kind: "plume_oracle"`, which is loaded as a
first-class target for the `plume` block. Skipping it would have thrown away
the only check in the set that says the physics is *right* rather than merely
*agreed on*.

## The `plume` block

The schema's `plume` block is what makes plume rise a reported, compared
result rather than an unexplained shift in `deathsK`. Each runner fills it
from `observed_field` on the document's own `plume_layer` and `stack_layer`
observeds and its `E_<pathway>_L<layer>` aggregates — no runner recomputes
plume rise, and none of them contains the word ASME. That is the claim being
made: the *engine* derived the assignment from the spec.

* `sr_layer` — digest and histogram of the per-record SR emission layer. THE
  quantity plume rise exists to produce; everything downstream follows from it.
* `stack_layer` — the same for the model layer the stack top sits in, which is
  what selects the meteorology the rise is computed from. It is there to
  localize a disagreement: a wrong stack layer means the met gather is wrong,
  a right stack layer with a wrong SR layer means the ASME expression is.
* `pathways.<SR name>.by_sr_layer` — emitted mass in each layer, per pathway.
  The physics, in tons.

Both layer assignments are **integer-valued**, so `compare_results.py`
compares them EXACTLY, the way it compares `ppl` — across every pair of live
records, and against `records/plume_oracle.json`. A float tolerance there
would hide a real disagreement about which layer a record emits into. The
per-layer masses are floats and get the ordinary `RTOL_FIELD = 1e-12`.

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

`isrm.esm` reproduces the InMAP source-receptor tutorial
(<https://inmap.run/blog/2019/04/20/sr/>), **plume rise included**: the
document states the ASME rise itself and charges each emission record to the
SR emission layer its plume reaches. Emitting everything at ground level
instead gives `sum(deathsK) = 7524.918845602511`; the blog gets `6928.959583`.
Plume rise changes exactly one intermediate quantity — the
`(source cell, SR layer)` pair each record is charged to — and every
downstream number follows from it.

That pair is checkable **without the 330 GB SR matrix**. It needs thirteen
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
horizontal grid, that the SR matrix's `layers` array is `[0, 1, 2]` — and
refuses to write the record if any of them drifts.

Two of its checks are free cross-validation against the document itself: the
five per-pathway emission totals and the `ppl` digest must equal those in
`records/ground-level-only/python.json`. Those two survived the plume-rise
re-baselining that retired that record — plume rise moves mass between layers,
never into or out of a pathway, and changes which layer a record emits into,
never which cell — so a live record must still carry them, and the comparator
checks that. Agreeing there proves the oracle ingests the FF10 zip
the way `isrm.esm` does, so a disagreement in the SR-layer assignment is about
plume rise and nothing else.

### The record

`records/plume_oracle.json`. `kind: "plume_oracle"`, so the comparator does not
read it as a results record — but it does not skip it either: it is the target
the `plume` block is checked against.

The headline field is `sr_layer.sha256` — sha256 over the per-record SR layer
as ASCII decimals joined by `,` in record order, the same wire format as
`ppl_sha256`. Alongside it: `stack_layer` and `plume_model_layer` histograms
and digests, `plume_height_m` (a `field_summary` plus `mean`), `branch_usage`
(which of ASME's four branches each record took), and the per-pathway emission
mass split across SR layers 0/1/2 in short tons/yr.

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
  VOC    total      33452.8   L0  6.3%  L1 13.1%  L2 80.5%
  NOx    total    1314462.9   L0  3.1%  L1  5.3%  L2 91.5%
  NH3    total      25012.5   L0 12.3%  L1  6.6%  L2 81.2%
  SOx    total    1571216.9   L0  0.5%  L1  4.3%  L2 95.2%
  PM25   total     140822.7   L0  3.5%  L1  5.8%  L2 90.7%

sr_layer sha256: 808e0971a2eda1de1ffc53e242f7ea3fd9bbbda85c3b61702a56f71dd12b434b
```

### One deliberate deviation from InMAP

654 records (13269.9 short tons/yr, **0.43% of emitted mass**) have plumes
above the top of model layer 7 — out of the 52411-cell ground grid and into
the 9324-cell high-altitude grid.

InMAP has a latent defect there. `sr.Reader.layerFracs` clamps such an emission
into SR layer 2, which is right (model layers 3–26 all clamp). But the
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
