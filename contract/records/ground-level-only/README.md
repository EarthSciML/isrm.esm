# `ground-level-only/` — the four pre-plume-rise baselines

These four full-scale records are the state of this repo *before* `isrm.esm`
stated plume rise. Every emission in them is charged to SR layer 0, because
that is what the document said at the time:
`data_loaders.ISRM_SR.metadata.x_esd.gated_select.axes[0]` was the literal
`{"fixed": [0]}`, and the four stack parameters were read, unit-converted and
coupled into the model but consumed by nothing.

| file | binding | mode | model |
| --- | --- | --- | --- |
| `julia.json` | Julia | `runtime_observed_graph` | `isrm_pushdown.esm` |
| `python.json` | Python | `runtime_observed_graph` | `isrm_pushdown.esm` |
| `rust.json` | Rust | `runtime_observed_graph` | `isrm_pushdown.esm` |
| `julia-oracle.json` | Julia | `oracle_step0` | `isrm.esm` |

All four report `sum(deathsK) = 7524.9188456…` and
`sum(deathsL) = 16979.6321714…`, agreeing to 3.1e-15.

## Why they are no longer compared against live runs

Plume rise moves most of the emitted mass into SR layers 1 and 2 — at full
scale SOx is 0.5 / 4.3 / 95.2 % across layers 0/1/2 — and mass in a higher
layer disperses differently, so every concentration and every death total
changes. A current run lands near the blog's `6928.959583`, about 8 % below
these. Comparing the two is not a tolerance question; they are answers to
different models, and `compare_results.py` will say so loudly (its float
tolerance is 1e-12, and one record carrying a `plume` block while the other
does not is itself a failure).

The live records are checked against `../plume_oracle.json` instead, which is
the reference for the physics that replaced them.

## Why they are kept

* **They are the evidence for the pushdown-era claim.** Three bindings and a
  hand-written STEP-0 oracle agreeing to 3.1e-15 on a 330 GB contraction is
  the result the previous phase of this repo was built to produce, and it is
  still true of the document those records were emitted from (tag
  `pushdown-era`).
* **Two of their numbers survive the re-baselining, and are still asserted
  on.** `ppl` — the 1,520-cell emission-bearing support set,
  `sha256 = 6f784d7e…` — and the five per-pathway `emis_sum` totals are
  invariant under plume rise, which changes which *layer* a record emits into
  and how much mass sits in each, never which *cell* it emits from nor how
  much a pathway emits in total. `contract/plume_oracle.py` asserts against
  exactly those values (see its `EXPECT` block), which is what proves its FF10
  ingest and its cell containment match the document's rather than merely
  matching themselves.
* **They are the before half of a measured change.** The gap between
  7524.92 and the current number is what plume rise is worth, in deaths.

## Re-running them

The document that produced them is at tag `pushdown-era`; `isrm.esm` on `main`
no longer has a ground-level-only mode, and should not grow one — plume rise
is physics, not an option.
