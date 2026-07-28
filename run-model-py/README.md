# `run-model-py` — ISRM through the Python binding's evaluation of the `.esm`

Computes the InMAP ISRM tutorial from `isrm_pushdown.esm` using **EarthSciAST's
Python binding** and **EarthSciIO's Python readers**, and emits a
[`contract/results_schema.json`](../contract/results_schema.json) record with
`mode: "runtime_observed_graph"`.

## The bar this meets

Every number in the record comes from the binding's own evaluation of the `.esm`
graph. There is no hand-written STEP-0 arithmetic in this directory. In
particular:

| quantity | where it comes from |
|---|---|
| `ppl` — the emission-bearing source cells | `value_invention.materialize_value_invention` running the model's own producer aggregate, gated by the spatial `join.overlap` broad phase |
| `E_VOC` … `E_PM25` — binned per-source emissions | `numpy_interpreter.eval_expr` on the `.esm` observed expressions |
| `conc_*`, `TotalPM25`, `deathsK`, `deathsL` | same |

`run_model.py` contributes **orchestration only**: which observed to evaluate
when (from each aggregate's own declared `args`), and where the bytes come from.

### Why `isrm_pushdown.esm` and not `isrm.esm`

The three model variants differ in where the source-cell support set comes from:

| model | `emis_src_cells` | `ppl` derived by |
|---|---|---|
| `isrm_clean.esm` | `src_cells` (52,411) | n/a — SR would be 52411² ≈ 330 GB |
| `isrm.esm` | `interval` of size `N_PPL` | **the runner**, supplied as `src_cell_of_ppl` |
| `isrm_pushdown.esm` | `derived` via `join.overlap` | **the graph** |

Only the pushdown variant keeps the spatial join inside the spec. Driving
`isrm.esm` would mean this runner computed the support set itself, which is
exactly the hand-written math the bar excludes — and the support set is one of
the numbers under test (`ppl` is compared **exactly** across bindings, per
CONFORMANCE_SPEC §5.5).

Running it required implementing the overlap gate in the Python binding
(CONFORMANCE_SPEC §5.5.6, which had named Python as not yet implementing it).

## What is *not* from the graph

`inputs.py` prepares exactly the arrays the `.esm` declares as `parameter`s, and
nothing beyond them:

* the EGU FF10 records, read with EarthSciIO's `FF10Reader`;
* their LCC projection into the grid plane — the `.esm` declares `X`/`Y` as
  parameters, so the projection is upstream of the spec. The Julia runner does
  the identical thing with the identical constants
  (`run-model-jl-pushdown/l3_common.jl`); the two MUST agree here or the derived
  support set diverges for a reason that is not the model;
* the grid `W/S/E/N` and `TotalPop`/`MortalityRate`, read from the ISRM zarr.

## The gated fetch

Value invention derives the members *first*; the SR fetch is then built from
them — `SR[layer 0, ppl, :]`, 1,520 rows of 52,411, never the whole matrix. The
selection is pushed down to the zarr reader, so only the intersecting chunks are
fetched, and pathways are fetched and evicted one at a time so peak disk is one
pathway (~2.7 GiB) rather than five.

## Running it

Needs Python ≥ 3.11 (zarr 3.x requires it; the repo's 3.9 `.venv` cannot read
zarr at all). Set up an environment with the two packages installed editable:

```bash
python3.12 -m venv /path/to/venv       # NOT under /tmp on a cluster with tmpfs
/path/to/venv/bin/pip install numpy 'zarr>=3' fsspec numcodecs requests
/path/to/venv/bin/pip install -e ../../EarthSciAST/pkg/earthsci-ast-py \
                              -e '../../EarthSciIO[zarr]'
```

```bash
ISRM_FIRSTN=200 python run_model.py    # reduced, ~1 min  -> results_reduced.json
python run_model.py                    # full             -> results.json
```

Environment overrides: `ISRM_MODEL`, `ISRM_RUNMODEL`, `EGU_ZIP`, `ISRM_ZARR_URL`,
`ISRM_N_SRC`, `ISRM_SCRATCH`, `ISRM_SR_DIR`, `ISRM_FIRSTN`. **`ISRM_SCRATCH` must
be disk-backed** — on this cluster `/tmp` is a tmpfs, so SR blobs written there
consume the same cgroup memory the model needs.

## Cross-checking against Julia

```bash
python3 ../contract/compare_results.py \
    ../run-model-jl-pushdown/results.json results.json
```

At reduced scale (200 records, 9 source cells) the two bindings agree:

```
ppl: count=9 sha256 MATCH
pathway PrimaryPM25  emis rel 0.00e+00   conc rel 1.72e-16
pathway SOA          emis rel 0.00e+00   conc rel 1.66e-16
pathway pNH4         emis rel 0.00e+00   conc rel 2.04e-16
pathway pNO3         emis rel 0.00e+00   conc rel 0.00e+00
pathway pSO4         emis rel 0.00e+00   conc rel 0.00e+00
total_pm25       max rel diff 0.000e+00   (bit-identical)
deaths.krewski   max rel diff 2.126e-14
deaths.lepeule   max rel diff 1.701e-14
99 checks, 0 failed
```

The `ppl` set is **byte-identical** (it is integer-valued, so §5.5 requires
that) and the binned emission sums are **bit-identical**. The ~1e-14 in the
deaths is float reduction order over a 1,520-term contraction, which is what the
1e-12 field tolerance exists for — asserting tighter would be asserting on
summation order rather than on the model.
