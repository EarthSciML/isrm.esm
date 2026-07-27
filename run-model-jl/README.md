# isrm.esm — Julia runner

Reproduces the [InMAP ISRM tutorial](https://inmap.run/blog/2022/12/15/tutorial/)
(minus the "Bonus! Visualization" section) from the language-agnostic
[`isrm.esm`](../isrm.esm) model, driven by the EarthSciAST Julia binding and the
EarthSciIO `zarr` + `ff10` readers.

## What it does

**`run-model.jl` — Part A: full-scale numerical reproduction.**
Fetches the ISRM source–receptor matrices from `s3://inmap-model/isrm_v1.2.1.zarr`
and the 2016fd EGU point emissions (FF10), performs the LCC spatial join and the
SR contraction, and reports attributable-death totals:

| quantity | this runner | tutorial target | rel. err |
|---|---|---|---|
| `sum(deathsK)` (Krewski, RR 1.06) | **7524.92** | 7524.84 | 0.001% |
| `sum(deathsL)` (Lepeule, RR 1.14) | **16979.63** | 16979.45 | 0.001% |
| `Σ TotalPM25` | 20619.02 | — | — |

It is **resumable and disk-safe**: each pathway's concentration vector is
checkpointed under `checkpoints/`, and its ~3 GB of SR chunk blobs are evicted
before the next pathway, so an interrupted run resumes without re-doing completed
pathways and peak disk stays ≈ one pathway (not all five).

**`part_b.jl` — Part B: end-to-end through the runtime.**
Loads the authored `isrm.esm` with `EA.load`, builds it with `EA.prepare`, and
evaluates its observeds with the framework's own state-free-observed evaluator
(`EA._observed_field`) on a small controlled input, checking every stage of the
graph (LCC projection → rectangle-containment emission binning → SR contraction →
`fact·Σconc` → exp/log health) against a plain-Julia oracle. All stages match to
machine precision — i.e. the runtime faithfully executes the authored `.esm`.

**`plot.jl` — the `.esm` example plot** (`isrm_demo` / `totalpm25_map`): a scatter
of `TotalPM25` and attributable deaths vs receptor cell-centroid easting across
the 52,411-cell grid (the unstructured mesh has no structured raster / polygon
choropleth in-format, so scatter is the authored alternative → `isrm_totalpm25.png`).

## Running

```bash
# One-time: dev-track the sibling EarthSciAST + EarthSciIO checkouts, then
# instantiate + precompile. (`Pkg.instantiate()` alone is NOT enough — the
# registered releases lack the zarr + ff10 readers this runner is built on.)
julia --project=. setup.jl

# the 72 MB EGU FF10 zip, expected at data/2016fd_inputs_point.zip:
mkdir -p data && curl -L -o data/2016fd_inputs_point.zip \
  https://gaftp.epa.gov/air/emismod/2016/alpha/2016fd/emissions/2016fd_inputs_point.zip

# Part A (fetches ~14 GB of SR chunks over the run; ~25 min; resumable)
ISRM_SR_ROOT=/scratch/$USER/isrm_sr julia --project=. run-model.jl

julia --project=. part_b.jl     # Part B (seconds; no network)
julia --project=. plot.jl       # example plot (needs Part A's result)
```

`setup.jl` finds the checkouts at `../../EarthSciAST/pkg/EarthSciAST.jl` and
`../../EarthSciIO/julia` by default; override with `EA_PATH` / `IO_PATH`. The
resulting `Manifest.toml` is gitignored because it pins machine-specific dev
paths — `setup.jl` is how a fresh checkout becomes runnable.

`ISRM_SR_ROOT` redirects the bulky, pathway-by-pathway-evicted SR chunk cache;
point it at fast **local** disk when the checkout lives on a network filesystem.
Network-robustness knobs (all optional): `EARTHSCIIO_HTTP_TIMEOUT`,
`EARTHSCIIO_HTTP_RETRIES`, `EARTHSCIIO_LOCK_STALE_AGE`.

## Output: the cross-language contract record

Part A also writes `results.json` in the shared shape defined by
[`../contract/results_schema.json`](../contract/results_schema.json), so this
runner's numbers can be diffed against the Python and Rust runners':

```bash
python3 ../contract/compare_results.py run-model-jl/results.json [...]
```

It is tagged `"mode": "oracle_step0"` — these numbers come from `run-model.jl`'s
own hand-written STEP-0 arithmetic over the loader outputs, **not** from
evaluating the `.esm` observed graph. It is the reference oracle the
runtime-driven runners are checked against, not itself evidence that a binding
executes the specification.
