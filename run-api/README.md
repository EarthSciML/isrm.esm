# `run-api` — the ISRM through EarthSciLab, not through a local engine

`run_isrm_demo.py` runs the **`isrm_demo`** analysis declared inside
`isrm_point.esm` on [EarthSciLab](https://earthscilab.com)'s API, and draws the
plot that analysis declares.

**This is not a fourth binding.** `run-rs`, `run-jl` and `run-py` drive
`earthsci-ast` in-process on this machine, and they need the engine, EarthSciIO,
a 69 MB EPA zip and a multi-GB chunk cache to do it. This drives the *same* core
through the product: EarthSciLab's hardware fetches the 14 GB of SR slabs, and
what arrives here is 52,411 numbers and a bill. The point of the totals it
prints is that the two must agree.

```
python3 run_isrm_demo.py --quote-only     # price it; no account needed
python3 run_isrm_demo.py --records 200    # reduced: minutes
python3 run_isrm_demo.py                  # FULL SCALE: ~50 min, ~$0.04
```

Only `matplotlib` is needed beyond the standard library, and only to draw —
`--quote-only` and signing in work from a bare interpreter.

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

## The notebook

[`isrm_esm.ipynb`](isrm_esm.ipynb) is not this script in notebook form. It is one
function:

```python
receptors = run_isrm(gdf)
```

`gdf` is a GeoDataFrame of emissions — points, polygons or lines, in any CRS,
with a column per pollutant in kg/yr (`PM25`/`NOx`/`NH3`/`SOx`/`VOC`) and, for
points, optional `STKHGT`/`STKDIAM`/`STKTEMP`/`STKVEL` in m/m/K/(m/s). Back comes
a GeoDataFrame of the 52,411 receptor cells with `TotalPM25`, `deathsK` and
`deathsL`, its geometry the cells' own rectangles in the grid's Lambert conformal
— so the variable-resolution mesh is drawn as it is, with `receptors.plot()`.
Sections 3, 4 and 5 run it on the EGU point inventory, on Illinois counties and
on Illinois interstates.

**It does not read these documents.** `isrm_point.esm` and its siblings expect an
inventory in the inventory's own units — FF10's long format, column indices and
77-entry pollutant-code map for the point document; a single `EMIS` column for
the other two. The notebook reads a generated sibling set instead —
`isrm_gdf_point.esm`, `isrm_gdf_point_flat.esm`, `isrm_gdf_polygon.esm`,
`isrm_gdf_line.esm` — written by
[`data/make_gdf_documents.py`](../data/make_gdf_documents.py) from the published
four: same physics, one uploaded shapefile with one column per pollutant in SI,
and all five pathways wired for every geometry rather than PM2.5 alone for the
area and line ones. Regenerate them with

```
python data/make_gdf_documents.py           # write
python data/make_gdf_documents.py --check   # fail if out of date
```

and note that the notebook fetches them from `raw.githubusercontent.com`, so a
regenerated document has to be **pushed** before a run sees it.

Three things `run_isrm` does that a document cannot do for itself, each because
the reader or the engine requires it:

- **One row per shapefile PART, emissions split between the parts.** The reader
  replicates a record's attributes onto each of its parts, so a county of
  mainland-plus-islands would otherwise carry its whole emission two or three
  times. `isrm_polygon.esm` says whose job that is: "the layer builder's job and
  not this document's."
- **Interior rings dropped.** A hole is written as another part and read as its
  own record, so it would add its emission share rather than subtract its area.
- **Polylines cut into two-vertex segments**, each road's emission apportioned by
  length — because the projection-pushdown rewrite needs the binning aggregate to
  declare exactly two ranges, so the segment has to be the record.

It also prunes the pathways your frame has no column for. The gated fetch is most
of what a run costs and it is per pathway, so a PM2.5-only frame fetches one
source-receptor slab instead of five.

Section 3 downloads the 72 MB FF10 zip to build its point frame; sections 4 and 5
fetch a Census boundary file and a TIGER road file. Nothing is written to disk
except one temporary directory, because GDAL writes a shapefile as four files
with a shared stem.

## A local file cannot be dispatched

`isrm_polygon.esm` reads `file://data/polygon_emissions_17.zip`, and that fails on
the server twice over — `file://` is deliberately not on the run allowlist (a
document is a fetch instruction to a dispatched run, so a loader pointing at
`file:///proc/self/environ` would be a worker reading its own secrets into a
dataset the caller can download), and the bytes are on your machine anyway.

The notebook's fix is to upload the layer as an EarthSciLab **dataset** and
re-point the loader at it. That works because the run allowlist admits our own
bucket by deriving the entry from the run's *own output target*, which is what
lets a run read back a dataset it did not write. It is three calls — `POST
/datasets`, `PUT /datasets/{id}/objects?key=…`, `POST /datasets/{id}/commit` —
and then one edited string in the document.

## Signing in

The first run that starts a job prints a code and opens a browser:

```
  Sign in to EarthSciLab. Your code is:  RRGQ-BJVS
  Open: https://beloved-melody-28.authkit.app/device?user_code=RRGQ-BJVS
```

That is the **OAuth 2.0 device authorization grant** — WorkOS's "CLI Auth",
which AuthKit serves with no extra configuration. It is the right flow here for
a reason that is specific rather than stylistic: a full-scale run takes about an
hour and a WorkOS *access* token is short-lived, so a hand-pasted token expires
long before the answer exists. The device grant hands back a **refresh** token,
which this script stores at `~/.earthscilab/credentials.json` (mode 0600),
rotates on every use, and spends to mint a fresh access token before each
request. The token it produces is an ordinary AuthKit user JWT — same JWKS, same
`sub` — so the API needed no change to accept it, and the run bills to your own
account.

`--login` signs in and stops; `--logout` forgets the credential.

## What it sends

| | |
|---|---|
| `POST /quote` | prices the document. No auth, no database — the pre-login preview. |
| `POST /runs` | `{"esm": …, "kind": "evaluate", "observeds": [...], "max_price": …}` |
| `GET /runs/{id}/events` | SSE progress, reconnecting if the socket drops |
| `GET /datasets/{id}/field` | one 1-D array per observed, out of the run's Zarr store |

Four things about that request are worth knowing before you change it. All four
are in the script's own header comment at more length; the short version:

1. **`kind: "evaluate"`.** `models.ISRM` has no `D(·)` anywhere, `isrm_demo`'s
   time span is `0 -> 0`, and the whole answer is the observed graph. Dispatched
   as a `simulate` run the engine refuses the document outright with
   `Invalid parameter 'src_E'`.

2. **The template library is inlined first.** `isrm_point.esm` reaches its shared
   body by `{"ref": "./isrm_base.esm"}`, and a document that arrives over the
   wire has no directory to be relative to — the raw document fails on the
   runner with `template-library file not found`. The script merges the
   library's `expression_templates` into the model and drops the import key.
   Deliberately *not* a full esm-spec §9.7 resolve, which would also close the
   metaparameters and bind `N_REC` to its declared default of 0 before the FF10
   loader has counted the records.

3. **`--records N` is a document edit**, not a request parameter — a loader-level
   `select` range on every source that discovers its own extent. A run's scale is
   a property of the request, and our requests carry documents.

4. **The observeds are read off the analysis**, not hardcoded: the plots' own
   `x`/`y` variables, plus `TotalPM25`, `deathsK`, `deathsL` for the totals. Ask
   for a name the build cannot produce and the run fails, which is the point —
   the "give me everything" mode silently skips an observed it cannot evaluate.

## Checking the answer

The totals are compared against whichever `run-*/results*.json` in this repo was
recorded at the same record count, so a run is checkable the moment it lands:

| records | `sum(deathsK)` | from |
|---|---|---|
| 200 | 49.11639491165982 | `run-rs/results_reduced.json` |
| 2,000 | 363.47671096747285 | `run-rs/results_isrm_point_reduced.json` |
| 43,650 (full) | 7022.724781368745 | `run-rs/results.json` |

Expect agreement to ~1e-15, not bit-identity: the shims sum with Kahan
compensation and this folds naively, which the repo README already records as a
~2.9e-13 difference on *bit-identical* fields. A disagreement bigger than that is
a real one.

The full-scale record labels itself `model: "isrm.esm"`; EarthSciLab's
`docs/isrm.md` records that its numbers are in fact the point document's against
the 1.2.2 store. The script prints that mismatch rather than hiding it.

## Two things that will look wrong and are not

**The quote is the same for 200 records and for 43,650.** The estimator floors a
static evaluation at a fixed duration and cannot yet tell the two apart, so both
price at $0.04. The difference between them is wall-clock and your own time, not
money.

**Progress crawls near the end.** The fraction is *phase*-weighted across eight
`PreparePhase`s that differ by four orders of magnitude in cost: `Rewrite` is
milliseconds and `GatedFetch` is tens of gigabytes off S3. 77% is the fetch 16%
in. A long crawl there is the I/O, not a hang.

## If the connection drops

Nothing is lost — the run is server-side. The event stream replays everything
already recorded before it streams, so the script reconnects and picks up what it
missed, falling back to polling `GET /runs/{id}` if the stream will not come
back. To walk away entirely, note the run id it prints and re-attach later:

```
python3 run_isrm_demo.py --run-id 490ace35-…
```

Pass `--records` again when re-attaching to a reduced run: the script has no way
to ask a finished run how many records it ingested, so it keys the reference
totals off what you tell it, and a reduced run compared against the full-scale
record will look wrong by two orders of magnitude.
