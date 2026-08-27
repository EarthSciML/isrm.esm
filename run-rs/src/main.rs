// =============================================================================
// run-rs — the RUST binding drives an `isrm_*.esm` end to end through the
// PUBLIC earthsci_ast surface. `ISRM_MODEL` picks which: `isrm_point.esm` (the
// EGU point inventory, with plume rise) or `isrm_polygon.esm` (an area
// inventory allocated by polygon/cell overlap). NOTHING MODEL-SHAPED LIVES
// HERE: this file names no pollutant, no column, no grid extent, no record
// count and no observed — the document's `metadata.x_esd.report` block names
// what a run reports, so the same shim drives either geometry.
//
//   * `esm_problem(doc, tspan, {providers, pushdown_rewrite: true})` — the automatic
//     projection-pushdown rewrite runs inside the engine; the SR provider
//     gates derive from the rewrite's own record
//     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
//     gate and hand-builds NO `Selection`;
//   * EVERY provider comes FROM THE DOCUMENT (`providers_from_document`:
//     format = `metadata.esio_format`, URL = `source.url_template`) — the SR
//     slabs, the grid, the population, AND the EGU FF10 table, whose ingest
//     the loader now declares in full (esm-spec §8.9): `reader_options` (the
//     zip member glob + header row), `codes` (POLID text -> the pathway enum,
//     an unrecognised code dropping the record), `record_filter` (no
//     coordinate / no annual total is not a record) and `extent` (the
//     surviving count binds N_REC). The src-cell rectangles are the
//     `select` range `W[0:N_SRC]` on their own loader variables;
//   * every reported number is the binding's evaluation of the document's
//     observed graph (`observed_field`) — NO hand-written STEP-0 math here,
//     and NO hand LCC projection (raw emis_lon/emis_lat are the parameters;
//     X/Y are in-model observeds the engine projects at build time).
//
//   FULL run  (default)          -> report sum(deathsK/L) against the tutorial
//   REDUCED   (ISRM_FIRSTN=n)    -> first n emission records, totals reported
//
// Emits the cross-language contract record (contract/results_schema.json) with
// model=<the driven .esm>, mode="runtime_observed_graph", binding="rust".
// =============================================================================

mod contract;
mod paths;

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use earthsci_ast::esio_provider::providers_from_document;
use earthsci_ast::{esm_problem, observed_field, PrepareProvider, ProblemOptions};
use serde_json::{Value, json};

/// One reported pathway, as the DOCUMENT names it.
struct PathwaySpec {
    sr_array: String,
    /// The per-source-cell emission observeds feeding this pathway. A LIST:
    /// `isrm_point.esm` splits a pathway across the three SR emission layers a
    /// record's plume falls between, `isrm_polygon.esm` has one because an area
    /// source emits at the ground.
    emissions: Vec<String>,
    concentration: String,
}

/// `metadata.x_esd.report` — what this document says a run should report.
///
/// The one place a pollutant, a pathway or an observed name enters this file,
/// and it enters FROM THE DOCUMENT. A document that declares no `plume` gets no
/// plume block; one that claims no published result gets no oracle check.
struct Report {
    pathways: Vec<PathwaySpec>,
    total_pm25: String,
    deaths_k: String,
    deaths_l: String,
    /// An observed over the RECORD axis, whose length IS N_REC. `extent` binds
    /// that count inside the engine, but no binding exposes the resolved
    /// metaparameter on its prepared model, so the count is read as the length
    /// of a field over that axis — named by the document, because which field
    /// it is depends on the geometry.
    record_field: String,
    /// `(sr_lower, stack_layer, weights)`, absent when the document states no
    /// plume rise — an area source has no stack, so there is nothing to rise.
    plume: Option<(String, String, Vec<String>)>,
    oracle: Option<Oracle>,
}

/// A published national total a document may claim to reproduce.
///
/// `published` is a REFERENCE POINT, not a target and not an assertion: for the
/// InMAP source-receptor tutorial the document declines BOTH of InMAP's
/// plume-rise defects — the high-plume source-index defect (a plume above model
/// layer 7 keeps an index built in the coarse 9324-cell grid and is then read
/// against the 52411-cell ground grid, misplacing 654 of 43650 records onto the
/// wrong source cell) and the inverted layerFracs interpolation (which puts
/// 6.25% of emitted mass on the wrong side of a split) — so the run lands ABOVE
/// it by about 1.35%, deliberately.
///
/// `corrected` is what THAT document computes at full scale with CORRECT
/// physics, measured 2026-08-20 against the repaired `isrm_v1.2.1.zarr`: not an
/// external oracle but a regression lock on the document's own output, whose
/// authority comes from the NumPy oracle agreeing on the weights and from the
/// InMAP-faithful configuration having matched the live service to 8.9e-9
/// before the physics was corrected. It is what is actually CHECKED.
///
/// `rel`: cross-binding spread on this document is ~4e-18 relative, so this is
/// loose by many orders of magnitude and catches a real change, not float noise.
#[derive(Clone, Copy)]
struct Oracle {
    published: (f64, f64),
    corrected: (f64, f64),
    rel: f64,
}

/// Keyed by the tag a document puts in `metadata.x_esd.report.oracle`, so which
/// published result a run is measured against is the DOCUMENT's claim and not
/// this file's guess. A document with no `oracle` tag — `isrm_polygon.esm`,
/// whose example emission layer this repository builds and nobody has published
/// a total for — is reported, not graded.
fn oracle_for(tag: &str) -> Option<Oracle> {
    match tag {
        "inmap_sr_tutorial" => Some(Oracle {
            published: (6928.959583, 15623.924632),
            corrected: (7022.724781368745, 15835.993595627131),
            rel: 1e-9,
        }),
        _ => None,
    }
}

fn report_block(doc: &Value) -> Result<Report, String> {
    let rep = doc
        .pointer("/metadata/x_esd/report")
        .ok_or("no metadata.x_esd.report block — this runner reads the reported \
                pathway and observed names from the document rather than carrying \
                a table of its own")?;
    let str_at = |v: &Value, k: &str| -> Result<String, String> {
        v[k].as_str()
            .map(str::to_string)
            .ok_or_else(|| format!("metadata.x_esd.report: `{k}` must be a string"))
    };
    let list_at = |v: &Value, k: &str| -> Result<Vec<String>, String> {
        v[k].as_array()
            .map(|a| a.iter().filter_map(|x| x.as_str().map(str::to_string)).collect())
            .ok_or_else(|| format!("metadata.x_esd.report: `{k}` must be a list of names"))
    };
    let pathways = rep["pathways"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or("metadata.x_esd.report.pathways must be a non-empty list")?
        .iter()
        .map(|e| {
            Ok(PathwaySpec {
                sr_array: str_at(e, "sr_array")?,
                emissions: list_at(e, "emissions")?,
                concentration: str_at(e, "concentration")?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let plume = match rep.get("plume") {
        Some(p) if p.is_object() => Some((
            str_at(p, "sr_lower")?,
            str_at(p, "stack_layer")?,
            list_at(p, "weights")?,
        )),
        _ => None,
    };
    Ok(Report {
        pathways,
        total_pm25: str_at(rep, "total_pm25")?,
        deaths_k: str_at(&rep["deaths"], "krewski")?,
        deaths_l: str_at(&rep["deaths"], "lepeule")?,
        record_field: str_at(rep, "record_field")?,
        plume,
        oracle: rep["oracle"].as_str().and_then(oracle_for),
    })
}

fn main() {
    if let Err(e) = run() {
        eprintln!("\nERROR: {e}");
        std::process::exit(1);
    }
}

/// VmHWM (peak resident set) from /proc, in bytes.
fn peak_rss_bytes() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find(|l| l.starts_with("VmHWM:")).and_then(|l| {
                l.split_whitespace()
                    .nth(1)
                    .and_then(|v| v.parse::<u64>().ok())
            })
        })
        .map(|kb| kb * 1024)
        .unwrap_or(0)
}

/// A metaparameter's declared default, read from the document (so no grid
/// extent is written down here).
fn metaparam(doc: &Value, name: &str) -> usize {
    doc["metaparameters"][name]["default"].as_u64().unwrap_or(0) as usize
}

/// The loaders that DISCOVER their own extent (`extent.metaparameter`) — the
/// record-bearing tables of the document, whatever they happen to be called.
/// The two knobs below are scale/locality concerns of a RUN, not of the model,
/// and both are expressed in the document's own vocabulary.
fn record_loaders(doc: &Value) -> Vec<String> {
    doc["data_sources"]
        .as_object()
        .map(|dls| {
            dls.iter()
                .filter(|(_, ld)| ld.pointer("/extent/metaparameter").is_some())
                .map(|(name, _)| name.clone())
                .collect()
        })
        .unwrap_or_default()
}

/// How many (source array, gate) pairs the DOCUMENT declares — the number of
/// model arrays the pushdown rewrite must end up gating, summed over the
/// sources that declare a `gated_select`. Derived from the document rather
/// than written down, so splitting or merging a gated source keeps the check
/// honest.
fn declared_gated_arrays(doc: &Value) -> usize {
    doc["data_sources"]
        .as_object()
        .map(|dls| {
            dls.values()
                .filter_map(|ld| ld.pointer("/metadata/x_esd/gated_select/applies_to"))
                .filter_map(Value::as_array)
                .map(|a| a.len())
                .sum()
        })
        .unwrap_or(0)
}

/// Every model PARAMETER fed by data source `src` — `(model, parameter)` pairs.
///
/// From esm 1.0.0 a source declares no variables of its own: the binding lives
/// on the consuming parameter as `update: {kind: "data", source: …}`
/// (esm-spec §6.3), so "what does this source deliver" is answered by walking
/// the models, not the source.
fn source_parameters(doc: &Value, src: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    if let Some(models) = doc["models"].as_object() {
        for (mname, m) in models {
            if let Some(vars) = m["variables"].as_object() {
                for (vname, v) in vars {
                    let up = &v["update"];
                    if up["kind"] == "data" && up["source"] == src {
                        out.push((mname.clone(), vname.clone()));
                    }
                }
            }
        }
    }
    out
}

/// REDUCED runs: truncate every record-discovering source to its first `n`
/// DELIVERED records with a `select` range (esm-spec §8.9.2). Because the
/// selection follows the source's own `record_filter`, this picks the same
/// records a post-filter `truncate(n)` would — and `extent` then re-discovers
/// the smaller N_REC by itself.
///
/// Written on each consuming PARAMETER rather than on the source, because
/// `select.axes` is one entry per NATIVE array dimension and a source may
/// deliver arrays of different ranks: the polygon layer's `geometry` is
/// `[record, vertex, xy]` while its emission column is `[record]`, so no single
/// source-level list is right for both. The rank comes from the parameter's own
/// declared `shape`, and the record axis is axis 0 by definition of a record
/// table.
fn truncate_records(doc: &mut Value, n: usize) {
    for name in record_loaders(doc) {
        for (mname, vname) in source_parameters(doc, &name) {
            let rank = doc["models"][&mname]["variables"][&vname]["shape"]
                .as_array()
                .map(|a| a.len())
                .unwrap_or(1);
            let mut axes = vec![json!({"range": {"start": 0, "stop": n}})];
            axes.extend(std::iter::repeat_n(json!("all"), rank.saturating_sub(1)));
            doc["models"][&mname]["variables"][&vname]["update"]["from"]["select"] =
                json!({ "axes": axes });
        }
    }
}

fn run() -> Result<(), String> {
    let t0 = Instant::now();
    let firstn: Option<usize> = std::env::var("ISRM_FIRSTN")
        .ok()
        .and_then(|s| s.parse().ok());
    let reduced = firstn.is_some();
    let model_path = paths::model();
    let mut doc: Value = serde_json::from_str(
        &std::fs::read_to_string(&model_path).map_err(|e| format!("read {model_path:?}: {e}"))?,
    )
    .map_err(|e| format!("parse {model_path:?}: {e}"))?;
    let report = report_block(&doc).map_err(|e| format!("{}: {e}", model_path.display()))?;
    if reduced {
        println!("REDUCED run — first {} emission records", firstn.unwrap());
    } else if let Some(o) = report.oracle {
        println!(
            "FULL run — whole domain (correct physics: deathsK≈{:.2}, deathsL≈{:.2}; the \
             tutorial's {:.2}/{:.2} is a reference, not a target)",
            o.corrected.0, o.corrected.1, o.published.0, o.published.1
        );
    } else {
        println!("FULL run — whole domain (the document claims no published total)");
    }
    println!("model:   {}", model_path.display());
    println!("scratch: {}", paths::scratch().display());
    println!("cache:   {}", paths::esio_cache().display());
    if let Some(n) = firstn {
        truncate_records(&mut doc, n);
    }

    // ---- providers FROM THE DOCUMENT — ALL of them --------------------------
    // Including the FF10 table: the loader declares its own reader options,
    // code map, record filter and extent, so there is nothing left here to
    // read, map, filter or count.
    println!("building providers from the document ...");
    let t = Instant::now();
    let cache_root = paths::esio_cache();
    // A local copy of a record source is a LOCALITY choice of this run
    // (gaftp.epa.gov is slow and flaky; the example polygon layer is built in
    // this repository and never needs fetching at all), so it is a url_override
    // rather than an edit to the document — and it is matched by the document's
    // OWN url basename, so this file names no source.
    let mut url_overrides: HashMap<String, String> = HashMap::new();
    for name in record_loaders(&doc) {
        let url = doc["data_sources"][&name]["source"]["url_template"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        if let Some(local) = paths::local_mirror(&url) {
            println!("  {name} mirrored from {}", local.display());
            url_overrides.insert(name, format!("file://{}", local.display()));
        }
    }
    let providers = providers_from_document(&doc, &cache_root, None, &url_overrides)
        .map_err(|e| e.to_string())?;
    println!(
        "  providers: {:?}",
        providers.iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>()
    );
    let t_providers = t.elapsed().as_secs_f64();

    // ---- PREPARE (extent -> rewrite -> coords -> VI -> gated fetch -> graph)-
    println!("esm_problem(pushdown_rewrite=true) — N_REC discovered by the loader ...");
    let t = Instant::now();
    let boxed: Vec<(String, Box<dyn PrepareProvider>)> = providers
        .into_iter()
        .map(|(k, p)| (k, Box::new(p) as Box<dyn PrepareProvider>))
        .collect();
    // EarthSciAST phase 4: `prepare` is replaced by problem CONSTRUCTION, which
    // absorbs the same pipeline. Providers move from a positional argument onto
    // `build_providers`, and `tspan` is the one genuinely new input -- this
    // driver never integrates, so it is nominal and nothing reads it. The
    // document declares no `D(...)` equation, so the default `Compile::Auto`
    // skips the right-hand-side compile, and with `default-features = false`
    // this crate never links the solver at all.
    let opts = ProblemOptions {
        base_path: model_path.parent().map(|p| p.to_path_buf()),
        pushdown_rewrite: true,
        verbose: true,
        build_providers: boxed,
        ..Default::default()
    };
    let prep = esm_problem(&doc, (0.0, 1.0), opts).map_err(|e| e.to_string())?;

    // ---- the gate covers EVERY declared SR array ---------------------------
    // A malformed `E_*` or `conc_*` body does not fail: the pathway simply
    // drops out of the rewrite's `applies_to` list, the rest of the rewrite
    // reports success, and the un-gated array is then fetched WHOLE — 330 GB,
    // which surfaces hours later as a memory failure rather than an error.
    // The document says how many arrays it declared for gating; anything less
    // here is that silent drop, so stop on it.
    let expect_gated = declared_gated_arrays(&doc);
    let gated: Vec<String> = prep.document()["metadata"]["x_esd"]["pushdown"]["gated_select"]
        ["applies_to"]
        .as_array()
        .map(|a| {
            a.iter()
                .filter_map(|v| v.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default();
    println!("gated arrays: {} of {expect_gated} declared", gated.len());
    if gated.len() != expect_gated {
        return Err(format!(
            "the pushdown rewrite gated {} arrays but the document declares {expect_gated} \
             ({gated:?}) — a pathway dropped out of the gate silently, and its SR array \
             would be fetched UNGATED. Check that each conc_*_L* body is a plain two-factor \
             SR*E product and that the containment ifelse is the FIRST ifelse in every E_* body.",
            gated.len()
        ));
    }

    let n_rec = observed_field(&prep, &report.record_field)
        .map(|a| a.len())
        .map_err(|e| format!("no field over the record axis to size N_REC from: {e}"))?;
    let t_prep = t.elapsed().as_secs_f64();
    println!(
        "PREPARE done in {t_prep:.1} s  (peak RSS so far: {:.2} GiB)",
        peak_rss_bytes() as f64 / (1u64 << 30) as f64
    );

    // ---- the engine-derived support set (for the contract record) -----------
    let producer_id = prep.document()["metadata"]["x_esd"]["pushdown"]["producer_id"]
        .as_str()
        .ok_or("no metadata.x_esd.pushdown record — did the rewrite fire?")?
        .to_string();
    let mut members: Vec<i64> = prep
        .members()
        .get(&producer_id)
        .ok_or_else(|| format!("no value-invention members for {producer_id}"))?
        .clone();
    members.sort_unstable();
    let n_ppl = members.len();
    let n_src = metaparam(&doc, "N_SRC");
    let n_rcv = metaparam(&doc, "N_RCV");
    println!("engine-derived support set: |members| = {n_ppl} of {n_src} source cells");
    if !reduced && report.oracle.is_some() && n_ppl != 1520 {
        println!("  WARNING: expected 1520 emission-bearing cells at full scale");
    }

    // ---- results through the prepared document's own graph ------------------
    let field = |n: &str| -> Result<Vec<f64>, String> {
        Ok(observed_field(&prep, n)
            .map_err(|e| e.to_string())?
            .iter()
            .copied()
            .collect())
    };
    let dk = field(&report.deaths_k)?;
    let dl = field(&report.deaths_l)?;
    let tp = field(&report.total_pm25)?;
    let mut pathways = BTreeMap::new();
    let mut emis_by_layer: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for spec in &report.pathways {
        // `emissions` is a LIST: isrm_point.esm splits a pathway across the
        // three SR emission layers a record's plume falls between,
        // isrm_polygon.esm has one because an area source emits at the ground.
        // `emis_sum` is the pathway TOTAL either way — plume rise moves mass
        // between layers, never into or out of a pathway — so it stays
        // comparable across both, and to the ground-level-only baselines.
        let mut e: Vec<f64> = Vec::new();
        let mut by_layer: Vec<f64> = Vec::new();
        for evar in &spec.emissions {
            let el = field(evar)?;
            by_layer.push(contract::compensated_sum(&el));
            e.extend(el);
        }
        let c = field(&spec.concentration)?;
        pathways.insert(
            spec.sr_array.clone(),
            contract::Pathway {
                emis_sum: contract::compensated_sum(&e),
                conc_sum: contract::compensated_sum(&c),
                conc_max: c.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            },
        );
        emis_by_layer.insert(spec.sr_array.clone(), by_layer);
    }

    // The layer assignment itself — a SPLIT, not a single layer: InMAP's
    // sr.Reader.layerFracs charges a record to two SR layers whenever its model
    // layer falls between two entries of `layers`. `sr_lower` is the lower index
    // (integer, compared exactly) and w_sr0/1/2 are the three shares. These are
    // the document's OWN observeds, read through the same `observed_field` path
    // as everything else — this runner does not know what ASME is, and must
    // not: the point of the contract's `plume` block is that the ENGINE produced
    // the assignment from the spec. contract/plume_oracle.py computes the same
    // quantity independently, from the meteorology arrays and without the SR
    // matrix, and compare_results.py checks the two against each other. A
    // document with no plume rise to state — an area source has no stack —
    // declares no `plume` key and emits no `plume` block.
    let plume: Option<Value> = match &report.plume {
        Some((lower, stack, weights)) => {
            let sr_lower = field(lower)?;
            let stack_layer = field(stack)?;
            let w: Vec<Vec<f64>> = weights
                .iter()
                .map(|n| field(n))
                .collect::<Result<_, String>>()?;
            if w.len() != 3 {
                return Err(format!(
                    "metadata.x_esd.report.plume.weights names {} fields; the SR matrix \
                     has three emission layers and the contract's plume block has three \
                     weight slots",
                    w.len()
                ));
            }
            Some(contract::plume_block(
                &sr_lower,
                &stack_layer,
                &[&w[0], &w[1], &w[2]],
                &emis_by_layer,
            )?)
        }
        None => None,
    };

    // Compensated throughout, so the reported totals are a property of the
    // data rather than of Rust's reduction order (contract::compensated_sum).
    let sk = contract::compensated_sum(&dk);
    let sl = contract::compensated_sum(&dl);
    println!("\n{}", "=".repeat(70));
    println!("  sum(deathsK) = {sk:?}");
    println!("  sum(deathsL) = {sl:?}");
    println!("  Σ TotalPM25  = {:?}", contract::compensated_sum(&tp));
    println!(
        "  Σ emitted    = {}",
        pathways
            .iter()
            .map(|(k, v)| format!("{k} {:?}", v.emis_sum))
            .collect::<Vec<_>>()
            .join(" / ")
    );
    if let (false, Some(o)) = (reduced, report.oracle) {
        let (ok, ol) = o.published;
        let (ck, cl) = o.corrected;
        println!(
            "  tutorial deathsK={ok}  deviation {:.6}%  (reference, not a target)",
            100.0 * (sk - ok) / ok
        );
        println!(
            "  tutorial deathsL={ol} deviation {:.6}%  (reference, not a target)",
            100.0 * (sl - ol) / ol
        );
        let rk = (sk - ck) / ck;
        let rl = (sl - cl) / cl;
        if rk.abs() > o.rel || rl.abs() > o.rel {
            println!(
                "  WARNING: {sk} / {sl} differs from the measured corrected-physics \
                 totals {ck} / {cl} by more than {:.0e} relative. That is a REGRESSION, \
                 not a tolerance: the two are the same document on the same store.",
                o.rel
            );
        } else {
            println!(
                "  matches the measured corrected-physics totals to {:.1e} / {:.1e}",
                rk.abs(),
                rl.abs()
            );
        }
    }
    if let Some(p) = &plume {
        println!(
            "  lower-SR-layer histogram (records per layer 0/1/2) = {}",
            p["sr_lower"]["histogram"]
        );
        println!(
            "  sr_lower sha256 = {}",
            p["sr_lower"]["sha256"].as_str().unwrap_or("?")
        );
        println!(
            "  Σ w_sr0/w_sr1/w_sr2 = {} / {} / {}   max|Σw - 1| = {}",
            p["weights"]["w_sr0"]["sum"],
            p["weights"]["w_sr1"]["sum"],
            p["weights"]["w_sr2"]["sum"],
            p["weights"]["max_sum_error"]
        );
        println!(
            "    (check it against `python3 contract/plume_oracle.py{}` — no SR matrix needed)",
            if reduced {
                format!(" --firstn {n_rec}")
            } else {
                String::new()
            }
        );
    }
    println!("{}", "=".repeat(70));

    // ---- contract record ----------------------------------------------------
    // Named after the MODEL, because two documents share this shim and must not
    // share one record file: isrm_point.esm and isrm_polygon.esm answer
    // different questions over the same grid.
    let stem = model_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "model".to_string());
    let out = paths::rs_dir().join(format!(
        "results_{stem}{}.json",
        if reduced { "_reduced" } else { "" }
    ));
    let mut timing = BTreeMap::new();
    timing.insert("wall_seconds".to_string(), t0.elapsed().as_secs_f64());
    timing.insert("providers_seconds".to_string(), t_providers);
    timing.insert("prepare_seconds".to_string(), t_prep);
    timing.insert("peak_rss_bytes".to_string(), peak_rss_bytes() as f64);
    contract::write_results(
        &out,
        &model_path.to_string_lossy(),
        "runtime_observed_graph",
        n_src,
        n_rcv,
        n_rec,
        &members,
        &pathways,
        &tp,
        &dk,
        &dl,
        &format!("rust / earthsci-ast {}", earthsci_ast::LIBRARY_VERSION),
        plume.as_ref(),
        &timing,
    )?;

    Ok(())
}
