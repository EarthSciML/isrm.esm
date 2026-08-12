// =============================================================================
// run-rs — the RUST binding drives the single clean `isrm.esm` end to end
// through the PUBLIC earthsci_ast surface. NOTHING MODEL-SHAPED LIVES HERE:
// this file names no pollutant, no column, no grid extent and no record count.
//
//   * `prepare(doc, {}, providers, pushdown_rewrite: true)` — the automatic
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
//   FULL run  (default)          -> assert sum(deathsK/L) ≈ 7524.92 / 16979.63
//   REDUCED   (ISRM_FIRSTN=n)    -> first n emission records, totals reported
//
// Emits the cross-language contract record (contract/results_schema.json) with
// model="isrm.esm", mode="runtime_observed_graph", binding="rust".
// =============================================================================

mod contract;
mod paths;

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use earthsci_ast::esio_provider::providers_from_document;
use earthsci_ast::prepare::{PrepareOptions, PrepareProvider, prepare};
use serde_json::{Value, json};

/// zarr array name -> (emissions observed, concentration observed)
const PW_OBS: [(&str, &str, &str); 5] = [
    ("SOA", "E_VOC", "conc_SOA"),
    ("pNO3", "E_NOx", "conc_pNO3"),
    ("pNH4", "E_NH3", "conc_pNH4"),
    ("pSO4", "E_SOx", "conc_pSO4"),
    ("PrimaryPM25", "E_PM25", "conc_PrimaryPM25"),
];

const ORACLE_K: f64 = 7524.918845602511;
const ORACLE_L: f64 = 16979.632171487083;

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
    doc["data_loaders"]
        .as_object()
        .map(|dls| {
            dls.iter()
                .filter(|(_, ld)| ld.pointer("/extent/metaparameter").is_some())
                .map(|(name, _)| name.clone())
                .collect()
        })
        .unwrap_or_default()
}

/// REDUCED runs: truncate every record-discovering loader to its first `n`
/// DELIVERED records with a loader-level `select` range (esm-spec §8.9.2).
/// Because the selection follows the loader's own `record_filter`, this picks
/// the same records the previous runners' post-filter `truncate(n)` did — and
/// `extent` then re-discovers the smaller N_REC by itself.
fn truncate_records(doc: &mut Value, n: usize) {
    for name in record_loaders(doc) {
        doc["data_loaders"][&name]["select"] =
            json!({"axes": [{"range": {"start": 0, "stop": n}}]});
    }
}

fn run() -> Result<(), String> {
    let t0 = Instant::now();
    let firstn: Option<usize> = std::env::var("ISRM_FIRSTN")
        .ok()
        .and_then(|s| s.parse().ok());
    let reduced = firstn.is_some();
    if reduced {
        println!("REDUCED run — first {} emission records", firstn.unwrap());
    } else {
        println!("FULL run — whole domain (target deathsK≈{ORACLE_K:.2}, deathsL≈{ORACLE_L:.2})");
    }
    let model_path = paths::model();
    println!("model:   {}", model_path.display());
    println!("scratch: {}", paths::scratch().display());
    println!("cache:   {}", paths::esio_cache().display());

    let mut doc: Value = serde_json::from_str(
        &std::fs::read_to_string(&model_path).map_err(|e| format!("read {model_path:?}: {e}"))?,
    )
    .map_err(|e| format!("parse {model_path:?}: {e}"))?;
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
    // A local copy of a record loader's source is a LOCALITY choice of this
    // run (gaftp.epa.gov is slow and flaky), so it is a url_override rather
    // than an edit to the document.
    let mut url_overrides: HashMap<String, String> = HashMap::new();
    let zip = paths::egu_zip();
    if zip.is_file() {
        for name in record_loaders(&doc) {
            url_overrides.insert(name, format!("file://{}", zip.display()));
        }
        println!("  record source mirrored from {}", zip.display());
    }
    let providers = providers_from_document(&doc, &cache_root, None, &url_overrides)
        .map_err(|e| e.to_string())?;
    println!(
        "  providers: {:?}",
        providers.iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>()
    );
    let t_providers = t.elapsed().as_secs_f64();

    // ---- PREPARE (extent -> rewrite -> coords -> VI -> gated fetch -> graph)-
    println!("prepare(pushdown_rewrite=true) — N_REC discovered by the loader ...");
    let t = Instant::now();
    let boxed: Vec<(String, Box<dyn PrepareProvider>)> = providers
        .into_iter()
        .map(|(k, p)| (k, Box::new(p) as Box<dyn PrepareProvider>))
        .collect();
    let opts = PrepareOptions {
        base_path: model_path.parent().map(|p| p.to_path_buf()),
        pushdown_rewrite: true,
        verbose: true,
        ..Default::default()
    };
    let prep = prepare(&doc, HashMap::new(), boxed, &opts).map_err(|e| e.to_string())?;
    let n_rec = prep
        .observed_field("X")
        .map(|a| a.len())
        .map_err(|e| format!("no projected emission coordinate to size N_REC from: {e}"))?;
    let t_prep = t.elapsed().as_secs_f64();
    println!(
        "PREPARE done in {t_prep:.1} s  (peak RSS so far: {:.2} GiB)",
        peak_rss_bytes() as f64 / (1u64 << 30) as f64
    );

    // ---- the engine-derived support set (for the contract record) -----------
    let producer_id = prep.doc["metadata"]["x_esd"]["pushdown"]["producer_id"]
        .as_str()
        .ok_or("no metadata.x_esd.pushdown record — did the rewrite fire?")?
        .to_string();
    let mut members: Vec<i64> = prep
        .members
        .get(&producer_id)
        .ok_or_else(|| format!("no value-invention members for {producer_id}"))?
        .clone();
    members.sort_unstable();
    let n_ppl = members.len();
    let n_src = metaparam(&doc, "N_SRC");
    let n_rcv = metaparam(&doc, "N_RCV");
    println!("engine-derived support set: |members| = {n_ppl} of {n_src} source cells");
    if !reduced && n_ppl != 1520 {
        println!("  WARNING: expected 1520 emission-bearing cells at full scale");
    }

    // ---- results through the prepared document's own graph ------------------
    let field = |n: &str| -> Result<Vec<f64>, String> {
        Ok(prep
            .observed_field(n)
            .map_err(|e| e.to_string())?
            .iter()
            .copied()
            .collect())
    };
    let dk = field("deathsK")?;
    let dl = field("deathsL")?;
    let tp = field("TotalPM25")?;
    let mut pathways = BTreeMap::new();
    for (arr, evar, cvar) in PW_OBS {
        let e = field(evar)?;
        let c = field(cvar)?;
        pathways.insert(
            arr.to_string(),
            contract::Pathway {
                emis_sum: contract::compensated_sum(&e),
                conc_sum: contract::compensated_sum(&c),
                conc_max: c.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            },
        );
    }

    // Compensated throughout, so the reported totals are a property of the
    // data rather than of Rust's reduction order (contract::compensated_sum).
    let sk = contract::compensated_sum(&dk);
    let sl = contract::compensated_sum(&dl);
    println!("\n{}", "=".repeat(70));
    println!("  sum(deathsK) = {sk:?}");
    println!("  sum(deathsL) = {sl:?}");
    println!("  Σ TotalPM25  = {:?}", contract::compensated_sum(&tp));
    let mut ok = true;
    if !reduced {
        ok = (sk - ORACLE_K).abs() <= 1e-4 * ORACLE_K && (sl - ORACLE_L).abs() <= 1e-4 * ORACLE_L;
        println!("  target deathsK={ORACLE_K}  rel.err {:.6}%", 100.0 * (sk - ORACLE_K) / ORACLE_K);
        println!("  target deathsL={ORACLE_L} rel.err {:.6}%", 100.0 * (sl - ORACLE_L) / ORACLE_L);
        println!("PHASE 4 FULL: {}", if ok { "PASS" } else { "FAIL" });
    }
    println!("{}", "=".repeat(70));

    // ---- contract record ----------------------------------------------------
    let out = paths::rs_dir().join(if reduced {
        "results_reduced.json"
    } else {
        "results.json"
    });
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
        &format!("rust / earthsci-ast {}", env!("CARGO_PKG_VERSION")),
        &timing,
    )?;

    if !ok {
        return Err("PHASE 4 FULL: oracle MISMATCH".to_string());
    }
    Ok(())
}
