// =============================================================================
// run-rs — Phase 4 (clean consolidation): the RUST binding drives the single
// clean `isrm.esm` end to end through the PUBLIC earthsci_ast surface.
//
//   * `prepare(doc, ca, providers, pushdown_rewrite: true)` — the automatic
//     projection-pushdown rewrite runs inside the engine; the SR provider
//     gates derive from the rewrite's own record
//     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
//     gate and hand-builds NO `Selection`;
//   * the SR / grid / pop / mortality providers come FROM THE DOCUMENT
//     (`providers_from_document`: format = `metadata.esio_format`, URL =
//     `source.url_template`);
//   * the EGU FF10 records are read through EarthSciIO's OWN ff10 reader with
//     the document-declared container options (`metadata.x_esd`: the zip
//     `member_filter` glob and `skip_header_row`) — no hand unzip, no hand
//     header-strip; the POLID→code map is the document's `pollutant_codes`;
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
use earthsciio::format::{ArrayData, Ff10Reader, Reader, Selection};
use ndarray::{ArrayD, IxDyn};
use serde_json::Value;

/// zarr array name -> (emissions observed, concentration observed)
const PW_OBS: [(&str, &str, &str); 5] = [
    ("SOA", "E_VOC", "conc_SOA"),
    ("pNO3", "E_NOx", "conc_pNO3"),
    ("pNH4", "E_NH3", "conc_pNH4"),
    ("pSO4", "E_SOx", "conc_pSO4"),
    ("PrimaryPM25", "E_PM25", "conc_PrimaryPM25"),
];

const N_SRC: usize = 52411; // metaparameter default; SR source/receptor cells
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

fn as_f64_vec(d: &ArrayData) -> Vec<f64> {
    match d {
        ArrayData::F64(v) => v.clone(),
        ArrayData::I64(v) => v.iter().map(|&x| x as f64).collect(),
        ArrayData::I32(v) => v.iter().map(|&x| x as f64).collect(),
        ArrayData::Bool(v) => v.iter().map(|&x| if x { 1.0 } else { 0.0 }).collect(),
        ArrayData::Str(_) => Vec::new(),
    }
}

fn arr1(v: Vec<f64>) -> ArrayD<f64> {
    let n = v.len();
    ArrayD::from_shape_vec(IxDyn(&[n]), v).expect("1-D shape always matches its own length")
}

/// EGU FF10 ingest — through EarthSciIO's ff10 reader, options FROM THE
/// DOCUMENT (`x_esd`: zip `member_filter` glob, `skip_header_row`, the
/// POLID→pathway-code integer map). Input PLUMBING only: read the table, map
/// the codes, drop the unrecognised/non-finite records (per the document's
/// note), truncate for a reduced run. The reader concatenates members in
/// sorted order — identical to the zip order for this archive, so ISRM_FIRSTN
/// selects the same records as every previous runner.
#[allow(clippy::type_complexity)]
fn read_egu(
    egu_meta: &Value,
    firstn: Option<usize>,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String> {
    let xe = &egu_meta["metadata"]["x_esd"];
    let mut zip = paths::egu_zip();
    if !zip.is_file() {
        // Fall back to the document's declared source, fetched once through
        // the EarthSciIO cache (content-addressed under the cache root).
        let url = egu_meta["source"]["url_template"]
            .as_str()
            .ok_or("EGU_Emis source.url_template missing")?;
        println!("  EGU zip not found at {} — fetching {url} via the cache ...", zip.display());
        let cache = earthsciio::Cache::builder()
            .data_dir(paths::esio_cache().join("EGU_Emis"))
            .build()
            .map_err(|e| format!("EGU cache: {e}"))?;
        let blob = cache
            .fetch(&earthsciio::FetchRequest::new(url).loader("EGU_Emis"))
            .map_err(|e| format!("EGU fetch {url}: {e}"))?;
        zip = blob.path;
    }
    let reader = Ff10Reader::new()
        .member_glob(xe["member_filter"].as_str().ok_or("EGU_Emis x_esd.member_filter missing")?)
        .skip_header_row(xe["skip_header_row"].as_bool().unwrap_or(false));
    let vars: Vec<String> = ["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let ds = reader
        .read_native(&zip, &vars, &Selection::All)
        .map_err(|e| format!("FF10 decode {zip:?}: {e}"))?;
    let get = |n: &str| {
        ds.variables
            .get(n)
            .ok_or_else(|| format!("FF10: no {n} column"))
    };
    let polid = match &get("POLID")?.data {
        ArrayData::Str(v) => v.clone(),
        other => return Err(format!("FF10: POLID decoded as {other:?}, expected strings")),
    };
    let ann = as_f64_vec(&get("ANN_VALUE")?.data);
    let lon = as_f64_vec(&get("LONGITUDE")?.data);
    let lat = as_f64_vec(&get("LATITUDE")?.data);

    let codes: HashMap<String, f64> = xe["pollutant_codes"]
        .as_object()
        .ok_or("EGU_Emis x_esd.pollutant_codes missing")?
        .iter()
        .map(|(k, v)| (k.to_uppercase(), v.as_f64().unwrap_or(0.0)))
        .collect();

    let n_all = polid.len();
    let (mut lo, mut la, mut an, mut co) = (vec![], vec![], vec![], vec![]);
    for i in 0..n_all {
        let c = codes
            .get(polid[i].trim().to_uppercase().as_str())
            .copied()
            .unwrap_or(0.0);
        if c > 0.0 && lon[i].is_finite() && lat[i].is_finite() && ann[i].is_finite() {
            lo.push(lon[i]);
            la.push(lat[i]);
            an.push(ann[i]);
            co.push(c);
        }
    }
    println!("  EGU FF10 records: {n_all} read, {} recognised-pathway kept", lo.len());
    if let Some(n) = firstn {
        lo.truncate(n);
        la.truncate(n);
        an.truncate(n);
        co.truncate(n);
    }
    Ok((lo, la, an, co))
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

    let doc: Value = serde_json::from_str(
        &std::fs::read_to_string(&model_path).map_err(|e| format!("read {model_path:?}: {e}"))?,
    )
    .map_err(|e| format!("parse {model_path:?}: {e}"))?;

    // ---- inputs: EGU records via the document-declared FF10 reader options --
    println!("reading EGU emissions (EarthSciIO ff10: member_glob + skip_header_row) ...");
    let t = Instant::now();
    let (lon, lat, ann, code) = read_egu(&doc["data_loaders"]["EGU_Emis"], firstn)?;
    let n_rec = lon.len();
    let t_egu = t.elapsed().as_secs_f64();
    println!("  N_REC = {n_rec}  in {t_egu:.1} s");

    // ---- providers FROM THE DOCUMENT ----------------------------------------
    println!("building providers from the document (ISRM_SR) ...");
    let cache_root = paths::esio_cache();
    let providers = providers_from_document(&doc, &cache_root, Some(&["ISRM_SR"]), &HashMap::new())
        .map_err(|e| e.to_string())?;
    println!(
        "  providers: {:?}",
        providers.iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>()
    );

    // ---- src-cell rectangles: the [0:N_SRC] prefix of the full grid ---------
    // The document declares src_W/src_S/src_E/src_N ([src_cells]) alongside the
    // full-grid W/S/E/N ([pop_cells]); the loader note pins the prefix
    // relationship. Slicing is input plumbing, not model math. (The rect
    // providers are sampled from a second document-built set; the shared cache
    // makes the re-read inside prepare a local hit.)
    println!("materializing grid prefix (src_W/S/E/N) ...");
    let t = Instant::now();
    let mut ca: HashMap<String, ArrayD<f64>> = HashMap::new();
    let mut rect_provs =
        providers_from_document(&doc, &cache_root, Some(&["ISRM_SR"]), &HashMap::new())
            .map_err(|e| e.to_string())?;
    for k in ["W", "S", "E", "N"] {
        let key = format!("ISRM_SR.{k}");
        let p = rect_provs
            .iter_mut()
            .find(|(pk, _)| pk == &key)
            .ok_or_else(|| format!("no document provider for {key}"))?;
        let full = p.1.sample().map_err(|e| format!("sample {key}: {e}"))?;
        let v: Vec<f64> = full.iter().take(N_SRC).copied().collect();
        ca.insert(format!("src_{k}"), arr1(v));
    }
    let t_grid = t.elapsed().as_secs_f64();
    println!("  grid prefix in {t_grid:.1} s");

    // ---- const arrays: the EGU loader's variables, src rects ----------------
    ca.insert("EGU_Emis.lon".into(), arr1(lon));
    ca.insert("EGU_Emis.lat".into(), arr1(lat));
    ca.insert("EGU_Emis.annual".into(), arr1(ann));
    ca.insert("EGU_Emis.pollutant".into(), arr1(code));
    for v in ["stkhgt", "stkdiam", "stktemp", "stkvel"] {
        ca.insert(format!("EGU_Emis.{v}"), arr1(vec![0.0; n_rec])); // declared, unused
    }

    // ---- PREPARE (rewrite -> coords -> VI -> gated SR fetch -> graph eval) --
    println!("prepare(pushdown_rewrite=true) — N_REC={n_rec} ...");
    let t = Instant::now();
    let boxed: Vec<(String, Box<dyn PrepareProvider>)> = providers
        .into_iter()
        .map(|(k, p)| (k, Box::new(p) as Box<dyn PrepareProvider>))
        .collect();
    let mut metaparameters = BTreeMap::new();
    metaparameters.insert("N_REC".to_string(), n_rec as i64);
    let opts = PrepareOptions {
        metaparameters,
        base_path: model_path.parent().map(|p| p.to_path_buf()),
        pushdown_rewrite: true,
        verbose: true,
        ..Default::default()
    };
    let prep = prepare(&doc, ca, boxed, &opts).map_err(|e| e.to_string())?;
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
    println!("engine-derived support set: |members| = {n_ppl} of {N_SRC} source cells");
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
    timing.insert("egu_seconds".to_string(), t_egu);
    timing.insert("grid_seconds".to_string(), t_grid);
    timing.insert("prepare_seconds".to_string(), t_prep);
    timing.insert("peak_rss_bytes".to_string(), peak_rss_bytes() as f64);
    contract::write_results(
        &out,
        &model_path.to_string_lossy(),
        "runtime_observed_graph",
        N_SRC,
        N_SRC,
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
