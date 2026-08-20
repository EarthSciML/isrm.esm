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
//   FULL run  (default)          -> report sum(deathsK/L) against the tutorial
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

/// zarr array name -> (the per-SR-layer emissions observeds, concentration
/// observed). The emissions side grew a LAYER dimension when the document
/// started stating plume rise: a record is charged to the SR layer its plume
/// reaches, so a pathway's emissions are three arrays, not one. The
/// concentration side did not — `conc_<p>` is the document's own sum over the
/// three layered contractions.
const PW_OBS: [(&str, [&str; 3], &str); 5] = [
    ("SOA", ["E_VOC_L0", "E_VOC_L1", "E_VOC_L2"], "conc_SOA"),
    ("pNO3", ["E_NOx_L0", "E_NOx_L1", "E_NOx_L2"], "conc_pNO3"),
    ("pNH4", ["E_NH3_L0", "E_NH3_L1", "E_NH3_L2"], "conc_pNH4"),
    ("pSO4", ["E_SOx_L0", "E_SOx_L1", "E_SOx_L2"], "conc_pSO4"),
    (
        "PrimaryPM25",
        ["E_PM25_L0", "E_PM25_L1", "E_PM25_L2"],
        "conc_PrimaryPM25",
    ),
];

/// The InMAP source-receptor tutorial's published national totals
/// (<https://inmap.run/blog/2019/04/20/sr/>), which account for plume rise.
///
/// These are a TARGET, not an assertion. The document deliberately does not
/// reproduce InMAP's high-plume source-index defect (a plume above model layer
/// 7 keeps an index built in the coarse 9324-cell grid and is then read against
/// the 52411-cell ground grid), which misplaces 654 of 43650 records — 0.43% of
/// emitted mass — onto the wrong source cell. So the run is expected to land
/// NEAR rather than ON these, and the deviation is printed rather than failed.
const ORACLE_K: f64 = 6928.959583;
const ORACLE_L: f64 = 15623.924632;
/// Beyond this the deviation is no longer explainable by the clean-physics
/// choice above and is worth investigating.
/// MEASURED at full scale on 2026-08-19 (Julia): deathsK 6983.9385617781645
/// (+0.79%) and deathsL 15752.315804140908 (+0.82%) against the published
/// totals. The misplaced group is 0.43% of emitted mass but buys about twice
/// that in deaths, because putting it back on the cells the emissions came
/// from puts it back over people. This threshold sits just above what was
/// measured.
///
/// STALE, and left in place only as a loose upper bound: it was fitted to a
/// full-scale run that predates sr.Reader.layerFracs and the [0, 3, 6] fix, and
/// no full-scale run of the current document has been made. At reduced scale
/// the document now reproduces the live inmap cloud service to 9e-9, so the
/// real threshold is expected to be far tighter — see SERVICE_DEATHS in
/// contract/compare_results.py.
const ORACLE_NOTABLE_REL: f64 = 8.3e-3;

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

/// How many (loader array, gate) pairs the DOCUMENT declares — the number of
/// model arrays the pushdown rewrite must end up gating, summed over the
/// loaders that declare a `gated_select`. Derived from the document rather
/// than written down, so splitting or merging a gated loader keeps the check
/// honest.
fn declared_gated_arrays(doc: &Value) -> usize {
    doc["data_loaders"]
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

    // ---- the gate covers EVERY declared SR array ---------------------------
    // A malformed `E_*` or `conc_*` body does not fail: the pathway simply
    // drops out of the rewrite's `applies_to` list, the rest of the rewrite
    // reports success, and the un-gated array is then fetched WHOLE — 330 GB,
    // which surfaces hours later as a memory failure rather than an error.
    // The document says how many arrays it declared for gating; anything less
    // here is that silent drop, so stop on it.
    let expect_gated = declared_gated_arrays(&doc);
    let gated: Vec<String> = prep.doc["metadata"]["x_esd"]["pushdown"]["gated_select"]
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
    let mut emis_by_layer: BTreeMap<String, [f64; 3]> = BTreeMap::new();
    for (arr, evars, cvar) in PW_OBS {
        // A pathway's emissions are now spread over the three SR layers; the
        // record's `emis_sum` is the pathway TOTAL, so it stays comparable to
        // the ground-level-only baselines (plume rise moves mass between
        // layers, never into or out of a pathway).
        let mut e: Vec<f64> = Vec::new();
        // How much mass plume rise put in each SR layer — the physics made
        // visible as tons, per pathway.
        let mut by_layer = [0.0f64; 3];
        for (layer, evar) in evars.iter().copied().enumerate() {
            let el = field(evar)?;
            by_layer[layer] = contract::compensated_sum(&el);
            e.extend(el);
        }
        let c = field(cvar)?;
        pathways.insert(
            arr.to_string(),
            contract::Pathway {
                emis_sum: contract::compensated_sum(&e),
                conc_sum: contract::compensated_sum(&c),
                conc_max: c.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            },
        );
        emis_by_layer.insert(arr.to_string(), by_layer);
    }

    // The layer assignment itself — now a SPLIT, not a single layer: InMAP's
    // sr.Reader.layerFracs charges a record to two SR layers whenever its model
    // layer falls between two entries of `layers`. `sr_lower` is the lower index
    // (integer, compared exactly) and w_sr0/1/2 are the three shares. These are
    // the document's OWN observeds, read through the same `observed_field` path
    // as everything else — this runner does not know what ASME is, and must
    // not: the point of the contract's `plume` block is that the ENGINE produced
    // the assignment from the spec. contract/plume_oracle.py computes the same
    // quantity independently, from the meteorology arrays and without the SR
    // matrix, and compare_results.py checks the two against each other.
    let sr_lower = field("sr_lower")?;
    let stack_layer = field("stack_layer")?;
    let w0 = field("w_sr0")?;
    let w1 = field("w_sr1")?;
    let w2 = field("w_sr2")?;
    let plume = contract::plume_block(
        &sr_lower,
        &stack_layer,
        &[&w0, &w1, &w2],
        &emis_by_layer,
    )?;

    // Compensated throughout, so the reported totals are a property of the
    // data rather than of Rust's reduction order (contract::compensated_sum).
    let sk = contract::compensated_sum(&dk);
    let sl = contract::compensated_sum(&dl);
    println!("\n{}", "=".repeat(70));
    println!("  sum(deathsK) = {sk:?}");
    println!("  sum(deathsL) = {sl:?}");
    println!("  Σ TotalPM25  = {:?}", contract::compensated_sum(&tp));
    if !reduced {
        let rk = (sk - ORACLE_K) / ORACLE_K;
        let rl = (sl - ORACLE_L) / ORACLE_L;
        println!("  tutorial deathsK={ORACLE_K}  deviation {:.6}%", 100.0 * rk);
        println!("  tutorial deathsL={ORACLE_L} deviation {:.6}%", 100.0 * rl);
        if rk.abs() > ORACLE_NOTABLE_REL || rl.abs() > ORACLE_NOTABLE_REL {
            println!(
                "  WARNING: deviation exceeds {:.2}% — more than the above-layer-7 \
                 group has been measured to be worth (0.43% of emitted mass, \
                 +0.79%/+0.82% of deaths at full scale), so something else differs.",
                100.0 * ORACLE_NOTABLE_REL
            );
        }
    }
    println!(
        "  lower-SR-layer histogram (records per layer 0/1/2) = {}",
        plume["sr_lower"]["histogram"]
    );
    println!(
        "  sr_lower sha256 = {}",
        plume["sr_lower"]["sha256"].as_str().unwrap_or("?")
    );
    println!(
        "  Σ w_sr0/w_sr1/w_sr2 = {} / {} / {}   max|Σw - 1| = {}",
        plume["weights"]["w_sr0"]["sum"],
        plume["weights"]["w_sr1"]["sum"],
        plume["weights"]["w_sr2"]["sum"],
        plume["weights"]["max_sum_error"]
    );
    println!(
        "    (check it against `python3 contract/plume_oracle.py{}` — no SR matrix needed)",
        if reduced {
            format!(" --firstn {}", n_rec)
        } else {
            String::new()
        }
    );
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
        Some(&plume),
        &timing,
    )?;

    Ok(())
}
