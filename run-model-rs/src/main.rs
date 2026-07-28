//! ISRM through the RUST binding's own evaluation of `isrm_pushdown.esm`.
//!
//! The strict bar: every number in the result record comes from the binding's
//! own evaluation of the `.esm` graph — no hand-written STEP-0 arithmetic here.
//!
//!   * the emission-bearing source-cell support set (`ppl`) is DERIVED by
//!     `earthsci_ast`'s value-invention engine running the model's own producer
//!     aggregate, gated by the spatial `join.overlap` broad phase
//!     (CONFORMANCE_SPEC §5.5.6);
//!   * every field — the per-source binned emissions `E_*`, the per-receptor
//!     concentrations `conc_*`, `TotalPM25`, `deathsK`/`deathsL` — is produced
//!     by `eval_expression_with_extents` on the `.esm`'s own observed
//!     expressions;
//!   * this file contributes ORCHESTRATION only: which observed to evaluate
//!     when (from each aggregate's own declared `args`), and where bytes come
//!     from.
//!
//! Architecturally identical to `run-model-py/run_model.py`, which is the point:
//! two engines, one spec, independently reproduced numbers.
//!
//! Usage:
//!   run-model-rs                  # full run  -> results.json
//!   ISRM_FIRSTN=200 run-model-rs  # reduced    -> results_reduced.json

mod contract;
mod inputs;
mod paths;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;
use std::time::Instant;

use earthsci_ast::aggregate::resolve_expr_ranges_with_extents;
use earthsci_ast::esio_provider::EsioProvider;
use earthsci_ast::provider::CadenceProvider;
use earthsci_ast::simulate_array::{eval_expression_with_extents, run_value_invention, Value as EvalValue};
use earthsci_ast::types::{EsmFile, Expr, IndexSet, Model, VariableType};
use earthsciio::format::{AxisSelect, Selection};
use earthsciio::{Cache, DataLoader};
use ndarray::{ArrayD, IxDyn};

/// zarr array name -> (model SR variable, emissions observed, concentration observed)
const PATHWAYS: [(&str, &str, &str, &str); 5] = [
    ("SOA", "SR_SOA", "E_VOC", "conc_SOA"),
    ("pNO3", "SR_pNO3", "E_NOx", "conc_pNO3"),
    ("pNH4", "SR_pNH4", "E_NH3", "conc_pNH4"),
    ("pSO4", "SR_pSO4", "E_SOx", "conc_pSO4"),
    ("PrimaryPM25", "SR_PrimaryPM25", "E_PM25", "conc_PrimaryPM25"),
];

const ORACLE_K: f64 = 7524.918845602511;
const ORACLE_L: f64 = 16979.632171487083;

fn main() {
    if let Err(e) = run() {
        eprintln!("\nERROR: {e}");
        std::process::exit(1);
    }
}

/// Every declared 0-D parameter's default — `fact`, the risk ratios, the LCC
/// constants. Load-time constants the graph reads directly.
fn scalar_params(model: &Model) -> (Vec<f64>, Vec<String>) {
    let mut names: Vec<String> = model
        .variables
        .iter()
        .filter(|(_, v)| {
            v.var_type == VariableType::Parameter && v.shape.is_none() && v.default.is_some()
        })
        .map(|(k, _)| k.clone())
        .collect();
    names.sort();
    let vals = names
        .iter()
        .map(|n| model.variables[n].default.unwrap_or(0.0))
        .collect();
    (vals, names)
}

/// name -> expression, for every `observed` variable.
fn observed_defs(model: &Model) -> HashMap<String, Expr> {
    model
        .variables
        .iter()
        .filter(|(_, v)| v.var_type == VariableType::Observed && v.expression.is_some())
        .map(|(k, v)| (k.clone(), v.expression.clone().unwrap()))
        .collect()
}

/// The factor names an expression declares it reads (`aggregate.args`), gathered
/// recursively. The dependency edges come from the MODEL, not from a
/// re-derivation here.
fn declared_args(expr: &Expr, out: &mut HashSet<String>) {
    match expr {
        // A bare symbol is a factor reference — this is how an aggregate's
        // `args` list arrives, and also how `index(<name>, …)` names its array.
        Expr::Variable(name) => {
            out.insert(name.clone());
        }
        Expr::Operator(node) => {
            for a in &node.args {
                declared_args(a, out);
            }
            node.for_each_child(&mut |c| declared_args(c, out));
        }
        _ => {}
    }
}

/// Dependency order over the observeds: an observed follows every observed it
/// names. A cycle is a malformed model and is a hard error.
fn observed_order(defs: &HashMap<String, Expr>) -> Result<Vec<String>, String> {
    let deps: HashMap<String, HashSet<String>> = defs
        .iter()
        .map(|(n, e)| {
            let mut a = HashSet::new();
            declared_args(e, &mut a);
            a.retain(|x| defs.contains_key(x) && x != n);
            (n.clone(), a)
        })
        .collect();

    let mut ordered = Vec::new();
    let mut done: HashSet<String> = HashSet::new();
    let mut pending: HashSet<String> = defs.keys().cloned().collect();
    while !pending.is_empty() {
        let mut ready: Vec<String> = pending
            .iter()
            .filter(|n| deps[*n].iter().all(|d| done.contains(d)))
            .cloned()
            .collect();
        if ready.is_empty() {
            let mut rest: Vec<_> = pending.into_iter().collect();
            rest.sort();
            return Err(format!("cyclic observed dependency among {rest:?}"));
        }
        ready.sort(); // deterministic tie-break
        for n in ready {
            ordered.push(n.clone());
            done.insert(n.clone());
            pending.remove(&n);
        }
    }
    Ok(ordered)
}

/// Fetch `SR[layer 0, members, :]` for each pathway through the opt-in
/// EarthSciIO bridge.
///
/// `members1` are the 1-based full-grid source-cell ids value invention just
/// derived; EarthSciIO indexes 0-based. The selection goes down to the zarr
/// reader so only the intersecting chunks are fetched — the 1,520 rows that
/// matter, never the 52,411-row matrix. Pathways are fetched and evicted one at
/// a time so peak disk is one pathway (~2.7 GiB), not five.
fn fetch_sr(members1: &[i64]) -> Result<HashMap<String, ArrayD<f64>>, String> {
    let root = paths::scratch().join(
        std::env::var("ISRM_SR_DIR").unwrap_or_else(|_| "rs_l3_cache_sr".to_string()),
    );
    let idx0: Vec<usize> = members1.iter().map(|&m| (m - 1) as usize).collect();
    let sel = Selection::Orthogonal(vec![
        AxisSelect::Indices(vec![0]),
        AxisSelect::Indices(idx0),
        AxisSelect::All,
    ]);

    let mut out = HashMap::new();
    for (zname, mname, _, _) in PATHWAYS {
        let blobs = root.join("v1").join("blobs");
        if blobs.is_dir() {
            let _ = std::fs::remove_dir_all(&blobs);
        }
        let cache = Arc::new(
            Cache::builder()
                .data_dir(&root)
                .offline(false)
                .build()
                .map_err(|e| format!("cache: {e}"))?,
        );
        let loader = DataLoader::new(format!("ISRM_SR[{zname}]"), "zarr", paths::zarr_url())
            .variables([zname.to_string()]);
        let mut p = EsioProvider::builder(loader, cache)
            .var(zname, mname)
            .select(sel.clone())
            .build()
            .map_err(|e| format!("esio provider {zname}: {e}"))?;

        let t = Instant::now();
        let fields = p
            .materialize()
            .map_err(|e| format!("gated fetch {zname}: {e}"))?;
        let f = fields
            .get(mname)
            .ok_or_else(|| format!("gated fetch {zname}: {mname} missing"))?;
        // (1, |ppl|, rcv) -> (|ppl|, rcv): the model's SR is
        // [emis_src_cells, rcv_cells]; emisLayer 0 is already sliced by the
        // selection.
        let shape = f.array.shape().to_vec();
        let arr = if shape.len() == 3 && shape[0] == 1 {
            f.array
                .clone()
                .into_shape_with_order(IxDyn(&[shape[1], shape[2]]))
                .map_err(|e| format!("reshape {zname}: {e}"))?
        } else {
            f.array.clone()
        };
        println!(
            "    [gated] fetched {zname} -> {:?}  in {:.1} s",
            arr.shape(),
            t.elapsed().as_secs_f64()
        );
        out.insert(mname.to_string(), arr);
    }
    let _ = std::fs::remove_dir_all(root.join("v1").join("blobs"));
    Ok(out)
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

    println!("building inputs ...");
    let inp = inputs::build_inputs(firstn)?;
    println!("  N_REC={}  N_SRC={}", inp.n_rec, inp.n_src);

    // ---- load the model with metaparameters bound ---------------------------
    // The document declares index-set sizes as metaparameter NAMES; binding them
    // through the loader API is the supported substitution (esm-spec §9.7.6),
    // and is what the Julia/Python runners do with their `resolve_sizes!`.
    let mut metaparameters = BTreeMap::new();
    metaparameters.insert("N_SRC".to_string(), inp.n_src as i64);
    metaparameters.insert("N_RCV".to_string(), inp.n_src as i64);
    metaparameters.insert("N_POP".to_string(), inp.n_src as i64);
    metaparameters.insert("N_LAYER".to_string(), 3i64);
    metaparameters.insert("N_REC".to_string(), inp.n_rec as i64);

    let model_path = paths::model();
    let json_text = std::fs::read_to_string(&model_path)
        .map_err(|e| format!("read {model_path:?}: {e}"))?;
    let opts = earthsci_ast::parse::LoadOptions {
        base_path: model_path.parent().map(|p| p.to_path_buf()),
        metaparameters,
    };
    let file: EsmFile = earthsci_ast::parse::load_with_options(&json_text, &opts)
        .map_err(|e| format!("load {model_path:?}: {e}"))?;
    let index_sets: HashMap<String, IndexSet> = file.index_sets.clone().unwrap_or_default();
    let model = file
        .models
        .as_ref()
        .and_then(|m| m.get("ISRM"))
        .ok_or("model 'ISRM' not found in the document")?
        .clone();

    let (param_vals, param_names) = scalar_params(&model);
    let mut arrays = inp.as_arrays();

    // ---- VALUE INVENTION: the graph derives its own support set -------------
    println!("value invention (spatial overlap gate -> emis_src_cells) ...");
    let t = Instant::now();
    let vi = run_value_invention(&model, &index_sets, Some(&arrays))
        .map_err(|e| format!("value invention: {e}"))?;
    let faq = "emis_src_cells_faq";
    let members: Vec<i64> = vi
        .members
        .get(faq)
        .ok_or_else(|| format!("value invention produced no members for {faq}"))?
        .iter()
        .map(|k| match k {
            earthsci_ast::relational::Key::Int(i) => Ok(*i),
            other => Err(format!("expected a scalar cell id member, got {other:?}")),
        })
        .collect::<Result<_, String>>()?;
    let n_ppl = members.len();
    println!(
        "  |emis_src_cells| = {n_ppl}   in {:.1} s",
        t.elapsed().as_secs_f64()
    );
    if !reduced && n_ppl != 1520 {
        println!("  WARNING: expected 1520 emission-bearing cells at full scale, got {n_ppl}");
    }

    // Hook 1: the derived set's `member_factor` is fed back as a const factor so
    // the downstream cell_W/cell_S/… gathers read the compact derived axis.
    if let Some(mf) = index_sets
        .get("emis_src_cells")
        .and_then(|is| is.member_factor.clone())
    {
        let v: Vec<f64> = members.iter().map(|&m| m as f64).collect();
        arrays.insert(
            mf.clone(),
            ArrayD::from_shape_vec(IxDyn(&[v.len()]), v).map_err(|e| e.to_string())?,
        );
    }

    // ---- Hook 2: the gated fetch, driven by the members just derived --------
    println!("gated SR fetch: layer 0, {n_ppl} of {} source cells, all receptors", inp.n_src);
    arrays.extend(fetch_sr(&members)?);

    // ---- evaluate the observed graph ---------------------------------------
    let defs = observed_defs(&model);
    let order = observed_order(&defs)?;
    let extents = vi.extents.clone();
    println!("evaluating {} observeds through the .esm graph ...", order.len());

    let mut fields: HashMap<String, Vec<f64>> = HashMap::new();
    for name in &order {
        let mut expr = defs[name].clone();
        // Ranges over the invented axis need its extent before they can resolve.
        resolve_expr_ranges_with_extents(&mut expr, &index_sets, &extents)
            .map_err(|e| format!("resolve ranges for {name}: {e}"))?;
        let t = Instant::now();
        let val = eval_expression_with_extents(
            &expr,
            &arrays,
            &param_vals,
            &param_names,
            0.0,
            &extents,
        )
        .map_err(|e| format!("evaluate {name}: {e}"))?;
        let arr: ArrayD<f64> = match val {
            EvalValue::Array(a) => *a,
            EvalValue::Scalar(s) => ArrayD::from_elem(IxDyn(&[1]), s),
        };
        println!(
            "  {name:<18} shape={:?}  {:>7.1} s",
            arr.shape(),
            t.elapsed().as_secs_f64()
        );
        let flat: Vec<f64> = arr.iter().cloned().collect();
        // Feed each result back so downstream observeds resolve it by name — the
        // same dependency-ordered substitution the array build performs.
        arrays.insert(name.clone(), arr);
        fields.insert(name.clone(), flat);
    }

    let get = |n: &str| -> Result<&Vec<f64>, String> {
        fields.get(n).ok_or_else(|| format!("observed {n} not evaluated"))
    };
    let dk = get("deathsK")?;
    let dl = get("deathsL")?;
    let tp = get("TotalPM25")?;
    // Compensated throughout, so the reported totals are a property of the data
    // rather than of Rust's reduction order (see contract::compensated_sum).
    let sk: f64 = contract::compensated_sum(dk);
    let sl: f64 = contract::compensated_sum(dl);
    println!("\n{}", "=".repeat(70));
    println!("  sum(deathsK) = {sk:?}");
    println!("  sum(deathsL) = {sl:?}");
    println!("  Σ TotalPM25  = {:?}", contract::compensated_sum(tp));

    let mut ok = true;
    if !reduced {
        ok = (sk - ORACLE_K).abs() <= 1e-4 * ORACLE_K && (sl - ORACLE_L).abs() <= 1e-4 * ORACLE_L;
        println!(
            "  target deathsK={ORACLE_K}  rel.err {:.4}%",
            100.0 * (sk - ORACLE_K) / ORACLE_K
        );
        println!(
            "  target deathsL={ORACLE_L} rel.err {:.4}%",
            100.0 * (sl - ORACLE_L) / ORACLE_L
        );
        println!("RUST FULL: {}", if ok { "PASS" } else { "FAIL" });
    }
    println!("{}", "=".repeat(70));

    let mut pathways = BTreeMap::new();
    for (zname, _, evar, cvar) in PATHWAYS {
        let e = get(evar)?;
        let c = get(cvar)?;
        pathways.insert(
            zname.to_string(),
            contract::Pathway {
                emis_sum: contract::compensated_sum(e),
                conc_sum: contract::compensated_sum(c),
                conc_max: c.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            },
        );
    }

    let out = paths::rs_dir().join(if reduced {
        "results_reduced.json"
    } else {
        "results.json"
    });
    let mut timing = BTreeMap::new();
    timing.insert("wall_seconds".to_string(), t0.elapsed().as_secs_f64());
    contract::write_results(
        &out,
        &model_path.to_string_lossy(),
        "runtime_observed_graph",
        inp.n_src,
        inp.n_src,
        inp.n_rec,
        &members,
        &pathways,
        tp,
        dk,
        dl,
        &format!("earthsci-ast {}", env!("CARGO_PKG_VERSION")),
        &timing,
    )?;

    if !ok {
        return Err("RUST FULL: oracle MISMATCH".to_string());
    }
    Ok(())
}
