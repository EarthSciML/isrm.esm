//! Emit a `contract/results_schema.json` record from the Rust runner.
//!
//! The Rust mirror of `contract/results.jl` and `contract/results.py`. The
//! hashing and sampling rules MUST match `contract/compare_results.py`, which is
//! the definition of record. All three `sample_indices` implementations are
//! deliberately pure INTEGER arithmetic so the index set cannot drift between
//! languages' float rounding.

use std::collections::BTreeMap;

use serde_json::{json, Map, Value};
use sha2::{Digest, Sha256};

const SAMPLE_N: usize = 25;

/// Neumaier (compensated) summation.
///
/// The contract's `sum` must be a property of the DATA, not of the summing
/// language. Julia's `sum` is pairwise and CPython 3.12's `sum()` is
/// Neumaier-compensated — both near-exact — while Rust's `Iterator::sum` is a
/// naive sequential fold. On a reduced ISRM run all three `total_pm25` fields
/// were BIT-IDENTICAL (same sha256, same samples) and yet the naive Rust total
/// differed from the other two by 2.9e-13, purely from accumulation error.
///
/// That gap grows with element count and magnitude, so at full scale a naive
/// total could exceed the 1e-12 field tolerance and report a cross-language
/// disagreement where the fields are provably identical. Compensating here
/// removes the reduction algorithm as a source of difference.
pub fn compensated_sum(xs: &[f64]) -> f64 {
    let mut s = 0.0f64;
    let mut c = 0.0f64; // running compensation for lost low-order bits
    for &x in xs {
        let t = s + x;
        c += if s.abs() >= x.abs() {
            (s - t) + x // s is larger: low-order bits of x were lost
        } else {
            (x - t) + s // x is larger: low-order bits of s were lost
        };
        s = t;
    }
    s + c
}

/// The fixed 1-based sample indices every runner reports.
pub fn sample_indices(n_rcv: usize) -> Vec<usize> {
    let d = SAMPLE_N - 1;
    (0..SAMPLE_N)
        .map(|k| 1 + (k * (n_rcv - 1) + d / 2) / d)
        .collect()
}

/// sha256 over sorted 1-based ids as ASCII decimals joined by "," (no spaces).
pub fn ppl_sha256(ids: &[i64]) -> String {
    let mut v = ids.to_vec();
    v.sort_unstable();
    let s = v
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("{:x}", Sha256::digest(s.as_bytes()))
}

/// sha256 over a float field as little-endian IEEE-754 float64 bytes.
pub fn field_sha256(v: &[f64]) -> String {
    let mut h = Sha256::new();
    for x in v {
        h.update(x.to_le_bytes());
    }
    format!("{:x}", h.finalize())
}

/// Summarize one length-n_rcv field into the schema's FieldSummary shape.
fn field_summary(v: &[f64]) -> Value {
    let sum: f64 = compensated_sum(v);
    let min = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sample: Vec<f64> = sample_indices(v.len()).iter().map(|&i| v[i - 1]).collect();
    json!({
        "sum": sum,
        "min": min,
        "max": max,
        "sample": sample,
        "sha256": field_sha256(v),
    })
}

/// One pathway's intermediates.
pub struct Pathway {
    pub emis_sum: f64,
    pub conc_sum: f64,
    pub conc_max: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn write_results(
    path: &std::path::Path,
    model: &str,
    mode: &str,
    n_src: usize,
    n_rcv: usize,
    n_rec: usize,
    ppl: &[i64],
    pathways: &BTreeMap<String, Pathway>,
    total_pm25: &[f64],
    deaths_k: &[f64],
    deaths_l: &[f64],
    binding_version: &str,
    timing: &BTreeMap<String, f64>,
) -> Result<(), String> {
    if mode != "runtime_observed_graph" && mode != "oracle_step0" {
        return Err(format!(
            "mode must be \"runtime_observed_graph\" or \"oracle_step0\", got {mode}"
        ));
    }
    let mut ids = ppl.to_vec();
    ids.sort_unstable();

    let mut pw = Map::new();
    for (k, v) in pathways {
        pw.insert(
            k.clone(),
            json!({ "emis_sum": v.emis_sum, "conc_sum": v.conc_sum, "conc_max": v.conc_max }),
        );
    }

    let mut timing_map = Map::new();
    for (k, v) in timing {
        timing_map.insert(k.clone(), json!(v));
    }

    let basename = std::path::Path::new(model)
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| model.to_string());

    let rec = json!({
        "binding": "rust",
        "binding_version": binding_version,
        "model": basename,
        "mode": mode,
        "grid": { "n_src": n_src, "n_rcv": n_rcv, "n_rec": n_rec },
        "ppl": { "count": ids.len(), "sha256": ppl_sha256(&ids), "ids": ids },
        "pathways": Value::Object(pw),
        "total_pm25": field_summary(total_pm25),
        "deaths": { "krewski": field_summary(deaths_k), "lepeule": field_summary(deaths_l) },
        "timing": Value::Object(timing_map),
    });

    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("mkdir {dir:?}: {e}"))?;
    }
    let text = serde_json::to_string_pretty(&rec).map_err(|e| format!("serialize: {e}"))?;
    std::fs::write(path, text).map_err(|e| format!("write {path:?}: {e}"))?;
    println!(
        "wrote {}  (mode={mode}, |ppl|={}, sum deathsK={:?})",
        path.display(),
        ids.len(),
        compensated_sum(deaths_k)
    );
    Ok(())
}
