//! Machine-specific paths — the Rust mirror of `run-model-jl-pushdown/paths.jl`
//! and `run-model-py/paths.py`. Every location is overridable by environment
//! variable so the runner is portable; nothing is hardcoded to one machine.
//!
//! IMPORTANT on this cluster: the scratch directory must be DISK-backed. The
//! root filesystem here is `tmpfs`, so SR chunk blobs written under `/tmp`
//! consume the same cgroup memory budget the model needs and OOM-kill the job.

use std::path::{Path, PathBuf};

fn env_or(key: &str, fallback: impl Into<String>) -> String {
    std::env::var(key).unwrap_or_else(|_| fallback.into())
}

/// This crate's directory (…/isrm.esm/run-model-rs).
pub fn rs_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The isrm.esm repository root.
pub fn repo() -> PathBuf {
    rs_dir().parent().map(Path::to_path_buf).unwrap_or_default()
}

/// The model this runner drives.
///
/// `isrm_pushdown.esm` is the only variant whose emission-bearing source-cell
/// set (`ppl`) is DERIVED by the graph. `isrm.esm` takes it as a
/// `src_cell_of_ppl` parameter, which would put the spatial join outside the
/// spec — and `ppl` is one of the numbers under test (compared EXACTLY across
/// bindings, CONFORMANCE_SPEC §5.5).
pub fn model() -> PathBuf {
    match std::env::var("ISRM_MODEL") {
        Ok(p) => PathBuf::from(p),
        Err(_) => repo().join("isrm_pushdown.esm"),
    }
}

/// The Julia oracle runner — the source of the shared FF10 input zip.
pub fn runmodel() -> PathBuf {
    match std::env::var("ISRM_RUNMODEL") {
        Ok(p) => PathBuf::from(p),
        Err(_) => repo().join("run-model-jl"),
    }
}

pub fn egu_zip() -> PathBuf {
    match std::env::var("EGU_ZIP") {
        Ok(p) => PathBuf::from(p),
        Err(_) => runmodel().join("data").join("2016fd_inputs_point.zip"),
    }
}

pub fn zarr_url() -> String {
    env_or("ISRM_ZARR_URL", "s3://inmap-model/isrm_v1.2.1.zarr/")
}

/// SR source cells on the SR axis (== receptor cells == population cells).
pub fn n_src() -> usize {
    env_or("ISRM_N_SRC", "52411").parse().unwrap_or(52411)
}

/// A DISK-backed scratch root. Mirrors the Julia/Python resolution order.
pub fn scratch() -> PathBuf {
    if let Ok(p) = std::env::var("ISRM_SCRATCH") {
        return PathBuf::from(p);
    }
    let user = std::env::var("USER")
        .or_else(|_| std::env::var("LOGNAME"))
        .unwrap_or_else(|_| "user".to_string());
    let local = PathBuf::from("/scratch.local").join(&user);
    if local.is_dir() {
        return local.join("isrm-esm");
    }
    if let Ok(s) = std::env::var("SCRATCH") {
        let p = PathBuf::from(s);
        if p.is_dir() {
            return p.join("isrm-esm");
        }
    }
    std::env::temp_dir().join("isrm-esm")
}
