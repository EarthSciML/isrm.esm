//! Machine-independent path resolution for the Phase-4 Rust shim — the Rust
//! mirror of `run-py/paths.py` / `run-jl/paths.jl`. Everything derives from
//! this crate's own location (hence the isrm.esm checkout) plus environment
//! overrides, so the shim is portable across machines.
//!
//!   ISRM_MODEL       the .esm to drive     (default: <repo>/isrm.esm)
//!   ISRM_SCRATCH     bulk scratch root     (default: first writable of
//!                    /scratch.local/$USER/isrm-esm, $SCRATCH/isrm-esm, tempdir)
//!   EGU_ZIP          FF10 point-source zip (default: <repo>/run-model-jl/data/
//!                    2016fd_inputs_point.zip)
//!   ISRM_ESIO_CACHE  EarthSciIO cache root (default: $ISRM_SCRATCH/
//!                    run-rs-esio-cache; point it at an existing cache — e.g.
//!                    run-jl's — to reuse already-fetched SR chunk blobs, the
//!                    cache format is cross-language)
//!
//! SR chunk blobs land under the cache root, so it wants GBs of fast local
//! disk — NOT a network filesystem and NEVER a tmpfs /tmp (it eats the memory
//! cgroup). The scratch resolution below prefers /scratch.local/$USER for
//! exactly that reason; tempdir is a last resort.

use std::path::{Path, PathBuf};

/// This crate's directory (…/isrm.esm/run-rs).
pub fn rs_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The isrm.esm repository root.
pub fn repo() -> PathBuf {
    rs_dir().parent().map(Path::to_path_buf).unwrap_or_default()
}

/// The model this runner drives — the single clean `isrm.esm`.
pub fn model() -> PathBuf {
    match std::env::var("ISRM_MODEL") {
        Ok(p) => PathBuf::from(p),
        Err(_) => repo().join("isrm.esm"),
    }
}

pub fn egu_zip() -> PathBuf {
    match std::env::var("EGU_ZIP") {
        Ok(p) => PathBuf::from(p),
        Err(_) => repo()
            .join("run-model-jl")
            .join("data")
            .join("2016fd_inputs_point.zip"),
    }
}

/// A DISK-backed scratch root (mirrors the Julia/Python resolution order).
pub fn scratch() -> PathBuf {
    if let Ok(p) = std::env::var("ISRM_SCRATCH") {
        let p = PathBuf::from(p);
        let _ = std::fs::create_dir_all(&p);
        return p;
    }
    let user = std::env::var("USER")
        .or_else(|_| std::env::var("LOGNAME"))
        .unwrap_or_else(|_| "user".to_string());
    for base in [
        Some(PathBuf::from("/scratch.local").join(&user)),
        std::env::var("SCRATCH").ok().map(PathBuf::from),
    ]
    .into_iter()
    .flatten()
    {
        if base.is_dir() {
            let c = base.join("isrm-esm");
            if std::fs::create_dir_all(&c).is_ok() {
                return c;
            }
        }
    }
    let c = std::env::temp_dir().join("isrm-esm");
    let _ = std::fs::create_dir_all(&c);
    c
}

/// The EarthSciIO cache root for the document-declared providers.
pub fn esio_cache() -> PathBuf {
    match std::env::var("ISRM_ESIO_CACHE") {
        Ok(p) => PathBuf::from(p),
        Err(_) => scratch().join("run-rs-esio-cache"),
    }
}
