//! HARNESS INPUT PREPARATION — the Rust mirror of `run-model-py/inputs.py`.
//!
//! Everything here is *outside* the model: it produces exactly the arrays
//! `isrm_pushdown.esm` declares as `parameter`s, and nothing more. Reading bytes
//! and projecting lon/lat into the grid's plane is the sanctioned impure
//! boundary on the input side — the `.esm` declares `X`/`Y` as parameters, so
//! the projection is upstream of the spec, and all three runners do the
//! identical thing with the identical constants. Everything DOWNSTREAM of these
//! arrays is evaluated from the graph.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

use earthsciio::format::{ArrayData, Ff10Reader, Reader, Selection};
use earthsciio::{Cache, DataLoader, Provider};
use ndarray::{ArrayD, IxDyn};

use crate::paths;

// --------------------------------------------------------------------------- //
// Snyder 1987 spherical Lambert Conformal Conic forward.
//
// Matches lambert_conformal.esm and the validated run-model.jl implementation
// constant for constant. The runners MUST agree here or their X/Y differ and the
// derived support set diverges for a reason that is not the model.
// --------------------------------------------------------------------------- //
const LAT_1: f64 = 33.0;
const LAT_2: f64 = 45.0;
const LAT_0: f64 = 40.0;
const LON_0: f64 = -97.0;
const LCC_R: f64 = 6_370_997.0;
const D2R: f64 = 0.017_453_292_519_943_295;

fn lcc_t(lat: f64) -> f64 {
    (0.785_398_163_397_448_3 + lat * 0.008_726_646_259_971_648).tan()
}

fn lcc_n() -> f64 {
    ((LAT_1 * D2R).cos() / (LAT_2 * D2R).cos()).ln() / (lcc_t(LAT_2) / lcc_t(LAT_1)).ln()
}

fn lcc_f() -> f64 {
    let n = lcc_n();
    (LAT_1 * D2R).cos() * lcc_t(LAT_1).powf(n) / n
}

fn lcc_rho(lat: f64) -> f64 {
    LCC_R * lcc_f() / lcc_t(lat).powf(lcc_n())
}

/// Forward projection -> (X, Y) in metres.
pub fn lcc_forward(lon: f64, lat: f64) -> (f64, f64) {
    let theta = lcc_n() * (lon - LON_0) * D2R;
    let rho = lcc_rho(lat);
    (rho * theta.sin(), lcc_rho(LAT_0) - rho * theta.cos())
}

// --------------------------------------------------------------------------- //
// POLID -> pathway CODE (STEP 0.6). Verbatim from run-model.jl so all three
// runners recognise exactly the same 43,650 records. The .esm's is_VOC/is_NOx/…
// observeds select by code BAND, so the code — not the pathway name — is the
// model's input.
// --------------------------------------------------------------------------- //
const VOC_SET: &[&str] = &[
    "VOC", "VOC_INV", "XYL", "TOL", "TERP", "PAR", "OLE", "NVOL", "MEOH", "ISOP", "IOLE", "FORM",
    "ETOH", "ETHA", "ETH", "ALD2", "ALDX", "CB05_ALD2", "CB05_ALDX", "CB05_BENZENE", "CB05_ETH",
    "CB05_ETHA", "CB05_ETOH", "CB05_FORM", "CB05_IOLE", "CB05_ISOP", "CB05_MEOH", "CB05_OLE",
    "CB05_PAR", "CB05_TERP", "CB05_TOL", "CB05_XYL", "ETHANOL", "NHTOG", "NMOG",
];
const PM25_SET: &[&str] = &[
    "PM25-PRI", "PM2_5", "DIESEL-PM25", "PAL", "PCA", "PCL", "PEC", "PFE", "PK", "PMG", "PMN",
    "PMOTHR", "PNH4", "PNO3", "POC", "PSI", "PSO4", "PTI",
];
const NOX_SET: &[&str] = &["NOX", "HONO", "NO", "NO2"];
const NH3_SET: &[&str] = &["NH3"];
const SOX_SET: &[&str] = &["SO2"];

pub fn pathway_code(polid: &str) -> f64 {
    let p = polid.trim().to_ascii_uppercase();
    if VOC_SET.contains(&p.as_str()) {
        1.0
    } else if PM25_SET.contains(&p.as_str()) {
        42.0
    } else if NOX_SET.contains(&p.as_str()) {
        36.0
    } else if NH3_SET.contains(&p.as_str()) {
        40.0
    } else if SOX_SET.contains(&p.as_str()) {
        41.0
    } else {
        0.0
    }
}

/// The `.esm`'s declared parameter arrays, ready to bind as factor arrays.
pub struct Inputs {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub emis_annual: Vec<f64>,
    pub pollutant: Vec<f64>,
    pub w: Vec<f64>,
    pub s: Vec<f64>,
    pub e: Vec<f64>,
    pub n: Vec<f64>,
    pub total_pop: Vec<f64>,
    pub mortality_rate: Vec<f64>,
    pub n_rec: usize,
    pub n_src: usize,
}

impl Inputs {
    /// As a name -> dense 1-D array map, the shape both the value-invention
    /// engine and the standalone evaluator consume.
    pub fn as_arrays(&self) -> HashMap<String, ArrayD<f64>> {
        let mut m = HashMap::new();
        let mut put = |k: &str, v: &Vec<f64>| {
            m.insert(
                k.to_string(),
                ArrayD::from_shape_vec(IxDyn(&[v.len()]), v.clone())
                    .expect("1-D shape always matches its own length"),
            );
        };
        put("X", &self.x);
        put("Y", &self.y);
        put("emis_annual", &self.emis_annual);
        put("pollutant", &self.pollutant);
        put("W", &self.w);
        put("S", &self.s);
        put("E", &self.e);
        put("N", &self.n);
        put("TotalPop", &self.total_pop);
        put("MortalityRate", &self.mortality_rate);
        m
    }
}

fn as_f64(d: &ArrayData) -> Vec<f64> {
    match d {
        ArrayData::F64(v) => v.clone(),
        ArrayData::I64(v) => v.iter().map(|&x| x as f64).collect(),
        ArrayData::I32(v) => v.iter().map(|&x| x as f64).collect(),
        ArrayData::Bool(v) => v.iter().map(|&x| if x { 1.0 } else { 0.0 }).collect(),
        ArrayData::Str(_) => Vec::new(),
    }
}

fn as_strings(d: &ArrayData) -> Vec<String> {
    match d {
        ArrayData::Str(v) => v.clone(),
        _ => Vec::new(),
    }
}

/// The EGU point CSV members of the 2016fd inputs zip — same selection rule as
/// run-model.jl (name contains "egu", ends in ".csv").
fn egu_members(zip_path: &Path) -> Result<Vec<String>, String> {
    let file = std::fs::File::open(zip_path).map_err(|e| format!("open {zip_path:?}: {e}"))?;
    let mut zip = zip::ZipArchive::new(file).map_err(|e| format!("read zip: {e}"))?;
    let mut out = Vec::new();
    for i in 0..zip.len() {
        let f = zip.by_index(i).map_err(|e| format!("zip entry {i}: {e}"))?;
        let name = f.name().to_string();
        let lower = name.to_ascii_lowercase();
        if lower.contains("egu") && lower.ends_with(".csv") {
            out.push(name);
        }
    }
    out.sort();
    Ok(out)
}

/// One EGU zip member -> a decoded dataset, with the column-name header row
/// stripped first.
///
/// These 2016fd members carry a `country_cd,region_cd,…` header row that is NOT
/// a `#` comment, and neither EarthSciIO's Rust nor its Python nor its Julia
/// FF10 reader skips one (the fixed 77-column schema comes from the spec, not
/// from the file). All three runners work around it the same way — filter the
/// header line, then hand the reader a clean CSV — so they see identical
/// records. Kept in the runner rather than patched into the reader so the
/// bindings' input paths stay equivalent.
fn read_ff10_member(
    zip_path: &Path,
    member: &str,
    variables: &[String],
    workdir: &Path,
) -> Result<HashMap<String, ArrayData>, String> {
    std::fs::create_dir_all(workdir).map_err(|e| format!("mkdir {workdir:?}: {e}"))?;
    let file = std::fs::File::open(zip_path).map_err(|e| format!("open {zip_path:?}: {e}"))?;
    let mut zip = zip::ZipArchive::new(file).map_err(|e| format!("read zip: {e}"))?;
    let mut text = String::new();
    zip.by_name(member)
        .map_err(|e| format!("zip member {member}: {e}"))?
        .read_to_string(&mut text)
        .map_err(|e| format!("read member {member}: {e}"))?;

    let cleaned: String = text
        .lines()
        .filter(|ln| !ln.trim_start().to_ascii_lowercase().starts_with("country_cd"))
        .collect::<Vec<_>>()
        .join("\n");

    let base = Path::new(member).file_name().unwrap_or_default();
    let tmp = workdir.join(base);
    std::fs::write(&tmp, cleaned).map_err(|e| format!("write {tmp:?}: {e}"))?;

    let reader = Ff10Reader::new();
    let ds = reader
        .read_native(&tmp, variables, &Selection::All)
        .map_err(|e| format!("FF10 decode {member}: {e}"))?;
    let _ = std::fs::remove_file(&tmp);

    Ok(ds
        .variables
        .into_iter()
        .map(|(k, f)| (k, f.data))
        .collect())
}

/// Read every EGU FF10 point record, keeping only the RECOGNISED-pathway, finite
/// ones — the same filter run-model.jl applies, so N_REC matches (43,650).
pub fn read_emissions() -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String> {
    let zip_path = paths::egu_zip();
    let workdir = paths::scratch().join("rs_ff10_tmp");
    let vars: Vec<String> = ["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"]
        .iter()
        .map(|s| s.to_string())
        .collect();

    let members = egu_members(&zip_path)?;
    if members.is_empty() {
        return Err(format!("no EGU csv members found in {zip_path:?}"));
    }

    let (mut lon, mut lat, mut ann, mut code) = (vec![], vec![], vec![], vec![]);
    for m in &members {
        let ds = read_ff10_member(&zip_path, m, &vars, &workdir)?;
        let polid = as_strings(ds.get("POLID").ok_or("FF10: no POLID column")?);
        let a = as_f64(ds.get("ANN_VALUE").ok_or("FF10: no ANN_VALUE column")?);
        let lo = as_f64(ds.get("LONGITUDE").ok_or("FF10: no LONGITUDE column")?);
        let la = as_f64(ds.get("LATITUDE").ok_or("FF10: no LATITUDE column")?);
        for i in 0..polid.len() {
            let c = pathway_code(&polid[i]);
            // Drop unrecognised pathways and any non-finite coordinate/value.
            if c > 0.0 && lo[i].is_finite() && la[i].is_finite() && a[i].is_finite() {
                lon.push(lo[i]);
                lat.push(la[i]);
                ann.push(a[i]);
                code.push(c);
            }
        }
    }
    Ok((lon, lat, ann, code))
}

/// W/S/E/N + TotalPop/MortalityRate for the first N_SRC cells of the ISRM zarr —
/// the container grid, exactly the prefix run-model.jl bins against.
pub fn read_grid() -> Result<HashMap<String, Vec<f64>>, String> {
    let root = paths::scratch().join("rs_cache_meta");
    let cache = Arc::new(
        Cache::builder()
            .data_dir(&root)
            .offline(false)
            .build()
            .map_err(|e| format!("cache: {e}"))?,
    );
    let n_src = paths::n_src();
    let mut out = HashMap::new();
    for group in [
        vec!["W", "S", "E", "N"],
        vec!["TotalPop", "MortalityRate"],
    ] {
        let loader = DataLoader::new("ISRM_grid", "zarr", paths::zarr_url())
            .variables(group.iter().map(|s| s.to_string()));
        let mut p = Provider::new(loader, cache.clone(), None)
            .map_err(|e| format!("grid provider: {e}"))?;
        let fields = p.materialize().map_err(|e| format!("grid read: {e}"))?;
        for v in &group {
            let f = fields
                .get(*v)
                .ok_or_else(|| format!("zarr: variable {v} missing"))?;
            let mut vals = as_f64(&f.data);
            vals.truncate(n_src);
            out.insert(v.to_string(), vals);
        }
    }
    Ok(out)
}

/// Assemble every declared parameter array. `firstn` truncates the emission
/// record list for a reduced, fast run (mirrors L3_FIRSTN / ISRM_FIRSTN).
pub fn build_inputs(firstn: Option<usize>) -> Result<Inputs, String> {
    let (mut lon, mut lat, mut ann, mut code) = read_emissions()?;
    if let Some(n) = firstn {
        let n = n.min(lon.len());
        lon.truncate(n);
        lat.truncate(n);
        ann.truncate(n);
        code.truncate(n);
    }
    let mut x = Vec::with_capacity(lon.len());
    let mut y = Vec::with_capacity(lon.len());
    for i in 0..lon.len() {
        let (xi, yi) = lcc_forward(lon[i], lat[i]);
        x.push(xi);
        y.push(yi);
    }
    let g = read_grid()?;
    let take = |k: &str| -> Result<Vec<f64>, String> {
        g.get(k).cloned().ok_or_else(|| format!("grid: {k} missing"))
    };
    Ok(Inputs {
        n_rec: x.len(),
        n_src: paths::n_src(),
        x,
        y,
        emis_annual: ann,
        pollutant: code,
        w: take("W")?,
        s: take("S")?,
        e: take("E")?,
        n: take("N")?,
        total_pop: take("TotalPop")?,
        mortality_rate: take("MortalityRate")?,
    })
}
