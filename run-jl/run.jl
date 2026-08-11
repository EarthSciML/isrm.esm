#!/usr/bin/env julia
# =============================================================================
# run.jl — Phase 1 (clean consolidation): the Julia binding drives the single
# clean `isrm.esm` end to end through the PUBLIC EarthSciAST surface.
#
#   * `prepare(file; providers, const_arrays, pushdown_rewrite=true)` — the
#     automatic projection-pushdown rewrite runs inside the engine; the SR
#     provider gates are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * the SR / grid / pop / mortality providers come FROM THE DOCUMENT
#     (`providers_from_document`: format = `metadata.esio_format`, URL =
#     `source.url_template`);
#   * every reported number is the binding's evaluation of the document's
#     observed graph (`observed_field`) — NO hand-written STEP-0 math here.
#
#   FULL run  (default)          → assert sum(deathsK/L) ≈ 7524.92 / 16979.63
#   REDUCED   (ISRM_FIRSTN=n)    → first n emission records, totals reported
#
# Emits the cross-language contract record (contract/results_schema.json) with
# model="isrm.esm", mode="runtime_observed_graph", binding="julia".
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc                     # activates EarthSciIO's zarr codec extension
using ZipFile
import GeometryOps, GeoInterface # STRtree broad-phase fast path for the join gate
import JSON
const EA = EarthSciAST

include(joinpath(@__DIR__, "paths.jl"))
# contract/results.jl must be included at TOP LEVEL (world age).
include(joinpath(@__DIR__, "..", "contract", "results.jl"))

const T0 = time()
const N_SRC = 52411             # metaparameter default; SR source/receptor cells
const ORACLE_K = 7524.918845602511
const ORACLE_L = 16979.632171487083
peak_rss_bytes() = parse(Int, split(read("/proc/self/statm", String))[2]) * 4096

# =============================================================================
# EGU FF10 ingest — TEMPORARY (Phase 2 closes this).
#
# LOUD CAVEAT: the document's EGU_Emis data_loader declares the zip URL and the
# ff10 esio_format, but the format stack cannot yet express (a) zip MEMBER
# selection (`*egu*.csv` inside 2016fd_inputs_point.zip) or (b) the
# `country_cd,...` column-header row these members carry, which the FF10
# reader's fixed 77-column schema cannot skip. Until Phase 2 lands those two
# capabilities, THIS BLOCK stands in for the EGU_Emis provider: it reads the
# zip at $EGU_ZIP, strips the header row, and feeds each member through
# EarthSciIO's own FF10Reader (`read_native`). The POLID→pathway-code map
# mirrors the document's G3 note (`is_*` observeds select by code band) and is
# copied verbatim from the validated run-model.jl so record recognition is
# identical. Everything DOWNSTREAM of these raw arrays — projection, spatial
# join, support-set invention, binning, contraction, deaths — is the engine
# evaluating the document.
# =============================================================================
const VOC_SET  = Set(["VOC","VOC_INV","XYL","TOL","TERP","PAR","OLE","NVOL","MEOH","ISOP","IOLE","FORM","ETOH","ETHA","ETH","ALD2","ALDX","CB05_ALD2","CB05_ALDX","CB05_BENZENE","CB05_ETH","CB05_ETHA","CB05_ETOH","CB05_FORM","CB05_IOLE","CB05_ISOP","CB05_MEOH","CB05_OLE","CB05_PAR","CB05_TERP","CB05_TOL","CB05_XYL","ETHANOL","NHTOG","NMOG"])
const PM25_SET = Set(["PM25-PRI","PM2_5","DIESEL-PM25","PAL","PCA","PCL","PEC","PFE","PK","PMG","PMN","PMOTHR","PNH4","PNO3","POC","PSI","PSO4","PTI"])
const NOX_SET  = Set(["NOX","HONO","NO","NO2"])
const NH3_SET  = Set(["NH3"])
const SOX_SET  = Set(["SO2"])
# pathway CODE bands, per the document's is_* observeds (G3): SOA 1-35,
# pNO3 36-39, pNH4 40, pSO4 41, PrimaryPM25 42-59; 0 = unrecognised (dropped).
pathcode(p) = (u = uppercase(strip(p));
    u in VOC_SET ? 1.0 : u in PM25_SET ? 42.0 : u in NOX_SET ? 36.0 :
    u in NH3_SET ? 40.0 : u in SOX_SET ? 41.0 : 0.0)

# One zip member → NativeDataset via EarthSciIO's FF10 reader, header stripped.
# (Minimal logic from the validated run-model.jl `read_ff10_egu`.)
function read_ff10_egu(zippath, member, vars)
    text = ""
    let r = ZipFile.Reader(zippath)
        for f in r.files; f.name == member && (text = String(read(f)); break); end
        close(r)
    end
    lines = split(text, '\n')
    kept = filter(ln -> !startswith(lowercase(strip(ln)), "country_cd"), lines)
    tmp = joinpath(mktempdir(SCRATCH), basename(member))
    write(tmp, join(kept, '\n'))
    nds = EarthSciIO.read_native(EarthSciIO.FF10Reader(), tmp; variables=vars)
    rm(dirname(tmp); recursive=true, force=true)
    return nds
end

"""Raw EGU emission-record arrays (lon, lat, annual, pathway code), recognised
+ finite records only — the same filter the validated runners apply. Members
are iterated in ZIP order (matching run-model.jl / the pushdown-era reduced
records, so ISRM_FIRSTN truncation selects the same records)."""
function read_egu(zippath; firstn=nothing)
    isfile(zippath) || error("EGU zip not found at $zippath — set EGU_ZIP")
    members = String[]
    let r = ZipFile.Reader(zippath)
        for f in r.files
            occursin("egu", lowercase(f.name)) && endswith(lowercase(f.name), ".csv") &&
                push!(members, f.name)
        end
        close(r)
    end
    isempty(members) && error("no *egu*.csv members in $zippath")
    println("  EGU members: ", members)
    lon = Float64[]; lat = Float64[]; ann = Float64[]; code = Float64[]
    for m in members
        nds = read_ff10_egu(zippath, m, ["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"])
        append!(code, Float64[pathcode(String(p)) for p in nds["POLID"].data])
        append!(ann,  Float64.(nds["ANN_VALUE"].data))
        append!(lon,  Float64.(nds["LONGITUDE"].data))
        append!(lat,  Float64.(nds["LATITUDE"].data))
    end
    n_all = length(code)
    keep = [code[i] > 0.0 && isfinite(lon[i]) && isfinite(lat[i]) && isfinite(ann[i])
            for i in 1:n_all]
    lon = lon[keep]; lat = lat[keep]; ann = ann[keep]; code = code[keep]
    println("  EGU FF10 records: $n_all read, $(length(lon)) recognised-pathway kept")
    if firstn !== nothing
        n = min(firstn, length(lon))
        lon = lon[1:n]; lat = lat[1:n]; ann = ann[1:n]; code = code[1:n]
    end
    return lon, lat, ann, code
end

# zarr workaround (unchanged from the validated runners): the SR arrays carry NO
# `.zattrs` object (live GET → 404, which `_fetch_bytes_optional` does not
# swallow), so seed an empty `{}` as a cache HIT before any SR fetch.
function seed_empty_zattrs(cache_root, base, arrays)
    cache = EarthSciIO.Cache(; root=cache_root)
    base = rstrip(base, '/')
    for arr in arrays
        url = "$base/$arr/.zattrs"
        key = EarthSciIO.cache_key(url)
        if EarthSciIO.get_blob(cache.store, key) === nothing
            tmp = EarthSciIO.staging_path(cache.store)
            write(tmp, "{}")
            EarthSciIO.put_blob!(cache.store, key, tmp)
        end
    end
end

# =============================================================================
function main()
    firstn = haskey(ENV, "ISRM_FIRSTN") ? parse(Int, ENV["ISRM_FIRSTN"]) : nothing
    reduced = firstn !== nothing
    println(reduced ? "REDUCED run — first $firstn emission records" :
            "FULL run — whole domain (target deathsK≈$(round(ORACLE_K, digits=2)), deathsL≈$(round(ORACLE_L, digits=2)))")
    println("model:   $MODEL")
    println("scratch: $SCRATCH")

    # ---- inputs: EGU records (TEMPORARY shim ingest, see above) -------------
    println("reading EGU emissions ..."); flush(stdout)
    lon, lat, ann, code = read_egu(EGU_ZIP; firstn=firstn)
    N_REC = length(lon)
    println("  N_REC = $N_REC"); flush(stdout)

    # ---- providers FROM THE DOCUMENT ----------------------------------------
    doc_raw = JSON.parsefile(MODEL)
    cache_root = joinpath(SCRATCH, "run-jl-esio-cache")
    sr_arrays = String[a for a in keys(doc_raw["data_loaders"]["ISRM_SR"]["metadata"]["x_esd"]["arrays"])
                       if !(a in ("TotalPop", "MortalityRate", "W", "S", "E", "N"))]
    seed_empty_zattrs(joinpath(cache_root, "ISRM_SR"),
                      doc_raw["data_loaders"]["ISRM_SR"]["source"]["url_template"], sr_arrays)
    println("building providers from the document (ISRM_SR; EGU_Emis is the Phase-2 gap) ...")
    providers = EA.providers_from_document(doc_raw; cache_root=cache_root,
                                           loaders=["ISRM_SR"])
    println("  providers: ", sort(collect(keys(providers)))); flush(stdout)

    # ---- src-cell rectangles: the [1:N_SRC] prefix of the full grid ---------
    # The document declares src_W/src_S/src_E/src_N ([src_cells]) alongside the
    # full-grid W/S/E/N ([pop_cells]); the loader note pins the prefix
    # relationship ("the first 52411 build the cell rectangles"). Slicing is
    # input plumbing, not model math.
    println("materializing grid prefix (src_W/S/E/N) ..."); flush(stdout)
    t_grid = @elapsed src_rects = Dict{String,Any}(
        "src_$k" => Float64.(EA.provider_sample(providers["ISRM_SR.$k"], 0.0)[1:N_SRC])
        for k in ("W", "S", "E", "N"))
    println("  grid prefix in $(round(t_grid, digits=1)) s"); flush(stdout)

    # ---- const arrays: the EGU loader's variables (temporary), src rects ----
    ca = Dict{String,Any}(src_rects)
    ca["EGU_Emis.lon"] = lon;  ca["EGU_Emis.lat"] = lat
    ca["EGU_Emis.annual"] = ann; ca["EGU_Emis.pollutant"] = code
    for v in ("stkhgt", "stkdiam", "stktemp", "stkvel")
        ca["EGU_Emis.$v"] = zeros(N_REC)     # declared, unused downstream
    end

    # ---- load + PREPARE (rewrite + record-derived gating inside the engine) --
    println("loading document (N_REC=$N_REC) ..."); flush(stdout)
    file = EA.load(MODEL; metaparameters=Dict("N_REC" => N_REC))
    insp = EA.BuildInspection()
    println("prepare(pushdown_rewrite=true) — rewrite → VI → gated SR fetch → build ...")
    flush(stdout)
    t_prep = @elapsed prep = EA.prepare(file; providers=providers, const_arrays=ca,
                                        inspect=insp, pushdown_rewrite=true)
    println("PREPARE done in $(round(t_prep, digits=1)) s  (peak RSS so far: ",
            round(peak_rss_bytes() / 2^30, digits=2), " GiB)"); flush(stdout)

    # ---- the engine-derived support set (for the contract record) -----------
    mf_keys = [k for k in keys(insp.const_arrays) if startswith(String(k), "pd_member_factor__")]
    isempty(mf_keys) && error("no pd_member_factor__* const array — did the rewrite fire?")
    members = sort!(Int.(insp.const_arrays[first(mf_keys)]))
    n_ppl = length(members)
    println("engine-derived support set: |members| = $n_ppl of $N_SRC source cells")
    !reduced && n_ppl != 1520 &&
        println("  WARNING: expected 1520 emission-bearing cells at full scale")

    # ---- evaluate the observed graph ----------------------------------------
    function rt(v)
        print("  evaluating observed $v ... "); flush(stdout)
        local fld
        t = @elapsed fld = EA.observed_field(prep, insp, v)
        println(round(t, digits=1), " s"); flush(stdout)
        return fld
    end
    t_eval = @elapsed begin
        dK = rt("deathsK"); dL = rt("deathsL"); tp = rt("TotalPM25")
        # per-pathway intermediates through the SAME runtime path, so a
        # disagreement localizes to one pathway instead of only the totals.
        PW_OBS = ["SOA"         => ("E_VOC",  "conc_SOA"),
                  "pNO3"        => ("E_NOx",  "conc_pNO3"),
                  "pNH4"        => ("E_NH3",  "conc_pNH4"),
                  "pSO4"        => ("E_SOx",  "conc_pSO4"),
                  "PrimaryPM25" => ("E_PM25", "conc_PrimaryPM25")]
        pathways = Dict{String,Any}()
        for (arr, (evar, cvar)) in PW_OBS
            Ep = rt(evar); cp = rt(cvar)
            pathways[arr] = (emis_sum = sum(Ep), conc_sum = sum(cp), conc_max = maximum(cp))
        end
    end
    println("EVAL done in $(round(t_eval, digits=1)) s")

    sK = sum(dK); sL = sum(dL)
    println("\n", "="^70)
    println("  sum(deathsK) = ", sK)
    println("  sum(deathsL) = ", sL)
    println("  Σ TotalPM25  = ", sum(tp))
    okK = isapprox(sK, ORACLE_K; rtol=1e-4)
    okL = isapprox(sL, ORACLE_L; rtol=1e-4)
    if !reduced
        println("  target deathsK=$ORACLE_K  rel.err ",
                round(100 * (sK - ORACLE_K) / ORACLE_K, digits=6), "%")
        println("  target deathsL=$ORACLE_L rel.err ",
                round(100 * (sL - ORACLE_L) / ORACLE_L, digits=6), "%")
        println("PHASE 1 FULL: ", (okK && okL) ? "PASS" : "FAIL")
    end
    println("="^70)

    # ---- contract record ----------------------------------------------------
    out = joinpath(@__DIR__, reduced ? "results_reduced.json" : "results.json")
    write_results(out;
        binding_version = "julia $(VERSION) / EarthSciAST $(pkgversion(EarthSciAST))",
        model = MODEL, mode = "runtime_observed_graph",
        n_src = N_SRC, n_rcv = N_SRC, n_rec = N_REC,
        ppl = members,
        pathways = pathways,
        total_pm25 = tp, deathsK = dK, deathsL = dL,
        timing = Dict("wall_seconds" => time() - T0,
                      "build_seconds" => t_prep,
                      "eval_seconds" => t_eval,
                      "peak_rss_bytes" => peak_rss_bytes()))
    if !reduced
        (okK && okL) || error("full-scale totals off oracle (see above)")
    end
    return 0
end

exit(main())
