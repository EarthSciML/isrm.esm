#!/usr/bin/env julia
# =============================================================================
# run.jl — the JULIA binding drives the single clean `isrm.esm` end to end
# through the PUBLIC EarthSciAST surface.
#
#   * `prepare(file; providers, const_arrays, pushdown_rewrite=true)` — the
#     automatic projection-pushdown rewrite runs inside the engine; the SR
#     provider gates are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * the SR / grid / pop / mortality providers come FROM THE DOCUMENT
#     (`providers_from_document`: format = `metadata.esio_format`, URL =
#     `source.url_template`);
#   * the EGU FF10 records are read through EarthSciIO's OWN ff10 reader with
#     the document-declared container options (`metadata.x_esd`: the zip
#     `member_filter` glob and `skip_header_row`) — no hand unzip, no hand
#     header-strip; the POLID→code map is the document's `pollutant_codes`;
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
# EGU FF10 ingest — through EarthSciIO's ff10 reader, options FROM THE DOCUMENT.
#
# The document's EGU_Emis loader declares the zip container (`x_esd.container`),
# the member glob (`x_esd.member_filter`, "*egu*"), the EPA column-header quirk
# (`x_esd.skip_header_row`), and the POLID→pathway-code integer map
# (`x_esd.pollutant_codes`) the is_* observeds compare against. This function
# contributes input PLUMBING only: read the table, map the codes, drop the
# unrecognised/non-finite records (per the document's note), truncate for a
# reduced run. The reader concatenates members in sorted (lexicographic) order —
# identical to the zip order for this archive, so ISRM_FIRSTN selects the same
# records as every previous runner. Everything DOWNSTREAM of these raw arrays —
# projection, spatial join, support-set invention, binning, contraction,
# deaths — is the engine evaluating the document.
# =============================================================================

"""The EGU zip to read: `\$EGU_ZIP` (default `<repo>/data/…`) when present,
otherwise the document's `source.url_template` fetched once through the
EarthSciIO cache (content-addressed under `cache_root`)."""
function resolve_egu_zip(egu_meta, cache_root)
    isfile(EGU_ZIP) && return EGU_ZIP
    url = String(egu_meta["source"]["url_template"])
    println("  EGU zip not found at $EGU_ZIP — fetching $url via the EarthSciIO cache ...")
    flush(stdout)
    entry = EarthSciIO.fetch_blob(EarthSciIO.Cache(; root=joinpath(cache_root, "EGU_Emis")),
                                  url; source_loader="EGU_Emis")
    return entry.path
end

"""Raw EGU emission-record arrays (lon, lat, annual, pathway code), recognised
+ finite records only — the filter the document's note declares."""
function read_egu(egu_meta, zippath; firstn=nothing)
    xe = egu_meta["metadata"]["x_esd"]
    nds = EarthSciIO.read_native(EarthSciIO.FF10Reader(), zippath;
        member_glob = String(xe["member_filter"]),
        skip_header_row = Bool(get(xe, "skip_header_row", false)),
        variables = ["POLID", "ANN_VALUE", "LONGITUDE", "LATITUDE"])
    codes = Dict{String,Float64}(uppercase(String(k)) => Float64(v)
                                 for (k, v) in xe["pollutant_codes"])
    code = Float64[get(codes, uppercase(strip(String(p))), 0.0) for p in nds["POLID"].data]
    ann  = Float64.(nds["ANN_VALUE"].data)
    lon  = Float64.(nds["LONGITUDE"].data)
    lat  = Float64.(nds["LATITUDE"].data)
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

    doc_raw = JSON.parsefile(MODEL)
    cache_root = joinpath(SCRATCH, "run-jl-esio-cache")

    # ---- inputs: EGU records via the document-declared FF10 reader options --
    println("reading EGU emissions (EarthSciIO ff10: member_glob + skip_header_row) ...")
    flush(stdout)
    egu_meta = doc_raw["data_loaders"]["EGU_Emis"]
    lon, lat, ann, code = read_egu(egu_meta, resolve_egu_zip(egu_meta, cache_root);
                                   firstn=firstn)
    N_REC = length(lon)
    println("  N_REC = $N_REC"); flush(stdout)

    # ---- providers FROM THE DOCUMENT ----------------------------------------
    sr_arrays = String[a for a in keys(doc_raw["data_loaders"]["ISRM_SR"]["metadata"]["x_esd"]["arrays"])
                       if !(a in ("TotalPop", "MortalityRate", "W", "S", "E", "N"))]
    seed_empty_zattrs(joinpath(cache_root, "ISRM_SR"),
                      doc_raw["data_loaders"]["ISRM_SR"]["source"]["url_template"], sr_arrays)
    println("building providers from the document (ISRM_SR) ...")
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

    # ---- const arrays: the EGU loader's variables, src rects ----------------
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
        println("FULL SCALE: ", (okK && okL) ? "PASS" : "FAIL")
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
