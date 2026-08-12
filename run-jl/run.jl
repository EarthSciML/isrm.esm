#!/usr/bin/env julia
# =============================================================================
# run.jl — the JULIA binding drives the single clean `isrm.esm` end to end
# through the PUBLIC EarthSciAST surface. NOTHING MODEL-SHAPED LIVES HERE: this
# file names no pollutant, no column, no grid extent and no record count.
#
#   * `prepare(doc; providers, pushdown_rewrite=true)` — the automatic
#     projection-pushdown rewrite runs inside the engine; the SR provider gates
#     are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * EVERY provider comes FROM THE DOCUMENT (`providers_from_document`:
#     format = `metadata.esio_format`, URL = `source.url_template`) — the SR
#     slabs, the grid, the population, AND the EGU FF10 table, whose ingest the
#     loader now declares in full (esm-spec §8.9): `reader_options` (the zip
#     member glob + header row), `codes` (POLID text -> the pathway enum, an
#     unrecognised code dropping the record), `record_filter` (no coordinate /
#     no annual total is not a record) and `extent` (the surviving count binds
#     N_REC). The src-cell rectangles are the `select` range `W[0:N_SRC]` on
#     their own loader variables;
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
const ORACLE_K = 7524.918845602511
const ORACLE_L = 16979.632171487083
# Peak resident set. `/proc/self/statm` field 2 is the CURRENT resident page
# count and is Linux-only; `Sys.maxrss()` is the high-water mark and is portable,
# so it is the fallback wherever /proc does not exist (macOS).
peak_rss_bytes() = isfile("/proc/self/statm") ?
    parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 : Int(Sys.maxrss())

"""A metaparameter's declared default, read from the document (so no grid
extent is written down here)."""
metaparam(doc, name) = Int(get(get(get(doc, "metaparameters", Dict()), name, Dict()),
                               "default", 0))

"""The loaders that DISCOVER their own extent (`extent.metaparameter`) — the
record-bearing tables of the document, whatever they happen to be called. The
two knobs below are scale/locality concerns of a RUN, not of the model, and both
are expressed in the document's own vocabulary."""
record_loaders(doc) = String[String(name) for (name, ld) in get(doc, "data_loaders", Dict())
                             if ld isa AbstractDict &&
                                get(ld, "extent", Dict()) isa AbstractDict &&
                                haskey(get(ld, "extent", Dict()), "metaparameter")]

"""The record count the loaders DISCOVERED — the delivered length of any one of
an extent-declaring loader's variables. They are aligned by construction (that
is exactly what `record_filter` guarantees), so the first one answers for all,
and this file still names no column of any particular model."""
function discovered_records(doc, insp)
    for name in record_loaders(doc), v in keys(doc["data_loaders"][name]["variables"])
        a = get(insp.const_arrays, "$name.$v", nothing)
        a === nothing || return length(a)
    end
    error("no record-discovering loader delivered an array to size N_REC from")
end

"""REDUCED runs: truncate every record-discovering loader to its first `n`
DELIVERED records with a loader-level `select` range (esm-spec §8.9.2). Because
the selection follows the loader's own `record_filter`, this picks the same
records the previous runners' post-filter truncation did — and `extent` then
re-discovers the smaller N_REC by itself."""
function truncate_records!(doc, n)
    for name in record_loaders(doc)
        doc["data_loaders"][name]["select"] =
            Dict("axes" => [Dict("range" => Dict("start" => 0, "stop" => n))])
    end
    return doc
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
    reduced && truncate_records!(doc_raw, firstn)

    # ---- providers FROM THE DOCUMENT — ALL of them --------------------------
    # Including the FF10 table: the loader declares its own reader options, code
    # map, record filter and extent, so there is nothing left here to read, map,
    # filter or count.
    sr_arrays = String[a for a in keys(doc_raw["data_loaders"]["ISRM_SR"]["metadata"]["x_esd"]["arrays"])
                       if !(a in ("TotalPop", "MortalityRate", "W", "S", "E", "N"))]
    seed_empty_zattrs(joinpath(cache_root, "ISRM_SR"),
                      doc_raw["data_loaders"]["ISRM_SR"]["source"]["url_template"], sr_arrays)
    println("building providers from the document ...")
    # A local copy of a record loader's source is a LOCALITY choice of this run
    # (gaftp.epa.gov is slow and flaky), so it is a url_override rather than an
    # edit to the document.
    url_overrides = Dict{String,String}()
    if isfile(EGU_ZIP)
        for name in record_loaders(doc_raw)
            url_overrides[name] = "file://" * EGU_ZIP
        end
        println("  record source mirrored from $EGU_ZIP")
    end
    t_providers = @elapsed providers =
        EA.providers_from_document(doc_raw; cache_root=cache_root,
                                   url_overrides=url_overrides)
    println("  providers: ", sort(collect(keys(providers)))); flush(stdout)

    # ---- PREPARE (extent → rewrite → coords → VI → gated fetch → graph) ------
    println("prepare(pushdown_rewrite=true) — N_REC discovered by the loader ...")
    flush(stdout)
    insp = EA.BuildInspection()
    t_prep = @elapsed prep = EA.prepare(doc_raw; providers=providers, base_path=ISRM_DIR,
                                        inspect=insp, pushdown_rewrite=true)
    N_REC = discovered_records(doc_raw, insp)
    println("PREPARE done in $(round(t_prep, digits=1)) s  (peak RSS so far: ",
            round(peak_rss_bytes() / 2^30, digits=2), " GiB)"); flush(stdout)

    # ---- the engine-derived support set (for the contract record) -----------
    mf_keys = [k for k in keys(insp.const_arrays) if startswith(String(k), "pd_member_factor__")]
    isempty(mf_keys) && error("no pd_member_factor__* const array — did the rewrite fire?")
    members = sort!(Int.(insp.const_arrays[first(mf_keys)]))
    n_ppl = length(members)
    n_src = metaparam(doc_raw, "N_SRC")
    n_rcv = metaparam(doc_raw, "N_RCV")
    println("engine-derived support set: |members| = $n_ppl of $n_src source cells")
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
        n_src = n_src, n_rcv = n_rcv, n_rec = N_REC,
        ppl = members,
        pathways = pathways,
        total_pm25 = tp, deathsK = dK, deathsL = dL,
        timing = Dict("wall_seconds" => time() - T0,
                      "providers_seconds" => t_providers,
                      "build_seconds" => t_prep,
                      "eval_seconds" => t_eval,
                      "peak_rss_bytes" => peak_rss_bytes()))
    if !reduced
        (okK && okL) || error("full-scale totals off oracle (see above)")
    end
    return 0
end

exit(main())
