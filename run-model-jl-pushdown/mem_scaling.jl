#!/usr/bin/env julia
# =============================================================================
# mem_scaling.jl — how does Phase 1's memory scale with |ppl|, and WHAT drives it?
#
# The full-scale run (|ppl|=1520) was OOM-killed at 25.1 GB anon-rss ~15 min into
# `_observed_field("deathsK")`, still climbing. The gated FETCH is no longer the
# problem (the zarr scatter fix cut it 9.17 -> 1.56 GB and the run reached
# "BUILD done in 590.8 s"); the growth is in the observed EVALUATION, which is
# genuinely live data and so cannot be shrunk by --heap-size-hint.
#
# This runs the SAME pipeline as L3_full.jl at a reduced |ppl| and reports, per
# stage, both RSS and TRUE LIVENESS.
#
# IMPORTANT — RSS IS NOT A LIVENESS SIGNAL. glibc retains freed buffers in its
# arenas, so RSS overstates live memory (measured elsewhere in this project: a
# Python case showed 313 MB RSS vs 29.6 MB live). We therefore report
# Base.gc_live_bytes() as the primary number and RSS only as context.
#
#   FIRSTN=<n>   emission records to keep. Maps to |ppl|:
#                200->9  1000->56  2500->152  5000->230  10000->441  20000->922
#                43650->1520 (full)
#   ALLOCPROF=1  attach the allocation profiler to the deathsK evaluation and
#                write a .pprof (needs PProf; see mem_scaling_setup.jl)
#
# Peak RSS is sampled by a background task, since the peak occurs BETWEEN the
# explicit probe points.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc
import GeometryOps, GeoInterface
import JSON
# MUST be top level: `using` inside the `if ALLOCPROF` block below is not in
# scope when `Profile.Allocs.@profile` on the following line is macro-expanded.
using Profile
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

const FIRSTN    = parse(Int, get(ENV, "FIRSTN", "2500"))
const ALLOCPROF = get(ENV, "ALLOCPROF", "0") == "1"
const CPUPROF   = get(ENV, "CPUPROF", "0") == "1"
const SR_ROOT   = joinpath(SCRATCH, get(ENV, "L3_SR_DIR", "scaling_sr"))
const SR_MAP = ["SOA"=>"SR_SOA", "pNO3"=>"SR_pNO3", "pNH4"=>"SR_pNH4",
                "pSO4"=>"SR_pSO4", "PrimaryPM25"=>"SR_PrimaryPM25"]

rss_gb()  = parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^30
live_gb() = Base.gc_live_bytes() / 2^30

# ---- background peak-RSS sampler (the peak is between probe points) ---------
const PEAK_RSS  = Ref(0.0)
const PEAK_LIVE = Ref(0.0)
const SAMPLING  = Ref(true)
function start_sampler()
    Threads.@spawn begin
        while SAMPLING[]
            PEAK_RSS[]  = max(PEAK_RSS[],  rss_gb())
            PEAK_LIVE[] = max(PEAK_LIVE[], live_gb())
            sleep(0.25)
        end
    end
end

const STAGES = Tuple{String,Float64,Float64,Float64}[]
function stage!(label)
    GC.gc(); GC.gc()
    push!(STAGES, (label, live_gb(), rss_gb(), PEAK_RSS[]))
    println(rpad(label, 34), " live=", rpad(round(live_gb(), digits=2), 7),
            " rss=", rpad(round(rss_gb(), digits=2), 7),
            " peakRSS=", round(PEAK_RSS[], digits=2)); flush(stdout)
end

seed_empty_zattrs(cache, base, arrays) = begin
    base = rstrip(base, '/')
    for arr in arrays
        url = "$base/$arr/.zattrs"; key = EarthSciIO.cache_key(url)
        if EarthSciIO.get_blob(cache.store, key) === nothing
            tmp = EarthSciIO.staging_path(cache.store); write(tmp, "{}")
            EarthSciIO.put_blob!(cache.store, key, tmp)
        end
    end
end
evict_blobs!(root) = (bd=joinpath(root,"v1","blobs"); isdir(bd) && rm(bd;recursive=true,force=true))

mutable struct GatedSR
    cache::Any; provs::Dict{String,Any}; gate::Dict{String,Any}; calls::Vector{Any}
end
EA.provider_gate_spec(g::GatedSR) = g.gate
EA.provider_is_gated(g::GatedSR) = true
EA.provider_supports_selection(g::GatedSR) = true
EA.provider_refresh_times(g::GatedSR) = Float64[]
function _native(selection)
    axes = Any[]
    for ax in selection
        if ax isa Colon;       push!(axes, "all")
        elseif ax isa Integer; push!(axes, Dict("indices"=>[Int(ax)-1]))
        else;                  push!(axes, Dict("indices"=>Int.(collect(ax)).-1))
        end
    end
    Dict("axes"=>axes)
end
function EA.provider_sample(g::GatedSR, t::Real; selection=nothing)
    selection === nothing && error("GatedSR must be sampled WITH a selection")
    push!(g.calls, deepcopy(selection))
    native = _native(selection)
    out = Dict{String,Any}()
    for (zname, mname) in SR_MAP
        evict_blobs!(SR_ROOT); seed_empty_zattrs(g.cache, ZARR_URL, [zname])
        nd = EarthSciIO.materialize(g.provs[zname], Float64(t); select=native)
        out[mname] = Array{Float64}(nd[zname].data)
    end
    evict_blobs!(SR_ROOT)
    return out
end

# =============================================================================
start_sampler()
println("="^76)
println("mem_scaling: FIRSTN=$FIRSTN  ALLOCPROF=$ALLOCPROF  threads=$(Threads.nthreads())")
println("="^76)
stage!("start")

inp = build_inputs(; firstn=FIRSTN)
stage!("inputs built (N_REC=$(inp.N_REC))")

evict_blobs!(SR_ROOT)
sr_cache = EarthSciIO.Cache(; root=SR_ROOT)
provs = Dict{String,Any}(z => EarthSciIO.const_provider(sr_cache, ZARR_URL; format="zarr", variables=[z])
                         for (z,_) in SR_MAP)
gate = Dict{String,Any}("axes"=>Any[Dict("fixed"=>[0]), Dict("gated_by"=>"emis_src_cells"), "all"],
                        "applies_to"=>[m for (_,m) in SR_MAP])
gated = GatedSR(sr_cache, provs, gate, Any[])

doc0 = JSON.parsefile(MODEL)
mp = Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,"N_LAYER"=>3,"N_REC"=>inp.N_REC)
doc = resolve_sizes!(deepcopy(doc0), mp)
f = EA.load(deepcopy(doc); base_path=ISRM_DIR)
ca = Dict{String,Any}(
    "X"=>inp.X, "Y"=>inp.Y, "emis_annual"=>inp.emis_annual, "pollutant"=>inp.pollutant,
    "emis_lon"=>zeros(inp.N_REC), "emis_lat"=>zeros(inp.N_REC),
    "stkhgt"=>zeros(inp.N_REC), "stkdiam"=>zeros(inp.N_REC),
    "stktemp"=>zeros(inp.N_REC), "stkvel"=>zeros(inp.N_REC),
    "W"=>inp.W, "S"=>inp.S, "E"=>inp.E, "N"=>inp.N,
    "TotalPop"=>inp.TotalPop, "MortalityRate"=>inp.MortalityRate)
stage!("const arrays ready")

insp = EA.BuildInspection()
t_build = @elapsed EA.build_evaluator(doc; model_name="ISRM", const_arrays=ca, inspect=insp,
    _gated_providers=Dict{String,Any}("ISRM_SR"=>gated), _sample_time=0.0)
n_ppl = length(gated.calls[1][2])
stage!("build_evaluator done (|ppl|=$n_ppl)")
println("  build took ", round(t_build, digits=1), " s"); flush(stdout)

# ---- the expensive stage ----------------------------------------------------
peak_before = PEAK_RSS[]
EA._CELLWISE_FASTPATH_HITS[] = 0; EA._CELLWISE_FASTPATH_MISS[] = 0
local dK, t_dK
if CPUPROF
    # Statistical CPU sampler. deathsK runs ~745 s at ppl=9, so the default
    # buffer would overflow long before the end — size it for the whole run.
    Profile.clear()
    Profile.init(n = 20_000_000, delay = 0.005)
    t_dK = @elapsed Profile.@profile begin
        dK = EA._observed_field(insp, f, "ISRM", "deathsK")[1]
    end
    println("\n", "="^76)
    println("CPU PROFILE — flat, top self-time frames")
    println("="^76)
    Profile.print(format = :flat, sortedby = :count, mincount = 200, maxdepth = 200)
    println("\n", "="^76)
    println("CPU PROFILE — tree, top hot path (>2% of samples)")
    println("="^76)
    Profile.print(format = :tree, mincount = 500, maxdepth = 45)
    try
        @eval using PProf
        PProf.pprof(out = joinpath(@__DIR__, "cpu_firstn$(FIRSTN).pb.gz"), web = false)
        println("\nwrote cpu_firstn$(FIRSTN).pb.gz (view: pprof -http=:8080 <file>)")
    catch e
        println("\n(PProf unavailable: ", sprint(showerror, e), ")")
    end
elseif ALLOCPROF
    Profile.Allocs.clear()
    # SAMPLE_RATE must stay LOW here. This stage allocates tens of GB; the
    # allocation profiler retains a record per sampled allocation, so
    # sample_rate=1.0 would itself exhaust the cgroup. Attribution stays
    # proportionally correct at a low rate — we want the SHAPE (which types /
    # sites dominate), not an exact byte count.
    rate = parse(Float64, get(ENV, "SAMPLE_RATE", "0.001"))
    println("  (alloc profiler sample_rate=$rate)"); flush(stdout)
    t_dK = @elapsed Profile.Allocs.@profile sample_rate=rate begin
        dK = EA._observed_field(insp, f, "ISRM", "deathsK")[1]
    end
    res = Profile.Allocs.fetch()
    total = sum(a.size for a in res.allocs)
    println("\nALLOC PROFILE: ", length(res.allocs), " allocations, ",
            round(total/2^30, digits=3), " GB total allocated")
    # top allocating types
    bytype = Dict{Any,Int}()
    for a in res.allocs; bytype[a.type] = get(bytype, a.type, 0) + a.size; end
    println("\ntop allocating TYPES:")
    for (ty, sz) in sort(collect(bytype), by=last, rev=true)[1:min(12,end)]
        println("  ", rpad(round(sz/2^20, digits=1), 10), " MB  ", ty)
    end
    # top allocation sites
    bysite = Dict{String,Int}()
    for a in res.allocs
        fr = isempty(a.stacktrace) ? nothing : a.stacktrace[1]
        key = fr === nothing ? "<unknown>" : string(fr.file, ":", fr.line, " ", fr.func)
        bysite[key] = get(bysite, key, 0) + a.size
    end
    println("\ntop allocation SITES:")
    for (site, sz) in sort(collect(bysite), by=last, rev=true)[1:min(15,end)]
        println("  ", rpad(round(sz/2^20, digits=1), 10), " MB  ", site)
    end
    try
        @eval using PProf
        PProf.Allocs.pprof(res; out=joinpath(@__DIR__, "allocs_firstn$(FIRSTN).pb.gz"), web=false)
        println("\nwrote allocs_firstn$(FIRSTN).pb.gz (view: pprof -http=:8080 <file>)")
    catch e
        println("\n(PProf unavailable: ", sprint(showerror, e), ")")
    end
else
    t_dK = @elapsed (dK = EA._observed_field(insp, f, "ISRM", "deathsK")[1])
end
stage!("deathsK evaluated")
println("  deathsK took ", round(t_dK, digits=1), " s   sum=", sum(dK),
        "   [fastpath hits=", EA._CELLWISE_FASTPATH_HITS[],
        " miss=", EA._CELLWISE_FASTPATH_MISS[], "]")
println("  peak RSS DURING deathsK = ", round(PEAK_RSS[] - peak_before, digits=2),
        " GB above pre-stage peak"); flush(stdout)

SAMPLING[] = false
println("\n", "="^76)
println("SUMMARY  FIRSTN=$FIRSTN  |ppl|=$n_ppl  N_REC=$(inp.N_REC)")
println("="^76)
println(rpad("stage",34), rpad("live GB",10), rpad("rss GB",10), "peakRSS GB")
for (l, lv, rs, pk) in STAGES
    println(rpad(l,34), rpad(round(lv,digits=2),10), rpad(round(rs,digits=2),10), round(pk,digits=2))
end
println("\nMACHINE-READABLE: ppl=$n_ppl nrec=$(inp.N_REC) build_s=$(round(t_build,digits=1)) ",
        "deathsK_s=$(round(t_dK,digits=1)) peak_rss_gb=$(round(PEAK_RSS[],digits=3)) ",
        "peak_live_gb=$(round(PEAK_LIVE[],digits=3)) sum_deathsK=$(sum(dK))")
