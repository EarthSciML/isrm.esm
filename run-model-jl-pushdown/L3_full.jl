#!/usr/bin/env julia
# =============================================================================
# L3 — full ISRM oracle through the projection-pushdown, on REAL S3 data.
#
# Drives the authored isrm_pushdown.esm through EA.prepare with a REAL gated SR
# provider wrapping the live s3://inmap-model/isrm_v1.2.1.zarr store.  Value-
# invention derives the emission-bearing source cells (ppl); the gated provider
# is DEFERRED and then fetches ONLY SR[layer0, ppl, :] — never the 330 GB whole
# matrix — and the downstream graph yields deaths.
#
#   FULL run  (default)          → assert sum(deathsK/L) ≈ 7524.92 / 16979.63.
#   REDUCED   (L3_REGION set)    → restrict emissions to a subregion (small ppl,
#                                  bounded disk) and cross-check against a direct
#                                  STEP-0 oracle computed on the same ppl rows.
#
# DISK-SAFE gated provider: wraps five per-array zarr providers sharing one cache
# and fetches+EVICTS pathway-by-pathway, so peak disk ≈ one pathway (~2 GB) not
# ~10 GB; the compact [ppl, rcv] Float64 slabs accumulate in RAM.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc
import GeometryOps, GeoInterface
import JSON
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

const SR_ROOT = joinpath(SCRATCH, get(ENV, "L3_SR_DIR", "l3_cache_sr"))

# zarr array name → model SR variable name (applies_to / const-factor key)
const SR_MAP = ["SOA"=>"SR_SOA", "pNO3"=>"SR_pNO3", "pNH4"=>"SR_pNH4",
                "pSO4"=>"SR_pSO4", "PrimaryPM25"=>"SR_PrimaryPM25"]
const FACT=28766.639; const POP_SCALE=1.0465819687408728; const MORT_SCALE=1.025229357798165
const RR_K=1.06; const RR_L=1.14

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
disk_used_gib(root) = (isdir(root) ? parse(Float64, split(readchomp(`du -sk $root`))[1])/2^20 : 0.0)

# ---- DISK-SAFE gated SR provider -------------------------------------------
mutable struct GatedSR
    cache::Any
    provs::Dict{String,Any}          # zarr name → EarthSciIO zarr provider (1 var each)
    gate::Dict{String,Any}
    calls::Vector{Any}
end
EA.provider_gate_spec(g::GatedSR) = g.gate
EA.provider_is_gated(g::GatedSR) = true
EA.provider_supports_selection(g::GatedSR) = true
EA.provider_refresh_times(g::GatedSR) = Float64[]

# translate the neutral per-axis selection into EarthSciIO's native select
function _native(selection)
    axes = Any[]
    for ax in selection
        if ax isa Colon;            push!(axes, "all")
        elseif ax isa Integer;      push!(axes, Dict("indices"=>[Int(ax)-1]))
        else;                       push!(axes, Dict("indices"=>Int.(collect(ax)).-1))
        end
    end
    return Dict("axes"=>axes)
end

function EA.provider_sample(g::GatedSR, t::Real; selection=nothing)
    selection === nothing && error("GatedSR must be sampled WITH a selection (never wholesale)")
    push!(g.calls, deepcopy(selection))
    native = _native(selection)
    out = Dict{String,Any}()
    for (zname, mname) in SR_MAP
        evict_blobs!(SR_ROOT)                    # fresh, empty blob store for this pathway
        seed_empty_zattrs(g.cache, ZARR_URL, [zname])
        nd = EarthSciIO.materialize(g.provs[zname], Float64(t); select=native)
        out[mname] = Array{Float64}(nd[zname].data)   # compact (1, |ppl|, rcv), kept in RAM
        du = disk_used_gib(SR_ROOT)
        println("    [gated] fetched $zname → $(size(out[mname]))  (pathway disk $(round(du,digits=2)) GiB)"); flush(stdout)
    end
    evict_blobs!(SR_ROOT)                         # free the last pathway's chunks
    return out
end

# ---- reduction knob (reduced run) -------------------------------------------
region = nothing; firstn = nothing
if haskey(ENV, "L3_REGION")
    region = Tuple(parse.(Float64, split(ENV["L3_REGION"], ",")))   # xlo,xhi,ylo,yhi
    println("REDUCED run — subregion X∈[$(region[1]),$(region[2])] Y∈[$(region[3]),$(region[4])]")
elseif haskey(ENV, "L3_FIRSTN")
    firstn = parse(Int, ENV["L3_FIRSTN"])
    println("REDUCED run — first $firstn emission records")
else
    println("FULL run — whole domain (target deathsK≈7524.92, deathsL≈16979.63)")
end
const REDUCED = region !== nothing || firstn !== nothing

println("building inputs ...")
inp = build_inputs(; region=region, firstn=firstn)
println("  N_REC=", inp.N_REC)

# ---- construct the gated provider ------------------------------------------
evict_blobs!(SR_ROOT)
sr_cache = EarthSciIO.Cache(; root=SR_ROOT)
provs = Dict{String,Any}(zname =>
    EarthSciIO.const_provider(sr_cache, ZARR_URL; format="zarr", variables=[zname])
    for (zname, _) in SR_MAP)
gate = Dict{String,Any}(
    "axes" => Any[Dict("fixed"=>[0]), Dict("gated_by"=>"emis_src_cells"), "all"],
    "applies_to" => [mname for (_, mname) in SR_MAP])
gated = GatedSR(sr_cache, provs, gate, Any[])

# ---- load + build via the DICT FRONT-DOOR (VI runs here; bare const arrays) --
doc0 = JSON.parsefile(MODEL)
mp = Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,"N_LAYER"=>3,"N_REC"=>inp.N_REC)
doc = resolve_sizes!(deepcopy(doc0), mp)
f = EA.load(deepcopy(doc); base_path=ISRM_DIR)
ca = Dict{String,Any}(
    "X"=>inp.X, "Y"=>inp.Y, "emis_annual"=>inp.emis_annual, "pollutant"=>inp.pollutant,
    "emis_lon"=>zeros(inp.N_REC), "emis_lat"=>zeros(inp.N_REC),   # consumer-less now (X/Y are params)
    "stkhgt"=>zeros(inp.N_REC), "stkdiam"=>zeros(inp.N_REC),
    "stktemp"=>zeros(inp.N_REC), "stkvel"=>zeros(inp.N_REC),
    "W"=>inp.W, "S"=>inp.S, "E"=>inp.E, "N"=>inp.N,
    "TotalPop"=>inp.TotalPop, "MortalityRate"=>inp.MortalityRate)

insp = EA.BuildInspection()
println("build_evaluator (value-invention → gated fetch → build) ...")
t_prep = @elapsed EA.build_evaluator(doc; model_name="ISRM", const_arrays=ca, inspect=insp,
    _gated_providers=Dict{String,Any}("ISRM_SR"=>gated), _sample_time=0.0)
println("BUILD done in ", round(t_prep, digits=1), " s")

# ---- selection assertions ---------------------------------------------------
@assert length(gated.calls) >= 1 "gated provider was never sampled"
members = gated.calls[1][2]           # selection = [layer_int, members_vec, Colon]
n_ppl = length(members)
println("gated selection: layer=", gated.calls[1][1], "  |members|=", n_ppl,
        "  rcv=", gated.calls[1][3])
if !REDUCED
    @assert n_ppl == length(inp.ppl0) "expected $(length(inp.ppl0)) ppl, got $n_ppl"
    @assert Int.(members) == sort(inp.ppl0 .+ 1) "members ≠ validated run-model ppl+1"
    println("  members == validated run-model ppl (", n_ppl, ") ✓  ONLY these rows fetched (not $(inp.N_SRC))")
end

# ---- evaluate deaths --------------------------------------------------------
rt(v)=(fld=EA._observed_field(insp, f, "ISRM", v); fld===nothing ? error("no $v") : fld[1])
dK = rt("deathsK"); dL = rt("deathsL"); tp = rt("TotalPM25")
sK = sum(dK); sL = sum(dL)
println("\n", "="^70)
println("  sum(deathsK) = ", sK)
println("  sum(deathsL) = ", sL)
println("  Σ TotalPM25  = ", sum(tp))

if !REDUCED
    okK = isapprox(sK, 7524.918845602511; rtol=1e-4)
    okL = isapprox(sL, 16979.632171487083; rtol=1e-4)
    println("  target deathsK=7524.918845602511  rel.err ", round(100*(sK-7524.918845602511)/7524.918845602511, digits=4), "%")
    println("  target deathsL=16979.632171487083 rel.err ", round(100*(sL-16979.632171487083)/16979.632171487083, digits=4), "%")
    println("L3 FULL: ", (okK && okL) ? "PASS" : "FAIL")
    println("="^70)
    (okK && okL) || error("L3 full oracle MISMATCH")
else
    # cross-check reduced deaths against a direct STEP-0 oracle on the same ppl rows
    println("  (reduced run — cross-check against direct oracle below)")
    println("="^70)
    oK, oL = reduced_oracle(inp, Int.(members))
    println("  direct-oracle sum(deathsK)=", oK, "   sum(deathsL)=", oL)
    mdK = abs(sK-oK); mdL = abs(sL-oL)
    println("  |Δ| deathsK=", mdK, "  deathsL=", mdL)
    pass = mdK < 1e-6*max(1.0,abs(oK)) && mdL < 1e-6*max(1.0,abs(oL))
    println("L3 REDUCED (real-S3 pushdown == direct oracle): ", pass ? "PASS" : "FAIL")
    pass || error("reduced pushdown ≠ direct oracle")
end
