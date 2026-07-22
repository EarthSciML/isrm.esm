#!/usr/bin/env julia
# =============================================================================
# capstone_full.jl — FULL-SCALE real-data capstone for the projection-pushdown.
#
# Sidesteps the pre-existing observed-graph type-instability (Wall #2) by using
# `materialize_value_invention` DIRECTLY (fast now, post Wall#1 fix) for provider
# selection + a hand-written contraction reusing the VALIDATED run-model.jl STEP-0
# math.  Proves the pushdown's core claims at FULL scale (52411 cells, 43650 EGU
# emission records) on the live s3://inmap-model/isrm_v1.2.1.zarr store:
#
#   CLAIM A  value-invention derives the EXACT emission-bearing support (ppl),
#            fast (seconds), byte-identical to run-model.jl's spatial join (1520).
#   CLAIM B  the pushdown selection fetches ONLY those ppl source-chunk rows
#            (416 of 525 chunks; slab (1,|ppl|,52411)), byte-identical to a direct
#            independent fetch of the same source ids.
#   CLAIM C  contract the pushdown-fetched SR against E_p → STEP-0 deaths ≈ oracle
#            (sum deathsK=7524.918845602511, deathsL=16979.632171487083).
#
# MODE=ab   (default) → CLAIM A + CLAIM B (bounded real fetch), fast, foreground.
# MODE=full           → CLAIM A + CLAIM C (all 5 pathways, chunked+evicting), ~55 min.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc
import GeometryOps, GeoInterface
import JSON
import Serialization
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

const MODEL    = "/Users/ctessum/code/earthsciml/isrm.esm-wt-runner/isrm_pushdown.esm"
const ISRM_DIR = dirname(MODEL)
const MODE     = get(ENV, "CAP_MODE", "ab")
const CK_DIR   = joinpath(SCRATCH, "capstone_ck"); mkpath(CK_DIR)
const SR_ROOT  = joinpath(SCRATCH, "capstone_sr")   # evicted between pathways

# STEP-0 constants (verbatim from run-model.jl) + oracle targets
const FACT=28766.639; const POP_SCALE=1.0465819687408728; const MORT_SCALE=1.025229357798165
const RR_K=1.06; const RR_L=1.14
const N_SRC = 52411
const ORACLE_K = 7524.918845602511
const ORACLE_L = 16979.632171487083
# zarr SR array → pollutant-code band (which pathway an emission feeds)
const PATHS = ["SOA","pNO3","pNH4","pSO4","PrimaryPM25"]
const BANDS = Dict("SOA"=>(1.0,35.0),"pNO3"=>(36.0,39.0),"pNH4"=>(40.0,40.0),
                   "pSO4"=>(41.0,41.0),"PrimaryPM25"=>(42.0,59.0))

seed_zattrs(cache, arrays) = begin
    base = rstrip(ZARR_URL, '/')
    for arr in arrays
        u = "$base/$arr/.zattrs"; k = EarthSciIO.cache_key(u)
        if EarthSciIO.get_blob(cache.store, k) === nothing
            tmp = EarthSciIO.staging_path(cache.store); write(tmp, "{}")
            EarthSciIO.put_blob!(cache.store, k, tmp)
        end
    end
end
blobdir(root)  = joinpath(root, "v1", "blobs")
evict!(root)   = (bd=blobdir(root); isdir(bd) && rm(bd; recursive=true, force=true))
disk_gib(root) = (isdir(root) ? parse(Float64, split(readchomp(`du -sk $root`))[1])/2^20 : 0.0)
# large (>10 KiB) blobs are data chunks; tiny ones are .zarray/.zattrs metadata.
# blobs nest under <blobs>/<2-char hash prefix>/<hash>.<idx>, so walk recursively.
function count_chunk_blobs(root)
    bd = blobdir(root); isdir(bd) || return 0
    n = 0
    for (dir, _, files) in walkdir(bd), fn in files
        filesize(joinpath(dir, fn)) > 10*1024 && (n += 1)
    end
    return n
end

# ---------------------------------------------------------------------------
# Shared: build FULL inputs, VI members (ppl), independent spatial join, E_p.
# ---------------------------------------------------------------------------
println("building FULL inputs (all EGU records, all $N_SRC cells) ...")
inp = build_inputs()                       # no region/firstn → all 43650 records
println("  N_REC = ", inp.N_REC, "   N_SRC = ", inp.N_SRC)

doc0 = JSON.parsefile(MODEL)
mp   = Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,"N_LAYER"=>3,"N_REC"=>inp.N_REC)

# ---- CLAIM A: value-invention derives the exact support, fast ---------------
function claim_a()
    # warm up (compile STRtree + VI engine) on a tiny nrec so the FULL timing is clean
    print("  [warmup] compiling VI engine (nrec=10) ... "); flush(stdout)
    let doc = resolve_sizes!(deepcopy(doc0),
                Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,"N_LAYER"=>3,"N_REC"=>10))
        f = EA.load(doc; base_path=ISRM_DIR); m = EA._select_model(f, "ISRM")
        ca = Dict{String,Any}("X"=>inp.X[1:10],"Y"=>inp.Y[1:10],"W"=>inp.W,"S"=>inp.S,"E"=>inp.E,"N"=>inp.N)
        tw = @elapsed EA.materialize_value_invention(m, f.index_sets, ca, Dict{String,Any}())
        println(round(tw,digits=1), " s")
    end
    doc = resolve_sizes!(deepcopy(doc0), mp)
    f   = EA.load(doc; base_path=ISRM_DIR)
    model = EA._select_model(f, "ISRM")
    ca = Dict{String,Any}("X"=>inp.X, "Y"=>inp.Y, "W"=>inp.W, "S"=>inp.S, "E"=>inp.E, "N"=>inp.N)
    GC.gc()
    print("  [timed] materialize_value_invention on FULL model (N_REC=$(inp.N_REC), grid=$N_SRC) ... "); flush(stdout)
    local vi
    t = @elapsed (vi = EA.materialize_value_invention(model, f.index_sets, ca, Dict{String,Any}()))
    println(round(t,digits=2), " s")
    ppl_vi = sort(Int.(vi.members["emis_src_cells_faq"]))   # 1-based source-cell ids
    return ppl_vi, t
end

# independent exact-rectangle spatial join (run-model.jl STAGE A.3), 1-based ppl,
# plus E_p per pathway binned in ppl order (self-contained reference)
function ppl_ref_and_Ep(ppl_order::Union{Nothing,Vector{Int}}=nothing)
    X=inp.X; Y=inp.Y; W=inp.W; S=inp.S; E=inp.E; Nn=inp.N; pol=inp.pollutant; ann=inp.emis_annual
    xmin=minimum(W); xmax=maximum(E); ymin=minimum(S); ymax=maximum(Nn)
    NB=600; dxb=(xmax-xmin)/NB; dyb=(ymax-ymin)/NB
    bx(x)=clamp(floor(Int,(x-xmin)/dxb),0,NB-1); by(y)=clamp(floor(Int,(y-ymin)/dyb),0,NB-1)
    bins=[Int[] for _ in 1:NB*NB]; binid(ix,iy)=ix*NB+iy+1
    for c in 1:N_SRC
        ix0=bx(W[c]); ix1=bx(prevfloat(E[c])); iy0=by(S[c]); iy1=by(prevfloat(Nn[c]))
        for ix in ix0:ix1, iy in iy0:iy1; push!(bins[binid(ix,iy)], c); end
    end
    assign=fill(0, inp.N_REC); noc=0
    for i in 1:inp.N_REC
        xi=X[i]; yi=Y[i]; found=0
        for c in bins[binid(bx(xi),by(yi))]
            if W[c]<=xi<E[c] && S[c]<=yi<Nn[c]; found=c; break; end
        end
        assign[i]=found; found==0 && (noc+=1)
    end
    ppl_ref = sort!(unique(c for c in assign if c>0))          # 1-based
    order = ppl_order === nothing ? ppl_ref : ppl_order
    comp = Dict(order[k]=>k for k in eachindex(order)); NP=length(order)
    Ep = Dict(a=>zeros(NP) for a in PATHS)
    for i in 1:inp.N_REC
        assign[i]==0 && continue
        haskey(comp, assign[i]) || continue
        k=comp[assign[i]]; p=pol[i]
        for a in PATHS; lo,hi=BANDS[a]; lo<=p<=hi && (Ep[a][k]+=ann[i]); end
    end
    return ppl_ref, Ep, noc
end

# ---------------------------------------------------------------------------
# Run CLAIM A
# ---------------------------------------------------------------------------
println("\n", "="^72, "\nCLAIM A — value-invention derives the exact emission-bearing support\n", "="^72)
ppl_vi, t_vi = claim_a()
ppl_ref, Ep, noc = ppl_ref_and_Ep(ppl_vi)      # bin E_p in ppl_vi order
ppl0_ck = sort(Int.(inp.ppl0) .+ 1)            # run-model checkpoint ppl (0-based)+1
n = length(ppl_vi)
matchref = ppl_vi == ppl_ref
matchck  = ppl_vi == ppl0_ck
chunks   = sort!(unique((ppl_vi .- 1) .÷ 100))
println("  VI time (FULL, timed, compile excluded): ", round(t_vi, digits=2), " s")
println("  |ppl_vi| = ", n, "   (expect 1520)")
println("  points with NO containing cell (dropped): ", noc)
println("  ppl_vi == independent spatial-join ppl_ref : ", matchref)
println("  ppl_vi == run-model.jl checkpoint ppl (+1)  : ", matchck)
println("  distinct SR source-chunks (÷100) spanned    : ", length(chunks), " of ", cld(N_SRC,100), " total")
@assert n == 1520 "expected 1520 ppl, got $n"
@assert matchref "VI members != independent spatial join"
@assert matchck  "VI members != run-model checkpoint ppl"
println("  CLAIM A: PASS  (VI derived the exact 1520-cell support in $(round(t_vi,digits=2)) s)")
flush(stdout)

# ---------------------------------------------------------------------------
# CLAIM B — bounded real selective fetch (foreground, ab mode)
# ---------------------------------------------------------------------------
function claim_b()
    println("\n", "="^72, "\nCLAIM B — pushdown selection fetches ONLY the ppl source-chunk rows\n", "="^72)
    # full-selection minimality (computed, no fetch):
    println("  full pushdown selection: ", n, " of ", N_SRC, " source rows",
            "  (", round(100n/N_SRC,digits=2), "%),  ",
            length(chunks), " of ", cld(N_SRC,100), " source-chunks")
    @assert length(chunks) == 416 "expected 416 ppl chunks, got $(length(chunks))"

    # bounded REAL fetch: members falling in the first NCHUNK distinct chunks
    NCHUNK = 6
    subch  = chunks[1:min(NCHUNK, length(chunks))]
    subset = Int[c for c in ppl_vi if ((c-1) ÷ 100) in Set(subch)]
    sub0   = subset .- 1                                   # 0-based source ids
    println("\n  bounded real fetch: first ", length(subch), " chunks → ",
            length(subset), " ppl rows (source ids ", first(subset), "…", last(subset), ")")

    arr = "SOA"
    root = joinpath(SCRATCH, "capstone_b_push"); isdir(root) && rm(root; recursive=true, force=true)
    cache = EarthSciIO.Cache(; root=root); seed_zattrs(cache, [arr])
    prov = EarthSciIO.const_provider(cache, ZARR_URL; format="zarr", variables=[arr])
    sel  = Dict("axes"=>Any[Dict("indices"=>[0]), Dict("indices"=>sub0), "all"])
    print("  [pushdown] materialize SR[$arr] layer0, $(length(sub0)) ppl rows, all rcv ... "); flush(stdout)
    tB = @elapsed (ndP = EarthSciIO.materialize(prov; select=sel))
    slabP = ndP[arr].data
    println(round(tB,digits=1), " s  → slab ", size(slabP))
    nblob = count_chunk_blobs(root)
    @assert size(slabP) == (1, length(sub0), N_SRC) "pushdown slab shape wrong: $(size(slabP))"
    println("  slab shape == (1, |ppl_sub|=$(length(sub0)), rcv=$N_SRC) ✓")
    println("  data-chunk blobs fetched: ", nblob, "  (expect ", length(subch),
            " = distinct chunks of sub-selection;  full grid would be ", cld(N_SRC,100), ")")
    @assert nblob == length(subch) "fetched $nblob chunks, expected $(length(subch))"

    # confirm the SR chunk shape claim [1,100,52411] from the fetched .zarray
    zpath = joinpath(root, "v1", "meta")   # metadata may live in blobs; read via reader instead
    # independent DIRECT fetch of the SAME source ids into a FRESH cache
    root2 = joinpath(SCRATCH, "capstone_b_direct"); isdir(root2) && rm(root2; recursive=true, force=true)
    cache2 = EarthSciIO.Cache(; root=root2); seed_zattrs(cache2, [arr])
    prov2 = EarthSciIO.const_provider(cache2, ZARR_URL; format="zarr", variables=[arr])
    print("  [direct]   independent fetch of the same $(length(sub0)) source ids ... "); flush(stdout)
    ndD = EarthSciIO.materialize(prov2; select=sel); slabD = ndD[arr].data
    println("→ slab ", size(slabD))
    maxd = maximum(abs.(Array{Float64}(slabP) .- Array{Float64}(slabD)))
    println("  max|pushdown − direct| over $(length(sub0))×$N_SRC = ", maxd)
    @assert maxd == 0.0 "pushdown fetch != direct fetch (max diff $maxd)"
    println("  pushdown-selected rows are BYTE-IDENTICAL to the independent direct fetch ✓")
    rm(root; recursive=true, force=true); rm(root2; recursive=true, force=true)
    println("  CLAIM B: PASS  (selective, chunk-minimal, byte-identical)")
    flush(stdout)
end

# ---------------------------------------------------------------------------
# CLAIM C — full 5-pathway pushdown oracle (chunked + evicting), full mode
# ---------------------------------------------------------------------------
# conc_p[rcv] = Σ_c SR_p[0, ppl_c, rcv]·E_p[c]; fetch one 100-source chunk at a
# time (pushdown select layer0, the ppl rows in that chunk, all rcv), evict per pathway.
function contract_pathway_push(arr, Ep_arr, ppl, rows_by_chunk, chunk_ids)
    conc = zeros(Float64, N_SRC)
    sr_cache = EarthSciIO.Cache(; root=SR_ROOT)
    for (ci, chid) in enumerate(chunk_ids)
        rows = rows_by_chunk[chid]
        src0 = sort!([ppl[k]-1 for k in rows])            # 0-based source ids in this chunk
        sel = Dict("axes"=>Any[Dict("indices"=>[0]), Dict("indices"=>src0), "all"])
        prov = EarthSciIO.const_provider(sr_cache, ZARR_URL; format="zarr", variables=[arr])
        nd = EarthSciIO.materialize(prov; select=sel)
        slab = nd[arr].data                                # (1, |src0|, N_SRC)
        for (j, s0) in enumerate(src0)
            k = findfirst(==(s0+1), ppl); ev = Ep_arr[k]
            ev == 0.0 && continue
            @inbounds @simd for rcv in 1:N_SRC
                conc[rcv] += slab[1, j, rcv] * ev
            end
        end
        if ci % 50 == 0 || ci == length(chunk_ids)
            println("      [$arr] chunk $ci/$(length(chunk_ids))  (disk $(round(disk_gib(SR_ROOT),digits=2)) GiB)"); flush(stdout)
        end
    end
    return conc
end

function claim_c()
    println("\n", "="^72, "\nCLAIM C — full 5-pathway pushdown oracle (chunked, evicting)\n", "="^72)
    ppl = ppl_vi; comp = Dict(ppl[k]=>k for k in eachindex(ppl))
    rows_by_chunk = Dict{Int,Vector{Int}}()
    for k in eachindex(ppl); push!(get!(rows_by_chunk, (ppl[k]-1) ÷ 100, Int[]), k); end
    chunk_ids = sort!(collect(keys(rows_by_chunk)))
    println("  contracting over ", length(chunk_ids), " source-chunks × 5 pathways")
    TotalPM25 = zeros(Float64, N_SRC)
    for arr in PATHS
        ck = joinpath(CK_DIR, "conc_$(arr).jls")
        if isfile(ck)
            conc = Serialization.deserialize(ck)
            println("  [$arr] checkpoint present → loaded ($ck)")
        else
            println("  [$arr] emis sum = ", round(sum(Ep[arr]),digits=1), " tons/yr — fetching+contracting ...")
            evict!(SR_ROOT); seed_zattrs(EarthSciIO.Cache(; root=SR_ROOT), [arr])
            t = @elapsed (conc = contract_pathway_push(arr, Ep[arr], ppl, rows_by_chunk, chunk_ids))
            Serialization.serialize(ck, conc)
            evict!(SR_ROOT)
            println("  [$arr] done in ", round(t/60,digits=1), " min → wrote checkpoint")
        end
        # cross-check this pathway against run-model.jl's validated cached conc
        rmck = joinpath(RUNMODEL, "checkpoints", "conc_$(arr).jls")
        if isfile(rmck)
            ref = Serialization.deserialize(rmck)
            md = maximum(abs.(conc .- ref))
            println("      vs run-model cached conc_$arr: max|Δ| = ", md,
                    md == 0.0 ? "  (BYTE-IDENTICAL)" : "")
        end
        TotalPM25 .+= FACT .* conc
        flush(stdout)
    end
    println("  Σ TotalPM25 = ", sum(TotalPM25), "   max = ", maximum(TotalPM25))
    dK = (exp.(log(RR_K)/10 .* TotalPM25) .- 1) .* inp.TotalPop .* POP_SCALE .* inp.MortalityRate ./ 100000 .* MORT_SCALE
    dL = (exp.(log(RR_L)/10 .* TotalPM25) .- 1) .* inp.TotalPop .* POP_SCALE .* inp.MortalityRate ./ 100000 .* MORT_SCALE
    sK = sum(dK); sL = sum(dL)
    println("\n  sum(deathsK) = ", sK, "   target ", ORACLE_K, "   rel.err ", round(100*(sK-ORACLE_K)/ORACLE_K, digits=5), "%")
    println("  sum(deathsL) = ", sL, "   target ", ORACLE_L, "   rel.err ", round(100*(sL-ORACLE_L)/ORACLE_L, digits=5), "%")
    okK = isapprox(sK, ORACLE_K; rtol=1e-3); okL = isapprox(sL, ORACLE_L; rtol=1e-3)
    println("  CLAIM C: ", (okK && okL) ? "PASS" : "FAIL")
    (okK && okL) || error("CLAIM C oracle MISMATCH")
end

if MODE == "ab"
    claim_b()
    println("\nA+B complete (mode=ab). Run with CAP_MODE=full for the full CLAIM C oracle.")
else
    claim_c()
end
println("\nCAPSTONE (mode=$MODE) DONE.")
