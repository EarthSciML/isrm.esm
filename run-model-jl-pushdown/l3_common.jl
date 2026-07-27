# Shared L3 inputs: real EGU emissions (from run-model.jl's stage1 checkpoint),
# the SAME LCC forward projection run-model.jl validated, and the full-grid
# geometry / pop / mortality read (fully) from the real ISRM zarr.  Used by the
# VI probe and the full L3.  NO SR chunks are fetched here (SR is the gated,
# deferred fetch).
import Serialization

# RUNMODEL / ZARR_URL / N_SRC / SCRATCH / MODEL / ISRM_DIR all resolve here.
include(joinpath(@__DIR__, "paths.jl"))

# ---- LCC constants + forward (verbatim from run-model.jl) -------------------
const LAT_1=33.0; const LAT_2=45.0; const LAT_0=40.0; const LON_0=-97.0; const LCC_R=6370997.0
const D2R=0.017453292519943295
lcc_t(lat)=tan(0.7853981633974483+lat*0.008726646259971648)
const LCC_N=log(cos(LAT_1*D2R)/cos(LAT_2*D2R))/log(lcc_t(LAT_2)/lcc_t(LAT_1))
const LCC_F=cos(LAT_1*D2R)*lcc_t(LAT_1)^LCC_N/LCC_N
lcc_rho(lat)=LCC_R*LCC_F/lcc_t(lat)^LCC_N
const LCC_RHO0=lcc_rho(LAT_0)
function lcc_forward(lon, lat)
    theta = LCC_N*(lon-LON_0)*D2R; rho = lcc_rho(lat)
    return (rho*sin(theta), LCC_RHO0-rho*cos(theta))
end

# ---- pollutant → pathway code (STEP 0.6). is_* masks pick the pathway; the
# .esm's is_VOC etc. select by code band, so we encode the same pathway as an
# integer pollutant code (1=VOC, 36=NOx, 40=NH3, 41=SOx, 42=PM25). ------------
pathcode(p) = p==:VOC ? 1.0 : p==:NOx ? 36.0 : p==:NH3 ? 40.0 : p==:SOx ? 41.0 :
              p==:PM25 ? 42.0 : 0.0

"""Build the full-scale inSRM inputs. Returns a NamedTuple. `subset` optionally
restricts to emission records whose projected X falls in [xlo,xhi] & Y in
[ylo,yhi] (a subregion, to shrink ppl for a disk-bounded real-S3 run)."""
function build_inputs(; region=nothing, firstn=nothing)
    st = Serialization.deserialize(joinpath(RUNMODEL, "checkpoints", "stage1.jls"))
    lon = Float64.(st["lon"]); lat = Float64.(st["lat"]); ann = Float64.(st["ann"])
    paths = st["paths"]
    # project every recognized record with the validated LCC
    X = similar(lon); Y = similar(lat)
    for i in eachindex(lon); X[i],Y[i] = lcc_forward(lon[i], lat[i]); end
    pollutant = Float64[pathcode(p) for p in paths]

    keep = trues(length(lon))
    if region !== nothing
        xlo,xhi,ylo,yhi = region
        keep = [xlo<=X[i]<xhi && ylo<=Y[i]<yhi for i in eachindex(X)]
    elseif firstn !== nothing
        keep = falses(length(lon)); keep[1:min(firstn,length(lon))] .= true
    end
    X=X[keep]; Y=Y[keep]; ann=ann[keep]; pollutant=pollutant[keep]
    N_REC = length(X)

    # full-grid geometry / pop / mortality (small; read fully, cached in cache_meta)
    meta_cache = EarthSciIO.Cache(; root=joinpath(RUNMODEL, "cache_meta"))
    pg = EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr",
             variables=["W","S","E","N"])
    geo = EarthSciIO.materialize(pg)
    ph = EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr",
             variables=["TotalPop","MortalityRate"])
    h = EarthSciIO.materialize(ph)
    # container grid = the first N_SRC cells (the SR source cells), exactly as the
    # validated run-model.jl bins (pop_cells == SR cells; the [0:52411] prefix).
    W  = Float64.(geo["W"].data[1:N_SRC]);  S  = Float64.(geo["S"].data[1:N_SRC])
    E  = Float64.(geo["E"].data[1:N_SRC]);  Nn = Float64.(geo["N"].data[1:N_SRC])
    TotalPop      = Float64.(h["TotalPop"].data[1:N_SRC])
    MortalityRate = Float64.(h["MortalityRate"].data[1:N_SRC])

    return (X=X, Y=Y, emis_annual=ann, pollutant=pollutant, N_REC=N_REC,
            W=W, S=S, E=E, N=Nn, TotalPop=TotalPop, MortalityRate=MortalityRate,
            N_SRC=N_SRC, ppl0=Int.(st["ppl"]))   # ppl0 = run-model's 0-based ppl
end

# ---- front-door build helpers ----------------------------------------------
# VALUE-INVENTION runs ONLY through the AbstractDict front-door of build_evaluator
# (build.jl: "runs only through the AbstractDict front-door"); `prepare` flattens
# the coupled esm and namespaces vars, which breaks join-env resolution for the
# loader-fed / observed producer inputs (W/S/E/N/X/Y).  So we drive the SAME VI
# mechanism the committed L1 test uses: build_evaluator(doc_dict; bare const
# arrays + _gated_providers), on a metaparameter-resolved copy of the doc.

"""Replace string index-set sizes ("N_SRC", …) with integers from `mp`."""
function resolve_sizes!(doc::AbstractDict, mp::AbstractDict)
    iss = get(doc, "index_sets", nothing); iss === nothing && return doc
    for (_, is) in iss
        is isa AbstractDict || continue
        sz = get(is, "size", nothing)
        if sz isa AbstractString && haskey(mp, sz)
            is["size"] = Int(mp[sz])
        end
    end
    return doc
end
# Fetches SR[0, members, :] DIRECTLY (a fresh cache, independent of the pushdown
# path), bins the emissions into `members` by rectangle containment, contracts,
# and returns (sum deathsK, sum deathsL).  Used to cross-check a reduced run.
const _PATHBANDS = Dict("SOA"=>(1.0,35.0),"pNO3"=>(36.0,39.0),"pNH4"=>(40.0,40.0),
                        "pSO4"=>(41.0,41.0),"PrimaryPM25"=>(42.0,59.0))
function reduced_oracle(inp, members::Vector{Int})
    FACT=28766.639; POP_SCALE=1.0465819687408728; MORT_SCALE=1.025229357798165
    RR_K=1.06; RR_L=1.14
    NP = length(members); comp = Dict(members[k]=>k for k in 1:NP)
    W=inp.W; E=inp.E; S=inp.S; N=inp.N; X=inp.X; Y=inp.Y
    # E_p per pathway, indexed in members order
    Ep = Dict(p=>zeros(NP) for p in keys(_PATHBANDS))
    for r in 1:inp.N_REC
        # find the member cell (1-based) containing emission r
        cell = 0
        for m in members
            if W[m] <= X[r] < E[m] && S[m] <= Y[r] < N[m]; cell = m; break; end
        end
        cell == 0 && continue
        for (p,(lo,hi)) in _PATHBANDS
            lo <= inp.pollutant[r] <= hi && (Ep[p][comp[cell]] += inp.emis_annual[r])
        end
    end
    # fetch SR[0, members, :] directly for each pathway, contract
    root = joinpath(SCRATCH, "l3_oracle_sr")
    isdir(root) && rm(root; recursive=true, force=true)
    cache = EarthSciIO.Cache(; root=root)
    base = rstrip(ZARR_URL,'/')
    for arr in keys(_PATHBANDS)   # seed .zattrs
        u="$base/$arr/.zattrs"; k=EarthSciIO.cache_key(u)
        if EarthSciIO.get_blob(cache.store,k)===nothing
            tmp=EarthSciIO.staging_path(cache.store); write(tmp,"{}"); EarthSciIO.put_blob!(cache.store,k,tmp)
        end
    end
    sel = Dict("axes"=>Any[Dict("indices"=>[0]), Dict("indices"=>members.-1), "all"])
    NRCV = inp.N_SRC
    TotalPM25 = zeros(NRCV)
    for (arr, p) in ("SOA"=>"SOA","pNO3"=>"pNO3","pNH4"=>"pNH4","pSO4"=>"pSO4","PrimaryPM25"=>"PrimaryPM25")
        prov = EarthSciIO.const_provider(cache, ZARR_URL; format="zarr", variables=[arr])
        nd = EarthSciIO.materialize(prov; select=sel)
        slab = nd[arr].data      # (1, NP, NRCV)
        conc = zeros(NRCV)
        for j in 1:NP
            ev = Ep[p][j]; ev == 0.0 && continue
            @inbounds @simd for rcv in 1:NRCV
                conc[rcv] += slab[1,j,rcv]*ev
            end
        end
        TotalPM25 .+= FACT .* conc
    end
    dK = (exp.(log(RR_K)/10 .* TotalPM25) .- 1) .* inp.TotalPop .* POP_SCALE .* inp.MortalityRate ./ 100000 .* MORT_SCALE
    dL = (exp.(log(RR_L)/10 .* TotalPM25) .- 1) .* inp.TotalPop .* POP_SCALE .* inp.MortalityRate ./ 100000 .* MORT_SCALE
    return sum(dK), sum(dL)
end
