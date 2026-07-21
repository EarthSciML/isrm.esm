#!/usr/bin/env julia
# =============================================================================
# run-model.jl — Part (A): direct numerical reproduction of the InMAP ISRM
# tutorial mortality totals, driven ENTIRELY by the new EarthSciIO `zarr` + `ff10`
# readers and the STEP-0 arithmetic of isrm.esm.
#
# Reproduces https://inmap.run/blog/2022/12/15/tutorial/ (minus the Bonus viz):
#   conc_p[rcv]     = Σ_c SR_p[0, ppl_c, rcv] · E_p[c]        (SR contraction)
#   TotalPM25[rcv]  = fact · Σ_p conc_p[rcv]
#   deaths[rcv]     = (exp(ln(RR)/10 · TotalPM25) − 1) · pop · pop_scale
#                       · mort · mort_scale / 100000
#   validate  sum(deathsK) ≈ 7524.84 (Krewski, RR 1.06)
#             sum(deathsL) ≈ 16979.45 (Lepeule, RR 1.14)
#
# RESUMABLE + DISK-SAFE (only ~21 GB free, SR cache would be ~14 GB otherwise):
#   * Stage A (ppl pre-pass) checkpoints to checkpoints/stage1.jls
#   * Stage B contracts one pathway at a time; each conc vector is checkpointed to
#     checkpoints/conc_<arr>.jls, then that pathway's SR chunk blobs are EVICTED so
#     peak disk stays ~one pathway (~3 GB), not all five.
#   * A kill/limit mid-run loses at most the in-progress pathway; re-run resumes.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciIO, Blosc, ZipFile
import Serialization

const HERE     = @__DIR__
const ZARR_URL = "s3://inmap-model/isrm_v1.2.1.zarr/"
const ZIP_LOCAL = get(ENV, "EGU_ZIP", joinpath(HERE, "data", "2016fd_inputs_point.zip"))
const N_SRC    = 52411
const CK_DIR   = joinpath(HERE, "checkpoints"); mkpath(CK_DIR)
const META_ROOT = joinpath(HERE, "cache_meta")   # persistent (geo/pop, tiny)
const SR_ROOT   = joinpath(HERE, "cache_sr")     # evicted between pathways
const OUT       = joinpath(HERE, "run-model_result.jls")

# ---- STEP-0 constants (from isrm.esm / tutorial notebook) ----
const FACT       = 28766.639
const POP_SCALE  = 1.0465819687408728
const MORT_SCALE = 1.025229357798165
const RR_K       = 1.06
const RR_L       = 1.14
const LAT_1 = 33.0; const LAT_2 = 45.0; const LAT_0 = 40.0; const LON_0 = -97.0; const LCC_R = 6370997.0
const SR_ARRAYS = ["SOA","pNO3","pNH4","pSO4","PrimaryPM25"]

# ---- caches: geo/pop persist; SR chunks live in an evictable store ----
meta_cache = EarthSciIO.Cache(; root=META_ROOT)
sr_cache   = EarthSciIO.Cache(; root=SR_ROOT)

# WORKAROUND (zarr reader): the SR arrays have NO `.zattrs` object (live GET → 404);
# _fetch_bytes_optional swallows only CacheMiss, not a 404, so seed an empty `{}`
# `.zattrs` (a legit absent-attrs value) as a cache HIT so the fetch never hits net.
function seed_empty_zattrs(cache, base, arrays)
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
# evict all SR chunk blobs (keep the store dir), then re-seed the tiny .zattrs
function evict_sr!()
    bd = joinpath(SR_ROOT, "v1", "blobs")
    isdir(bd) && rm(bd; recursive=true, force=true)
    seed_empty_zattrs(sr_cache, ZARR_URL, SR_ARRAYS)
end

# ---- Snyder 1987 spherical Lambert Conformal Conic forward (matches lambert_conformal.esm) ----
const D2R = 0.017453292519943295
lcc_t(lat) = tan(0.7853981633974483 + lat*0.008726646259971648)
const LCC_N = log(cos(LAT_1*D2R)/cos(LAT_2*D2R)) / log(lcc_t(LAT_2)/lcc_t(LAT_1))
const LCC_F = cos(LAT_1*D2R) * lcc_t(LAT_1)^LCC_N / LCC_N
lcc_rho(lat) = LCC_R*LCC_F / lcc_t(lat)^LCC_N
const LCC_RHO0 = lcc_rho(LAT_0)
function lcc_forward(lon, lat)
    theta = LCC_N * (lon - LON_0) * D2R
    rho = lcc_rho(lat)
    return (rho*sin(theta), LCC_RHO0 - rho*cos(theta))
end

# ---- pollutant -> pathway map (STEP 0.6); keys are FF10 POLID strings ----
const VOC_SET  = Set(["VOC","VOC_INV","XYL","TOL","TERP","PAR","OLE","NVOL","MEOH","ISOP","IOLE","FORM","ETOH","ETHA","ETH","ALD2","ALDX","CB05_ALD2","CB05_ALDX","CB05_BENZENE","CB05_ETH","CB05_ETHA","CB05_ETOH","CB05_FORM","CB05_IOLE","CB05_ISOP","CB05_MEOH","CB05_OLE","CB05_PAR","CB05_TERP","CB05_TOL","CB05_XYL","ETHANOL","NHTOG","NMOG"])
const PM25_SET = Set(["PM25-PRI","PM2_5","DIESEL-PM25","PAL","PCA","PCL","PEC","PFE","PK","PMG","PMN","PMOTHR","PNH4","PNO3","POC","PSI","PSO4","PTI"])
const NOX_SET  = Set(["NOX","HONO","NO","NO2"])
const NH3_SET  = Set(["NH3"])
const SOX_SET  = Set(["SO2"])
pathway(pol) = pol in VOC_SET ? :VOC : pol in PM25_SET ? :PM25 : pol in NOX_SET ? :NOx :
               pol in NH3_SET ? :NH3 : pol in SOX_SET ? :SOx : :none

# =============================================================================
# STAGE A — ppl pre-pass (geo read + FF10 read + LCC join + emission binning)
# =============================================================================
const STAGE1 = joinpath(CK_DIR, "stage1.jls")

function read_ff10_egu(zippath, member, vars)
    text = ""
    let r = ZipFile.Reader(zippath)
        for f in r.files; f.name == member && (text = String(read(f)); break); end
        close(r)
    end
    lines = split(text, '\n')
    kept = filter(ln -> !startswith(lowercase(strip(ln)), "country_cd"), lines)
    tmp = tempname() * ".csv"
    write(tmp, join(kept, '\n'))
    nds = EarthSciIO.read_native(EarthSciIO.FF10Reader(), tmp; variables=vars)
    rm(tmp; force=true)
    return nds
end

function stage_a()
    if isfile(STAGE1)
        println("== STAGE A: loading checkpoint $STAGE1 ==")
        return Serialization.deserialize(STAGE1)
    end
    println("== STAGE A.1: read W/S/E/N (first $N_SRC) from zarr ==")
    p = EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr", variables=["W","S","E","N"])
    geo = EarthSciIO.materialize(p)
    W  = Float64.(@view geo["W"].data[1:N_SRC]); S = Float64.(@view geo["S"].data[1:N_SRC])
    E  = Float64.(@view geo["E"].data[1:N_SRC]); Nn = Float64.(@view geo["N"].data[1:N_SRC])
    println("  domain X[$(minimum(W)),$(maximum(E))]  Y[$(minimum(S)),$(maximum(Nn))]")

    println("== STAGE A.2: read EGU emissions from FF10 (egu members) ==")
    lon = Float64[]; lat = Float64[]; ann = Float64[]; pol = String[]
    members = String[]
    let r = ZipFile.Reader(ZIP_LOCAL)
        for f in r.files; occursin("egu", lowercase(f.name)) && endswith(f.name, ".csv") && push!(members, f.name); end
        close(r)
    end
    println("  EGU members: ", members)
    for m in members
        nds = read_ff10_egu(ZIP_LOCAL, m, ["POLID","ANN_VALUE","LONGITUDE","LATITUDE"])
        append!(pol, String.(nds["POLID"].data)); append!(ann, Float64.(nds["ANN_VALUE"].data))
        append!(lon, Float64.(nds["LONGITUDE"].data)); append!(lat, Float64.(nds["LATITUDE"].data))
    end
    n_all = length(pol); println("  total EGU FF10 records read: ", n_all)

    paths = pathway.(pol)
    keep = [paths[i] != :none && isfinite(lon[i]) && isfinite(lat[i]) && isfinite(ann[i]) for i in 1:n_all]
    lon = lon[keep]; lat = lat[keep]; ann = ann[keep]; paths = paths[keep]
    N_REC = length(lon)
    println("  recognized-pathway records N_REC = ", N_REC, " (dropped ", n_all-N_REC, ")")
    for pth in (:VOC,:NOx,:NH3,:SOx,:PM25); println("     pathway $pth : ", count(==(pth), paths), " records"); end

    println("== STAGE A.3: LCC project + spatial join (point -> source cell) ==")
    X = similar(lon); Y = similar(lat)
    for i in eachindex(lon); X[i], Y[i] = lcc_forward(lon[i], lat[i]); end
    xmin = minimum(W); xmax = maximum(E); ymin = minimum(S); ymax = maximum(Nn)
    NB = 600; dxb = (xmax - xmin)/NB; dyb = (ymax - ymin)/NB
    bx(x) = clamp(floor(Int,(x - xmin)/dxb), 0, NB-1)
    by(y) = clamp(floor(Int,(y - ymin)/dyb), 0, NB-1)
    bins = [Int[] for _ in 1:NB*NB]; binid(ix,iy) = ix*NB + iy + 1
    for c in 1:N_SRC
        ix0=bx(W[c]); ix1=bx(prevfloat(E[c])); iy0=by(S[c]); iy1=by(prevfloat(Nn[c]))
        for ix in ix0:ix1, iy in iy0:iy1; push!(bins[binid(ix,iy)], c); end
    end
    assign = fill(0, N_REC); noc = 0
    for i in 1:N_REC
        xi=X[i]; yi=Y[i]; found=-1
        for c in bins[binid(bx(xi),by(yi))]
            if W[c] <= xi < E[c] && S[c] <= yi < Nn[c]; found = c-1; break; end
        end
        assign[i] = found; found == -1 && (noc += 1)
    end
    println("  points with NO containing cell: ", noc)
    ppl = sort!(unique(c for c in assign if c >= 0)); N_PPL = length(ppl)
    comp = Dict(ppl[k] => k for k in 1:N_PPL)
    println("  N_PPL (unique emission-bearing source cells) = ", N_PPL)
    chunks = sort!(unique(c ÷ 100 for c in ppl))
    println("  distinct source-chunks (size 100) spanned by ppl: ", length(chunks))

    E_VOC=zeros(N_PPL); E_NOx=zeros(N_PPL); E_NH3=zeros(N_PPL); E_SOx=zeros(N_PPL); E_PM25=zeros(N_PPL)
    for i in 1:N_REC
        assign[i] < 0 && continue
        row = comp[assign[i]]; pth = paths[i]; v = ann[i]
        pth==:VOC && (E_VOC[row]+=v); pth==:NOx && (E_NOx[row]+=v); pth==:NH3 && (E_NH3[row]+=v)
        pth==:SOx && (E_SOx[row]+=v); pth==:PM25 && (E_PM25[row]+=v)
    end
    println("  binned emission sums (tons/yr): VOC=", round(sum(E_VOC)), " NOx=", round(sum(E_NOx)),
            " NH3=", round(sum(E_NH3)), " SOx=", round(sum(E_SOx)), " PM25=", round(sum(E_PM25)))

    st1 = Dict("ppl"=>ppl,"N_PPL"=>N_PPL,"N_REC"=>N_REC,"lon"=>lon,"lat"=>lat,"ann"=>ann,
               "paths"=>paths,"assign"=>assign,
               "E_VOC"=>E_VOC,"E_NOx"=>E_NOx,"E_NH3"=>E_NH3,"E_SOx"=>E_SOx,"E_PM25"=>E_PM25)
    Serialization.serialize(STAGE1, st1); println("  wrote ", STAGE1)
    return st1
end

# =============================================================================
# STAGE B — SR contraction, one pathway at a time, evicting between pathways
# =============================================================================
# conc_p[rcv] = Σ_c SR_p[0, ppl_c, rcv] · E_p[c]  — fetch one 100-source chunk at a
# time (select layer=[0], the ppl cells in that chunk, all rcv).
function contract_pathway(arrname, Ep, ppl, comp, rows_by_chunk, chunk_ids)
    conc = zeros(Float64, N_SRC)
    for (ci, chid) in enumerate(chunk_ids)
        rows = rows_by_chunk[chid]
        src0 = sort!([ppl[k] for k in rows])
        sel = Dict("axes" => Any[ Any[0], src0, "all" ])
        pp = EarthSciIO.const_provider(sr_cache, ZARR_URL; format="zarr",
                 variables=[arrname], reader_kwargs=(select=sel,))
        nd = EarthSciIO.materialize(pp)
        slab = nd[arrname].data
        for (j, s) in enumerate(src0)
            ev = Ep[comp[s]]; ev == 0.0 && continue
            @inbounds @simd for rcv in 1:N_SRC
                conc[rcv] += slab[1, j, rcv] * ev
            end
        end
        if ci % 25 == 0 || ci == length(chunk_ids)
            println("    [$arrname] chunk $ci/$(length(chunk_ids))"); flush(stdout)
        end
    end
    return conc
end

function stage_b(st1)
    ppl = st1["ppl"]; N_PPL = st1["N_PPL"]; comp = Dict(ppl[k] => k for k in 1:N_PPL)
    rows_by_chunk = Dict{Int,Vector{Int}}()
    for k in 1:N_PPL; push!(get!(rows_by_chunk, ppl[k] ÷ 100, Int[]), k); end
    chunk_ids = sort!(collect(keys(rows_by_chunk)))
    PATHS = [("SOA",st1["E_VOC"]),("pNO3",st1["E_NOx"]),("pNH4",st1["E_NH3"]),
             ("pSO4",st1["E_SOx"]),("PrimaryPM25",st1["E_PM25"])]

    println("== STAGE B: SR contraction over $(length(chunk_ids)) chunks x 5 pathways ==")
    concs = Dict{String,Vector{Float64}}()
    for (arrname, Ep) in PATHS
        ck = joinpath(CK_DIR, "conc_$(arrname).jls")
        if isfile(ck)
            println("  [$arrname] checkpoint present -> load $ck")
            concs[arrname] = Serialization.deserialize(ck); continue
        end
        println("  pathway array: ", arrname, "  (emis sum=", round(sum(Ep)), " tons/yr)")
        evict_sr!()   # fresh, empty SR blob store + re-seeded .zattrs for this array
        @time conc = contract_pathway(arrname, Ep, ppl, comp, rows_by_chunk, chunk_ids)
        Serialization.serialize(ck, conc); println("  wrote ", ck)
        concs[arrname] = conc
        evict_sr!()   # free ~3 GB before the next pathway
        flush(stdout)
    end
    return concs
end

# =============================================================================
# STAGE C — health (deaths), STEP 0.5, with [0:52411] positional prefix slice
# =============================================================================
function stage_c(concs)
    TotalPM25 = zeros(Float64, N_SRC)
    pm_components = Dict{String,Float64}()
    for arr in SR_ARRAYS
        c = FACT .* concs[arr]; TotalPM25 .+= c; pm_components[arr] = sum(c)
    end
    println("== STAGE C: health ==")
    println("  per-pathway Σ(pm) over receptors: ", pm_components)
    println("  Σ TotalPM25 = ", sum(TotalPM25), "   max = ", maximum(TotalPM25))

    p2 = EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr", variables=["TotalPop","MortalityRate"])
    h = EarthSciIO.materialize(p2)
    TotalPop      = Float64.(@view h["TotalPop"].data[1:N_SRC])
    MortalityRate = Float64.(@view h["MortalityRate"].data[1:N_SRC])
    deathsK = (exp.(log(RR_K)/10 .* TotalPM25) .- 1) .* TotalPop .* POP_SCALE .* MortalityRate ./ 100000 .* MORT_SCALE
    deathsL = (exp.(log(RR_L)/10 .* TotalPM25) .- 1) .* TotalPop .* POP_SCALE .* MortalityRate ./ 100000 .* MORT_SCALE
    sK = sum(deathsK); sL = sum(deathsL)
    println("  sum(deathsK) = ", sK, "   (target ~7524.84 Krewski)")
    println("  sum(deathsL) = ", sL, "   (target ~16979.45 Lepeule)")
    return (TotalPM25=TotalPM25, deathsK=deathsK, deathsL=deathsL, sum_deathsK=sK, sum_deathsL=sL)
end

# =============================================================================
function main()
    t0 = time()
    st1   = stage_a()
    concs = stage_b(st1)
    res   = stage_c(concs)
    Serialization.serialize(OUT, Dict(
        "ppl"=>st1["ppl"], "N_PPL"=>st1["N_PPL"], "N_REC"=>st1["N_REC"],
        "E_VOC"=>st1["E_VOC"],"E_NOx"=>st1["E_NOx"],"E_NH3"=>st1["E_NH3"],"E_SOx"=>st1["E_SOx"],"E_PM25"=>st1["E_PM25"],
        "TotalPM25"=>res.TotalPM25, "deathsK"=>res.deathsK, "deathsL"=>res.deathsL,
        "sum_deathsK"=>res.sum_deathsK, "sum_deathsL"=>res.sum_deathsL,
    ))
    println("wrote ", OUT)
    println(repeat("=", 70))
    println("PART (A) COMPLETE in $(round((time()-t0)/60, digits=1)) min")
    println("  sum(deathsK) = ", round(res.sum_deathsK, digits=2), "   target 7524.84   ",
            "rel.err ", round(100*(res.sum_deathsK-7524.84)/7524.84, digits=3), "%")
    println("  sum(deathsL) = ", round(res.sum_deathsL, digits=2), "   target 16979.45  ",
            "rel.err ", round(100*(res.sum_deathsL-16979.45)/16979.45, digits=3), "%")
    println(repeat("=", 70))
end
main()
