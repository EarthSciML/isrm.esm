#!/usr/bin/env julia
# L3 pre-flight #2 (offline, fast): drive the REAL authored isrm_pushdown.esm
# through EA.prepare with a MOCK gated SR provider (no S3), on a tiny controlled
# problem, and check deaths against a plain-Julia STEP-0 oracle.  De-risks
# everything EXCEPT the live-S3 fetch: the X/Y-observed-vs-const-supply question,
# the gated deferral + selective fetch wiring, the members-fed factor, the
# applies_to→model-var (SR_SOA) mapping, and the full downstream graph.
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
import GeometryOps, GeoInterface
import JSON
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

# ---- MOCK gated SR provider: records calls, slices synthetic SR per selection.
# Keyed by MODEL variable names (SR_SOA…) so applies_to resolves via
# _const_factor_aliases to the SR_* parameters.
mutable struct MockSR
    full::Dict{String,Array{Float64,3}}
    gate::Dict{String,Any}
    calls::Vector{Any}
end
EA.provider_gate_spec(m::MockSR) = m.gate
EA.provider_supports_selection(m::MockSR) = true
EA.provider_refresh_times(m::MockSR) = Float64[]
function EA.provider_sample(m::MockSR, ::Real; selection=nothing)
    if selection === nothing
        push!(m.calls, (:wholesale,))
        return Dict{String,Any}(k => v for (k, v) in m.full)
    end
    push!(m.calls, (:selection, deepcopy(selection)))
    lay, src, rcv = selection[1], selection[2], selection[3]
    return Dict{String,Any}(k => v[lay:lay, src, rcv] for (k, v) in m.full)
end

# ---- tiny controlled problem -----------------------------------------------
NSRC = 3; N_RCV = 3; N_POP = 3; N_REC = 5; N_LAYER = 3
# 5 emissions, one per pathway (codes 1 VOC, 36 NOx, 40 NH3, 41 SOx, 42 PM25)
emis_lon = [-90.0, -88.0, -95.0, -85.0, -92.0]
emis_lat = [ 40.0,  41.0,  38.0,  43.0,  39.0]
pollutant = [1.0, 36.0, 40.0, 41.0, 42.0]
emis_annual = [10.0, 20.0, 30.0, 40.0, 50.0]
XY = [lcc_forward(emis_lon[r], emis_lat[r]) for r in 1:N_REC]
X = [p[1] for p in XY]; Y = [p[2] for p in XY]
# two source cells split the points by easting; cell 3 unused (ppl = {1,2})
xmid = (minimum(X)+maximum(X))/2; PAD = 5.0e5
W = [minimum(X)-PAD, xmid,          9.9e9]
E = [xmid,           maximum(X)+PAD, 9.9e9+1]
S = [minimum(Y)-PAD, minimum(Y)-PAD, 0.0]
Nn= [maximum(Y)+PAD, maximum(Y)+PAD, 1.0]
TotalPop      = [100.0, 200.0, 300.0]
MortalityRate = [500.0, 600.0, 700.0]
HAND_MEMBERS = [1, 2]

FACT=28766.639; POP_SCALE=1.0465819687408728; MORT_SCALE=1.025229357798165
RR_K=1.06; RR_L=1.14

# synthetic full SR per pathway [layer, src, rcv]; layer-0 (1-based 1) used
SRV = Dict("SR_SOA"=>1.0,"SR_pNO3"=>2.0,"SR_pNH4"=>3.0,"SR_pSO4"=>4.0,"SR_PrimaryPM25"=>5.0)
fullSR = Dict{String,Array{Float64,3}}()
for (nm, b) in SRV
    A = Array{Float64}(undef, N_LAYER, NSRC, N_RCV)
    for l in 1:N_LAYER, s in 1:NSRC, r in 1:N_RCV
        A[l,s,r] = ((l-1)*1e6 + b*1000 + s*10 + r) * 1e-9   # realistic ~1e-6 SR magnitude
    end
    fullSR[nm] = A
end

# ---- plain-Julia STEP-0 ORACLE ---------------------------------------------
is_VOC(p)=1<=p<=35; is_NOx(p)=36<=p<=39; is_NH3(p)=p==40; is_SOx(p)=p==41; is_PM25(p)=42<=p<=59
maskfun = Dict("SR_SOA"=>is_VOC,"SR_pNO3"=>is_NOx,"SR_pNH4"=>is_NH3,"SR_pSO4"=>is_SOx,"SR_PrimaryPM25"=>is_PM25)
NP = length(HAND_MEMBERS)
cell_W=[W[HAND_MEMBERS[c]] for c in 1:NP]; cell_E=[E[HAND_MEMBERS[c]] for c in 1:NP]
cell_S=[S[HAND_MEMBERS[c]] for c in 1:NP]; cell_N=[Nn[HAND_MEMBERS[c]] for c in 1:NP]
contains(c,r)= cell_W[c]<=X[r]<cell_E[c] && cell_S[c]<=Y[r]<cell_N[c]
oracle_E(nm)=(m=maskfun[nm]; [sum((contains(c,r) ? 1.0 : 0.0)*emis_annual[r]*(m(pollutant[r]) ? 1.0 : 0.0) for r in 1:N_REC) for c in 1:NP])
srC(nm)=[fullSR[nm][1,HAND_MEMBERS[c],rcv] for c in 1:NP, rcv in 1:N_RCV]
oracle_conc(nm)=(Ep=oracle_E(nm); [sum(srC(nm)[c,rcv]*Ep[c] for c in 1:NP) for rcv in 1:N_RCV])
oracle_TotalPM25=[FACT*sum(oracle_conc(nm)[rcv] for nm in keys(SRV)) for rcv in 1:N_RCV]
oracle_deaths(rr)=[(exp(log(rr)/10*oracle_TotalPM25[rcv])-1)*TotalPop[rcv]*POP_SCALE*MortalityRate[rcv]*MORT_SCALE/100000 for rcv in 1:N_RCV]

# ---- load + build the REAL isrm_pushdown.esm via the DICT FRONT-DOOR --------
# VI runs only through the AbstractDict front-door; `prepare` flattens the coupled
# esm and breaks join-env resolution.  So we drive build_evaluator on a
# metaparameter-resolved doc with BARE const arrays + _gated_providers — exactly
# the committed L1 milestone's primary build path.
doc0 = JSON.parsefile(MODEL)
mp = Dict("N_SRC"=>NSRC,"N_RCV"=>N_RCV,"N_POP"=>N_POP,"N_LAYER"=>N_LAYER,"N_REC"=>N_REC)
doc = resolve_sizes!(deepcopy(doc0), mp)
f = EA.load(deepcopy(doc); base_path=ISRM_DIR)

gate = Dict{String,Any}(
    "axes" => Any[Dict("fixed"=>[0]), Dict("gated_by"=>"emis_src_cells"), "all"],
    "applies_to" => ["SR_SOA","SR_pNO3","SR_pNH4","SR_pSO4","SR_PrimaryPM25"])
mock = MockSR(fullSR, gate, Any[])

ca = Dict{String,Any}(
    "X"=>X, "Y"=>Y, "emis_lon"=>emis_lon, "emis_lat"=>emis_lat,
    "emis_annual"=>emis_annual, "pollutant"=>pollutant,
    "stkhgt"=>zeros(N_REC), "stkdiam"=>zeros(N_REC),
    "stktemp"=>zeros(N_REC), "stkvel"=>zeros(N_REC),
    "W"=>W, "S"=>S, "E"=>E, "N"=>Nn,
    "TotalPop"=>TotalPop, "MortalityRate"=>MortalityRate)

insp = EA.BuildInspection()
EA.build_evaluator(doc; model_name="ISRM", const_arrays=ca, inspect=insp,
    _gated_providers=Dict{String,Any}("ISRM_SR"=>mock), _sample_time=0.0)
println("BUILD OK")

# ---- deferral + selection assertions ---------------------------------------
whole = count(c->c[1]==:wholesale, mock.calls)
selc  = [c for c in mock.calls if c[1]==:selection]
println("wholesale calls (must be 0): ", whole)
println("selection calls: ", length(selc))
sel = selc[1][2]
println("  pushed selection: layer=", sel[1], "  members=", sel[2], "  rcv=", typeof(sel[3]))
mem_ok = sel[2] == HAND_MEMBERS
println("  members == hand {1,2}: ", mem_ok)
sc = get(insp.const_arrays, "ISRM.src_cell_of_ppl", get(insp.const_arrays,"src_cell_of_ppl",nothing))
println("  src_cell_of_ppl fed back: ", sc === nothing ? "MISSING" : Vector{Float64}(sc))

# ---- evaluate deaths through the runtime, compare to oracle ----------------
rt(v)=(fld=EA._observed_field(insp, f, "ISRM", v); fld===nothing ? error("no $v") : fld[1])
dK = rt("deathsK"); dL = rt("deathsL"); tp = rt("TotalPM25")
mdK = maximum(abs.(dK .- oracle_deaths(RR_K)))
mdL = maximum(abs.(dL .- oracle_deaths(RR_L)))
mdT = maximum(abs.(tp .- oracle_TotalPM25))
println("\nmax|Δ| TotalPM25=", mdT, "  deathsK=", mdK, "  deathsL=", mdL)
tol = 1e-6*max(1.0, maximum(abs.(oracle_deaths(RR_L))))
pass = whole==0 && mem_ok && sc!==nothing && mdK<tol && mdL<tol && mdT<1e-6*max(1.0,maximum(abs.(oracle_TotalPM25)))
println("\n", "="^60)
println("L3 OFFLINE GRAPH: ", pass ? "PASS" : "FAIL")
println("  sum(deathsK)=", sum(dK), "  sum(deathsL)=", sum(dL))
println("="^60)
pass || error("offline graph test FAILED")
