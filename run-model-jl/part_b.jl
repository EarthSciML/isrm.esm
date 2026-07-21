#!/usr/bin/env julia
# =============================================================================
# part_b.jl — END-TO-END through the runtime: prove that EarthSciAST.jl LOADS and
# EXECUTES the authored isrm.esm (not a hand-port of its math) and evaluates its
# observeds correctly.
#
# The authored observed graph is:
#     X,Y            = LCC(emis_lon, emis_lat)                    [emis_records]
#     is_p           = pollutant-code masks                       [emis_records]
#     cell_{W,S,E,N} = {W,S,E,N}[ppl]                             [emis_src_cells]
#     E_p[c]         = Σ_r [cell in rect] · emis_annual · is_p    [emis_src_cells]
#     conc_p[rcv]    = Σ_c SR_p[c,rcv] · E_p[c]                    [rcv_cells]
#     TotalPM25[rcv] = fact · Σ_p conc_p[rcv]                     [rcv_cells]
#     deathsK[rcv]   = (exp(ln rr_K/10·PM)-1)·pop·pop_scale·mort·mort_scale/1e5
#
# We drive it with tiny CONTROLLED inputs through EA.load → EA.prepare →
# EA._observed_field (the framework's own state-free-observed evaluator, the same
# machinery §6.6.5 inline assertions use), then check every stage against a plain-
# Julia ORACLE of the identical STEP-0 formulas (the validated run-model.jl math).
# Match ⇒ the runtime faithfully executes the .esm.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
import JSON
const EA = EarthSciAST
const MODEL = joinpath(dirname(@__DIR__), "isrm.esm")
const ISRM_DIR = dirname(@__DIR__)

# ---- tiny controlled problem -------------------------------------------------
const N_PPL = 2      # emission source cells (compact axis)
const N_RCV = 3      # receptor cells
const N_REC = 5      # emission records (one per pathway)
const N_POP = 3      # pop grid length (== N_RCV so the [1:N_RCV] prefix is all)
const N_SRC = 3

# STEP-0 scalars (== isrm.esm defaults)
const FACT=28766.639; const POP_SCALE=1.0465819687408728; const MORT_SCALE=1.025229357798165
const RR_K=1.06; const RR_L=1.14
const LAT_1=33.0; const LAT_2=45.0; const LAT_0=40.0; const LON_0=-97.0; const LCC_R=6370997.0

# Snyder spherical LCC forward (the oracle's copy of lambert_conformal_forward_x/y)
const D2R=0.017453292519943295
lcc_t(lat)=tan(0.7853981633974483+lat*0.008726646259971648)
const LCC_N=log(cos(LAT_1*D2R)/cos(LAT_2*D2R))/log(lcc_t(LAT_2)/lcc_t(LAT_1))
const LCC_F=cos(LAT_1*D2R)*lcc_t(LAT_1)^LCC_N/LCC_N
lcc_rho(lat)=LCC_R*LCC_F/lcc_t(lat)^LCC_N
const LCC_RHO0=lcc_rho(LAT_0)
lcc_fwd(lon,lat)=(lcc_rho(lat)*sin(LCC_N*(lon-LON_0)*D2R),
                  LCC_RHO0-lcc_rho(lat)*cos(LCC_N*(lon-LON_0)*D2R))

# emission points (lon,lat), one per pathway; pollutant code picks the pathway
emis_lon = [-90.0, -88.0, -95.0, -85.0, -92.0]
emis_lat = [ 40.0,  41.0,  38.0,  43.0,  39.0]
pollutant = Float64[1, 36, 40, 41, 42]      # VOC, NOx, NH3, SOx, PM25
emis_annual = [10.0, 20.0, 30.0, 40.0, 50.0]  # kg/yr

# project the points, then build 2 non-overlapping source-cell rects that split
# them by easting (points well inside → robust to any LCC rounding)
XY = [lcc_fwd(emis_lon[r], emis_lat[r]) for r in 1:N_REC]
X = [p[1] for p in XY]; Y = [p[2] for p in XY]
xmid = (minimum(X)+maximum(X))/2; PAD = 5.0e5
# pop cells 1,2 are the two source cells (ppl=[1,2]); cell 3 is an unused extra
W = [minimum(X)-PAD, xmid,          9.9e9]
E = [xmid,           maximum(X)+PAD, 9.9e9+1]
S = [minimum(Y)-PAD, minimum(Y)-PAD, 0.0]
Nn= [maximum(Y)+PAD, maximum(Y)+PAD, 1.0]
ppl = Float64[1, 2]                          # 1-based indices into W/S/E/N

# SR matrices (emis_src_cells × rcv_cells), distinct per pathway
SR = Dict(
  "SOA"        => [1e-6 2e-6 3e-6; 4e-6 5e-6 6e-6],
  "pNO3"       => [2e-6 1e-6 4e-6; 3e-6 6e-6 5e-6],
  "pNH4"       => [5e-6 4e-6 1e-6; 2e-6 3e-6 6e-6],
  "pSO4"       => [3e-6 6e-6 2e-6; 1e-6 4e-6 5e-6],
  "PrimaryPM25"=> [6e-6 5e-6 4e-6; 3e-6 2e-6 1e-6])
TotalPop      = [100.0, 200.0, 300.0]        # pop_cells (prefix [1:3] = receptors)
MortalityRate = [500.0, 600.0, 700.0]

# =============================================================================
# ORACLE — the validated run-model.jl STEP-0 math on these same inputs
# =============================================================================
is_VOC(p)  = 1  <= p <= 35
is_NOx(p)  = 36 <= p <= 39
is_NH3(p)  = p == 40
is_SOx(p)  = p == 41
is_PM25(p) = 42 <= p <= 59
maskfun = Dict("SOA"=>is_VOC,"pNO3"=>is_NOx,"pNH4"=>is_NH3,"pSO4"=>is_SOx,"PrimaryPM25"=>is_PM25)

cell_W = [W[Int(ppl[c])] for c in 1:N_PPL]; cell_E = [E[Int(ppl[c])] for c in 1:N_PPL]
cell_S = [S[Int(ppl[c])] for c in 1:N_PPL]; cell_N = [Nn[Int(ppl[c])] for c in 1:N_PPL]
contains(c,r) = cell_W[c] <= X[r] < cell_E[c] && cell_S[c] <= Y[r] < cell_N[c]

function oracle_E(pathway)
    m = maskfun[pathway]
    [sum((contains(c,r) ? 1.0 : 0.0) * emis_annual[r] * (m(pollutant[r]) ? 1.0 : 0.0)
         for r in 1:N_REC) for c in 1:N_PPL]
end
oracle_conc(pathway) = [sum(SR[pathway][c,rcv]*oracle_E(pathway)[c] for c in 1:N_PPL) for rcv in 1:N_RCV]
oracle_TotalPM25 = [FACT*sum(oracle_conc(p)[rcv] for p in keys(SR)) for rcv in 1:N_RCV]
oracle_deaths(rr) = [(exp(log(rr)/10*oracle_TotalPM25[rcv])-1)*TotalPop[rcv]*POP_SCALE*
                     MortalityRate[rcv]*MORT_SCALE/100000 for rcv in 1:N_RCV]

# =============================================================================
# RUNTIME — load + prepare the authored .esm, evaluate its observeds
# =============================================================================
doc = JSON.parsefile(MODEL)
f = EA.load(doc; base_path=ISRM_DIR,
            metaparameters=Dict("N_SRC"=>N_SRC,"N_RCV"=>N_RCV,"N_POP"=>N_POP,
                                "N_REC"=>N_REC,"N_PPL"=>N_PPL))
ca = Dict{String,Any}(
  "ISRM_SR.SOA"=>SR["SOA"], "ISRM_SR.pNO3"=>SR["pNO3"], "ISRM_SR.pNH4"=>SR["pNH4"],
  "ISRM_SR.pSO4"=>SR["pSO4"], "ISRM_SR.PrimaryPM25"=>SR["PrimaryPM25"],
  "ISRM_SR.W"=>W, "ISRM_SR.S"=>S, "ISRM_SR.E"=>E, "ISRM_SR.N"=>Nn,
  "ISRM_SR.TotalPop"=>TotalPop, "ISRM_SR.MortalityRate"=>MortalityRate,
  "EGU_Emis.lon"=>emis_lon, "EGU_Emis.lat"=>emis_lat, "EGU_Emis.annual"=>emis_annual,
  "EGU_Emis.pollutant"=>pollutant,
  "EGU_Emis.stkhgt"=>zeros(N_REC), "EGU_Emis.stkdiam"=>zeros(N_REC),
  "EGU_Emis.stktemp"=>zeros(N_REC), "EGU_Emis.stkvel"=>zeros(N_REC),
  "ISRM.ppl"=>ppl)

insp = EA.BuildInspection()
prep = EA.prepare(f; const_arrays=ca, inspect=insp)
println("PREPARE OK: ", prep)

rt(var) = (fld = EA._observed_field(insp, f, "ISRM", var);
           fld === nothing ? error("observed $var not evaluable") : fld[1])

checks = [
  ("E_VOC",     rt("E_VOC"),     oracle_E("SOA")),
  ("E_NOx",     rt("E_NOx"),     oracle_E("pNO3")),
  ("E_PM25",    rt("E_PM25"),    oracle_E("PrimaryPM25")),
  ("conc_SOA",  rt("conc_SOA"),  oracle_conc("SOA")),
  ("conc_pSO4", rt("conc_pSO4"), oracle_conc("pSO4")),
  ("TotalPM25", rt("TotalPM25"), oracle_TotalPM25),
  ("deathsK",   rt("deathsK"),   oracle_deaths(RR_K)),
  ("deathsL",   rt("deathsL"),   oracle_deaths(RR_L)),
]
function report(checks)
    println("\n", rpad("observed",12), rpad("runtime (EA.esm)",30), rpad("oracle (STEP-0)",30), "max|Δ|")
    allok = true
    for (name, got, exp) in checks
        d = maximum(abs.(got .- exp))
        ok = d < 1e-9 * max(1.0, maximum(abs.(exp)))
        allok &= ok
        println(rpad(name,12), rpad(string(round.(got,sigdigits=6)),30),
                rpad(string(round.(exp,sigdigits=6)),30), d, ok ? "" : "  <-- MISMATCH")
    end
    return allok
end
allok = report(checks)
println("\n", allok ? "PART (B) PASS: runtime executes isrm.esm == STEP-0 oracle on all stages" :
                      "PART (B) FAIL: see mismatches above")
