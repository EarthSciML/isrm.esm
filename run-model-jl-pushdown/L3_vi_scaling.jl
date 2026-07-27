#!/usr/bin/env julia
# Characterize the value-invention producer scaling: the producer ranges
# r∈emis_records × c∈pop_cells, and _vi_enumerate walks the FULL Cartesian
# product (value_invention.jl rec(1)), gating each tuple by the broad-phase
# candidate set — so cost ≈ N_REC × N_POP.  Time VI at increasing N_REC (first-N
# emissions, full 52411 cells) to confirm linear-in-N_REC scaling and extrapolate
# the full-scale (43650) cost.  Also verify members are correct + STRtree loaded.
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc
import GeometryOps, GeoInterface
import JSON
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

println("STRtree fast path available: ",
        hasmethod(EA.build_spatial_index, Tuple{AbstractVector}))

inp = build_inputs()   # full 43650 emissions + full-grid 52411 cells
# run-model's 0-based assignment of each emission → containing cell (or -1)
st = Serialization.deserialize(joinpath(RUNMODEL, "checkpoints", "stage1.jls"))
assign = st["assign"]   # length 43650, 0-based cell id or -1

doc0 = JSON.parsefile(MODEL)

function time_vi(nrec)
    doc = deepcopy(doc0)
    resolve_sizes!(doc, Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,
                             "N_LAYER"=>3,"N_REC"=>nrec))
    f = EA.load(doc; base_path=dirname(MODEL))
    model = EA._select_model(f, "ISRM")
    ca = Dict{String,Any}("X"=>inp.X[1:nrec], "Y"=>inp.Y[1:nrec],
                          "W"=>inp.W, "S"=>inp.S, "E"=>inp.E, "N"=>inp.N)
    GC.gc()
    t = @elapsed vi = EA.materialize_value_invention(model, f.index_sets, ca, Dict{String,Any}())
    mem = sort(collect(Int, vi.members["emis_src_cells_faq"]))
    expect = sort(unique(assign[i]+1 for i in 1:nrec if assign[i] >= 0))
    return t, mem, expect
end

# warm up compilation on a tiny slice first (so the first timing isn't dominated
# by JIT), then measure the scaling.
print("warmup (nrec=20) ... "); flush(stdout)
tw, _, _ = time_vi(20); println(round(tw, digits=1), " s (incl. compile)")

println("\n nrec |  VI time (s) | pairs=nrec*52411 | |members| | members==oracle")
println("-"^72)
prev_t = 0.0; prev_n = 0
for nrec in (100, 400, 1600, 6400)
    t, mem, expect = time_vi(nrec)
    ok = mem == expect
    rate = nrec*52411 / t / 1e6
    println(lpad(nrec,5), " | ", lpad(round(t,digits=2),11), " | ",
            lpad(nrec*52411, 16), " | ", lpad(length(mem),8), " | ", ok,
            "   (", round(rate,digits=1), " Mpair/s)")
    global prev_t = t; global prev_n = nrec
end
# extrapolate to the full 43650 from the last measured rate
full_est = prev_t * 43650 / prev_n
println("-"^72)
println("linear extrapolation to full N_REC=43650: ≈ ", round(full_est/60, digits=1),
        " min for VI producer alone (enumerates 43650×52411 = ",
        round(43650*52411/1e9, digits=2), "e9 tuples)")
