#!/usr/bin/env julia
# L3 pre-flight (cheap, NO SR fetch): does the pushdown value-invention producer
# derive the SAME 1520 ppl run-model.jl validated, at FULL scale
# (43650 emissions × 52411 pop cells), and how long does the overlap-gated
# producer take?  De-risks the enumeration-scaling question before any S3 SR
# fetch / disk commitment.
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc                          # activate the blosc decode extension
import GeometryOps, GeoInterface     # activate the STRtree broad-phase fast path
import JSON
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))

const MODEL = "/Users/ctessum/code/earthsciml/isrm.esm-wt-runner/isrm_pushdown.esm"
const ISRM_DIR = dirname(MODEL)

println("building full-scale inputs (43650 emissions, LCC-projected) ...")
inp = build_inputs()
println("  N_REC=", inp.N_REC, "  N_SRC=", inp.N_SRC,
        "  X∈[", round(minimum(inp.X)), ",", round(maximum(inp.X)), "]",
        "  Y∈[", round(minimum(inp.Y)), ",", round(maximum(inp.Y)), "]")

doc = JSON.parsefile(MODEL)
f = EA.load(doc; base_path=ISRM_DIR,
    metaparameters=Dict("N_SRC"=>inp.N_SRC, "N_RCV"=>inp.N_SRC,
                        "N_POP"=>inp.N_SRC, "N_LAYER"=>3, "N_REC"=>inp.N_REC))
model = EA._select_model(f, "ISRM")

ca = Dict{String,Any}("X"=>inp.X, "Y"=>inp.Y,
                      "W"=>inp.W, "S"=>inp.S, "E"=>inp.E, "N"=>inp.N)

println("running materialize_value_invention (overlap-gated producer) ...")
GC.gc()
t = @elapsed vi = EA.materialize_value_invention(model, f.index_sets, ca,
                                                 Dict{String,Any}())
faq = "emis_src_cells_faq"
mem = sort(collect(Int, vi.members[faq]))
ext = vi.extents[faq]
expected = sort(inp.ppl0 .+ 1)          # run-model's 0-based ppl → 1-based members

println("\n", "="^70)
println("VI PROBE RESULT")
println("  producer materialize time : ", round(t, digits=1), " s")
println("  invented |members|        : ", length(mem), "   (extent=", ext, ")")
println("  run-model.jl ppl count    : ", length(expected))
println("  members == run-model ppl+1: ", mem == expected)
if mem != expected
    only_vi  = setdiff(mem, expected); only_rm = setdiff(expected, mem)
    println("  in VI not run-model (", length(only_vi), "): ", first(only_vi, 10))
    println("  in run-model not VI (", length(only_rm), "): ", first(only_rm, 10))
end
println("  member range              : [", minimum(mem), ",", maximum(mem), "]")
println("="^70)
mem == expected || error("VI members DIFFER from validated run-model ppl")
println("VI PROBE PASS")
