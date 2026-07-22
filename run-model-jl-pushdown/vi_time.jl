#!/usr/bin/env julia
# Clean VI-ONLY timing (no S3 SR fetch, no build): time the value-invention
# PRODUCER materialization at increasing N_REC over the full 52411-cell grid, to
# quantify the emis_records × pop_cells full-product enumeration cost.
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc
import GeometryOps, GeoInterface
import JSON
const EA = EarthSciAST
include(joinpath(@__DIR__, "l3_common.jl"))
const MODEL = "/Users/ctessum/code/earthsciml/isrm.esm-wt-runner/isrm_pushdown.esm"

inp = build_inputs()
doc0 = JSON.parsefile(MODEL)
function vi_time(nrec)
    doc = resolve_sizes!(deepcopy(doc0),
        Dict("N_SRC"=>inp.N_SRC,"N_RCV"=>inp.N_SRC,"N_POP"=>inp.N_SRC,"N_LAYER"=>3,"N_REC"=>nrec))
    f = EA.load(doc; base_path=dirname(MODEL))
    model = EA._select_model(f, "ISRM")
    ca = Dict{String,Any}("X"=>inp.X[1:nrec],"Y"=>inp.Y[1:nrec],
                          "W"=>inp.W,"S"=>inp.S,"E"=>inp.E,"N"=>inp.N)
    GC.gc()
    t = @elapsed vi = EA.materialize_value_invention(model, f.index_sets, ca, Dict{String,Any}())
    return t, length(vi.members["emis_src_cells_faq"])
end
print("warmup nrec=10 (compiles STRtree+VI) ... "); flush(stdout)
tw,_ = vi_time(10); println(round(tw,digits=1), " s (one-time compile incl.)"); flush(stdout)
println("\n nrec | VI time (s) | pairs(nrec*52411) | |members|"); flush(stdout)
for nrec in (20, 100, 300)
    t, m = vi_time(nrec)
    println(lpad(nrec,5), " | ", lpad(round(t,digits=3),11), " | ", lpad(nrec*52411,17), " | ", m)
    flush(stdout)
end
# extrapolate from nrec=300 to full 43650 (enumeration is linear in nrec)
t300,_ = vi_time(300)
println("\nextrapolation: VI(43650) ≈ ", round(t300*43650/300/60,digits=1),
        " min  (full product = ", round(43650*52411/1e9,digits=2), "e9 tuples)"); flush(stdout)
