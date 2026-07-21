#!/usr/bin/env julia
# plot.jl — renders the isrm.esm example plot (`isrm_demo` → `totalpm25_map`):
# a scatter of TotalPM25 (ug/m3) vs receptor cell-centroid easting (m, LCC).
# The 52411-cell ISRM grid is an unstructured variable-resolution polygon mesh
# with a single flat index, so a structured field_snapshot / true polygon
# choropleth is out of format scope (gap G5); the authored example uses scatter.
# A second panel shows attributable deaths (Krewski) the same way.
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciIO, Blosc, Serialization, Plots

const HERE = @__DIR__
const N_SRC = 52411
const ZARR_URL = "s3://inmap-model/isrm_v1.2.1.zarr/"
res = Serialization.deserialize(joinpath(HERE, "run-model_result.jls"))
TotalPM25 = res["TotalPM25"]::Vector{Float64}
deathsK   = res["deathsK"]::Vector{Float64}

# receptor centroid easting rcv_cx = (W+E)/2 for the 52411 receptor cells (== the
# source-cell prefix). Read W/E from the persistent meta cache (already fetched).
meta_cache = EarthSciIO.Cache(; root=joinpath(HERE, "cache_meta"))
geo = EarthSciIO.materialize(EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr", variables=["W","E"]))
W = Float64.(@view geo["W"].data[1:N_SRC]); E = Float64.(@view geo["E"].data[1:N_SRC])
rcv_cx = (W .+ E) ./ 2

gr()
p1 = scatter(rcv_cx, TotalPM25; ms=1, msw=0, alpha=0.25, legend=false,
    xlabel="receptor centroid x (m, LCC)", ylabel="PM2.5 (ug/m3)",
    title="isrm_demo / totalpm25_map  (ΣTotalPM25=$(round(sum(TotalPM25),digits=1)))")
p2 = scatter(rcv_cx, deathsK; ms=1, msw=0, alpha=0.25, legend=false, color=:firebrick,
    xlabel="receptor centroid x (m, LCC)", ylabel="attributable deaths (Krewski)",
    title="deathsK by cell  (Σ=$(round(sum(deathsK),digits=1)))")
plt = plot(p1, p2; layout=(2,1), size=(1000,800), dpi=130)
out = joinpath(HERE, "isrm_totalpm25.png")
savefig(plt, out)
println("wrote ", out)
