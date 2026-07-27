#!/usr/bin/env julia
# =============================================================================
# mem_probe_fetch.jl — where does the memory go during the gated SR fetch?
#
# The full L3 run was OOM-killed at ~24 GB anon-RSS by a 40 GB cgroup, DURING
# the fetch/build phase (the log never reached "BUILD done"). The five SR slabs
# are only ~637 MB each (3.2 GB total), so the footprint is elsewhere. This
# isolates ONE pathway's fetch and reports RSS around it.
#
# Hypothesis under test: EarthSciIO's zarr reader holds all 416 decompressed
# chunk buffers simultaneously while assembling the (1, |ppl|, 52411) output,
# rather than scattering each chunk and freeing it — ~8.7 GB per pathway.
#
# PPL=<n> limits how many ppl rows to request (default 1520 = full).
# Prints RSS at each step so the growth curve is visible even if killed.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciIO, Blosc
import Serialization
include(joinpath(@__DIR__, "paths.jl"))

rss_gb() = parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^30
step(label) = (GC.gc(); println(rpad(label, 42), " RSS = ", round(rss_gb(), digits=2), " GB"); flush(stdout))

const NPPL = parse(Int, get(ENV, "PPL", "1520"))
const ARR  = get(ENV, "ARR", "SOA")
const ROOT = joinpath(SCRATCH, "memprobe_sr")

step("start")

# the validated ppl support from run-model.jl's checkpoint (0-based -> 0-based ids)
st  = Serialization.deserialize(joinpath(RUNMODEL, "checkpoints", "stage1.jls"))
ppl0 = sort(Int.(st["ppl"]))[1:min(NPPL, length(st["ppl"]))]
nchunk = length(unique(ppl0 .÷ 100))
println("requesting $(length(ppl0)) ppl rows spanning $nchunk source-chunks of array $ARR")
step("after loading ppl checkpoint")

isdir(ROOT) && rm(ROOT; recursive=true, force=true)
cache = EarthSciIO.Cache(; root=ROOT)
let base = rstrip(ZARR_URL, '/')                     # seed the absent .zattrs
    u = "$base/$ARR/.zattrs"; k = EarthSciIO.cache_key(u)
    if EarthSciIO.get_blob(cache.store, k) === nothing
        tmp = EarthSciIO.staging_path(cache.store); write(tmp, "{}")
        EarthSciIO.put_blob!(cache.store, k, tmp)
    end
end
prov = EarthSciIO.const_provider(cache, ZARR_URL; format="zarr", variables=[ARR])
sel  = Dict("axes" => Any[Dict("indices" => [0]), Dict("indices" => ppl0), "all"])
step("after provider construction")

println("materializing SR[$ARR][0, $(length(ppl0)) rows, :] ..."); flush(stdout)
t = @elapsed (nd = EarthSciIO.materialize(prov; select=sel))
raw = nd[ARR].data
step("after materialize (raw, $(eltype(raw)))")
println("  materialize took ", round(t, digits=1), " s   size=", size(raw),
        "   nominal = ", round(prod(size(raw)) * sizeof(eltype(raw)) / 2^30, digits=2), " GB")

slab = Array{Float64}(raw)
step("after Array{Float64} conversion")
println("  Float64 slab nominal = ", round(prod(size(slab)) * 8 / 2^30, digits=2), " GB")

nd = nothing; raw = nothing
step("after dropping the raw NativeDataset")

println("\nSUMMARY: peak-ish RSS above vs 3.2 GB of Float64 slabs for all 5 pathways.")
println("If RSS >> nominal, the reader is retaining per-chunk buffers.")
rm(ROOT; recursive=true, force=true)

# -----------------------------------------------------------------------------
# MEASURED (2026-07-27, 40 GB SLURM-step cgroup, ISRM SR on s3://inmap-model):
#
#   after provider construction              RSS = 0.43 GB
#   after materialize (1520 rows, 416 chunks) RSS = 9.17 GB   <- 15x the result
#   nominal Float64 slab                            0.59 GB
#   after GC + dropping the NativeDataset    RSS = 8.68 GB
#
#   Two pathways in ONE process, GC.gc() between them:
#     after SOA   RSS = 9.12 GB
#     after pNO3  RSS = 9.53 GB              <- PLATEAUS: reclaimable, not a leak
#
# CONCLUSION: the churn is real (~8.7 GB of decompressed chunk buffers per
# pathway) but fully reclaimable. The full run's OOM was NOT retention — Julia
# sizes its GC heap from TOTAL SYSTEM RAM (188 GB here) and cannot see the 40 GB
# cgroup, so it grows past the cap before collecting. Run the real thing with
# `julia --heap-size-hint=<N>G`, sized to the cgroup MINUS whatever else shares
# it, not to `free -g`.
# -----------------------------------------------------------------------------
