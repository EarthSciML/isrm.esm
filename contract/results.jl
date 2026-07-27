# =============================================================================
# results.jl — emit a contract/results_schema.json record from a Julia runner.
#
#   include(joinpath(@__DIR__, "..", "contract", "results.jl"))
#   write_results("run-model-jl/results.json";
#                 binding_version = "...", model = "isrm.esm",
#                 mode = "runtime_observed_graph",
#                 n_src = 52411, n_rcv = 52411, n_rec = 43650,
#                 ppl = ppl_ids,                       # 1-based, any order
#                 pathways = Dict("SOA" => (emis_sum=..., conc_sum=..., conc_max=...), ...),
#                 total_pm25 = TotalPM25, deathsK = dK, deathsL = dL,
#                 timing = Dict("wall_seconds" => t))
#
# The hashing and sampling rules here MUST match contract/compare_results.py.
# =============================================================================
import JSON
import SHA

const SAMPLE_N = 25

"""Fixed 1-based sample indices — pure integer arithmetic, mirroring
`sample_indices()` in compare_results.py."""
function sample_indices(n_rcv::Integer)
    d = SAMPLE_N - 1
    return [1 + div(k * (n_rcv - 1) + div(d, 2), d) for k in 0:d]
end

"""sha256 over sorted 1-based ids as ASCII decimals joined by "," (no spaces)."""
function ppl_sha256(ids)
    s = join(string.(sort(collect(Int, ids))), ",")
    return bytes2hex(SHA.sha256(s))
end

"""sha256 over a float field as little-endian IEEE-754 float64 bytes."""
function field_sha256(v::AbstractVector{<:Real})
    buf = IOBuffer()
    for x in v
        write(buf, htol(Float64(x)))
    end
    return bytes2hex(SHA.sha256(take!(buf)))
end

"""Summarize one length-n_rcv field into the schema's FieldSummary shape."""
function field_summary(v::AbstractVector{<:Real})
    vv = Float64.(v)
    idx = sample_indices(length(vv))
    return Dict{String,Any}(
        "sum"    => sum(vv),
        "min"    => minimum(vv),
        "max"    => maximum(vv),
        "sample" => [vv[i] for i in idx],
        "sha256" => field_sha256(vv),
    )
end

function write_results(path::AbstractString;
                       binding_version::AbstractString = string(VERSION),
                       model::AbstractString,
                       mode::AbstractString,
                       n_src::Integer, n_rcv::Integer, n_rec::Integer,
                       ppl,
                       pathways::AbstractDict,
                       total_pm25::AbstractVector,
                       deathsK::AbstractVector,
                       deathsL::AbstractVector,
                       include_ppl_ids::Bool = true,
                       timing::Union{Nothing,AbstractDict} = nothing)
    mode in ("runtime_observed_graph", "oracle_step0") ||
        error("mode must be \"runtime_observed_graph\" or \"oracle_step0\", got $mode")
    ids = sort(collect(Int, ppl))

    pw = Dict{String,Any}()
    for (k, v) in pathways
        entry = Dict{String,Any}("emis_sum" => Float64(v.emis_sum),
                                 "conc_sum" => Float64(v.conc_sum))
        hasproperty(v, :conc_max) && (entry["conc_max"] = Float64(v.conc_max))
        pw[String(k)] = entry
    end

    rec = Dict{String,Any}(
        "binding"         => "julia",
        "binding_version" => binding_version,
        "model"           => basename(model),
        "mode"            => mode,
        "grid"    => Dict("n_src" => Int(n_src), "n_rcv" => Int(n_rcv), "n_rec" => Int(n_rec)),
        "ppl"     => merge(Dict{String,Any}("count"  => length(ids),
                                            "sha256" => ppl_sha256(ids)),
                           include_ppl_ids ? Dict{String,Any}("ids" => ids) : Dict{String,Any}()),
        "pathways"   => pw,
        "total_pm25" => field_summary(total_pm25),
        "deaths"     => Dict("krewski" => field_summary(deathsK),
                             "lepeule" => field_summary(deathsL)),
    )
    timing === nothing || (rec["timing"] = Dict{String,Any}(String(k) => v for (k, v) in timing))

    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.print(io, rec, 2)
    end
    println("wrote $path  (mode=$mode, |ppl|=$(length(ids)), ",
            "sum deathsK=", rec["deaths"]["krewski"]["sum"], ")")
    return rec
end
