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
#                 plume = plume_block(sr_lower = ..., stack_layer = ...,
#                                     weights = Dict("w_sr0" => ..., "w_sr1" => ..., "w_sr2" => ...),
#                                     emis_by_sr_layer = Dict("SOA" => [e0, e1, e2], ...)),
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

"""sha256 over an ORDERED integer sequence as ASCII decimals joined by ",".

The same wire format as `ppl_sha256` — ASCII decimals, "," separator, no spaces
— but order-preserving. `ppl` is a member *set*, so it is sorted before hashing;
a per-record integer field (e.g. the plume-rise SR-layer assignment, one value
per emission record) is a *sequence* whose order is part of the value, so it
must not be sorted. Same convention, two uses; mirrors `int_seq_sha256` in
results.py and contract.rs."""
function int_seq_sha256(values)
    s = join((string(Int(v)) for v in values), ",")
    return bytes2hex(SHA.sha256(s))
end

"""sha256 over sorted 1-based ids as ASCII decimals joined by "," (no spaces)."""
ppl_sha256(ids) = int_seq_sha256(sort(collect(Int, ids)))

"""An integer-valued observed, read back off the graph as Float64, as Ints.

The document's `sr_lower` / `stack_layer` are sums of 0.0/1.0 indicators, so
every value is an exact integer in Float64 and this conversion is lossless. A
value that is NOT integral means the observed is no longer the indicator sum it
is supposed to be — a real disagreement about the physics — so it errors rather
than rounding it away."""
function as_int_seq(values, label::AbstractString)
    out = Vector{Int}(undef, length(values))
    for (i, x) in enumerate(values)
        f = Float64(x)
        isinteger(f) || error("$label[$i] = $f is not integral; an integer-valued " *
                              "observed came back fractional, so the layer assignment " *
                              "is not what the document states")
        out[i] = Int(f)
    end
    return out
end

"""Counts of 0, 1, 2, ... over a non-negative integer sequence.

At least `min_bins` bins so a reduced run, where a layer may simply be empty,
emits the same shape a full run does; more if the data needs them, because a
value the schema does not expect must be VISIBLE rather than dropped off the end
of a fixed-width histogram."""
function histogram(values, min_bins::Integer)
    any(v -> v < 0, values) && error("negative value in a layer assignment")
    n = max(Int(min_bins), isempty(values) ? 0 : maximum(values) + 1)
    bins = zeros(Int, n)
    for v in values
        bins[v + 1] += 1
    end
    return bins
end

"""max over records of |w0 + w1 + w2 - 1|.

InMAP's `layerFracs` conserves mass: whether a record lands wholly in one SR
layer or is split across two, its weights sum to exactly 1. So this is float
noise or a broken document, and nothing in between."""
function weight_sum_error(w0, w1, w2)
    length(w0) == length(w1) == length(w2) ||
        error("the three SR-layer weight fields differ in length")
    isempty(w0) && return 0.0
    return maximum(abs(Float64(a) + Float64(b) + Float64(c) - 1.0)
                   for (a, b, c) in zip(w0, w1, w2))
end

"""The schema's `plume` block, from the document's OWN observeds.

`sr_lower`, `stack_layer` and the three `weights` fields are per-record
observeds read straight off the graph — nothing here recomputes plume rise, and
nothing here knows what ASME is. `weights` maps "w_sr0"/"w_sr1"/"w_sr2" to those
observeds; `emis_by_sr_layer` maps each SR array name to the three
`sum(E_<pathway>_L<layer>)` totals, in layer order.

Two digests are integer sequences in record order, hashed the way `ppl` is, so
they can be compared EXACTLY — against the other bindings and against
`contract/records/plume_oracle.json`, which computes the same assignment from
the meteorology arrays without touching the SR matrix. The weights themselves
are floats (`sr.Reader.layerFracs` interpolates a plume between two SR layers)
and get the FieldSummary treatment every other float here gets."""
function plume_block(; sr_lower, stack_layer, weights::AbstractDict,
                     emis_by_sr_layer::AbstractDict)
    sr = as_int_seq(sr_lower, "sr_lower")
    sl = as_int_seq(stack_layer, "stack_layer")
    w = Dict(k => Float64.(collect(weights[k])) for k in ("w_sr0", "w_sr1", "w_sr2"))
    wblock = Dict{String,Any}(
        "count"         => length(w["w_sr0"]),
        "max_sum_error" => weight_sum_error(w["w_sr0"], w["w_sr1"], w["w_sr2"]))
    for (k, v) in w
        wblock[k] = field_summary(v)
    end
    return Dict{String,Any}(
        "sr_lower" => Dict{String,Any}("count"     => length(sr),
                                       "histogram" => histogram(sr, 3),
                                       "sha256"    => int_seq_sha256(sr)),
        "stack_layer" => Dict{String,Any}("count"     => length(sl),
                                          "histogram" => histogram(sl, 4),
                                          "sha256"    => int_seq_sha256(sl)),
        "weights" => wblock,
        "pathways" => Dict{String,Any}(
            String(k) => Dict{String,Any}("by_sr_layer" => Float64.(collect(v)))
            for (k, v) in emis_by_sr_layer),
    )
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
                       plume::Union{Nothing,AbstractDict} = nothing,
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
    plume === nothing || (rec["plume"] = plume)
    timing === nothing || (rec["timing"] = Dict{String,Any}(String(k) => v for (k, v) in timing))

    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.print(io, rec, 2)
    end
    println("wrote $path  (mode=$mode, |ppl|=$(length(ids)), ",
            "sum deathsK=", rec["deaths"]["krewski"]["sum"], ")")
    return rec
end
