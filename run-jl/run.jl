#!/usr/bin/env julia
# =============================================================================
# run.jl — the JULIA binding drives the single clean `isrm.esm` end to end
# through the PUBLIC EarthSciAST surface. NOTHING MODEL-SHAPED LIVES HERE: this
# file names no pollutant, no column, no grid extent and no record count.
#
#   * `prepare(doc; providers, pushdown_rewrite=true)` — the automatic
#     projection-pushdown rewrite runs inside the engine; the SR provider gates
#     are derived from the rewrite's own record
#     (`metadata.x_esd.pushdown.gated_select`), so this file hand-authors NO
#     gate dict and implements NO provider protocol;
#   * EVERY provider comes FROM THE DOCUMENT (`providers_from_document`:
#     format = `metadata.esio_format`, URL = `source.url_template`) — the SR
#     slabs, the grid, the population, AND the EGU FF10 table, whose ingest the
#     loader now declares in full (esm-spec §8.9): `reader_options` (the zip
#     member glob + header row), `codes` (POLID text -> the pathway enum, an
#     unrecognised code dropping the record), `record_filter` (no coordinate /
#     no annual total is not a record) and `extent` (the surviving count binds
#     N_REC). The src-cell rectangles are the `select` range `W[0:N_SRC]` on
#     their own loader variables;
#   * every reported number is the binding's evaluation of the document's
#     observed graph (`observed_field`) — NO hand-written STEP-0 math here.
#
#   FULL run  (default)          → report sum(deathsK/L) against the tutorial
#   REDUCED   (ISRM_FIRSTN=n)    → first n emission records, totals reported
#
# Emits the cross-language contract record (contract/results_schema.json) with
# model="isrm.esm", mode="runtime_observed_graph", binding="julia".
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc                     # activates EarthSciIO's zarr codec extension
import GeometryOps, GeoInterface # STRtree broad-phase fast path for the join gate
import JSON
const EA = EarthSciAST

include(joinpath(@__DIR__, "paths.jl"))
# contract/results.jl must be included at TOP LEVEL (world age).
include(joinpath(@__DIR__, "..", "contract", "results.jl"))

const T0 = time()
# The InMAP source-receptor tutorial's published national totals
# (https://inmap.run/blog/2019/04/20/sr/), which account for plume rise.
# A REFERENCE POINT, not a target: the document declines BOTH of InMAP's
# plume-rise defects — the high-plume source-index defect (a plume above model
# layer 7 keeps an index built in the coarse 9324-cell grid, then read against
# the 52411-cell ground grid, misplacing 654 of 43650 records) and the inverted
# layerFracs interpolation (6.25% of emitted mass on the wrong side of a
# split). So a run lands ABOVE these by about 1.35%, deliberately. What is
# checked is CORRECTED_K / CORRECTED_L below.
const ORACLE_K = 6928.959583
const ORACLE_L = 15623.924632
# What this document computes at full scale with CORRECT physics — measured
# 2026-08-20 against the repaired store. Unlike ORACLE_K/L this is not an
# external oracle but a regression lock on the document's own output; its
# authority comes from the NumPy oracle agreeing on the weights and from the
# InMAP-faithful configuration having matched the live service to 8.9e-9 before
# the physics was corrected.
const CORRECTED_K = 7022.724781368745
const CORRECTED_L = 15835.993595627131
# Cross-binding spread on this document is ~4e-18 relative, so this is loose by
# many orders of magnitude and catches a real change rather than float noise.
const CORRECTED_REL = 1e-9
# Peak resident set. `/proc/self/statm` field 2 is the CURRENT resident page
# count and is Linux-only; `Sys.maxrss()` is the high-water mark and is portable,
# so it is the fallback wherever /proc does not exist (macOS).
peak_rss_bytes() = isfile("/proc/self/statm") ?
    parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 : Int(Sys.maxrss())

"""A metaparameter's declared default, read from the document (so no grid
extent is written down here)."""
metaparam(doc, name) = Int(get(get(get(doc, "metaparameters", Dict()), name, Dict()),
                               "default", 0))

"""The data sources that DISCOVER their own extent (`extent.metaparameter`) —
the record-bearing tables of the document, whatever they happen to be called. The
two knobs below are scale/locality concerns of a RUN, not of the model, and both
are expressed in the document's own vocabulary."""
record_loaders(doc) = String[String(name) for (name, ld) in get(doc, "data_sources", Dict())
                             if ld isa AbstractDict &&
                                get(ld, "extent", Dict()) isa AbstractDict &&
                                haskey(get(ld, "extent", Dict()), "metaparameter")]

"""Every model PARAMETER fed by data source `src` — `(model, parameter)` pairs.

From esm 1.0.0 a source declares no variables of its own: the binding lives on
the consuming parameter as `update: {kind: "data", source: …, from: …}`
(esm-spec §6.3). So "which columns does this source deliver" is answered by
walking the models, not the source."""
function source_parameters(doc, src)
    out = Tuple{String,String}[]
    for (mname, m) in get(doc, "models", Dict())
        m isa AbstractDict || continue
        for (vname, v) in get(m, "variables", Dict())
            v isa AbstractDict || continue
            up = get(v, "update", nothing); up isa AbstractDict || continue
            get(up, "kind", "") == "data" && get(up, "source", "") == src &&
                push!(out, (String(mname), String(vname)))
        end
    end
    return out
end

"""The record count the sources DISCOVERED — the delivered length of any one
parameter an extent-declaring source feeds. They are aligned by construction
(that is exactly what `record_filter` guarantees), so the first one answers for
all, and this file still names no column of any particular model."""
function discovered_records(doc, insp)
    for name in record_loaders(doc), (mname, vname) in source_parameters(doc, name)
        # 1.0.0 keys the delivered column by the CONSUMING parameter, so try the
        # flattened name and the bare one rather than the old "source.var".
        for key in ("$mname.$vname", vname)
            a = get(insp.const_arrays, key, nothing)
            a === nothing || return length(a)
        end
    end
    error("no record-discovering source delivered an array to size N_REC from")
end

"""The sources that declare a `gated_select`, and the arrays each gates — the
document's own statement of how many model arrays the pushdown rewrite must end
up gating. Derived rather than written down, so splitting or merging a gated
source keeps the check honest."""
function gated_loader_arrays(doc)
    out = Dict{String,Vector{String}}()
    for (name, ld) in get(doc, "data_sources", Dict())
        ld isa AbstractDict || continue
        md = get(ld, "metadata", Dict()); md isa AbstractDict || continue
        xe = get(md, "x_esd", Dict()); xe isa AbstractDict || continue
        gs = get(xe, "gated_select", nothing); gs isa AbstractDict || continue
        out[String(name)] = String[String(a) for a in get(gs, "applies_to", [])]
    end
    return out
end

"""REDUCED runs: truncate every record-discovering source to its first `n`
DELIVERED records with a source-level `select` range (esm-spec §8.9.2). Because
the selection follows the source's own `record_filter`, this picks the same
records the previous runners' post-filter truncation did — and `extent` then
re-discovers the smaller N_REC by itself."""
function truncate_records!(doc, n)
    for name in record_loaders(doc)
        doc["data_sources"][name]["select"] =
            Dict("axes" => [Dict("range" => Dict("start" => 0, "stop" => n))])
    end
    return doc
end

# zarr workaround (unchanged from the validated runners): the SR arrays carry NO
# `.zattrs` object (live GET → 404, which `_fetch_bytes_optional` does not
# swallow), so seed an empty `{}` as a cache HIT before any SR fetch.
function seed_empty_zattrs(cache_root, base, arrays)
    cache = EarthSciIO.Cache(; root=cache_root)
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

# =============================================================================
function main()
    firstn = haskey(ENV, "ISRM_FIRSTN") ? parse(Int, ENV["ISRM_FIRSTN"]) : nothing
    reduced = firstn !== nothing
    println(reduced ? "REDUCED run — first $firstn emission records" :
            "FULL run — whole domain (target deathsK≈$(round(ORACLE_K, digits=2)), deathsL≈$(round(ORACLE_L, digits=2)))")
    println("model:   $MODEL")
    println("scratch: $SCRATCH")

    doc_raw = JSON.parsefile(MODEL)
    cache_root = joinpath(SCRATCH, "run-jl-esio-cache")
    reduced && truncate_records!(doc_raw, firstn)

    # ---- providers FROM THE DOCUMENT — ALL of them --------------------------
    # Including the FF10 table: the loader declares its own reader options, code
    # map, record filter and extent, so there is nothing left here to read, map,
    # filter or count.
    # Only the GATED loaders' arrays: those are the SR slabs, and they are the
    # ones the store has no `.zattrs` for. The 1-D grid arrays are read whole
    # through the ordinary path and must keep whatever attrs the store has.
    gated = gated_loader_arrays(doc_raw)
    for (lname, arrs) in gated
        seed_empty_zattrs(joinpath(cache_root, lname),
                          doc_raw["data_sources"][lname]["source"]["url_template"], arrs)
    end
    println("building providers from the document ...")
    # A local copy of a record loader's source is a LOCALITY choice of this run
    # (gaftp.epa.gov is slow and flaky), so it is a url_override rather than an
    # edit to the document.
    url_overrides = Dict{String,String}()
    if isfile(EGU_ZIP)
        for name in record_loaders(doc_raw)
            url_overrides[name] = "file://" * EGU_ZIP
        end
        println("  record source mirrored from $EGU_ZIP")
    end
    t_providers = @elapsed providers =
        EA.providers_from_document(doc_raw; cache_root=cache_root,
                                   url_overrides=url_overrides)
    println("  providers: ", sort(collect(keys(providers)))); flush(stdout)

    # ---- the gate must cover EVERY declared SR array -----------------------
    # A malformed E_* or conc_* body does not fail: the pathway simply drops out
    # of the rewrite's `applies_to` list, the rest of the rewrite reports
    # success, and the un-gated array is then fetched WHOLE — 330 GB, which
    # surfaces hours later as a memory failure rather than an error. The
    # document says how many arrays it declared for gating; anything less is
    # that silent drop.
    #
    # This runs BEFORE `prepare` rather than after, because in this binding the
    # gated fetch happens INSIDE `prepare` and `PreparedModel.run_doc` is the
    # FLATTENED document, which no longer carries `metadata.x_esd`. So the
    # record is read from the rewrite itself — the same call `prepare` makes
    # internally, on the same document, and idempotent — which also puts the
    # stop before the fetch instead of after it.
    expect_gated = sum(length(v) for v in values(gated); init=0)
    rewritten = EA.desugar_pushdown(EA.serialize_esm_file(
        EA.load(doc_raw; base_path=ISRM_DIR)))
    applies = get(get(get(get(get(rewritten, "metadata", Dict()), "x_esd", Dict()),
                          "pushdown", Dict()), "gated_select", Dict()), "applies_to", [])
    gated_now = String[String(a) for a in applies]
    println("gated arrays: $(length(gated_now)) of $expect_gated declared")
    length(gated_now) == expect_gated || error(
        "the pushdown rewrite gated $(length(gated_now)) arrays but the document " *
        "declares $expect_gated ($gated_now) — a pathway dropped out of the gate " *
        "silently, and its SR array would be fetched UNGATED. Check that each " *
        "conc_*_L* body is a plain two-factor SR*E product and that the " *
        "containment ifelse is the FIRST ifelse in every E_* body.")

    # ---- PREPARE (extent → rewrite → coords → VI → gated fetch → graph) ------
    println("prepare(pushdown_rewrite=true) — N_REC discovered by the loader ...")
    flush(stdout)
    insp = EA.BuildInspection()
    t_prep = @elapsed prep = EA.prepare(doc_raw; providers=providers, base_path=ISRM_DIR,
                                        inspect=insp, pushdown_rewrite=true)
    N_REC = discovered_records(doc_raw, insp)
    println("PREPARE done in $(round(t_prep, digits=1)) s  (peak RSS so far: ",
            round(peak_rss_bytes() / 2^30, digits=2), " GiB)"); flush(stdout)

    # ---- the engine-derived support set (for the contract record) -----------
    mf_keys = [k for k in keys(insp.const_arrays) if startswith(String(k), "pd_member_factor__")]
    isempty(mf_keys) && error("no pd_member_factor__* const array — did the rewrite fire?")
    members = sort!(Int.(insp.const_arrays[first(mf_keys)]))
    n_ppl = length(members)
    n_src = metaparam(doc_raw, "N_SRC")
    n_rcv = metaparam(doc_raw, "N_RCV")
    println("engine-derived support set: |members| = $n_ppl of $n_src source cells")
    !reduced && n_ppl != 1520 &&
        println("  WARNING: expected 1520 emission-bearing cells at full scale")

    # ---- evaluate the observed graph ----------------------------------------
    function rt(v)
        print("  evaluating observed $v ... "); flush(stdout)
        local fld
        t = @elapsed fld = EA.observed_field(prep, insp, v)
        println(round(t, digits=1), " s"); flush(stdout)
        return fld
    end
    t_eval = @elapsed begin
        dK = rt("deathsK"); dL = rt("deathsL"); tp = rt("TotalPM25")
        # per-pathway intermediates through the SAME runtime path, so a
        # disagreement localizes to one pathway instead of only the totals.
        # The emissions side grew a LAYER dimension when the document started
        # stating plume rise: a record is charged to the SR layer its plume
        # reaches, so a pathway's emissions are three arrays. The record's
        # `emis_sum` stays the pathway TOTAL — plume rise moves mass between
        # layers, never into or out of a pathway — so it remains comparable to
        # the ground-level-only baselines. The concentration side did not grow:
        # `conc_<p>` is the document's own sum over the three contractions.
        PW_OBS = ["SOA"         => (["E_VOC_L0",  "E_VOC_L1",  "E_VOC_L2"],  "conc_SOA"),
                  "pNO3"        => (["E_NOx_L0",  "E_NOx_L1",  "E_NOx_L2"],  "conc_pNO3"),
                  "pNH4"        => (["E_NH3_L0",  "E_NH3_L1",  "E_NH3_L2"],  "conc_pNH4"),
                  "pSO4"        => (["E_SOx_L0",  "E_SOx_L1",  "E_SOx_L2"],  "conc_pSO4"),
                  "PrimaryPM25" => (["E_PM25_L0", "E_PM25_L1", "E_PM25_L2"], "conc_PrimaryPM25")]
        pathways = Dict{String,Any}()
        emis_by_layer = Dict{String,Vector{Float64}}()
        for (arr, (evars, cvar)) in PW_OBS
            Es = [vec(rt(v)) for v in evars]
            Ep = reduce(vcat, Es)
            cp = rt(cvar)
            pathways[arr] = (emis_sum = sum(Ep), conc_sum = sum(cp), conc_max = maximum(cp))
            # How much mass plume rise put in each SR layer — the physics made
            # visible as tons, per pathway.
            emis_by_layer[arr] = [sum(e) for e in Es]
        end
        # The layer assignment itself — now a SPLIT, not a single layer:
        # InMAP's sr.Reader.layerFracs charges a record to two SR layers
        # whenever its model layer falls between two entries of `layers`.
        # `sr_lower` is the lower index (integer, compared exactly) and
        # w_sr0/1/2 are the three shares. These are the document's OWN observeds,
        # read through the same `observed_field` path as everything else — this
        # runner does not know what ASME is, and must not: the point of the
        # contract's `plume` block is that the ENGINE produced the assignment
        # from the spec. contract/plume_oracle.py computes the same quantity
        # independently, from the meteorology arrays and without the SR matrix,
        # and compare_results.py checks the two against each other.
        sr_lower = rt("sr_lower")
        stack_layer = rt("stack_layer")
        weights = Dict(w => vec(rt(w)) for w in ("w_sr0", "w_sr1", "w_sr2"))
    end
    println("EVAL done in $(round(t_eval, digits=1)) s")

    plume = plume_block(sr_lower = sr_lower, stack_layer = stack_layer,
                        weights = weights, emis_by_sr_layer = emis_by_layer)

    sK = sum(dK); sL = sum(dL)
    println("\n", "="^70)
    println("  sum(deathsK) = ", sK)
    println("  sum(deathsL) = ", sL)
    println("  Σ TotalPM25  = ", sum(tp))
    if !reduced
        println("  tutorial deathsK=$ORACLE_K  deviation ",
                round(100 * (sK - ORACLE_K) / ORACLE_K, digits=6),
                "%  (reference, not a target)")
        println("  tutorial deathsL=$ORACLE_L deviation ",
                round(100 * (sL - ORACLE_L) / ORACLE_L, digits=6),
                "%  (reference, not a target)")
        rK = (sK - CORRECTED_K) / CORRECTED_K
        rL = (sL - CORRECTED_L) / CORRECTED_L
        if abs(rK) > CORRECTED_REL || abs(rL) > CORRECTED_REL
            println("  WARNING: $sK / $sL differs from the measured ",
                    "corrected-physics totals $CORRECTED_K / $CORRECTED_L by more ",
                    "than $CORRECTED_REL relative. That is a REGRESSION, not a ",
                    "tolerance: the two are the same document on the same store.")
        else
            println("  matches the measured corrected-physics totals to ",
                    abs(rK), " / ", abs(rL))
        end
    end
    println("  lower-SR-layer histogram (records per layer 0/1/2) = ",
            plume["sr_lower"]["histogram"])
    println("  sr_lower sha256 = ", plume["sr_lower"]["sha256"])
    println("  Σ w_sr0/w_sr1/w_sr2 = ",
            join((plume["weights"][w]["sum"] for w in ("w_sr0", "w_sr1", "w_sr2")), " / "),
            "   max|Σw - 1| = ", plume["weights"]["max_sum_error"])
    println("    (check it against `python3 contract/plume_oracle.py",
            reduced ? " --firstn $firstn`" : "`", " — no SR matrix needed)")
    println("="^70)

    # ---- contract record ----------------------------------------------------
    out = joinpath(@__DIR__, reduced ? "results_reduced.json" : "results.json")
    write_results(out;
        binding_version = "julia $(VERSION) / EarthSciAST $(pkgversion(EarthSciAST))",
        model = MODEL, mode = "runtime_observed_graph",
        n_src = n_src, n_rcv = n_rcv, n_rec = N_REC,
        ppl = members,
        pathways = pathways,
        total_pm25 = tp, deathsK = dK, deathsL = dL,
        plume = plume,
        timing = Dict("wall_seconds" => time() - T0,
                      "providers_seconds" => t_providers,
                      "build_seconds" => t_prep,
                      "eval_seconds" => t_eval,
                      "peak_rss_bytes" => peak_rss_bytes()))
    return 0
end

exit(main())
