#!/usr/bin/env julia
# =============================================================================
# run.jl — the JULIA binding drives an `isrm_*.esm` end to end through the
# PUBLIC EarthSciAST surface. `ISRM_MODEL` picks which: `isrm_point.esm` (the
# EGU point inventory, with plume rise) or `isrm_polygon.esm` (an area
# inventory allocated by polygon/cell overlap). NOTHING MODEL-SHAPED LIVES
# HERE: this file names no pollutant, no column, no grid extent, no record
# count and no observed — the document's `metadata.x_esd.report` block names
# what a run reports, so the same shim drives either geometry.
#
#   * `esm_problem(doc, tspan; providers, pushdown_rewrite=true)` — the automatic
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
# model="isrm_point.esm", mode="runtime_observed_graph", binding="julia".
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciAST
using EarthSciIO
using Blosc                     # activates EarthSciIO's zarr codec extension
using Shapefile                 # activates EarthSciIO's shapefile reader extension
import GeometryOps, GeoInterface # STRtree broad-phase fast path for the join gate
import JSON
const EA = EarthSciAST

include(joinpath(@__DIR__, "paths.jl"))
# contract/results.jl must be included at TOP LEVEL (world age).
include(joinpath(@__DIR__, "..", "contract", "results.jl"))

const T0 = time()
# The published national totals a document may claim to reproduce. For
# `inmap_sr_tutorial` — https://inmap.run/blog/2019/04/20/sr/, which accounts
# for plume rise — `published` is a REFERENCE POINT, not a target: the document
# declines BOTH of InMAP's plume-rise defects, the high-plume source-index
# defect (a plume above model layer 7 keeps an index built in the coarse
# 9324-cell grid, then read against the 52411-cell ground grid, misplacing 654
# of 43650 records) and the inverted layerFracs interpolation (6.25% of emitted
# mass on the wrong side of a split). So a run lands ABOVE it by about 1.35%,
# deliberately, and what is CHECKED is `corrected`.
# Keyed by the tag a document puts in `metadata.x_esd.report.oracle`, so which
# published result a run is measured against is the DOCUMENT's claim and not
# this file's guess. A document with no `oracle` tag — isrm_polygon.esm, whose
# example emission layer this repository builds and nobody has published a
# total for — is reported, not graded.
#
#   published: the tutorial's own two numbers.
#   corrected: what THIS document computes at full scale with correct physics,
#              measured 2026-08-20 against the repaired store. Not an external
#              oracle but a regression lock on the document's own output; its
#              authority comes from the NumPy oracle agreeing on the weights and
#              from the InMAP-faithful configuration having matched the live
#              service to 8.9e-9 before the physics was corrected.
#   rel:       cross-binding spread on this document is ~4e-18 relative, so this
#              is loose by many orders of magnitude and catches a real change
#              rather than float noise.
const ORACLES = Dict(
    "inmap_sr_tutorial" => (published = (6928.959583, 15623.924632),
                            corrected = (7022.724781368745, 15835.993595627131),
                            rel = 1e-9))
# Peak resident set. `/proc/self/statm` field 2 is the CURRENT resident page
# count and is Linux-only; `Sys.maxrss()` is the high-water mark and is portable,
# so it is the fallback wherever /proc does not exist (macOS).
peak_rss_bytes() = isfile("/proc/self/statm") ?
    parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 : Int(Sys.maxrss())

"""`metadata.x_esd.report` — what this document says a run should report.

The one place a pollutant, a pathway or an observed name enters this file, and
it enters FROM THE DOCUMENT. A document that declares no `plume` gets no plume
block; one that claims no published result gets no oracle check."""
function report_block(doc)
    md = get(doc, "metadata", Dict()); md isa AbstractDict || (md = Dict())
    xe = get(md, "x_esd", Dict()); xe isa AbstractDict || (xe = Dict())
    rep = get(xe, "report", nothing)
    (rep isa AbstractDict && !isempty(get(rep, "pathways", []))) || error(
        "$MODEL: no metadata.x_esd.report block — this runner reads the " *
        "reported pathway and observed names from the document rather than " *
        "carrying a table of its own")
    return rep
end

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
DELIVERED records with a `select` range (esm-spec §8.9.2). Because the selection
follows the source's own `record_filter`, this picks the same records a
post-filter truncation would — and `extent` then re-discovers the smaller N_REC
by itself.

Written on each consuming PARAMETER rather than on the source, because
`select.axes` is one entry per NATIVE array dimension and a source may deliver
arrays of different ranks: the polygon layer's `geometry` is
`[record, vertex, xy]` while its emission column is `[record]`, so no single
source-level list is right for both. The rank comes from the parameter's own
declared `shape`, and the record axis is axis 0 by definition of a record
table."""
function truncate_records!(doc, n)
    for name in record_loaders(doc), (mname, vname) in source_parameters(doc, name)
        v = doc["models"][mname]["variables"][vname]
        rank = length(get(v, "shape", []))
        v["update"]["from"]["select"] =
            Dict("axes" => vcat(Any[Dict("range" => Dict("start" => 0, "stop" => n))],
                                fill("all", max(rank - 1, 0))))
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
    doc_raw = JSON.parsefile(MODEL)
    report = report_block(doc_raw)
    oracle = get(ORACLES, String(get(report, "oracle", "")), nothing)
    println(reduced ? "REDUCED run — first $firstn emission records" :
            "FULL run — whole domain " * (oracle === nothing ?
                "(the document claims no published total)" :
                "(target deathsK≈$(round(oracle.published[1], digits=2)), " *
                "deathsL≈$(round(oracle.published[2], digits=2)))"))
    println("model:   $MODEL")
    println("scratch: $SCRATCH")

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
    # A local copy of a record source is a LOCALITY choice of this run
    # (gaftp.epa.gov is slow and flaky; the example polygon layer is built in
    # this repository and never needs fetching at all), so it is a url_override
    # rather than an edit to the document — and it is matched by the document's
    # OWN url basename, so this file names no source.
    url_overrides = Dict{String,String}()
    for name in record_loaders(doc_raw)
        url = get(get(doc_raw["data_sources"][name], "source", Dict()), "url_template", "")
        local_path = local_mirror(url)
        if !isempty(local_path)
            url_overrides[name] = "file://" * abspath(local_path)
            println("  $name mirrored from $local_path")
        end
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
    # This runs BEFORE construction rather than after, because in this binding
    # the gated fetch happens INSIDE it and `EsmProblem.run_doc` is the
    # FLATTENED document, which no longer carries `metadata.x_esd`. So the
    # record is read from the rewrite itself — the same call construction makes
    # internally, on the same document, and idempotent — which also puts the
    # stop before the fetch instead of after it.
    expect_gated = sum(length(v) for v in values(gated); init=0)
    rewritten = EA.desugar_pushdown(EA.serialize_esm_file(
        EA.load_document(doc_raw; base_path=ISRM_DIR)))
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
    println("esm_problem(pushdown_rewrite=true) — N_REC discovered by the loader ...")
    flush(stdout)
    insp = EA.BuildInspection()
    # EarthSciAST phase 4: `prepare` is replaced by problem CONSTRUCTION, which
    # absorbs the same pipeline. The only new argument is `tspan`, which
    # `prepare` did not take; this driver never integrates, so the interval is
    # nominal. `observed_field` is now TWO arguments — the problem owns its own
    # BuildInspection, so the caller no longer threads the same one through both
    # calls and hopes they match.
    t_prep = @elapsed prep = EA.esm_problem(doc_raw, (0.0, 1.0);
                                            providers=providers, base_path=ISRM_DIR,
                                            inspect=insp, pushdown_rewrite=true)
    N_REC = length(EA.observed_field(prep, String(report["record_field"])))
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
    !reduced && oracle !== nothing && n_ppl != 1520 &&
        println("  WARNING: expected 1520 emission-bearing cells at full scale")

    # ---- evaluate the observed graph ----------------------------------------
    function rt(v)
        print("  evaluating observed $v ... "); flush(stdout)
        local fld
        t = @elapsed fld = EA.observed_field(prep, v)
        println(round(t, digits=1), " s"); flush(stdout)
        return fld
    end
    sr_lower = nothing; stack_layer = nothing; weights = Dict{String,Vector{Float64}}()
    t_eval = @elapsed begin
        dK = rt(report["deaths"]["krewski"]); dL = rt(report["deaths"]["lepeule"])
        tp = rt(report["total_pm25"])
        # per-pathway intermediates through the SAME runtime path, so a
        # disagreement localizes to one pathway instead of only the totals.
        # `emissions` is a LIST: isrm_point.esm splits a pathway across the
        # three SR emission layers a record's plume falls between,
        # isrm_polygon.esm has one because an area source emits at the ground.
        # `emis_sum` is the pathway TOTAL either way — plume rise moves mass
        # between layers, never into or out of a pathway — so it stays
        # comparable across both, and to the ground-level-only baselines.
        pathways = Dict{String,Any}()
        emis_by_layer = Dict{String,Vector{Float64}}()
        for entry in report["pathways"]
            arr = String(entry["sr_array"])
            Es = [vec(rt(String(v))) for v in entry["emissions"]]
            cp = rt(String(entry["concentration"]))
            pathways[arr] = (emis_sum = sum(reduce(vcat, Es)),
                             conc_sum = sum(cp), conc_max = maximum(cp))
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
        # A document with no plume rise to state — an area source has no stack —
        # declares no `plume` key and emits no `plume` block.
        if haskey(report, "plume")
            pl = report["plume"]
            sr_lower = rt(String(pl["sr_lower"]))
            stack_layer = rt(String(pl["stack_layer"]))
            weights = Dict(String(w) => vec(rt(String(w))) for w in pl["weights"])
        end
    end
    println("EVAL done in $(round(t_eval, digits=1)) s")

    plume = haskey(report, "plume") ?
        plume_block(sr_lower = sr_lower, stack_layer = stack_layer,
                    weights = weights, emis_by_sr_layer = emis_by_layer) : nothing

    sK = sum(dK); sL = sum(dL)
    println("\n", "="^70)
    println("  sum(deathsK) = ", sK)
    println("  sum(deathsL) = ", sL)
    println("  Σ TotalPM25  = ", sum(tp))
    println("  Σ emitted    = ",
            join(("$k $(v.emis_sum)" for (k, v) in sort(collect(pathways), by=first)), " / "))
    if !reduced && oracle !== nothing
        oK, oL = oracle.published
        cK, cL = oracle.corrected
        println("  tutorial deathsK=$oK  deviation ",
                round(100 * (sK - oK) / oK, digits=6), "%  (reference, not a target)")
        println("  tutorial deathsL=$oL deviation ",
                round(100 * (sL - oL) / oL, digits=6), "%  (reference, not a target)")
        rK = (sK - cK) / cK
        rL = (sL - cL) / cL
        if abs(rK) > oracle.rel || abs(rL) > oracle.rel
            println("  WARNING: $sK / $sL differs from the measured ",
                    "corrected-physics totals $cK / $cL by more than $(oracle.rel) ",
                    "relative. That is a REGRESSION, not a tolerance: the two are ",
                    "the same document on the same store.")
        else
            println("  matches the measured corrected-physics totals to ",
                    abs(rK), " / ", abs(rL))
        end
    end
    if plume !== nothing
        println("  lower-SR-layer histogram (records per layer 0/1/2) = ",
                plume["sr_lower"]["histogram"])
        println("  sr_lower sha256 = ", plume["sr_lower"]["sha256"])
        println("  Σ w_sr0/w_sr1/w_sr2 = ",
                join((plume["weights"][String(w)]["sum"] for w in report["plume"]["weights"]), " / "),
                "   max|Σw - 1| = ", plume["weights"]["max_sum_error"])
        println("    (check it against `python3 contract/plume_oracle.py",
                reduced ? " --firstn $firstn`" : "`", " — no SR matrix needed)")
    end
    println("="^70)

    # ---- contract record ----------------------------------------------------
    # Named after the MODEL, because two documents share this shim and must not
    # share one record file: isrm_point.esm and isrm_polygon.esm answer
    # different questions over the same grid.
    stem = splitext(basename(MODEL))[1]
    out = joinpath(@__DIR__, "results_$(stem)$(reduced ? "_reduced" : "").json")
    write_results(out;
        binding_version = "julia $(VERSION) / EarthSciAST $(EarthSciAST.LIBRARY_VERSION)",
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
