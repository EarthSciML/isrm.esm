#!/usr/bin/env julia
# =============================================================================
# L2 — LIVE-S3 selective fetch through the pushdown path.
#
# Constructs an EarthSciIO zarr const_provider for the REAL ISRM SR
# (s3://inmap-model/isrm_v1.2.1.zarr), one SR variable (SOA), and fetches ~8
# arbitrary source-cell rows two ways, asserting byte/tol equality and that the
# pushdown fetched ONLY the intersecting chunks:
#
#   (i)  PUSHDOWN batched multi-index select
#        materialize(prov; select=[{indices:[0]}, {indices: ids-1}, "all"])
#   (ii) INDEPENDENT reference:
#        (a) each id fetched on its own (single-index select) → stack & compare
#        (b) a whole source-CHUNK fetched via a :slice selector (a different
#            reader code path) → index rows in it & compare (proves the VALUES,
#            not just self-consistency between two index-selects)
#
# Chunk accounting is measured on a DEDICATED fresh cache that saw ONLY the
# batched fetch: exactly the distinct intersecting chunk objects are present and
# no non-intersecting chunk is. Also measures on-disk (compressed) chunk size to
# project the L3 disk budget.
# =============================================================================
import Pkg; Pkg.activate(@__DIR__; io=devnull)
using EarthSciIO
using Blosc   # activates the EarthSciIOBloscExt decode path for blosc chunks
import JSON

const ZARR_URL = "s3://inmap-model/isrm_v1.2.1.zarr/"
const ARR      = "SOA"
const HERE     = @__DIR__
const SCRATCH  = "/private/tmp/claude-501/-Users-ctessum-code-earthsciml-isrm-esm/dc49fdd1-a55e-425c-a6d0-6e168223d411/scratchpad"

# 8 arbitrary 1-based source ids (2 share a chunk → 7 distinct chunk-rows)
const IDS = [1, 100, 5000, 20000, 30000, 40000, 50000, 52411]

# ---- .zattrs 404 workaround: the SR arrays carry no .zattrs object (live GET
# 404s, which the reader raises as a plain error, not a swallowed CacheMiss);
# seed an empty {} as a cache HIT so the optional fetch never hits the net. ----
function seed_empty_zattrs(cache, base, arrays)
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

# enumerate every blob file under <root>/v1/blobs
function all_blob_files(root)
    bd = joinpath(root, "v1", "blobs")
    isdir(bd) || return String[]
    out = String[]
    for (dir, _, files) in walkdir(bd), f in files
        push!(out, joinpath(dir, f))
    end
    return out
end

blob_present(cache, url) =
    EarthSciIO.get_blob(cache.store, EarthSciIO.cache_key(url)) !== nothing
blob_path(cache, url) =
    EarthSciIO.get_blob(cache.store, EarthSciIO.cache_key(url))

sel(layer0, srcs0, rcv) = Dict("axes" => Any[
    Dict("indices" => layer0), Dict("indices" => srcs0), rcv])

println("="^74)
println("L2 — live-S3 selective fetch  (URL=$ZARR_URL  array=$ARR)")
println("  1-based source ids: ", IDS)
println("="^74)

base = rstrip(ZARR_URL, '/')

# ---------------------------------------------------------------------------
# STEP 0 — pin the .zarray facts (fetches ONLY metadata, never a chunk)
# ---------------------------------------------------------------------------
meta_cache = EarthSciIO.Cache(; root=joinpath(SCRATCH, "l2_cache_meta"))
seed_empty_zattrs(meta_cache, ZARR_URL, [ARR])
provM = EarthSciIO.const_provider(meta_cache, ZARR_URL; format="zarr", variables=[ARR])
shp   = EarthSciIO.array_shape(provM, ARR)
println("array_shape($ARR) = ", shp, "   (reads only .zarray, no chunk)")
@assert length(shp) == 3 "expected a 3-D SR array"
NSRC = shp[2]; NRCV = shp[3]
@assert all(1 .<= IDS .<= NSRC) "some id out of range 1..$NSRC"
# read the pinned chunking straight from the .zarray
zurl = "$base/$ARR/.zarray"
EarthSciIO.fetch_blob(meta_cache, zurl)   # ensure cached
zarray = JSON.parse(String(read(blob_path(meta_cache, zurl))))
CHUNKS  = Int.(zarray["chunks"])
DIM_SEP = let s = get(zarray, "dimension_separator", nothing); (s === nothing || s == "") ? "." : String(s) end
println("chunks  = ", CHUNKS, "   dim_sep = '", DIM_SEP, "'")
@assert CHUNKS[1] == 1 "layer chunk expected 1"
SRC_CHUNK = CHUNKS[2]

# distinct intersecting source chunk indices (0-based) + their chunk-object keys
src0        = IDS .- 1
src_chunks  = sort(unique(g ÷ SRC_CHUNK for g in src0))
chunkkey(c) = "0" * DIM_SEP * string(c) * DIM_SEP * "0"      # layer0 . srcchunk . rcv0
chunk_urls  = ["$base/$ARR/$(chunkkey(c))" for c in src_chunks]
println("distinct intersecting source-chunks (size $SRC_CHUNK): ", src_chunks,
        "   → ", length(src_chunks), " chunk objects")

# ---------------------------------------------------------------------------
# STEP 1 — (i) PUSHDOWN batched fetch into a DEDICATED fresh cache
# ---------------------------------------------------------------------------
rootB = joinpath(SCRATCH, "l2_cache_batched")
isdir(rootB) && rm(rootB; recursive=true, force=true)
cacheB = EarthSciIO.Cache(; root=rootB)
seed_empty_zattrs(cacheB, ZARR_URL, [ARR])
provB = EarthSciIO.const_provider(cacheB, ZARR_URL; format="zarr", variables=[ARR])
println("\n[i] pushdown batched materialize(select=[{indices:[0]}, {indices: ids-1}, all]) ...")
tB = @elapsed ndB = EarthSciIO.materialize(provB; select=sel([0], src0, "all"))
B = ndB[ARR].data
println("    got array size ", size(B), "   in ", round(tB, digits=1), " s")
@assert size(B) == (1, length(IDS), NRCV) "unexpected pushdown result shape $(size(B))"

# ---- chunk accounting on cacheB (saw ONLY the batched fetch) ----
present_chunks = [u for u in chunk_urls if blob_present(cacheB, u)]
# a non-intersecting chunk (source ids 100..199 0-based = chunk 1) must be ABSENT
noniz = "$base/$ARR/$(chunkkey(1))"
absent_ok = (1 in src_chunks) ? true : !blob_present(cacheB, noniz)
# total chunk blobs = all blobs minus the 2 metadata blobs (.zarray fetched, .zattrs seeded)
nblobs      = length(all_blob_files(rootB))
zarray_url  = "$base/$ARR/.zarray"
zattrs_url  = "$base/$ARR/.zattrs"
n_meta      = count(u -> blob_present(cacheB, u), [zarray_url, zattrs_url])
n_chunk_blobs = nblobs - n_meta
chunk_bytes = sum(filesize(blob_path(cacheB, u)) for u in present_chunks)
println("    chunk objects present : ", length(present_chunks), " / ",
        length(chunk_urls), " expected")
println("    total blobs in store  : ", nblobs, "  (", n_meta, " meta + ",
        n_chunk_blobs, " chunk)")
println("    non-intersecting chunk (0.1.0) absent: ", absent_ok)
println("    fetched compressed bytes: ", chunk_bytes, " (",
        round(chunk_bytes/2^20, digits=1), " MiB)  → ",
        round(chunk_bytes/length(present_chunks)/2^20, digits=2), " MiB/chunk")

# ---------------------------------------------------------------------------
# STEP 2 — (ii) INDEPENDENT reference fetches into a separate cache
# ---------------------------------------------------------------------------
rootR = joinpath(SCRATCH, "l2_cache_ref")
isdir(rootR) && rm(rootR; recursive=true, force=true)
cacheR = EarthSciIO.Cache(; root=rootR)
seed_empty_zattrs(cacheR, ZARR_URL, [ARR])
provR = EarthSciIO.const_provider(cacheR, ZARR_URL; format="zarr", variables=[ARR])

println("\n[ii-a] per-id single-index fetches (independent of the batched gather) ...")
R = Array{Float64}(undef, 1, length(IDS), NRCV)
for (j, id) in enumerate(IDS)
    nd = EarthSciIO.materialize(provR; select=sel([0], [id - 1], "all"))
    r  = nd[ARR].data
    @assert size(r) == (1, 1, NRCV)
    R[1, j, :] = r[1, 1, :]
end
maxdiff_a = maximum(abs.(B .- R))
println("       max|batched − per_id| = ", maxdiff_a)

println("[ii-b] whole-CHUNK :slice fetch (different reader path) as value oracle ...")
# fetch source-chunk 0 (0-based rows 0:99) via a slice selector, index ids 1 & 100
ndS = EarthSciIO.materialize(provR;
        select=Dict("axes"=>Any[Dict("indices"=>[0]),
                                Dict("slice"=>[0, SRC_CHUNK]), "all"]))
S = ndS[ARR].data                       # (1, SRC_CHUNK, NRCV)
@assert size(S) == (1, SRC_CHUNK, NRCV)
# batched positions of ids 1 (0-based 0) and 100 (0-based 99), both in chunk 0
j1   = findfirst(==(1), IDS);   j100 = findfirst(==(100), IDS)
d1   = maximum(abs.(B[1, j1,   :] .- S[1, 1,   :]))   # id 1   → within-chunk row 0
d100 = maximum(abs.(B[1, j100, :] .- S[1, 100, :]))   # id 100 → within-chunk row 99
println("       max|batched(id=1)   − slice_row0 |  = ", d1)
println("       max|batched(id=100) − slice_row99|  = ", d100)

# ---------------------------------------------------------------------------
# VERDICT
# ---------------------------------------------------------------------------
vals_ok   = maxdiff_a == 0.0 && d1 == 0.0 && d100 == 0.0
chunks_ok = length(present_chunks) == length(chunk_urls) &&
            n_chunk_blobs == length(chunk_urls) && absent_ok
pass = vals_ok && chunks_ok && size(B) == (1, length(IDS), NRCV)

println("\n", "="^74)
println("L2 RESULT: ", pass ? "PASS" : "FAIL")
println("  ids (1-based)         : ", IDS)
println("  distinct chunk-rows   : ", length(src_chunks), "  (chunks ", src_chunks, ")")
println("  chunk objects fetched : ", n_chunk_blobs, "  (== ", length(chunk_urls),
        " intersecting; ", nblobs, " total blobs incl. ", n_meta, " meta)")
println("  value agreement       : per-id max|Δ|=", maxdiff_a,
        "  slice-oracle max|Δ|=", max(d1, d100), "  (byte-exact: ", vals_ok, ")")
println("  compressed fetch      : ", round(chunk_bytes/2^20, digits=1),
        " MiB for ", length(present_chunks), " chunks (",
        round(chunk_bytes/length(present_chunks)/2^20, digits=2), " MiB/chunk)")
println("="^74)
pass || error("L2 FAILED")
