# =============================================================================
# paths.jl — machine-independent path resolution for the Phase-1 Julia shim.
#
# Everything derives from this file's own location (hence the isrm.esm checkout)
# plus optional environment overrides, so the shim is portable across machines.
#
#   ISRM_MODEL    the .esm to drive          (default: <repo>/isrm_point.esm;
#                 <repo>/isrm_polygon.esm is the area-source sibling. Both
#                 import the shared templates from isrm_base.esm)
#   ISRM_SCRATCH  bulk scratch root          (default: first writable of
#                 /scratch.local/$USER, $SCRATCH, tempdir())
#   EA_PATH       EarthSciAST.jl checkout    (default: <sibling>/EarthSciAST/pkg/EarthSciAST.jl)
#   IO_PATH       EarthSciIO julia checkout  (default: <sibling>/EarthSciIO/julia)
#   DATA_DIR      local mirrors of fetchable sources (default: <repo>/data).
#                 A data source whose url_template's basename names a file in
#                 here is read from disk instead of the network — the FF10 point
#                 inventory, the example polygon layer. Absent file = fetch the
#                 document's url_template through the EarthSciIO cache.
#   EGU_ZIP       DEPRECATED single-file alias, kept so an existing EGU_ZIP=…
#                 invocation still works; folded in as one more mirror, matched
#                 by its own basename.
#
# SR chunk blobs land under ISRM_SCRATCH, so it wants GBs of fast local disk —
# NOT a network filesystem and NEVER a tmpfs /tmp (it eats the memory cgroup).
# =============================================================================
if !@isdefined(ISRM_PATHS_LOADED)

const ISRM_PATHS_LOADED = true

const RUNJL_DIR  = @__DIR__                      # run-jl/
const REPO       = dirname(RUNJL_DIR)            # the isrm.esm checkout
const CODE_ROOT  = dirname(REPO)                 # where sibling checkouts live

const MODEL    = get(ENV, "ISRM_MODEL", joinpath(REPO, "isrm_point.esm"))
const ISRM_DIR = dirname(MODEL)

const EA_PATH = get(ENV, "EA_PATH", joinpath(CODE_ROOT, "EarthSciAST", "pkg", "EarthSciAST.jl"))
const IO_PATH = get(ENV, "IO_PATH", joinpath(CODE_ROOT, "EarthSciIO", "julia"))

const DATA_DIR = get(ENV, "DATA_DIR", joinpath(REPO, "data"))
const EGU_ZIP  = get(ENV, "EGU_ZIP", "")   # deprecated; see local_mirror

"""
    local_mirror(url) -> String

A local copy of `url`, or `""`. Mirroring is a LOCALITY choice of a run, never a
property of the model, so it is resolved by matching the document's own
`url_template` basename against `DATA_DIR` rather than by naming any particular
source here. That is what lets one shim mirror the FF10 point inventory for
isrm_point.esm and the example polygon layer for isrm_polygon.esm without
knowing either name.
"""
function local_mirror(url::AbstractString)
    base = basename(rstrip(String(url), '/'))
    isempty(base) && return ""
    (!isempty(EGU_ZIP) && basename(EGU_ZIP) == base && isfile(EGU_ZIP)) && return EGU_ZIP
    cand = joinpath(DATA_DIR, base)
    return isfile(cand) ? cand : ""
end

"""First writable scratch root among the candidates, as `<root>/isrm-esm`."""
function _resolve_scratch()
    env = get(ENV, "ISRM_SCRATCH", "")
    isempty(env) || (mkpath(env); return env)
    user = get(ENV, "USER", get(ENV, "LOGNAME", "user"))
    for base in ("/scratch.local/$user", get(ENV, "SCRATCH", ""), tempdir())
        isempty(base) && continue
        candidate = joinpath(base, "isrm-esm")
        try
            mkpath(candidate)
            return candidate
        catch
            continue
        end
    end
    error("no writable scratch root found; set ISRM_SCRATCH")
end

const SCRATCH = _resolve_scratch()

const ZARR_URL = get(ENV, "ISRM_ZARR_URL", "s3://inmap-model/isrm_v1.2.1.zarr/")

isfile(MODEL) || @warn "model file not found (set ISRM_MODEL)" MODEL

end # ISRM_PATHS_LOADED
