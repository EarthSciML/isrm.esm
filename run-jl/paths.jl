# =============================================================================
# paths.jl — machine-independent path resolution for the Phase-1 Julia shim.
#
# Everything derives from this file's own location (hence the isrm.esm checkout)
# plus optional environment overrides, so the shim is portable across machines.
#
#   ISRM_MODEL    the .esm to drive          (default: <repo>/isrm.esm)
#   ISRM_SCRATCH  bulk scratch root          (default: first writable of
#                 /scratch.local/$USER, $SCRATCH, tempdir())
#   EA_PATH       EarthSciAST.jl checkout    (default: <sibling>/EarthSciAST/pkg/EarthSciAST.jl)
#   IO_PATH       EarthSciIO julia checkout  (default: <sibling>/EarthSciIO/julia)
#   EGU_ZIP       FF10 point-source zip      (default: <repo>/run-model-jl/data/2016fd_inputs_point.zip)
#
# SR chunk blobs land under ISRM_SCRATCH, so it wants GBs of fast local disk —
# NOT a network filesystem and NEVER a tmpfs /tmp (it eats the memory cgroup).
# =============================================================================
if !@isdefined(ISRM_PATHS_LOADED)

const ISRM_PATHS_LOADED = true

const RUNJL_DIR  = @__DIR__                      # run-jl/
const REPO       = dirname(RUNJL_DIR)            # the isrm.esm checkout
const CODE_ROOT  = dirname(REPO)                 # where sibling checkouts live

const MODEL    = get(ENV, "ISRM_MODEL", joinpath(REPO, "isrm.esm"))
const ISRM_DIR = dirname(MODEL)

const EA_PATH = get(ENV, "EA_PATH", joinpath(CODE_ROOT, "EarthSciAST", "pkg", "EarthSciAST.jl"))
const IO_PATH = get(ENV, "IO_PATH", joinpath(CODE_ROOT, "EarthSciIO", "julia"))

const EGU_ZIP = get(ENV, "EGU_ZIP",
                    joinpath(REPO, "run-model-jl", "data", "2016fd_inputs_point.zip"))

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
