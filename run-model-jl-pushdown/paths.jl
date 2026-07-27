# =============================================================================
# paths.jl — machine-independent path resolution for the pushdown runners.
#
# Everything derives from this file's own location (hence the isrm.esm checkout)
# plus optional environment overrides, so the runners are portable across
# machines.  Include this (or l3_common.jl, which includes it) before using any
# of the constants below.
#
#   ISRM_MODEL     the .esm to drive          (default: <repo>/isrm_pushdown.esm)
#   ISRM_RUNMODEL  run-model.jl's directory   (default: <repo>/run-model-jl)
#                  — supplies checkpoints/stage1.jls and cache_meta/
#   ISRM_SCRATCH   bulk scratch root          (default: first writable of
#                  /scratch.local/$USER, $SCRATCH, tempdir())
#   EA_PATH        EarthSciAST.jl checkout    (default: <sibling>/EarthSciAST/pkg/EarthSciAST.jl)
#   IO_PATH        EarthSciIO julia checkout  (default: <sibling>/EarthSciIO/julia)
#
# SR chunk blobs land under ISRM_SCRATCH and are evicted between pathways, so it
# wants tens of GB of fast local disk — NOT a network filesystem.
# =============================================================================
if !@isdefined(ISRM_PATHS_LOADED)

const ISRM_PATHS_LOADED = true

const PUSHDOWN_DIR = @__DIR__                     # run-model-jl-pushdown/
const REPO         = dirname(PUSHDOWN_DIR)        # the isrm.esm checkout
const CODE_ROOT    = dirname(REPO)                # where sibling checkouts live

const MODEL    = get(ENV, "ISRM_MODEL", joinpath(REPO, "isrm_pushdown.esm"))
const ISRM_DIR = dirname(MODEL)
const RUNMODEL = get(ENV, "ISRM_RUNMODEL", joinpath(REPO, "run-model-jl"))

const EA_PATH = get(ENV, "EA_PATH", joinpath(CODE_ROOT, "EarthSciAST", "pkg", "EarthSciAST.jl"))
const IO_PATH = get(ENV, "IO_PATH", joinpath(CODE_ROOT, "EarthSciIO", "julia"))

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
const N_SRC    = 52411

isfile(MODEL) || @warn "model file not found (set ISRM_MODEL)" MODEL

end # ISRM_PATHS_LOADED
