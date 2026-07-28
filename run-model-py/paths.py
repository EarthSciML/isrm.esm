"""Machine-specific paths for the Python ISRM runner — the Python mirror of
``run-model-jl-pushdown/paths.jl``. Every location is overridable by environment
variable so the runner is portable; nothing here is hardcoded to one machine.

IMPORTANT (this cluster): the scratch directory must be DISK-backed. ``/tmp`` is
a tmpfs here, so putting SR chunk blobs there consumes the same cgroup memory
budget the model needs and OOM-kills the job.
"""

from __future__ import annotations

import os

PY_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(PY_DIR)
CODE_ROOT = os.path.dirname(REPO)

#: The model driven by this runner. isrm_pushdown.esm is the only variant whose
#: emission-bearing source-cell set (`ppl`) is DERIVED by the graph — isrm.esm
#: takes it as a `src_cell_of_ppl` parameter, which would put the spatial join
#: outside the spec and break the strict "every number from the graph" bar.
MODEL = os.environ.get("ISRM_MODEL", os.path.join(REPO, "isrm_pushdown.esm"))

#: The Julia oracle runner — the source of the shared FF10 input zip.
RUNMODEL = os.environ.get("ISRM_RUNMODEL", os.path.join(REPO, "run-model-jl"))

EGU_ZIP = os.environ.get("EGU_ZIP", os.path.join(RUNMODEL, "data", "2016fd_inputs_point.zip"))

ZARR_URL = os.environ.get("ISRM_ZARR_URL", "s3://inmap-model/isrm_v1.2.1.zarr/")

#: SR source cells on the SR axis (== receptor cells == population cells).
N_SRC = int(os.environ.get("ISRM_N_SRC", "52411"))


def _resolve_scratch() -> str:
    """A DISK-backed scratch root. Mirrors paths.jl's resolution order."""
    env = os.environ.get("ISRM_SCRATCH")
    if env:
        return env
    user = os.environ.get("USER") or os.environ.get("LOGNAME") or "user"
    local = os.path.join("/scratch.local", user)
    if os.path.isdir(local):
        return os.path.join(local, "isrm-esm")
    scratch = os.environ.get("SCRATCH")
    if scratch and os.path.isdir(scratch):
        return os.path.join(scratch, "isrm-esm")
    import tempfile

    return os.path.join(tempfile.gettempdir(), "isrm-esm")


SCRATCH = _resolve_scratch()
