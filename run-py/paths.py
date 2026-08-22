"""Machine-independent path resolution for the Phase-3 Python shim.

Everything derives from this file's own location (hence the isrm.esm checkout)
plus optional environment overrides, so the shim is portable across machines —
the Python mirror of run-jl/paths.jl.

  ISRM_MODEL       the .esm to drive       (default: <repo>/isrm_point.esm;
                   <repo>/isrm_polygon.esm is the area-source sibling. Both
                   import the shared templates from isrm_base.esm)
  ISRM_SCRATCH     bulk scratch root       (default: first writable of
                   /scratch.local/$USER, $SCRATCH, tempdir)
  EA_PATH          earthsci-ast-py checkout  (default: <sibling>/EarthSciAST/pkg/earthsci-ast-py)
  IO_PATH          EarthSciIO checkout       (default: <sibling>/EarthSciIO)
  DATA_DIR         local mirrors of fetchable sources (default: <repo>/data).
                   A data source whose url_template's basename names a file in
                   here is read from disk instead of the network — the FF10
                   point inventory, the example polygon layer. Absent file =
                   fetch the document's url_template through the cache.
  EGU_ZIP          DEPRECATED single-file alias, kept so an existing
                   EGU_ZIP=... invocation still works; it is folded in as one
                   more mirror, matched by its own basename.
  ISRM_ESIO_CACHE  EarthSciIO cache root     (default: $ISRM_SCRATCH/run-py-esio-cache;
                   point it at an existing cache — e.g. run-jl's — to reuse
                   already-fetched SR chunk blobs, the cache format is
                   cross-language)

SR chunk blobs land under the cache root, so it wants GBs of fast local
disk — NOT a network filesystem and NEVER a tmpfs /tmp (it eats the memory
cgroup). The scratch resolution below prefers /scratch.local/$USER for
exactly that reason; tempdir is a last resort.

EA_PATH / IO_PATH are only consulted when `earthsci_ast` / `earthsciio` are
not already importable (i.e. when no venv with editable installs is active):
their `src`/package dirs are appended to sys.path.
"""

from __future__ import annotations

import os
import sys
import tempfile

RUNPY_DIR = os.path.dirname(os.path.abspath(__file__))  # run-py/
REPO = os.path.dirname(RUNPY_DIR)  # the isrm.esm checkout
CODE_ROOT = os.path.dirname(REPO)  # where sibling checkouts live

MODEL = os.environ.get("ISRM_MODEL", os.path.join(REPO, "isrm_point.esm"))

EA_PATH = os.environ.get(
    "EA_PATH", os.path.join(CODE_ROOT, "EarthSciAST", "pkg", "earthsci-ast-py")
)
IO_PATH = os.environ.get("IO_PATH", os.path.join(CODE_ROOT, "EarthSciIO"))

DATA_DIR = os.environ.get("DATA_DIR", os.path.join(REPO, "data"))
EGU_ZIP = os.environ.get("EGU_ZIP", "")  # deprecated; see local_mirror


def local_mirror(url: str) -> str:
    """A local copy of ``url``, or ``""``.

    Mirroring is a LOCALITY choice of a run, never a property of the model, so
    it is resolved by matching the document's own ``url_template`` basename
    against ``DATA_DIR`` rather than by naming any particular source here. That
    is what lets one shim mirror the FF10 point inventory for isrm_point.esm and
    the example polygon layer for isrm_polygon.esm without knowing either name.
    """
    base = os.path.basename(url.rstrip("/"))
    if not base:
        return ""
    if EGU_ZIP and os.path.basename(EGU_ZIP) == base and os.path.isfile(EGU_ZIP):
        return EGU_ZIP
    cand = os.path.join(DATA_DIR, base)
    return cand if os.path.isfile(cand) else ""


def _resolve_scratch() -> str:
    """First writable scratch root among the candidates, as <root>/isrm-esm."""
    env = os.environ.get("ISRM_SCRATCH", "")
    if env:
        os.makedirs(env, exist_ok=True)
        return env
    user = os.environ.get("USER", os.environ.get("LOGNAME", "user"))
    for base in (f"/scratch.local/{user}", os.environ.get("SCRATCH", ""), tempfile.gettempdir()):
        if not base:
            continue
        candidate = os.path.join(base, "isrm-esm")
        try:
            os.makedirs(candidate, exist_ok=True)
            probe = os.path.join(candidate, ".probe")
            with open(probe, "w") as fh:
                fh.write("x")
            os.remove(probe)
            return candidate
        except OSError:
            continue
    raise RuntimeError("no writable scratch root found; set ISRM_SCRATCH")


SCRATCH = _resolve_scratch()

ESIO_CACHE = os.environ.get("ISRM_ESIO_CACHE", os.path.join(SCRATCH, "run-py-esio-cache"))


def ensure_importable() -> None:
    """Make earthsci_ast / earthsciio importable, preferring an active venv's
    installs; fall back to the EA_PATH / IO_PATH checkouts."""
    try:
        import earthsci_ast  # noqa: F401
    except ImportError:
        sys.path.insert(0, os.path.join(EA_PATH, "src"))
    try:
        import earthsciio  # noqa: F401
    except ImportError:
        sys.path.insert(0, IO_PATH)


if not os.path.isfile(MODEL):
    print(f"WARNING: model file not found at {MODEL} (set ISRM_MODEL)", file=sys.stderr)
