"""Deprecated shim for Step4.

Implementation moved to ``electrofit.pipeline.steps.step4``.
This module will be removed in a future release.
"""
from __future__ import annotations
import warnings as _warnings
from electrofit.pipeline.steps.step4 import main as _new_main  # noqa: F401
from electrofit.domain.sampling import select_frame_indices as _select_indices_impl  # new import

_warnings.warn(
    "electrofit.workflows.step4_extract_conforms is deprecated; use electrofit.pipeline.steps.step4",
    DeprecationWarning,
    stacklevel=2,
)

main = _new_main  # re-export

# Backwards compatibility for unit tests importing internal helper
def _select_indices(traj, n: int, method: str, seed):  # pragma: no cover - thin wrapper
    return _select_indices_impl(traj, n, method, seed)

if __name__ == "__main__":  # pragma: no cover
    _new_main()
