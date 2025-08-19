"""Deprecated shim for Step4.

Implementation moved to ``electrofit.pipeline.steps.step4``.
This module will be removed in a future release.
"""
from __future__ import annotations
import warnings as _warnings
from electrofit.pipeline.steps.step4 import main as _new_main  # noqa: F401

_warnings.warn(
    "electrofit.workflows.step4_extract_conforms is deprecated; use electrofit.pipeline.steps.step4",
    DeprecationWarning,
    stacklevel=2,
)

main = _new_main  # re-export

if __name__ == "__main__":  # pragma: no cover
    _new_main()
