"""Deprecated shim for Step6 legacy workflow.

Use ``electrofit.pipeline.steps.step6``. This module will be removed in a future release.
"""
from __future__ import annotations

import warnings
from electrofit.pipeline.steps.step6 import main, run_step6  # type: ignore F401

__all__ = ["main", "run_step6"]

warnings.warn(
    "electrofit.workflows.step6_extract_average_charges is deprecated; use electrofit.pipeline.steps.step6",
    DeprecationWarning,
    stacklevel=2,
)
