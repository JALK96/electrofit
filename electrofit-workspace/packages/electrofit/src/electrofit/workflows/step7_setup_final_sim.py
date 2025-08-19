"""Deprecated shim for Step7 legacy workflow.

Use ``electrofit.pipeline.steps.step7``.
"""
from __future__ import annotations

import warnings
from electrofit.pipeline.steps.step7 import main, run_step7  # type: ignore F401

__all__ = ["main", "run_step7"]

warnings.warn(
    "electrofit.workflows.step7_setup_final_sim is deprecated; use electrofit.pipeline.steps.step7",
    DeprecationWarning,
    stacklevel=2,
)
