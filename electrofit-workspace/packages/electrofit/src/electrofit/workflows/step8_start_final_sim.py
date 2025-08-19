"""Deprecated shim for Step8 legacy workflow. Use ``electrofit.pipeline.steps.step8``."""
from __future__ import annotations

import warnings
from electrofit.pipeline.steps.step8 import main, run_step8  # type: ignore F401

__all__ = ["main", "run_step8"]

warnings.warn(
    "electrofit.workflows.step8_start_final_sim is deprecated; use electrofit.pipeline.steps.step8",
    DeprecationWarning,
    stacklevel=2,
)
