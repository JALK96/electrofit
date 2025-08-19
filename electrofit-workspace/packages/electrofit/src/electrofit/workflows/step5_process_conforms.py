"""Deprecated Step5 orchestrator module.

Implementation relocated to ``electrofit.pipeline.steps.step5``. This shim
retains the historical import path for external users and emits a one-time
DeprecationWarning.
"""
from __future__ import annotations

import warnings
from electrofit.pipeline.steps.step5 import run_step5, main as _real_main

__all__ = ["run_step5", "main"]

_WARNED = False


def main(argv=None):  # pragma: no cover
    global _WARNED
    if not _WARNED:
        warnings.warn(
            "electrofit.workflows.step5_process_conforms is deprecated; use electrofit.pipeline.steps.step5",
            DeprecationWarning,
            stacklevel=2,
        )
        _WARNED = True
    return _real_main(argv)
