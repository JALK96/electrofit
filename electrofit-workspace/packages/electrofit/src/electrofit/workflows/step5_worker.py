"""Deprecated worker shim for Step5.

Process logic moved to ``electrofit.domain.charges.conformer_batch``. This
module preserves the import path & function name ``process_one`` by forwarding
to the new implementation while issuing a deprecation warning once.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from electrofit.domain.charges.conformer_batch import process_conformer_dir

__all__ = ["process_one"]

_WARNED = False


def process_one(
    conf_dir_str: str,
    project_root_str: str,
    override_cfg_path: str | None,
    multi_mol: bool,
    mock: bool,
    verbose: bool,
) -> tuple[str, bool, str]:  # pragma: no cover - thin wrapper
    global _WARNED
    if not _WARNED:
        warnings.warn(
            "electrofit.workflows.step5_worker.process_one is deprecated; use electrofit.domain.charges.conformer_batch.process_conformer_dir",
            DeprecationWarning,
            stacklevel=2,
        )
        _WARNED = True
    return process_conformer_dir(
        Path(conf_dir_str),
        Path(project_root_str),
        Path(override_cfg_path) if override_cfg_path else None,
        multi_mol,
        mock,
        verbose,
    )
