"""DEPRECATED shim for snapshot utilities.

This module is kept temporarily to avoid breaking imports. Logic moved to
`electrofit.infra.config_snapshot`.

Will be removed after the reorganisation deprecation window.
"""
from __future__ import annotations

import warnings
from electrofit.infra.config_snapshot import (
    compose_snapshot as build_snapshot_with_layers,  # backward alias
    CONFIG_ARG_HELP,
)

warnings.warn(
    "Import von 'electrofit.workflows.snapshot' ist veraltet; bitte 'electrofit.infra.config_snapshot' verwenden.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["build_snapshot_with_layers", "CONFIG_ARG_HELP"]
