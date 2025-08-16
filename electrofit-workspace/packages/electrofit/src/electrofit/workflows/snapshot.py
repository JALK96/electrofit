"""Snapshot utility helpers for electrofit workflow steps.

Provides common functions to (a) ensure an electrofit.toml snapshot exists in a
working directory, and (b) merge an override config into that snapshot with
consistent multi-molecule aware logging.

Exposes a shared ``CONFIG_ARG_HELP`` constant (German) so all workflow CLIs
describe ``--config`` identically.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Callable

from electrofit.config.merger import merge_into_snapshot, fill_missing_into_snapshot

# Shared help text for --config flag (imported by individual step scripts)
CONFIG_ARG_HELP = (
    "Pfad zu einer zusätzlichen electrofit.toml: Molekül- und Prozess-Dateien überschreiben Projekt-Defaults; "
    "diese Datei wirkt als zusätzliche stärkste Override-Ebene (ändert bestehende Werte, füllt nichts nur auf)."
)


def _ensure_snapshot(target_dir: Path, precedence: Iterable[Path]) -> Path | None:
    snap = target_dir / "electrofit.toml"
    if snap.is_file():
        return snap
    for cand in precedence:
        if cand.is_file():
            try:
                snap.write_text(cand.read_text())
            except Exception:
                pass
            break
    return snap if snap.exists() else None


def build_snapshot_with_layers(
    run_dir: Path,
    project_root: Path,
    molecule: str | None,
    multi_molecule: bool,
    log_fn: Callable[[str], None],
    upstream: Path | None = None,
    process_cfg: Path | None = None,
    molecule_input: Path | None = None,
    project_defaults: Path | None = None,
    extra_override: Path | None = None,
) -> Path | None:
    """Construct or update a snapshot in run_dir using unified precedence.

    Order:
      Base: existing snapshot OR first existing among [upstream, molecule_input, process_cfg, project_defaults].
      Strong overrides (replace values): molecule_input -> process_cfg -> extra_override.
      Fill-only: project_defaults (adds missing keys only).

    Returns snapshot path or None if creation failed.
    """
    # 1. Ensure snapshot exists (using broad precedence list)
    precedence: list[Path] = []
    snap_existing = run_dir / "electrofit.toml"
    precedence.append(snap_existing)
    if upstream: precedence.append(upstream)
    if molecule_input: precedence.append(molecule_input)
    if process_cfg: precedence.append(process_cfg)
    if project_defaults: precedence.append(project_defaults)
    snap = _ensure_snapshot(run_dir, precedence)
    if not snap or not snap.is_file():
        return None

    # 2. Strong overrides
    for layer in [molecule_input, process_cfg, extra_override]:
        if layer and layer.is_file():
            try:
                merge_into_snapshot(snap, layer, multi_molecule=multi_molecule, log_fn=log_fn)
            except Exception:
                pass

    # 3. Fill missing from project defaults
    if project_defaults and project_defaults.is_file():
        try:
            fill_missing_into_snapshot(snap, project_defaults, log_fn=log_fn)
        except Exception:
            pass
    return snap
