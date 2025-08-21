"""Utilities for unified optional termination of an enclosing GNU screen session.

Public entry point: maybe_terminate_screen(step, cfg)

Design:
- Adapter code must not terminate screen sessions.
- Pipeline steps decide post-success whether to terminate based on config.
- Logic combines global flag (gmx.runtime.exit_screen) and per-step screen.stepX flag
  via resolve_exit_screen().
- No config loading here: caller already has cfg.
- Silent by default (no info log); termination is a deliberate, opt-in action,
  we avoid log noise. Debug logging can be added later if needed.
"""
from __future__ import annotations

import os
from typing import Any
from electrofit.infra.screen_logic import resolve_exit_screen
from electrofit.infra.screen import terminate_session_if_active

__all__ = ["maybe_terminate_screen"]


def maybe_terminate_screen(step: str, cfg: Any) -> None:
    """Terminate active screen session if configuration requests it.

    Parameters
    ----------
    step : str
        Step name (e.g. "step3").
    cfg : Any
        Loaded configuration object for current run context.
    """
    if not cfg:
        return
    try:
        if resolve_exit_screen(step, cfg):
            session = os.environ.get("STY", "")
            if session:  # only attempt if inside screen
                terminate_session_if_active(session)
    except Exception:
        # Fail silent; termination is optional and should not break pipeline
        return
