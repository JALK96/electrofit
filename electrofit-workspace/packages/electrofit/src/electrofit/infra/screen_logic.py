"""Helper logic to resolve per-step screen termination decisions.

Centralizes combination of the global runtime flag (gmx.runtime.exit_screen)
with fine grained per-step overrides in Config.screen.stepX.

Rationale:
- Adapter-Ebene (z.B. GROMACS) soll keine Beendigungen mehr auslösen.
- Pipeline-Step Skripte rufen diesen Helper auf und führen ggf.
  terminate_session_if_active() aus.

Policy:
- Alle step*-Flags default False (opt-in) -> verhindert versehentliches
  Schließen von vom Benutzer manuell gestarteten screen Sitzungen.
- Effektives Ergebnis = global_flag AND step_flag.
- Kein Logging hier; der aufrufende Step kann auf Wunsch loggen.
"""
from __future__ import annotations

from typing import Any


def resolve_exit_screen(step: str, cfg: Any) -> bool:
    """Return True if screen session should be terminated after *step*.

    Parameters
    ----------
    step : str
        Step name wie "step3".
    cfg : Config
        Vollständig geladene Config Instanz.
    """
    try:
        runtime = getattr(getattr(cfg, 'gmx', None), 'runtime', None)
        global_flag = bool(getattr(runtime, 'exit_screen', True))
        screen_cfg = getattr(cfg, 'screen', None)
        step_flag = bool(getattr(screen_cfg, step, False))
        return global_flag and step_flag
    except Exception:  # defensive fallback
        return False

__all__ = ["resolve_exit_screen"]
