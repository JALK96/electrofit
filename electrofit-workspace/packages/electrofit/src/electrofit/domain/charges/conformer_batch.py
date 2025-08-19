"""Domain helpers for batch processing of conformer directories (Step5).

This module extracts the discovery + per‑conformer processing logic from the
previous ``workflows.step5_*`` modules so the new orchestrator
``pipeline.steps.step5`` can remain thin. The legacy public symbol name
``process_one`` is intentionally preserved via a shim in
``workflows.step5_worker`` for backwards compatibility with any external
scripts/tests that still import it.

Responsibilities:
  * Discover conformer directories (heuristic: directories containing a PDB).
  * Execute per‑conformer processing (config snapshot, logging, mock mode,
    domain charge pipeline).

Notes:
  * Heavy Gaussian/RESP work is delegated to the existing core shim
    ``core.process_conform`` which already manages scratch setup and finalise
    semantics. A future refactor can inline or adapt that to use
    ``domain.charges.process_conformer`` directly once scratch handling is
    extracted cleanly.
  * Advanced diagnostic monkey‑patching from the legacy worker is retained but
    guarded so it runs only once per process.
"""
from __future__ import annotations

from pathlib import Path
import os
import logging
import traceback
from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.infra.logging import log_run_header, reset_logging, setup_logging
from electrofit.core.process_conform import process_conform  # legacy shim -> domain
from electrofit.io.files import find_file_with_extension, strip_extension

__all__ = [
    "discover_conformer_dirs",
    "process_conformer_dir",
]


def discover_conformer_dirs(project_root: Path) -> list[Path]:
    """Return conformer directories containing at least one PDB file.

    Mirrors historical logic from ``step5_process_conforms._discover_conformer_dirs``.
    Empty list if project/process layout not present.
    """
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        return []
    dirs: list[Path] = []
    for mol_dir in sorted(process_dir.iterdir()):
        ec = mol_dir / "extracted_conforms"
        if not ec.is_dir():
            continue
        for conf in sorted(ec.iterdir()):
            if not conf.is_dir():
                continue
            if any(conf.glob("*.pdb")):
                dirs.append(conf)
    return dirs


_ADV_DIAG_INITIALISED = False


def _init_advanced_diagnostics():  # pragma: no cover - diagnostic path
    global _ADV_DIAG_INITIALISED
    if _ADV_DIAG_INITIALISED:
        return
    if os.environ.get("ELECTROFIT_DISABLE_ADV_DIAG") == "1":  # opt-out switch
        _ADV_DIAG_INITIALISED = True
        return
    try:
        import sys as _sys_d, faulthandler as _fh_d, signal as _sig_d, atexit as _ax_d, gc as _gc_d
        try:
            _fh_d.enable(file=_sys_d.stderr)
        except Exception:
            pass

        def _ax():
            try:
                print('[debug-atexit] worker normal shutdown', flush=True)
            except Exception:
                pass

        try:
            _ax_d.register(_ax)
        except Exception:
            pass

        def _sig_handler(sig, frame):  # noqa: ARG001
            try:
                print(f'[debug-signal] sig={sig}', flush=True)
                _fh_d.dump_traceback(file=_sys_d.stderr)
            except Exception:
                pass

        for _s in [getattr(_sig_d, n, None) for n in ('SIGTERM','SIGINT','SIGSEGV')]:
            if _s is not None:
                try:
                    _sig_d.signal(_s, _sig_handler)
                except Exception:
                    pass
        if hasattr(_gc_d, 'callbacks'):
            try:
                _gc_d.callbacks.append(lambda phase, info: phase=='stop' and print('[debug-gc] cycle', flush=True))  # type: ignore
            except Exception:
                pass
    except Exception:
        pass
    _ADV_DIAG_INITIALISED = True


def process_conformer_dir(
    conf_dir: Path,
    project_root: Path,
    override_cfg_path: Path | None,
    multi_mol: bool,
    mock: bool,
    verbose: bool,
) -> tuple[str, bool, str]:
    """Process a single conformer directory.

    Returns (relative_path, ok_flag, message).
    This is a near‑lift of ``workflows.step5_worker.process_one`` with minor
    cleanups and a direct call to ``core.process_conform`` for parity.
    """
    prev = os.getcwd()
    try:
        os.chdir(conf_dir)
        _init_advanced_diagnostics()
        # Monkeypatch sys/os exit for diagnostics (kept verbatim)
        import sys as _sys, os as _os  # local
        _orig_sys_exit = _sys.exit
        _orig_os_exit = _os._exit

        def _dbg_sys_exit(code=0):  # noqa: D401
            try:
                print(f"[debug-sys-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_sys_exit(code)

        def _dbg_os_exit(code=0):  # noqa: D401
            try:
                print(f"[debug-os-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_os_exit(code)

        _sys.exit = _dbg_sys_exit  # type: ignore
        _os._exit = _dbg_os_exit   # type: ignore

        # per‑conformer logging
        reset_logging()
        log_path = Path(os.getcwd()) / "process.log"
        setup_logging(str(log_path), also_console=verbose)
        log_run_header("step5")

        pdb_file = find_file_with_extension("pdb")
        if not pdb_file:
            raise FileNotFoundError("No PDB file found in conformer directory")
        mol_name = conf_dir.parent.parent.name if len(conf_dir.parts) >= 2 else conf_dir.parent.name
        snap_path = compose_snapshot(
            conf_dir,
            project_root,
            mol_name,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=conf_dir.parent / "electrofit.toml",  # parent is extracted_conforms root
            process_cfg=conf_dir.parent.parent / "electrofit.toml" if len(conf_dir.parts) >= 3 else None,
            molecule_input=project_root / "data" / "input" / mol_name / "electrofit.toml",
            project_defaults=project_root / "electrofit.toml",
            extra_override=override_cfg_path,
        )
        if not snap_path:
            logging.warning(f"[step5][warn] snapshot could not be created in {conf_dir}")
        cfg = load_config(project_root, context_dir=conf_dir)
        dump_config(cfg, header=True, log_fn=logging.info)
        proj = cfg.project
        molecule_name = proj.molecule_name or strip_extension(pdb_file)
        if mock:
            with open("executed.txt", "w") as f:
                f.write(f"run{molecule_name}")
        else:
            # Delegate heavy pipeline
            process_conform(
                molecule_name=molecule_name,
                pdb_file=pdb_file,
                base_scratch_dir=cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch",
                net_charge=proj.charge or 0,
                residue_name=proj.residue_name or "LIG",
                adjust_sym=getattr(proj, "adjust_symmetry", False),
                protocol=proj.protocol or "bcc",
                ignore_sym=getattr(proj, "ignore_symmetry", False),
            )
        rel = conf_dir.relative_to(project_root)
        print(f"[worker-return] {rel} ok", flush=True)
        return (str(rel), True, "ok")
    except Exception as e:  # pragma: no cover - defensive
        try:
            tb = traceback.format_exc()
            logging.error(f"Worker failure {conf_dir}: {e}\n{tb}")
        except Exception:
            pass
        try:
            rel = conf_dir.relative_to(project_root)
            rel_str = str(rel)
        except Exception:
            rel_str = str(conf_dir)
        print(f"[worker-exc] {rel_str}: {e}", flush=True)
        return (rel_str, False, str(e))
    finally:
        try:  # restore monkeypatches
            import sys as _sys, os as _os
            if '_orig_sys_exit' in locals():
                _sys.exit = _orig_sys_exit  # type: ignore
            if '_orig_os_exit' in locals():
                _os._exit = _orig_os_exit   # type: ignore
        except Exception:
            pass
        try:
            os.chdir(prev)
        except Exception:
            pass
