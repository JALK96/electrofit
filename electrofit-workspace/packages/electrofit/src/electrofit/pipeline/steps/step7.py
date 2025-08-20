"""Pipeline Step 7: Prepare final GROMACS production simulation directories.

Orchestrator (legacy workflows layer removed).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

from electrofit.domain.final_sim import prepare_final_sim_directory
from electrofit.infra.logging import setup_logging, reset_logging
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP

__all__ = ["main", "run_step7"]


def run_step7(project: Path, override_cfg: Path | None) -> int:
    process_dir = project / "process"
    if not process_dir.is_dir():
        print("[step7] No process directory.")
        return 0
    mol_dirs = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    multi_mol = len(mol_dirs) > 1
    done = 0
    for mol_dir in mol_dirs:
        dest_dir = mol_dir / "run_final_gmx_simulation"
        log_path = dest_dir / "process.log"
        dest_dir.mkdir(exist_ok=True)
        try:
            reset_logging()
            setup_logging(str(log_path), also_console=False)
        except Exception:  # pragma: no cover
            pass
        ok, msg = prepare_final_sim_directory(mol_dir, project, override_cfg, multi_mol)
        status = "[step7]" if ok else "[step7][skip]"
        print(f"{status} {mol_dir.name}: {msg}")
        if ok:
            done += 1
    print(f"[step7] Prepared final sim dirs for {done}/{len(mol_dirs)} molecules.")
    return 0


def main(argv: list[str] | None = None):  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step7: Prepare final simulation directories (copy topology, MDP, scripts).")
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    args = ap.parse_args(argv)
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    rc = run_step7(project_root, override_cfg)
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()
