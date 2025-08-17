import argparse
import fnmatch
import os
import shutil
import logging
from pathlib import Path

from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, reset_logging


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    args = ap.parse_args()
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        print("[step7] No process directory.")
        return
    mdp_source_dir = project_root / "data" / "MDP"
    bash_script_source = project_root / "scripts" / "gmx.sh"
    file_patterns = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]
    mol_dirs = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    multi_mol = len(mol_dirs) > 1
    done = 0
    for mol_dir in mol_dirs:
        results_dir = mol_dir / "results"
        if not results_dir.is_dir():
            print(f"[step7][skip] {mol_dir.name}: no results dir")
            continue
        dest_dir = mol_dir / "run_final_gmx_simulation"
        dest_dir.mkdir(exist_ok=True)
        # Set up per-molecule logging (file only, mirror step6 convention)
        log_path = dest_dir / "process.log"
        try:
            setup_logging(str(log_path), also_console=False)
        except Exception:
            pass
    # Snapshot propagation with unified precedence (mirrors step6):
    #   1. results/electrofit.toml
    #   2. process/<mol>/electrofit.toml
    #   3. data/input/<mol>/electrofit.toml
    #   4. project_root/electrofit.toml
        snap_dest = compose_snapshot(
            dest_dir,
            project_root,
            mol_dir.name,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=results_dir / "electrofit.toml",
            process_cfg=mol_dir / "electrofit.toml",
            molecule_input=project_root / "data" / "input" / mol_dir.name / "electrofit.toml",
            project_defaults=project_root / "electrofit.toml",
            extra_override=override_cfg,
        )
        cfg = load_config(project_root, context_dir=dest_dir, molecule_name=mol_dir.name)
        try:
            dump_config(cfg, log_fn=logging.info)
            for h in logging.getLogger().handlers:
                try:
                    h.flush()
                except Exception:
                    pass
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step7] config dump failed: {e}")
        # Locate acpype in results
        acpype_dir = None
        for sub in results_dir.iterdir():
            if sub.is_dir() and sub.name.endswith(".acpype"):
                acpype_dir = sub; break
        if not acpype_dir:
            print(f"[step7][skip] {mol_dir.name}: no .acpype dir in results")
            continue
        for fn in os.listdir(acpype_dir):
            for pat in file_patterns:
                if fnmatch.fnmatch(fn, pat):
                    shutil.copy(acpype_dir / fn, dest_dir)
                    logging.info(f"[step7] Copied {fn} -> {dest_dir}")
                    break
        md_dest = dest_dir / "MDP"
        if mdp_source_dir.is_dir():
            if md_dest.exists():
                shutil.rmtree(md_dest)
            shutil.copytree(mdp_source_dir, md_dest)
            logging.info(f"[step7] Copied MDP -> {md_dest}")
        else:
            print(f"[step7][warn] no MDP source {mdp_source_dir}")
        if bash_script_source.is_file():
            shutil.copy(bash_script_source, dest_dir / "gmx.sh")
        reset_logging()
        done += 1
    print(f"[step7] Prepared final sim dirs for {done}/{len(mol_dirs)} molecules.")


if __name__ == "__main__":  # pragma: no cover
    main()
