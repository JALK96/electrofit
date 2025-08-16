import argparse
import os
import logging
from pathlib import Path

from electrofit.config.loader import load_config, dump_config
from electrofit.external.gromacs import set_up_production
from electrofit.workflows.snapshot import build_snapshot_with_layers, CONFIG_ARG_HELP
from electrofit.logging import setup_logging, reset_logging, log_run_header


def _iter_final_dirs(project_root: Path):
    process = project_root / "process"
    if not process.is_dir():
        return
    for mol in process.iterdir():
        run_dir = mol / "run_final_gmx_simulation"
        if run_dir.is_dir():
            yield run_dir


def main():
    ap = argparse.ArgumentParser(description="Step8: Launch final GROMACS production runs from run_final_gmx_simulation dirs")
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    args = ap.parse_args()
    project_root = Path(args.project).resolve()
    setup_logging(str(project_root / "step8.log"), also_console=args.log_console)
    log_run_header("step8")
    # Multi-molecule context (affects logging annotations in mergers)
    process_root = project_root / "process"
    multi_mol = False
    if process_root.is_dir():
        pmols = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(pmols) > 1
        extra_override = Path(args.config).resolve() if getattr(args, "config", None) else None
    ran = 0
    for run_dir in _iter_final_dirs(project_root):
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        log_run_header("step8")
        mol = run_dir.parent.name
        build_snapshot_with_layers(
            run_dir,
            project_root,
            mol,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=run_dir.parent / "results" / "electrofit.toml",
            process_cfg=run_dir.parent / "electrofit.toml",
            molecule_input=project_root / "data" / "input" / mol / "electrofit.toml",
            project_defaults=project_root / "electrofit.toml",
            extra_override=extra_override,
        )
        cfg = load_config(project_root, context_dir=run_dir, molecule_name=mol)
        dump_config(cfg, log_fn=logging.info)
        sim = cfg.simulation
        box_type = sim.box.type
        d = sim.box.edge_nm
        cat = sim.ions.cation
        an = sim.ions.anion
        conc = sim.ions.concentration
        molecule = cfg.project.molecule_name or run_dir.parent.name
        base_scratch = cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch"
        # Forcefield now canonical under simulation section
        ff = getattr(getattr(cfg, 'simulation', None), 'forcefield', None) or "amber14sb.ff"
        runtime = cfg.gmx.runtime
        threads = runtime.threads
        pin = runtime.pin
        # Expect gro/top/itp present
        gro = next((f.name for f in run_dir.iterdir() if f.suffix == ".gro"), None)
        if not gro:
            logging.info(f"[step8][skip] {run_dir}: no .gro file")
            continue
        prev = os.getcwd()
        try:
            os.chdir(run_dir)
            set_up_production(
                m_gro=gro,
                MDP_dir="MDP",
                base_scratch_dir=base_scratch,
                molecule_name=molecule,
                box_type=box_type,
                cation=cat,
                anion=an,
                d=str(d),
                conc=str(conc),
                exit_screen=True,
                ff=ff,
                threads=threads,
                pin=pin,
            )
            ran += 1
        finally:
            os.chdir(prev)
    summary = f"[step8] Completed {ran} final simulation run(s)."
    print(summary)
    reset_logging()
    setup_logging(str(project_root / "step8.log"), also_console=args.log_console)
    logging.info(summary)


if __name__ == "__main__":  # pragma: no cover
    main()
