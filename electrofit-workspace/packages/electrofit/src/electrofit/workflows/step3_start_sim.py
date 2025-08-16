import argparse
import json
import os
import logging
from pathlib import Path

from electrofit.config.loader import load_config, dump_config
from electrofit.external.gromacs import set_up_production
from electrofit.workflows.snapshot import build_snapshot_with_layers, CONFIG_ARG_HELP
from electrofit.logging import setup_logging, reset_logging, log_run_header


def _iter_run_dirs(project_root: Path):
    """Yield all run_gmx_simulation directories under project/process/*/."""
    process = project_root / "process"
    if not process.is_dir():
        return
    for mol_dir in process.iterdir():
        run_dir = mol_dir / "run_gmx_simulation"
        if run_dir.is_dir():
            yield run_dir


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Step3: Setup GROMACS production directories based on run.json & electrofit.toml. "
            "Automates box/neutralization/ion concentration and writes required inputs."
        )
    )
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()), help="Path to the project root (the folder that contains 'process/')")
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    args = ap.parse_args()

    project_root = Path(args.project).resolve()
    # Determine multi-molecule context
    process_root = project_root / "process"
    multi_mol = False
    if process_root.is_dir():
        pmols = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(pmols) > 1

    # Base project-level config (for defaults only)
    # Global logging (summary only). Individual run dirs get their own logs.
    setup_logging(str(project_root / "step3.log"), also_console=True)
    try:
        log_run_header("step3")
    except Exception:  # pragma: no cover
        logging.info("electrofit unknown | step=step3")
    cfg = load_config(project_root)
    dump_config(cfg, log_fn=logging.info)

    # Pull simulation knobs from config with safe defaults
    sim = getattr(cfg, "simulation", None)
    ff       = (getattr(sim, "forcefield", None)          or "amber14sb.ff")
    box_type = (getattr(sim, "box_type", None)            or "dodecahedron")
    d        = (getattr(sim, "box_edge_distance", None)   or 1.2)
    cation   = (getattr(sim, "cation", None)              or "NA")
    anion    = (getattr(sim, "anion", None)               or "CL")
    conc     = (getattr(sim, "ion_concentration", None)   or 0.15)

    # Resolve & expand scratch base (allow ${USER} etc.)
    base_scratch = cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch"
    base_scratch = os.path.expanduser(os.path.expandvars(base_scratch))

    # Runtime flags for gmx (optional)
    runtime = getattr(getattr(cfg, "gmx", None), "runtime", None)
    threads = getattr(runtime, "threads", None) if runtime else None
    pin     = getattr(runtime, "pin", None) if runtime else None

    ran = 0
    for run_dir in _iter_run_dirs(project_root):
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        try:
            log_run_header("step3")
        except Exception:  # pragma: no cover
            logging.info("electrofit unknown | step=step3")
        mol = run_dir.parent.name
        upstream = run_dir.parent / "run_gau_create_gmx_in" / "electrofit.toml"
        molecule_input = project_root / "data" / "input" / mol / "electrofit.toml"
        process_cfg = run_dir.parent / "electrofit.toml"
        project_defaults = project_root / "electrofit.toml"
        snap = build_snapshot_with_layers(
            run_dir,
            project_root,
            mol,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=upstream,
            process_cfg=process_cfg,
            molecule_input=molecule_input,
            project_defaults=project_defaults,
            extra_override=Path(args.config) if args.config else None,
        )
        manifest = run_dir / "run.json"
        if not manifest.is_file():
            msg = f"[step3][skip] no run.json in {run_dir}"
            logging.info(msg)
            print(msg)
            continue

        with manifest.open("r") as f:
            meta = json.load(f)

        # Pull required fields with small fallbacks
        molecule = meta.get("molecule")
        gro_file = meta.get("gro")
        mdp_dir  = meta.get("mdp_dir") or (getattr(cfg.paths, "mdp_dir", None) or "MDP")

        if not molecule or not gro_file:
            msg = f"[step3][skip] Skip: incomplete manifest ({manifest})"
            logging.info(msg)
            print(msg)
            continue

        # Ensure files exist in the run directory
        gro_path = run_dir / gro_file
        mdp_path = run_dir / mdp_dir
        if not gro_path.is_file():
            msg = f"[step3][skip] Skip: GRO file missing {gro_path}"
            logging.info(msg)
            print(msg)
            continue
        # Require a matching TOP file (same stem) to proceed; skip placeholder runs
        top_path = gro_path.with_suffix('.top')
        if not top_path.is_file():
            msg = f"[step3][skip] Skip: TOP file missing {top_path}"
            logging.info(msg)
            print(msg)
            continue
        if not mdp_path.is_dir():
            msg = f"[step3][skip] MDP dir missing: {mdp_path}"
            logging.info(msg)
            print(msg)
            continue

        logging.info(f"[step3] Starting GROMACS production for {molecule} in {run_dir}")

        # Change working directory to the run directory so set_up_production finds inputs.
        prev_cwd = os.getcwd()
        try:
            os.chdir(run_dir)
            set_up_production(
                m_gro=gro_file,
                MDP_dir=mdp_dir,
                base_scratch_dir=base_scratch,
                molecule_name=molecule,
                box_type=box_type,
                cation=cation,
                anion=anion,
                d=str(d),
                conc=str(conc),
                exit_screen=True,
                ff=ff,
                threads=threads,
                pin=pin,
            )
            ran += 1
        finally:
            os.chdir(prev_cwd)

    summary = f"[step3] Completed {ran} run(s)." if ran else "[step3] No runs executed (no manifests found)."
    # Emit summary to stdout and append to global log
    print(summary)
    reset_logging()
    setup_logging(str(project_root / "step3.log"), also_console=True)
    logging.info(summary)


if __name__ == "__main__":
    main()
