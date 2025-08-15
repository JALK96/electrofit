import argparse
import json
import os
from pathlib import Path

from electrofit.config.loader import load_config
from electrofit.external.gromacs import set_up_production


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
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--project",
        default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()),
        help="Path to the project root (the folder that contains 'process/')",
    )
    ap.add_argument(
        "--config",
        help="Optional explicit path to electrofit.toml; defaults to project/electrofit.toml",
    )
    args = ap.parse_args()

    project_root = Path(args.project).resolve()
    # Prefer explicit config, else project-level TOML, else let loader decide
    cfg_path = Path(args.config).resolve() if args.config else None
    cfg = load_config(project_root, config_path=cfg_path)

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
        manifest = run_dir / "run.json"
        if not manifest.is_file():
            print(f"[step3] Skip: no run.json in {run_dir}")
            continue

        with manifest.open("r") as f:
            meta = json.load(f)

        # Pull required fields with small fallbacks
        molecule = meta.get("molecule")
        gro_file = meta.get("gro")
        mdp_dir  = meta.get("mdp_dir") or (getattr(cfg.paths, "mdp_dir", None) or "MDP")

        if not molecule or not gro_file:
            print(f"[step3] Skip: incomplete manifest in {manifest}")
            continue

        # Ensure files exist in the run directory
        gro_path = run_dir / gro_file
        mdp_path = run_dir / mdp_dir
        if not gro_path.is_file():
            print(f"[step3] Skip: GRO file missing: {gro_path}")
            continue
        if not mdp_path.is_dir():
            print(f"[step3] Skip: MDP dir missing: {mdp_path}")
            continue

        print(f"[step3] Starting GROMACS production for {molecule} in {run_dir}")

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

    if ran == 0:
        print("[step3] No runs executed (no manifests found).")
    else:
        print(f"[step3] Completed {ran} run(s).")


if __name__ == "__main__":
    main()
