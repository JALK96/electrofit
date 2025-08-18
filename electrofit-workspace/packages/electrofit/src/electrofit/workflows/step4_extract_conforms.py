import os
import argparse
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import random
import sys
from typing import List

import mdtraj as md
import numpy as np
from tqdm import tqdm
import logging

try:  # Python 3.11+ has tomllib
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml  # type: ignore
from electrofit.config.loader import load_config, dump_config
from electrofit.io.files import mol2_to_pdb_with_bonds
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_sampling_decision


def _select_indices(traj: md.Trajectory, n: int, method: str, seed: int | None) -> List[int]:
    """Select frame indices according to a sampling method.

    Methods:
      linear  : Evenly spaced frames (deterministic)
      random  : Uniform random without replacement (seeded)
      maxmin  : Farthest point (max-min RMSD) selection (deterministic given seed for start)
    """
    total = len(traj)
    if total == 0 or n <= 0:
        return []
    if n >= total:
        return list(range(total))

    method = method.lower()
    if method in {"linear", "even", "linspace"}:
        return list(np.linspace(0, total - 1, num=n, dtype=int))

    if method in {"random", "rand"}:
        rng = random.Random(seed)
        return rng.sample(range(total), n)

    if method in {"maxmin", "minmax", "farthest"}:
        # Farthest point sampling in RMSD space
        # Start frame choice (seeded): either 0 or seeded random
        if seed is not None:
            rng = random.Random(seed)
            first = rng.randrange(total)
        else:
            first = 0
        selected = [first]
        # Distances to current set: start with distances to first
        # md.rmsd takes reference index; we compute distances to 'first'
        min_d = md.rmsd(traj, traj, first)  # shape (total,)
        # Ensure selected distance set to -inf so it's never re-picked
        min_d[first] = -np.inf
        for _ in range(1, n):
            # Pick frame with maximal current min distance
            next_idx = int(np.argmax(min_d))
            selected.append(next_idx)
            if len(selected) == n:
                break
            # Update min distances with distances to new point
            d_new = md.rmsd(traj, traj, next_idx)
            # new minimal distance to selected set
            min_d = np.minimum(min_d, d_new)
            min_d[next_idx] = -np.inf
        return selected

    # Fallback
    return list(np.linspace(0, total - 1, num=n, dtype=int))


def _extract_for_molecule(mol_proc_dir: Path, project_root: Path, sample: int, method: str, seed: int | None, override_cfg: Path | None, multi_mol: bool, verbose: bool):
    """Extract conformers for one molecule directory.

    Legacy .ef parsing removed: configuration is resolved from layered TOMLs
    via load_config (project-level + per-molecule input + per-run context).
    """
    sim_dir = mol_proc_dir / "run_gmx_simulation"
    pis_dir = mol_proc_dir / "run_gau_create_gmx_in"
    if not sim_dir.is_dir():
        return False, "no sim dir"
    # Resolve config (context=sim_dir) to pick up molecule-specific overrides
    cfg = load_config(project_root, context_dir=sim_dir, molecule_name=mol_proc_dir.name)
    proj = cfg.project
    molecule_name = proj.molecule_name or mol_proc_dir.name
    residue_name = proj.residue_name or "LIG"
    adjust_sym = getattr(proj, "adjust_symmetry", False)
    protocol = getattr(proj, "protocol", "bcc")

    # Optional resp/symmetry inputs
    respin1_file = respin2_file = equiv_groups_file = None
    if protocol == "opt":
        respin1_file = pis_dir / ("ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN")
        respin2_file = pis_dir / "ANTECHAMBER_RESP2.IN"
    elif protocol == "bcc" and adjust_sym:
        # Look explicitly inside the per-molecule preparation dir (pis_dir) for symmetry groups
        json_candidates = sorted(pis_dir.glob("*.json"))
        if json_candidates:
            equiv_groups_file = json_candidates[0]
        else:
            # Warn only into per-molecule log (will set up logging shortly)
            pass

    input_mol2_file = project_root / "data" / "input" / mol_proc_dir.name / f"{molecule_name}.mol2"
    if not input_mol2_file.is_file():
        logging.warning(f"[step4][{mol_proc_dir.name}] missing input mol2 ({input_mol2_file.name}); continuing (bond insertion skipped)")
    extracted_conforms_dir = mol_proc_dir / "extracted_conforms"
    extracted_conforms_dir.mkdir(exist_ok=True)
    # Initialise logging early (per molecule) into extracted_conforms/process.log WITHOUT console handler
    log_path = extracted_conforms_dir / "process.log"
    setup_logging(str(log_path), also_console=False)
    # Now that logging is configured only to file, emit warnings / config snapshot info
    if protocol == "bcc" and adjust_sym and not equiv_groups_file:
        logging.warning(f"adjust_symmetry requested but no *.json in {pis_dir}")
    existing_snapshot = extracted_conforms_dir / "electrofit.toml"
    if existing_snapshot.is_file():
        try:
            cfg_existing = load_config(project_root, context_dir=extracted_conforms_dir, molecule_name=mol_proc_dir.name)
            logging.info(f"[config] existing extracted_conforms snapshot for {mol_proc_dir.name} detected -> dump below")
            dump_config(cfg_existing, log_fn=logging.info)
        except Exception as e:  # pragma: no cover
            logging.warning(f"failed to dump existing extracted_conforms snapshot: {e}")

    # Structured decision (sampling) + minimal relevant config logging
    try:
        symmetry_json_present = any(pis_dir.glob('*.json'))
        build_sampling_decision(
            protocol=protocol,
            adjust_sym=adjust_sym,
            ignore_sym=getattr(proj, 'ignore_symmetry', False),
            sampling_method=method,
            sample_count=sample,
            seed=seed,
            symmetry_json_present=symmetry_json_present,
        ).log('step4')
        log_relevant_config('step4', proj, ['molecule_name','residue_name','protocol','adjust_symmetry','ignore_symmetry'])
    except Exception:
        logging.debug('[step4][decisions] logging failed', exc_info=True)

    traj_path = sim_dir / "md_center.xtc"
    gro_path = sim_dir / "md.gro"
    if not traj_path.is_file() or not gro_path.is_file():
        return False, "missing md_center.xtc or md.gro"

    # Ensure snapshot at extracted_conforms root
    parent_cfg_target = compose_snapshot(
        extracted_conforms_dir,
        project_root,
        mol_proc_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        upstream=sim_dir / "electrofit.toml",
        process_cfg=mol_proc_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_proc_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    # (Legacy behaviour copied the .ef here; no longer needed with TOML-based config)

    # Log start & resolved layered config (now already have logging configured)
    logging.info(f"Step4: start extraction for molecule_dir={mol_proc_dir.name} method={method} sample={sample}")
    dump_config(cfg, log_fn=logging.info)
    raw_traj = md.load(str(traj_path), top=str(gro_path))
    # Residue inventory logging (diagnostics for mismatches)
    try:
        res_counts: dict[str, int] = {}
        for res in raw_traj.topology.residues:
            # Count atoms belonging to residue
            res_counts[res.name] = res_counts.get(res.name, 0) + len(res.atoms)
        inv_str = ", ".join(f"{k}:{v}" for k, v in sorted(res_counts.items()))
        logging.info(f"[step4][{mol_proc_dir.name}] residue inventory -> {inv_str}")
    except Exception as e:  # pragma: no cover
        logging.debug(f"[step4][{mol_proc_dir.name}] residue inventory logging failed: {e}")
    ipl = raw_traj.top.select(f"resname {residue_name}")
    if len(ipl) == 0:
        # Strikter: statt stiller Fallback jetzt harter Abbruch, damit Residueninkonsistenz früh sichtbar wird
        logging.error(f"[step4][{mol_proc_dir.name}] residue '{residue_name}' not found in trajectory topology; abort extraction (fix residue_name or input files).")
        return False, f"residue '{residue_name}' not in topology"
    traj = raw_traj.atom_slice(ipl)
    if traj.n_atoms == 0:
        logging.warning(f"[step4][{mol_proc_dir.name}] zero atoms after selection; skipping")
        return False, "no atoms after selection"

    total = len(traj)
    n = min(sample, total)
    indices = _select_indices(traj, n, method, seed)
    configs = [traj[i] for i in indices]
    logging.info(f"Selected indices (n={len(indices)}) -> {indices}")

    # Legacy pc.sh distribution removed; we no longer copy execution shell scripts.
    bash_script_path = project_root / "electrofit" / "bash" / "pc.sh"
    for i, c in enumerate(configs):
        conform_dir = extracted_conforms_dir / f"{molecule_name}c{i}"
        conform_dir.mkdir(exist_ok=True)
        logging.info(f"Creating conformer dir {conform_dir.name}")
        # Ensure each conformer has an up-to-date snapshot: simple copy/refresh from parent snapshot
        snap_local = conform_dir / "electrofit.toml"
        if parent_cfg_target.is_file():
            try:
                refresh = False
                if not snap_local.is_file():
                    refresh = True
                else:
                    try:
                        if parent_cfg_target.stat().st_mtime > snap_local.stat().st_mtime:
                            refresh = True
                        elif override_cfg:
                            # If an override was supplied this run, force refresh to reflect new layered state
                            refresh = True
                    except Exception:
                        refresh = True
                if refresh:
                    shutil.copy2(parent_cfg_target, snap_local)
                    logging.info(f"[snapshot] {'created' if not snap_local.exists() else 'refreshed'} {snap_local.relative_to(conform_dir)} from parent")
            except Exception as e:
                logging.debug(f"[snapshot] copy failure for {snap_local}: {e}")
        if bash_script_path.is_file():
            shutil.copy(str(bash_script_path), conform_dir)
        # Protocol-specific auxiliary files
        if protocol == "opt":
            if respin1_file and respin1_file.is_file():
                shutil.copy(str(respin1_file), conform_dir)
            if respin2_file and respin2_file.is_file():
                shutil.copy(str(respin2_file), conform_dir)
        # Symmetry groups JSON (for bcc + adjust_sym OR if already present), always copy so step5 can decide
        if not equiv_groups_file and (pis_dir / "equiv_groups.json").is_file():
            equiv_groups_file = pis_dir / "equiv_groups.json"
        if equiv_groups_file and equiv_groups_file.is_file():
            target_json = conform_dir / "equiv_groups.json"
            if not target_json.exists():
                try:
                    shutil.copy(str(equiv_groups_file), target_json)
                except Exception:
                    pass

    # Ensure extracted_conforms root snapshot remains consistent (already built once)

        conform_name = f"{molecule_name}c{i}.pdb"
        conform_path = conform_dir / conform_name
        c.save_pdb(str(conform_path))
    if input_mol2_file.is_file():
        mol2_to_pdb_with_bonds(input_file=str(input_mol2_file), existing_pdb_file=str(conform_path), verbose=verbose)
    logging.info(f"Wrote conformer {conform_name}")
    logging.info(f"Completed extraction: {len(configs)} conformers (method={method})")
    return True, f"Extracted {len(configs)} conformers (method={method}) to {extracted_conforms_dir}"


def _worker(args_tuple):
    mol_dir_str, project_root_str, sample, method, seed, override_cfg, multi_mol, verbose = args_tuple
    try:
        ok, msg = _extract_for_molecule(Path(mol_dir_str), Path(project_root_str), sample, method, seed, override_cfg, multi_mol, verbose)
        return (mol_dir_str, ok, msg, None)
    except Exception as e:
        return (mol_dir_str, False, "exception", str(e))


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Step4: Extract conformers from GROMACS trajectories (md_center.xtc + md.gro).\n"
            "Sampling methods: linear (evenly spaced), random (uniform without replacement), maxmin (farthest point RMSD).\n"
            "Configuration: sampling.count/method/seed in electrofit.toml (override per run via --config)."
        )
    )
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--sample", type=int, default=None, help="Number of conformers to sample (overrides config)")
    ap.add_argument("--sampling-method", default=None, help="Sampling method: linear|random|maxmin (overrides config)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for sampling where applicable")
    ap.add_argument("--workers", type=int, default=0, help="Parallel workers (0=auto, 1=sequential)")
    ap.add_argument("--no-progress", action="store_true", help="Disable progress bar output")
    ap.add_argument("--clean", action="store_true", help="Delete existing extracted_conforms directories and re-create fresh snapshots before extraction")
    ap.add_argument("--verbose", action="store_true", help="Emit per-conformer bond insertion info (quiet by default)")
    args = ap.parse_args()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    project_root = Path(args.project).resolve()
    # Load optional global config for sampling defaults
    sampling_cfg = {}
    candidate_cfgs: list[str] = []
    if getattr(args, "config", None):
        candidate_cfgs.append(args.config)
    candidate_cfgs.append(str(project_root / "electrofit.toml"))
    for c in candidate_cfgs:
        if not c:
            continue
        p = Path(c)
        if p.is_file():
            try:
                with p.open("rb") as f:
                    data = _toml.load(f)
                sampling_cfg = data.get("sampling", {}) or {}
            except Exception:  # pragma: no cover - robust fallback
                sampling_cfg = {}
            break

    method = args.sampling_method or sampling_cfg.get("method") or "linear"
    sample_count = args.sample or sampling_cfg.get("count") or 20
    seed = args.seed if args.seed is not None else sampling_cfg.get("seed")
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        print("No process directory found; nothing to do.")
        return
    mol_dirs = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    multi_mol = len(mol_dirs) > 1
    total = len(mol_dirs)
    if total == 0:
        print("[step4] No molecule directories found.")
        return

    # Optional destructive cleanup: remove existing extracted_conforms to guarantee fresh sampling
    if args.clean:
        for mdir in mol_dirs:
            ec_dir = mdir / "extracted_conforms"
            if ec_dir.is_dir():
                try:
                    shutil.rmtree(ec_dir)
                    print(f"[step4][clean] removed {ec_dir}")
                except Exception as e:  # pragma: no cover
                    print(f"[step4][clean][warn] failed to remove {ec_dir}: {e}")

    # Determine worker count
    if args.workers <= 0:
        # Leave 1 core free if possible
        cpu_count = max(1, multiprocessing.cpu_count() - 1)
        workers = min(cpu_count, total)
    else:
        workers = min(args.workers, total)

    extracted = 0
    results = []
    if workers == 1:
        iterator = mol_dirs
        if not args.no_progress:
            iterator = tqdm(iterator, desc="step4", unit="mol")
        for sub in iterator:
            ok, msg = _extract_for_molecule(sub, project_root, sample_count, method, seed, override_cfg, multi_mol, args.verbose)
            results.append((sub.name, ok, msg, None))
            if ok:
                extracted += 1
    else:
        tasks = [
            (str(p), str(project_root), sample_count, method, seed, override_cfg, multi_mol, args.verbose)
            for p in mol_dirs
        ]
        if not args.no_progress:
            pbar = tqdm(total=total, desc="step4", unit="mol")
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_map = {ex.submit(_worker, t): t[0] for t in tasks}
            for fut in as_completed(future_map):
                name, ok, msg, err = fut.result()
                results.append((Path(name).name, ok, msg, err))
                if ok:
                    extracted += 1
                if not args.no_progress:
                    pbar.update(1)
        if not args.no_progress:
            pbar.close()

    # Ordered output by molecule name
    for name, ok, msg, err in sorted(results, key=lambda x: x[0]):
        prefix = "[step4]" if ok else "[step4][skip]"
        detail = msg if err is None else f"{msg}: {err}"
        print(f"{prefix} {name}: {detail}")
    if extracted == 0:
        print("[step4] No conformers extracted.")
    else:
        print(f"[step4] Conformer extraction complete for {extracted}/{total} molecules.")


if __name__ == "__main__":
    main()
