"""Pipeline Step 4: Extract conformers from production GROMACS trajectories.

Core responsibilities (orchestration only):
  * Parse CLI arguments (sampling overrides, workers, clean mode).
  * Derive sampling defaults from layered project / override TOMLs.
  * Iterate molecule process directories and invoke per‑molecule extraction.
  * Manage optional parallel execution and progress reporting.
  * Emit per‑molecule logging into that molecule's ``extracted_conforms/process.log``.
  * Summarise results to stdout and project-level ``step.log``.

Heavy domain logic (frame selection, snapshot layering, file copies) stays local
for now but may later move into a dedicated domain service module to minimise
side‑effects during import. The legacy implementation lived in
``electrofit.workflows.step4_extract_conforms`` and is now deprecated; this
module supersedes it with a thinner orchestration focus.
"""
from __future__ import annotations
import argparse, os, shutil, multiprocessing, logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

import mdtraj as md  # noqa: F401
from tqdm import tqdm

try:  # Python 3.11+
    import tomllib as _toml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml  # type: ignore

from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_sampling_decision
from electrofit.domain.sampling import select_frame_indices, prepare_conformer_directory

__all__ = ["main"]


def _extract_for_molecule(
    mol_proc_dir: Path,
    project_root: Path,
    sample: int,
    method: str,
    seed: int | None,
    override_cfg: Path | None,
    multi_mol: bool,
    verbose: bool,
) -> Tuple[bool, str]:
    sim_dir = mol_proc_dir / "run_gmx_simulation"
    pis_dir = mol_proc_dir / "run_gau_create_gmx_in"
    if not sim_dir.is_dir():
        return False, "no sim dir"
    cfg = load_config(project_root, context_dir=sim_dir, molecule_name=mol_proc_dir.name)
    proj = cfg.project
    molecule_name = proj.molecule_name or mol_proc_dir.name
    residue_name = proj.residue_name or "LIG"
    adjust_sym = getattr(proj, "adjust_symmetry", False)
    protocol = getattr(proj, "protocol", "bcc")

    respin1_file = respin2_file = equiv_groups_file = None
    if protocol == "opt":
        respin1_file = pis_dir / ("ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN")
        respin2_file = pis_dir / "ANTECHAMBER_RESP2.IN"
    elif protocol == "bcc" and adjust_sym:
        json_candidates = sorted(pis_dir.glob("*.json"))
        if json_candidates:
            equiv_groups_file = json_candidates[0]

    input_mol2_file = project_root / "data" / "input" / mol_proc_dir.name / f"{molecule_name}.mol2"
    if not input_mol2_file.is_file():
        logging.warning(f"[step4][{mol_proc_dir.name}] missing input mol2 ({input_mol2_file.name}); continuing (bond insertion skipped)")
    extracted_conforms_dir = mol_proc_dir / "extracted_conforms"
    extracted_conforms_dir.mkdir(exist_ok=True)
    reset_logging()
    setup_logging(str(extracted_conforms_dir / "process.log"), also_console=False)
    existing_snapshot = extracted_conforms_dir / "electrofit.toml"
    if existing_snapshot.is_file():
        try:
            cfg_existing = load_config(project_root, context_dir=extracted_conforms_dir, molecule_name=mol_proc_dir.name)
            logging.info(f"[config] existing extracted_conforms snapshot for {mol_proc_dir.name} detected -> dump below")
            dump_config(cfg_existing, log_fn=logging.info)
        except Exception:
            logging.debug("[step4] existing snapshot dump failed", exc_info=True)

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

    logging.info(f"Step4: start extraction for molecule_dir={mol_proc_dir.name} method={method} sample={sample}")
    try:
        dump_config(cfg, log_fn=logging.info)
    except Exception:
        logging.debug("[step4] dump per-molecule config failed", exc_info=True)

    raw_traj = md.load(str(traj_path), top=str(gro_path))
    try:
        res_counts: dict[str, int] = {}
        for res in raw_traj.topology.residues:
            res_counts[res.name] = res_counts.get(res.name, 0) + len(res.atoms)
        inv_str = ", ".join(f"{k}:{v}" for k, v in sorted(res_counts.items()))
        logging.info(f"[step4][{mol_proc_dir.name}] residue inventory -> {inv_str}")
    except Exception:
        logging.debug("[step4] residue inventory logging failed", exc_info=True)
    ipl = raw_traj.top.select(f"resname {residue_name}")
    if len(ipl) == 0:
        logging.error(f"[step4][{mol_proc_dir.name}] residue '{residue_name}' not found in trajectory topology; abort extraction.")
        return False, f"residue '{residue_name}' not in topology"
    traj = raw_traj.atom_slice(ipl)
    if traj.n_atoms == 0:
        logging.warning(f"[step4][{mol_proc_dir.name}] zero atoms after selection; skipping")
        return False, "no atoms after selection"

    total = len(traj)
    n = min(sample, total)
    indices = select_frame_indices(traj, n, method, seed)
    configs = [traj[i] for i in indices]
    logging.info(f"Selected indices (n={len(indices)}) -> {indices}")

    for i, c in enumerate(configs):
        prepare_conformer_directory(
            conform_index=i,
            molecule_name=molecule_name,
            parent_cfg_target=parent_cfg_target,
            override_cfg=override_cfg,
            protocol=protocol,
            respin1_file=respin1_file,
            respin2_file=respin2_file,
            equiv_groups_file=equiv_groups_file,
            pis_dir=pis_dir,
            extracted_conforms_dir=extracted_conforms_dir,
            input_mol2_file=input_mol2_file,
            traj_frame_save_fn=c.save_pdb,
            verbose=verbose,
        )

    logging.info(f"Completed extraction: {len(configs)} conformers (method={method}) for {mol_proc_dir.name}")
    return True, f"Extracted {len(configs)} conformers (method={method}) to {extracted_conforms_dir}"


def _worker(args_tuple):  # pragma: no cover
    mol_dir_str, project_root_str, sample, method, seed, override_cfg, multi_mol, verbose = args_tuple
    try:
        ok, msg = _extract_for_molecule(Path(mol_dir_str), Path(project_root_str), sample, method, seed, override_cfg, multi_mol, verbose)
        return (mol_dir_str, ok, msg, None)
    except Exception as e:  # pragma: no cover
        return (mol_dir_str, False, "exception", str(e))


def main():  # pragma: no cover
    ap = argparse.ArgumentParser(
        description=(
            "Step4: Extract conformers from GROMACS trajectories (md_center.xtc + md.gro).\n"
            "Sampling methods: linear, random, maxmin. Defaults read from layered electrofit.toml unless overridden."
        )
    )
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--sample", type=int, default=None, help="Number of conformers to sample (override)")
    ap.add_argument("--sampling-method", default=None, help="Sampling method: linear|random|maxmin (override)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for sampling (override)")
    ap.add_argument("--workers", type=int, default=0, help="Parallel workers (0=auto, 1=sequential)")
    ap.add_argument("--no-progress", action="store_true", help="Disable progress bar output")
    ap.add_argument("--clean", action="store_true", help="Remove existing extracted_conforms dirs before extraction")
    ap.add_argument("--verbose", action="store_true", help="Verbose bond insertion logs")
    args = ap.parse_args()

    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    project_root = Path(args.project).resolve()
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        print("No process directory found; nothing to do.")
        return
    mol_dirs = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    multi_mol = len(mol_dirs) > 1
    if not mol_dirs:
        print("[step4] No molecule directories found.")
        return

    sampling_cfg = {}
    candidate_cfgs: list[str] = []
    if getattr(args, "config", None):
        candidate_cfgs.append(args.config)
    candidate_cfgs.append(str(project_root / "electrofit.toml"))
    for c in candidate_cfgs:
        p = Path(c)
        if p.is_file():
            try:
                with p.open("rb") as fh:
                    data = _toml.load(fh)
                sampling_cfg = data.get("sampling", {}) or {}
            except Exception:
                sampling_cfg = {}
            break
    method = args.sampling_method or sampling_cfg.get("method") or "linear"
    sample_count = args.sample or sampling_cfg.get("count") or 20
    seed = args.seed if args.seed is not None else sampling_cfg.get("seed")

    if args.clean:
        for mdir in mol_dirs:
            ec_dir = mdir / "extracted_conforms"
            if ec_dir.is_dir():
                try:
                    shutil.rmtree(ec_dir)
                    print(f"[step4][clean] removed {ec_dir}")
                except Exception as e:
                    print(f"[step4][clean][warn] failed to remove {ec_dir}: {e}")

    if args.workers <= 0:
        cpu_count = max(1, multiprocessing.cpu_count() - 1)
        workers = min(cpu_count, len(mol_dirs))
    else:
        workers = min(args.workers, len(mol_dirs))

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
        pbar = None
        if not args.no_progress:
            pbar = tqdm(total=len(mol_dirs), desc="step4", unit="mol")
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_map = {ex.submit(_worker, t): t[0] for t in tasks}
            for fut in as_completed(future_map):
                name, ok, msg, err = fut.result()
                results.append((Path(name).name, ok, msg, err))
                if ok:
                    extracted += 1
                if pbar:
                    pbar.update(1)
        if pbar:
            pbar.close()

    for name, ok, msg, err in sorted(results, key=lambda x: x[0]):
        prefix = "[step4]" if ok else "[step4][skip]"
        detail = msg if err is None else f"{msg}: {err}"
        print(f"{prefix} {name}: {detail}")
    summary = f"[step4] Conformer extraction complete for {extracted}/{len(mol_dirs)} molecules." if extracted else "[step4] No conformers extracted."
    print(summary)

    reset_logging()
    setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
    try:
        log_run_header("step4")
    except Exception:
        logging.info("electrofit unknown | step=step4")
    logging.info(summary)


if __name__ == "__main__":  # pragma: no cover
    main()
