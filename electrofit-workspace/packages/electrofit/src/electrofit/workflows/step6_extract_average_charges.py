import argparse
import fnmatch
import glob
import json
import os
import shutil
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from electrofit.io.mol2 import update_mol2_charges
from electrofit.cli.run_commands import run_acpype
from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.io.files import (
    adjust_atom_names,
    extract_charges_from_subdirectories,
    find_file_with_extension,
    load_symmetry_groups,
    parse_charges_from_mol2,
)
from electrofit.infra.logging import reset_logging, setup_logging
from electrofit.plotting.helpers import plot_charges_by_symmetry, plot_charges_by_atom

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH


def calculate_symmetric_group_averages(charges_dict_file, equivalent_groups_file):
    """
    Calculate the mean of the average charges for symmetric atoms.

    Parameters:
    - charges_dict_file: Path to the JSON file containing charges data for atoms.
        Example JSON:
        {
            "H1": {"average_charge": 0.12, "charges": [...]},
            "H3": {"average_charge": 0.11, "charges": [...]},
            "H6": {"average_charge": 0.13, "charges": [...]},
            ...
        }
    - equivalent_groups_file: Path to the JSON file containing symmetric atom groups.
        Example JSON:
        {
            "H1": ["H3", "H6"],
            "C1": ["C2"]
        }

    Returns:
    - updated_charges_dict: Dictionary with updated average charges for symmetric groups.
        Example:
        {
            "H1": {"average_charge": 0.12, "charges": [...]},
            "H3": {"average_charge": 0.12, "charges": [...]},
            "H6": {"average_charge": 0.12, "charges": [...]},
            ...
        }
    """
    # Load equivalent groups from the JSON file
    with open(equivalent_groups_file, "r") as f:
        equivalent_groups = json.load(f)

    # Load charges from the JSON file
    with open(charges_dict_file, "r") as f:
        charges_dict = json.load(f)

    updated_charges_dict = charges_dict.copy()  # Create a copy to modify

    # Iterate over each group in the equivalent groups
    for representative, group in equivalent_groups.items():
        # Add the representative atom to the group
        full_group = [representative] + group

        # Collect the average charges for all atoms in the group
        group_average_charges = [
            charges_dict[atom]["average_charge"]
            for atom in full_group
            if atom in charges_dict
        ]

        # Calculate the mean of the average charges
        if group_average_charges:
            group_average = sum(group_average_charges) / len(group_average_charges)

            # Update the average charge for all atoms in the group
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]["average_charge"] = group_average

    return updated_charges_dict


def plot_histograms(
    df, title, filename, adjusted_average_charges=None, color="darkred"
):
    """
    Plot histograms of the DataFrame and save the plot.

    Parameters:
        df (pd.DataFrame): DataFrame containing the charges data.
        title (str): Title of the plot.
        filename (str): Filename to save the plot.
        adjusted_average_charges (dict, optional): Dictionary of adjusted average charges.
        color (str): Color of the histogram bars.
    """
    axes = df.hist(bins=20, figsize=(15, 10), color=color, alpha=0.9, grid=False)

    for ax, col in zip(axes.flatten(), df.columns):
        ax.set_title("")
        ax.set_xlabel(col)
        ax.set_ylabel("")

        if adjusted_average_charges is not None:
            # Get the adjusted average charge from the dictionary
            mean_value = adjusted_average_charges.get(col, None)
        else:
            # Calculate the mean of the data
            mean_value = df[col].mean()

        if mean_value is not None:
            # Plot a vertical red dashed line at the mean
            ax.axvline(mean_value, color="black", linestyle="dashed", linewidth=1)

            # Add a red text label with the mean value
            ax.text(
                0.95,
                0.95,
                f"{mean_value:.2f}",
                color="black",
                fontsize=10,
                ha="right",
                va="top",
                transform=ax.transAxes,
            )

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename)
    plt.close()


def _process_one(mol_dir: Path, project_root: Path, override_cfg: Path | None, multi_mol: bool, remove_outlier: bool):
    """Process one molecule: aggregate charges across conformers and optionally symmetry-average."""
    extracted = mol_dir / "extracted_conforms"
    pis_dir = mol_dir / "run_gau_create_gmx_in"
    results_dir = mol_dir / "results"
    if not extracted.is_dir():
        return False, "no extracted_conforms"
    acpype_glob = list(pis_dir.glob("*.acpype"))
    if not acpype_glob:
        return False, "no acpype dir"
    ac_dir = acpype_glob[0]
    results_dir.mkdir(exist_ok=True)
    # Set up per-molecule logging (align with step4/5 style, no console spam)
    log_path = results_dir / "process.log"
    try:
        setup_logging(str(log_path), also_console=False)
    except Exception:
        pass  # fall back silently
    # Ensure snapshot in results dir.
    # Precedence mirrors step4 logic (context layering):
    #   1. extracted_conforms/electrofit.toml
    #   2. process/<mol>/electrofit.toml
    #   3. data/input/<mol>/electrofit.toml
    #   4. project_root/electrofit.toml
    snapshot = compose_snapshot(
        results_dir,
        project_root,
        mol_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        upstream=extracted / "electrofit.toml",
        process_cfg=mol_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    # Load config with context=results so last-level wins
    cfg = load_config(project_root, context_dir=results_dir, molecule_name=mol_dir.name)
    # Dump layered config into log (avoid direct stdout noise)
    try:
        dump_config(cfg, log_fn=logging.info)
        # Force flush of all file handlers so user immediately sees config in process.log
        for _h in logging.getLogger().handlers:
            try:
                _h.flush()
            except Exception:
                pass
    except Exception as e:  # pragma: no cover
        logging.debug(f"[step6] config dump failed: {e}")
    proj = cfg.project
    molecule_name = proj.molecule_name or mol_dir.name
    charge = proj.charge or 0
    atom_type = proj.atom_type or "gaff2"
    adjust_sym = getattr(proj, "adjust_symmetry", False)
    ignore_sym = getattr(proj, "ignore_symmetry", False)
    calc_group_average = getattr(proj, "calculate_group_average", False)

    # Find GAFF mol2 inside acpype
    mol2_source_file_path = None
    pattern = f"*{atom_type}.mol2"
    for fn in os.listdir(ac_dir):
        if fnmatch.fnmatch(fn, pattern):
            mol2_source_file_path = os.path.join(ac_dir, fn)
            break
    if not mol2_source_file_path:
        return False, "no matching mol2"

    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))
    atoms_dict = extract_charges_from_subdirectories(str(extracted), str(results_dir))
    with (results_dir / "charges_dict.json").open("w") as f:
        json.dump(atoms_dict, f, indent=2)
    with (results_dir / "initial_charges_dict.json").open("w") as f:
        json.dump(initial_charges_dict, f, indent=2)

    # Always produce per-atom distribution plot & average_charges.chg (legacy behaviour)
    try:
        plot_charges_by_atom(atoms_dict, initial_charges_dict, str(results_dir))
    except Exception as e:  # pragma: no cover - plotting robustness
        print(f"[step6][warn] plotting per-atom charges failed: {e}")

    if adjust_sym and not ignore_sym:
        # symmetry JSON expected already copied earlier; fallback search in pis_dir
        sym_json = extracted / "equiv_groups.json"
        if not sym_json.is_file():
            cand = list(pis_dir.glob("*.json"))
            if cand:
                sym_json = cand[0]
        if sym_json.is_file():
            equiv_group = load_symmetry_groups(str(sym_json))
            plot_charges_by_symmetry(atoms_dict, initial_charges_dict, str(results_dir), equiv_group)
            if calc_group_average and not remove_outlier:
                updated_charges_dict = calculate_symmetric_group_averages(results_dir / "charges_dict.json", sym_json)
                with (results_dir / "group_average_charges_dict.json").open("w") as f:
                    json.dump(updated_charges_dict, f, indent=2)
                # Write plain text & .chg
                group_txt = results_dir / "group_average_charges.txt"
                with group_txt.open("w") as f:
                    f.write("#Atom_Name\tAverage_Charge\n")
                    for atom, rec in updated_charges_dict.items():
                        f.write(f"{atom}\t{rec['average_charge']:.4f}\n")
                lines = [f"{rec['average_charge']:.4f}" for rec in updated_charges_dict.values()]
                (results_dir / "group_average_charges.chg").write_text("\n".join(lines) + "\n")
                # Prepare updated mol2 & run acpype with user charges
                updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                chg_file = results_dir / "group_average_charges.chg"
                logging.info("[step6] Updating MOL2 with group average charges and running acpype (user charges mode)...")
                update_mol2_charges(mol2_source_file_path, str(chg_file), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), atom_type, charges="user")
            else:
                # No group averaging requested -> still proceed with user charges update using average_charges.chg
                avg_chg = results_dir / "average_charges.chg"
                if avg_chg.is_file():
                    updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                    logging.info("[step6] Updating MOL2 with average charges (symmetry, no group average) and running acpype...")
                    update_mol2_charges(mol2_source_file_path, str(avg_chg), str(updated_mol2_out))
                    run_acpype(str(updated_mol2_out), charge, str(results_dir), atom_type, charges="user")
                else:
                    print(f"[step6][warn] (symmetry,no-group-average) missing average_charges.chg in {results_dir}")
        else:
            print(f"[step6][warn] no symmetry JSON for {molecule_name}; skipping symmetry averaging")
    else:
        # If we are NOT symmetry averaging or user opted out, still allow user-charges path using average_charges.chg
        # average_charges.chg already written by plot_charges_by_atom
        if not calc_group_average:  # mimic old condition (calc_group_average == False)
            avg_chg = results_dir / "average_charges.chg"
            if avg_chg.is_file():
                updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                logging.info("[step6] Updating MOL2 with average charges (no symmetry averaging) and running acpype...")
                update_mol2_charges(mol2_source_file_path, str(avg_chg), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), atom_type, charges="user")
            else:
                print(f"[step6][warn] expected average_charges.chg missing in {results_dir}")
    return True, "ok"


def main():  # pragma: no cover (CLI wrapper)
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--remove-outlier", action="store_true", help="Enable experimental outlier removal pipeline (legacy)")
    args = ap.parse_args()
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        print("[step6] No process directory.")
        return
    mol_dirs = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    multi_mol = len(mol_dirs) > 1
    done = 0
    for m in mol_dirs:
        ok, msg = _process_one(m, project_root, override_cfg, multi_mol, args.remove_outlier)
        status = "[step6]" if ok else "[step6][skip]"
        print(f"{status} {m.name}: {msg}")
        if ok:
            done += 1
    print(f"[step6] Completed {done}/{len(mol_dirs)} molecules.")


if __name__ == "__main__":  # pragma: no cover
    main()
