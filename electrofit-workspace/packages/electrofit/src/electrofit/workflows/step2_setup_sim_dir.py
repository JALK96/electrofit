import fnmatch
import os
import shutil
import json

from electrofit.config.loader import load_config

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
CONFIG_PATH = os.environ.get("ELECTROFIT_CONFIG_PATH")  # optional
project_path = PROJECT_PATH

# Define the base process directory
process_dir = os.path.join(project_path, "process")

# Path to the MD directory
mdp_source_dir = os.path.join(project_path, "data/MDP")

cfg = load_config(project_path, CONFIG_PATH)

# Helper to write a manifest for step3
def _write_manifest(dest_dir: str, files: dict[str, str], mdp_subdir: str = "MDP") -> None:
    """Write a small run.json manifest so step3 can run deterministically."""
    # molecule: derive from .gro without extension, strip a trailing "_GMX" if present
    gro = files.get("gro")
    molecule = None
    if gro:
        base = os.path.splitext(os.path.basename(gro))[0]
        molecule = base[:-4] if base.endswith("_GMX") else base
    manifest = {
        "molecule": molecule,
        "gro": os.path.basename(files.get("gro", "")),
        "top": os.path.basename(files.get("top", "")),
        "itp": os.path.basename(files.get("itp", "")),
        "posres": os.path.basename(files.get("posres", "")),
        "mdp_dir": mdp_subdir,
    }
    with open(os.path.join(dest_dir, "run.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote manifest: {os.path.join(dest_dir, 'run.json')}")

# File patterns to search for
file_patterns = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]

for folder_name in os.listdir(process_dir):
    folder_path = os.path.join(process_dir, folder_name)
    if not os.path.isdir(folder_path):
        continue

    run_gau_dir = os.path.join(folder_path, "run_gau_create_gmx_in")
    if not os.path.isdir(run_gau_dir):
        print(f"'run_gau_create_gmx_in' does not exist in {folder_path}")
        continue

    # Destination: run_gmx_simulation
    dest_dir = os.path.join(folder_path, "run_gmx_simulation")
    os.makedirs(dest_dir, exist_ok=True)

    # copy .ef file(s) from run_gau_* (your existing logic)
    for file in os.listdir(run_gau_dir):
        if file.endswith(".ef"):
            src = os.path.join(run_gau_dir, file)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(dest_dir, os.path.basename(src)))
                print(f"Copied {file} to {dest_dir}")
            else:
                print(f"Configuration input file (.ef) does not exist: {src}")

    # locate the acpype directory
    acpype_folder_path = None
    for subfolder_name in os.listdir(run_gau_dir):
        if subfolder_name.endswith(".acpype"):
            acpype_folder_path = os.path.join(run_gau_dir, subfolder_name)
            break
    if not acpype_folder_path:
        print(f"No acpype directory found in {run_gau_dir}")
        continue

    # Copy matching files from acpype
    for pattern in file_patterns:
        for file_name in os.listdir(acpype_folder_path):
            if fnmatch.fnmatch(file_name, pattern):
                shutil.copy(
                    os.path.join(acpype_folder_path, file_name),
                    dest_dir,
                )
                print(f"Copied {file_name} to {dest_dir}")

    # Copy the MDP directory
    md_dest_dir = os.path.join(dest_dir, "MDP")
    if os.path.exists(mdp_source_dir):
        # handle re-runs gracefully
        if os.path.exists(md_dest_dir):
            shutil.rmtree(md_dest_dir)
        shutil.copytree(mdp_source_dir, md_dest_dir)
        print(f"Copied MDP directory to {md_dest_dir}")
    else:
        print(f"MDP source directory does not exist: {mdp_source_dir}")


    # Build and write manifest for step3
    selected = {"gro": None, "itp": None, "top": None, "posres": None}
    for name in os.listdir(dest_dir):
        if name.endswith(".gro") and name.endswith("GMX.gro"):
            selected["gro"] = os.path.join(dest_dir, name)
        elif name.endswith(".itp") and name.endswith("GMX.itp"):
            selected["itp"] = os.path.join(dest_dir, name)
        elif name.endswith(".top"):
            selected["top"] = os.path.join(dest_dir, name)
        elif name.startswith("posre_") and name.endswith(".itp"):
            selected["posres"] = os.path.join(dest_dir, name)
    # Fallback scan if patterns above didn't match
    if not selected["gro"]:
        for name in os.listdir(dest_dir):
            if name.endswith(".gro"):
                selected["gro"] = os.path.join(dest_dir, name); break
    if not selected["itp"]:
        for name in os.listdir(dest_dir):
            if name.endswith(".itp") and not name.startswith("posre_"):
                selected["itp"] = os.path.join(dest_dir, name); break
    if not selected["top"]:
        for name in os.listdir(dest_dir):
            if name.endswith(".top"):
                selected["top"] = os.path.join(dest_dir, name); break
    if not selected["posres"]:
        for name in os.listdir(dest_dir):
            if name.startswith("posre_") and name.endswith(".itp"):
                selected["posres"] = os.path.join(dest_dir, name); break

    _write_manifest(dest_dir, selected, mdp_subdir="MDP")

print("Done!")