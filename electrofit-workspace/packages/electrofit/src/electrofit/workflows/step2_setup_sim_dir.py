import fnmatch
import os
import shutil

from electrofit.config import load_config, write_hpc_env_file

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
CONFIG_PATH = os.environ.get("ELECTROFIT_CONFIG_PATH")  # optional
project_path = PROJECT_PATH

# Define the base process directory
process_dir = os.path.join(project_path, "process")

# Path to the MD directory
mdp_source_dir = os.path.join(project_path, "data/MDP")

# Path to gmx.sh executable (in workspace scripts/)
bash_script_source = os.path.join(project_path, "scripts/gmx.sh")

cfg = load_config(project_path, CONFIG_PATH)

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

    # Copy gmx.sh
    if os.path.exists(bash_script_source):
        bash_dest_path = os.path.join(dest_dir, os.path.basename(bash_script_source))
        shutil.copy2(bash_script_source, bash_dest_path)
        print(f"Copied {os.path.basename(bash_script_source)} to {dest_dir}")
    else:
        print(f"Bash script does not exist: {bash_script_source}")

    # NEW: write the HPC env file for gmx.sh
    env_path = write_hpc_env_file(cfg, dest_dir)
    print(f"Wrote {env_path}")

print("Done!")