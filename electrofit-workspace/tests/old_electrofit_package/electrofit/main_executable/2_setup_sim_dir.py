import os
import sys
import shutil
import fnmatch

def find_project_root(current_dir, project_name="electrofit"):
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir  # Set root to the current project_name directory
        if parent_dir == current_dir:
            # We've reached the filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root  # Return the outermost match found
        current_dir = parent_dir
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)

# Define the base process directory
process_dir = os.path.join(project_path, "process")

# Path to the MD directory
mdp_source_dir = os.path.join(project_path, "data/MDP")

# Path to gmx.sh executable
bash_script_source = os.path.join(project_path, "electrofit/bash/gmx.sh") 


# File patterns to search for
file_patterns = [
    '*GMX.gro',
    '*GMX.itp',
    '*GMX.top',
    'posre_*.itp'
]

# Loop through each subdirectory in the process directory
for folder_name in os.listdir(process_dir):
    folder_path = os.path.join(process_dir, folder_name)
    
    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Define the 'run_gau_create_gmx_in' directory within this folder
        run_gau_dir = os.path.join(folder_path, 'run_gau_create_gmx_in')
        
        # Check if 'run_gau_create_gmx_in' exists
        if os.path.isdir(run_gau_dir):
            # Define the destination directory 'run_gmx_simulation' within the same working directory
            dest_dir = os.path.join(folder_path, 'run_gmx_simulation')
            os.makedirs(dest_dir, exist_ok=True)
            for file in os.listdir(run_gau_dir):
                if file.endswith('.ef'):
                    input_file_source = os.path.join(run_gau_dir, file)
                    # Proceed with processing


                    # Copy the bash script into dest_dir
                    if os.path.exists(input_file_source):
                        input_file_name = os.path.basename(input_file_source)
                        input_dest_path = os.path.join(dest_dir, input_file_name)
                        shutil.copy2(input_file_source, input_dest_path)
                        print(f"Copied {input_file_name} to {dest_dir}")
                    else:
                        print(f"Configuration input file (.ef) does not exist: {input_file_source}")

            # Find the acpype subdirectory within 'run_gau_create_gmx_in'
            for subfolder_name in os.listdir(run_gau_dir):
                if subfolder_name.endswith('.acpype'):
                    acpype_folder_path = os.path.join(run_gau_dir, subfolder_name)

                    # Copy files that match the patterns from the acpype folder
                    for pattern in file_patterns:
                        for file_name in os.listdir(acpype_folder_path):
                            if fnmatch.fnmatch(file_name, pattern):
                                source_file_path = os.path.join(acpype_folder_path, file_name)
                                # Copy the file to the destination directory
                                shutil.copy(source_file_path, dest_dir)
                                print(f"Copied {file_name} to {dest_dir}")

                    # Copy the MDP directory into dest_dir
                    md_dest_dir = os.path.join(dest_dir, 'MDP')
                    if os.path.exists(mdp_source_dir):
                        shutil.copytree(mdp_source_dir, md_dest_dir)
                        print(f"Copied MDP directory to {md_dest_dir}")
                    else:
                        print(f"MDP source directory does not exist: {mdp_source_dir}")

                    # Copy the bash script into dest_dir
                    if os.path.exists(bash_script_source):
                        bash_script_name = os.path.basename(bash_script_source)
                        bash_dest_path = os.path.join(dest_dir, bash_script_name)
                        shutil.copy2(bash_script_source, bash_dest_path)
                        print(f"Copied {bash_script_name} to {dest_dir}")
                    else:
                        print(f"Bash script does not exist: {bash_script_source}")

                    # Break after processing the acpype directory
                    break
            else:
                print(f"No acpype directory found in {run_gau_dir}")
        else:
            print(f"'run_gau_create_gmx_in' does not exist in {folder_path}")

print("Done!")