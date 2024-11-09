#!/usr/bin/env python
# coding: utf-8

import numpy as np
import mdtraj as md
import os
import sys
import shutil

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

from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.helper.config_parser import ConfigParser



process_dir = os.path.join(project_path, "old/process")

# Path to the bash script you want to copy
bash_script_path = os.path.join(project_path,"electrofit/bash/pc.sh")

# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):

    # Define paths for the run_gmx_simulation directory and extracted_conforms directory
    sim_dir = os.path.join(process_dir, sub_dir, "run_gmx_simulation")
    input_file_path = os.path.join(sim_dir, find_file_with_extension("ef"))

    config=ConfigParser(input_file_path)
    molecule_name = config.MoleculeName
    residue_name = config.ResidueName

    extracted_conforms_dir = os.path.join(process_dir, sub_dir, "extracted_conforms")

    # Check if required files are present
    traj_path = os.path.join(sim_dir, "md_center.xtc")
    gro_path = os.path.join(sim_dir, "md.gro")
    if not os.path.exists(traj_path) or not os.path.exists(gro_path):
        print(f"Skipping {sub_dir}: Required trajectory or structure file missing.")
        continue

    # Create extracted_conforms directory if it doesn't exist
    os.makedirs(extracted_conforms_dir, exist_ok=True)

    # Load the trajectory
    raw_traj = md.load(traj_path, top=gro_path)
    
    # Select atoms with the residue name derived from the binary string
    ipl = raw_traj.top.select(f"resname {residue_name}")  # Use "MOL" for legacy data
    traj = raw_traj.atom_slice(ipl)

    # Each frame is 100 ps, select frames at intervals of 1000 ps
    configs = [c for c, t in zip(traj, traj.time) if t % 1000 == 0]

    # Save each conformation in its own directory 
    for i, c in enumerate(configs):
        # Create a unique directory for each conformation
        conform_dir = os.path.join(extracted_conforms_dir, f"{molecule_name}c{i}")
        os.makedirs(conform_dir, exist_ok=True)

        # Copy the bash script into the conform_dir
        shutil.copy(bash_script_path, conform_dir)

        # Define the path for the PDB file
        conform_name = f"{molecule_name}c{i}.pdb"
        conform_path = os.path.join(conform_dir, conform_name)
        
        # Save the PDB file in its own directory
        c.save_pdb(conform_path)

    print(f"Extracted conformations and bash script saved for {sub_dir} in {extracted_conforms_dir}")