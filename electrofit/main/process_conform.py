import os
import logging
import sys

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

from electrofit.commands.run_commands import (
    run_command, 
    create_gaussian_input, 
    run_gaussian_calculation, 
    run_espgen, 
    run_resp
)
from electrofit.helper.file_manipulation import pdb_to_mol2, modify_gaussian_input
from electrofit.helper.setup_finalize_scratch import (
    finalize_scratch_directory,
    setup_scratch_directory,
)

def process_conform(molecule_name, pdb_file, base_scratch_dir, net_charge, residue_name, adjust_symm=False, exit_screen=False): 
    """
    Processes the conformers from the simulation output: Antechamber, Gaussian calculation, espgen, RESP fitting.

    Parameters:
    - molecule_name (str): Name of the molecule.
    - pdb_file (str): Path to the input pdb file.
    - base_scratch_dir (str): Base directory for scratch space.
    - net_charge (int): Net charge of the molecule.
    """
    # Define RESP input
    if adjust_symm:
        respin1_file = "ANTECHAMBER_RESP1_MOD.IN"
    else:
        respin1_file = "ANTECHAMBER_RESP1.IN"

    respin2_file = "ANTECHAMBER_RESP2.IN"

    # Setup scratch directory with RESP input files
    input_files = [pdb_file, respin1_file, respin2_file]  # Only include existing input files
    
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)

    try:

        # Step 0: Resolve Structure
        mol2_file = f"{molecule_name}.mol2"
        pdb_to_mol2(pdb_file, mol2_file, residue_name, cwd=scratch_dir)


        # Step 1: Generate Gaussian input file
        create_gaussian_input(mol2_file, molecule_name, net_charge, scratch_dir)

        # Modify gaussian input: replace "opt" with ""
        gaussian_input = f"{molecule_name}.gcrt"
        modify_gaussian_input(gaussian_input)

        # Step 2: Run Gaussian calculation
        run_gaussian_calculation(gaussian_input, molecule_name, scratch_dir)

        # Step 3: Run espgen to generate .esp file
        gesp_file = f"{molecule_name}.gesp"
        esp_file = f"{molecule_name}.esp"
        run_espgen(gesp_file, esp_file, scratch_dir)

        # Step 4: Prepare/Define RESP input files
        if adjust_symm:
            respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN")
        else:
            respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1.IN")

        respin2_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
        
        # Verify that RESP input files are present
        resp_files = [respin1_file, respin2_file]
        missing_resp_files = [file for file in resp_files if not os.path.isfile(file)]
        if missing_resp_files:
            logging.error(f"Missing RESP input files: {', '.join(missing_resp_files)}")
            raise FileNotFoundError(
                "RESP input files were not generated by antechamber."
            )
        else:
            logging.info("RESP input files generated successfully.")
        
        # Write symmetry output to check
        write_symmetry = os.path.join(project_path, "electrofit/execution/write_symmetry.py")
        
        if adjust_symm:
            # Write symmetry output of alterd RESP1_MOD.IN to check/compare
            run_command(f"python {write_symmetry} {respin1_file} symmetry_resp_MOD.txt", cwd=scratch_dir)

        else:
            # Write symmetry output of RESP1.IN to check/compare
            run_command(f"python {write_symmetry} {respin1_file} symmetry_resp.txt", cwd=scratch_dir)

        # Step 5: Run RESP fitting stages
        resp_output1 = f"{molecule_name}-resp1.out"
        resp_pch1 = f"{molecule_name}-resp1.pch"
        resp_chg1 = f"{molecule_name}-resp1.chg"
        run_resp(
            respin1_file, esp_file, resp_output1, resp_pch1, resp_chg1, scratch_dir
        )
        logging.info("RESP fitting stage 1 completed.")

        resp_output2 = f"{molecule_name}-resp2.out"
        resp_pch2 = f"{molecule_name}-resp2.pch"
        resp_chg2 = f"{molecule_name}-resp2.chg"
        run_resp(
            respin2_file, esp_file, resp_output2, resp_pch2, resp_chg2, scratch_dir, resp_chg1
        )
        logging.info("RESP fitting stage 2 completed.")

        # Step 6: Generate mol2 file with RESP charges
        mol2_resp = f"{molecule_name}_resp.mol2"

        update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
        run_command(f"python {update_mol2} {mol2_file} {resp_chg2} {mol2_resp}", cwd=scratch_dir)

        # Finalize scratch directory
        finalize_scratch_directory(original_dir, scratch_dir, input_files)

        # Exit screen session (close detached session)
        if exit_screen:
            run_command("exit")

    except Exception as e:
        logging.error(f"Error processing conform: {e}")
        finalize_scratch_directory(original_dir, scratch_dir, input_files)
        exit(1)  # Exit the script with a non-zero status

