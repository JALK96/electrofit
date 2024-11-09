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

from electrofit.commands.run_commands import run_command, create_gaussian_input, run_gaussian_calculation, run_espgen, gaussian_out_to_prepi, run_resp, run_acpype, generate_mol2_with_resp_charges
from electrofit.helper.file_manipulation import find_file_with_extension, strip_extension, extract_charge_from_folder_name, mol2_to_pdb_and_back
from electrofit.helper.setup_finalize_scratch import (
    finalize_scratch_directory,
    setup_scratch_directory,
)

def process_initial_structure(molecule_name, mol2_file, base_scratch_dir, additional_input, net_charge, residue_name, adjust_sym=False, atom_type="gaff2", exit_screen=False): #delet later additional_input
    """
    Processes the initial optimized structure: Antechamber, Gaussian optimization, espgen, RESP fitting.

    Parameters:
    - molecule_name (str): Name of the molecule.
    - mol2_file (str): Path to the input mol2 file.
    - base_scratch_dir (str): Base directory for scratch space.
    - net_charge (int): Net charge of the molecule.
    """
    # Setup scratch directory without RESP input files
    input_files = [mol2_file]  # Only include existing input files

    input_files.extend(additional_input) # delet later additional_input
    
    # Step 0: Resolve Structure
    mol2_to_pdb_and_back(mol2_file, mol2_file, residue_name)
    
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)

    try:

        # Step 1: Generate Gaussian input file
        #create_gaussian_input(mol2_file, molecule_name, net_charge, scratch_dir, atom_type)

        # Step 2: Run Gaussian optimization
        gaussian_input = f"{molecule_name}.gcrt"
        #run_gaussian_calculation(gaussian_input, molecule_name, scratch_dir)

        # Step 3: Run espgen to generate .esp file
        gesp_file = f"{molecule_name}.gesp"
        esp_file = f"{molecule_name}.esp"
        run_espgen(gesp_file, esp_file, scratch_dir)

        # Step 4: Prepare RESP input files
        g_out = f"{gaussian_input}.log"
        gaussian_out_to_prepi(g_out, scratch_dir)

        # At this point, antechamber should have generated RESP input files in scratch_dir
        respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1.IN")
        respin2_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
        
        # Verify that RESP input files are created
        resp_files = [respin1_file, respin2_file]
        missing_resp_files = [file for file in resp_files if not os.path.isfile(file)]
        if missing_resp_files:
            logging.error(f"Missing RESP input files: {', '.join(missing_resp_files)}")
            raise FileNotFoundError(
                "RESP input files were not generated by antechamber."
            )
        else:
            logging.info("RESP input files generated successfully.")

        # Create json file containing symmetry information for phosphate groups
        json_symmetry_file = find_file_with_extension("json")
        
        # Write symmetry output to check
        write_symmetry = os.path.join(project_path, "electrofit/execution/write_symmetry.py")
        run_command(f"python {write_symmetry} {respin1_file} symmetry_resp.txt", cwd=scratch_dir)

        if adjust_sym:
            # Edit symmetry in "ANTECHAMBER_RESP1.IN"
            edit_resp_script = os.path.join(project_path, "electrofit/execution/edit_resp.py")
            run_command(f"python {edit_resp_script} ANTECHAMBER_RESP1.IN {json_symmetry_file} ANTECHAMBER_RESP1_MOD.IN", cwd=scratch_dir)

            # Re-define the respin1-file
            respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN")

            # Write symmetry output of alterd RESP1_MOD.IN to check/compare
            run_command(f"python {write_symmetry} {respin1_file} symmetry_resp_MOD.txt", cwd=scratch_dir)


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

        # Step 6: Generate mol2 file with RESP charges (comment: the above step 5 is no longer necesary with this step)
        mol2_resp = f"{molecule_name}_resp.mol2"

        update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
        run_command(f"python {update_mol2} {mol2_file} {resp_chg2} {mol2_resp}", cwd=scratch_dir)

        #generate_mol2_with_resp_charges(g_out, mol2_resp, scratch_dir, atom_type)

        # Step 7: Run acpype to generate GROMACS input (.itp/.gro/.top)
        run_acpype(mol2_resp, net_charge, scratch_dir, atom_type)

        # Finalize scratch directory
        finalize_scratch_directory(original_dir, scratch_dir, input_files)

        # Exit screen session (close detached session)
        if exit_screen:
            run_command("exit")


    except Exception as e:
        logging.error(f"Error processing initial structure: {e}")
        finalize_scratch_directory(original_dir, scratch_dir, input_files)
        exit(1)  # Exit the script with a non-zero status


