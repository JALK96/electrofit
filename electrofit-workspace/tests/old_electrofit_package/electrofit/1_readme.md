
# Summary of Executed Operations Performed by Processing the Initial Structure

This document provides a comprehensive overview of the operations performed during the execution of the `process_initial_structure` script in the **Electrofit** package. The steps are generalized to avoid explicit file naming and include all functionalities, options, and conditions based on protocol and symmetry adjustments.

---

## Preliminary Steps

1. **Resolve Structure**
   - **Function:** `mol2_to_pdb_and_back()`
   - Converts the input MOL2 file to a PDB format and back to MOL2 using Open Babel to ensure correct file formatting and atom labeling.
   - Addresses potential issues with atom ordering and file inconsistencies.

2. **Find Project Root Directory**
   - **Function:** `find_project_root()`
   - Locates the root directory of the Electrofit project to ensure relative paths and imports work correctly.
   - Essential for accessing modules and scripts in different directories.

3. **Parse Configuration and Set Up Logging**
   - **Imports:**
     - `from electrofit.helper.config_parser import ConfigParser`
     - `from electrofit.helper.set_logging import setup_logging`
   - Parses the `input.ef` configuration file to retrieve parameters like `Protocol`, `AdjustSymmetry`, `IgnoreSymmetry`, etc.
   - Sets up logging for process tracking and debugging.

4. **Create Scratch Directory**
   - **Functions:** `setup_scratch_directory()`, `finalize_scratch_directory()`
   - Initializes a temporary working directory for intermediate file processing.
   - Prevents clutter in the main directory and isolates the workflow.

---

## Conditional Workflow Based on Protocol

The workflow diverges based on the `Protocol` specified in `input.ef`. The two supported protocols are:

- **`bcc`**: Uses AM1-BCC charges assigned by Antechamber.
- **`opt`**: Performs full Gaussian optimization followed by RESP fitting.

### Protocol: `bcc`

1. **Run ACPYPE to Generate GROMACS Input with AM1-BCC Charges**
   - **Function:** `run_acpype()`
   - **Command:**

     ```bash
     acpype -i input_file.mol2 -n <charge> -a <atom_type> -c bcc
     ```

   - Generates GROMACS-compatible topology and coordinate files using AM1-BCC charges.
   - Specifies the net charge and atom types (e.g., gaff, gaff2).

2. **Finalize Scratch Directory**
   - Copies output files back to the original directory.
   - Cleans up the scratch directory.

### Protocol: `opt`

1. **Run Antechamber to Generate Gaussian Input**
   - **Command:**

     ```bash
     antechamber -i input_file.mol2 -fi mol2 -nc <charge> -at <atom_type> -o gau_input.gcrt -fo gcrt -gv 1 -ge output.gesp
     ```

   - Prepares the molecule for Gaussian optimization without assigning charges.

2. **Run Gaussian Optimization**
   - **Function:** `run_gaussian_calculation()`
   - **Command:**

     ```bash
     rung16 gau_input.gcrt
     ```

   - Performs full geometry optimization using Gaussian.

3. **Execute ESPGEN**
   - **Function:** `run_espgen()`
   - **Command:**

     ```bash
     espgen -i input_file.gesp -o output_file.esp
     ```

   - Converts the Gaussian ESP output into a format suitable for RESP fitting.

4. **Generate RESP Input Files**
   - **Function:** `gaussian_out_to_prepi()`
   - **Command:**

     ```bash
     antechamber -i gau_input.gcrt.log -fi gout -o output_file.prepi -fo prepi -c resp -s 2
     ```

   - Processes Gaussian output to generate RESP input files (ANTECHAMBER_RESP1.IN, ANTECHAMBER_RESP2.IN).

5. **Adjust Symmetry (If Applicable)**
   - **Parameters:**
     - `AdjustSymmetry`
     - `IgnoreSymmetry`
   - **Scripts:** `edit_resp.py`, `write_symmetry.py`
   - Commands vary based on parameter values:
     - **AdjustSymmetry = True, IgnoreSymmetry = False:**

       ```bash
       python edit_resp.py ANTECHAMBER_RESP1.IN equiv_groups.json ANTECHAMBER_RESP1_MOD.IN
       ```

     - **AdjustSymmetry = True, IgnoreSymmetry = True:**

       ```bash
       python edit_resp.py ANTECHAMBER_RESP1.IN equiv_groups.json ANTECHAMBER_RESP1_MOD.IN --ignore_sym
       ```

     - **Verification:**  

       ```bash
       python write_symmetry.py ANTECHAMBER_RESP1_MOD.IN symmetry_resp_MOD.txt
       ```

6. **Run RESP Charge Fitting**
   - **Functions:** `run_resp()`
   - **Commands:**
     - **Stage 1:**

       ```bash
       resp -O -i ANTECHAMBER_RESP1.IN -o resp_output1.OUT -p resp_pch1.pch -t resp_chg1.chg -e output_file.esp
       ```

     - **Stage 2:**

       ```bash
       resp -O -i ANTECHAMBER_RESP2.IN -o resp_output2.OUT -p resp_pch2.pch -t resp_chg2.chg -e output_file.esp -q resp_chg1.chg
       ```

7. **Update MOL2 File with New RESP Charges**
   - **Function:** `update_mol2_charges()`
   - **Script:** `update_mol2.py`
   - **Command:**

     ```bash
     python update_mol2.py input_file.mol2 resp_chg2.chg output_file_resp.mol2
     ```

8. **Run ACPYPE for Topology and Parameter Generation**
   - **Function:** `run_acpype()`
   - **Command:**

     ```bash
     acpype -i output_file_resp.mol2 -n <charge> -a <atom_type> -c user
     ```

9. **Finalize Scratch Directory**
   - Copies output files back to the original directory.
   - Cleans up the scratch directory.

---

## Common Steps (Regardless of Protocol)

1. **Handle Screen Session (If Applicable)**
   - **Parameter:** `exit_screen`
   - Closes the detached screen session after processing if `exit_screen` is `True`.

2. **Error Handling and Logging**
   - Logs errors using Python’s logging module.
   - Ensures scratch directory is finalized to prevent data loss in case of errors.

---
