# Summary of Executed Operations Performed to Setup the Simulation Directory

This document provides an abstracted overview of the operations performed during the execution process of the `setup_sim_dir` script located in the main_executable package. This step is conducted after processing the initial structure and before conducting the simulation. Each step is generalized to avoid explicit file naming.

---

1. **Set Up Environment and Define Paths**
   - *Determines the script's directory and identifies the project root directory named `electrofit`.*
   - *Appends the project root directory to the system path for module imports.*
   - *Defines paths to key directories and files:*
     - `process_dir`: The base directory containing processed data.
     - `mdp_source_dir`: The source directory for MDP (Molecular Dynamics Parameter) files.
     - `bash_script_source`: The source path to the GROMACS bash script (`gmx.sh`).

2. **Specify File Patterns to Search For**
   - *Defines patterns of files to locate and copy for simulation setup:*
     - `*GMX.gro`
     - `*GMX.itp`
     - `*GMX.top`
     - `posre_*_resp.itp`

3. **Iterate Over Processed Directories**
   - *Loops through each subdirectory within the `process_dir` to locate necessary files and directories for simulation setup.*

4. **Check for `run_gau_create_gmx_in` Directory**
   - *For each subdirectory, checks if a `run_gau_create_gmx_in` directory exists.*
   - *If it exists:*
     - **Create Destination Directory**
       - *Defines and creates a `run_gmx_simulation` directory within the current working directory for simulation files.*
     - **Copy Configuration Files**
       - *Locates files ending with `.ef` within `run_gau_create_gmx_in` and copies them to `run_gmx_simulation`.*
       - **Command:**

         ```python
         shutil.copy2(run_gau_create_gmx_in/*.ef, run_gmx_simulation/) 
         ```

5. **Locate and Process ACPYPE Output**
   - *Within `run_gau_create_gmx_in`, searches for a subdirectory ending with `_resp.acpype`.*
   - *If found:*
     - **Copy Relevant Files from ACPYPE Output**
       - *Copies files matching the specified patterns from the ACPYPE directory to `run_gmx_simulation`.*
       - **Commands:**

         ```python
         shutil.copy2(acpype_output_dir/*GMX.gro, run_gmx_simulation/)
         shutil.copy2(acpype_output_dir/*GMX.itp, run_gmx_simulation/)
         shutil.copy2(acpype_output_dir/*GMX.top, run_gmx_simulation/)
         shutil.copy2(acpype_output_dir/posre_*_resp.itp, run_gmx_simulation/)
         ```

     - **Copy MDP Directory**
       - *Copies the MDP files directory from its source directory `mdp_source_dir` to `run_gmx_simulation/MDP` for simulation parameters.*
       - **Command:**

         ```python
         shutil.copy2(mdp_source_dir/MDP, run_gmx_simulation/MDP)
         ```

     - **Copy GROMACS Bash Script**
       - *Copies the `gmx.sh` bash script to `run_gmx_simulation` for executing GROMACS commands.*
       - **Command:**

         ```python
         shutil.copy2(electrofit/bash/gmx.sh, run_gmx_simulation/)
         ```

   - *If no ACPYPE directory is found:*
     - *Logs a message indicating the absence of the ACPYPE directory.*

6. **Handle Directories Without Necessary Files**
   - *If `run_gau_create_gmx_in` does not exist within a subdirectory, logs a message indicating its absence.*

7. **Finalize the Setup**
   - *After processing all subdirectories, confirms the completion of the setup process.*

---
