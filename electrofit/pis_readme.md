# Summary of Executed Operations Performed by Processing the Initial Structure

This document provides an abstracted overview of the operations performed during the execution process of "run_process_initial_structure.py". Each step is generalized to avoid explicit file naming.

1. **Create Scratch Directory**
   - *Initializes a temporary working directory for processing intermediate files, without interfering with regular update protocol.*

2. **Copy Input Files to Scratch Directory**
   - *Transfers necessary input files (i.e., molecular structure file: "input_file.mol2" and symmetry information file "equiv_groups.json") to the scratch directory for isolated processing.*

3. **Verify Scratch Directory Contents**
   - *Lists and logs all files and directories in the scratch directory to confirm successful copying of input files.*

4. **Run `antechamber` to generate gaussian input**
   - **Command:** `antechamber -i input_file.mol2 -fi mol2 -nc <charge> -at <atom_type> -o gau_input.gcrt -fo gcrt -gv 1 -ge output.gesp`
   - *Generates gau_input.gcrt and specifies in the input file that Gaussian should generate a .gesp file as output along with other output files.*

5. **Run `gaussian` Optimization**
   - **Command:** `rung16 gau_input.gcrt`
   - *Generates input_log_file.gcrt.log, molecule.chk, input_file.gesp, providing the necessary input to derive the electrostatic potential data for subsequent RESP fitting.*
  
6. **Execute `espgen`**
   - **Command:** `espgen -i input_file.gesp -o output_file.esp`
   - *Generates electrostatic potential data from Gaussian output files, essential for charge fitting.*

7. **Report `espgen` Output**
   - *Logs the creation of new electrostatic potential files and confirms the completion of `espgen` processing.*

8. **Execute `antechamber` to generate RESP input files for subsequent manual fitting** 
    - **Command:** `antechamber -i input_log_file.gcrt.log -fi gout -o output_file.prepi -fo prepi -c resp -s 2`
    - *Processes Gaussian output to prepare RESP charge assignments and preliminary parameter files. In this step the input files for symmetry aware RESP fitting are generated, i.e "resp_input1.IN" and "resp_input2.IN".*
    - *All other files generated in this step are not used for further processing.*
    - ----
    - **Substeps Executed by Antechamber:**

        8.1 **Assign Bond Types with `bondtype`**
        - **Command:** `bondtype -j full -i bond_type_input.ac0 -o bond_type_output.ac -f ac`
        - *Determines bond types based on input data to facilitate accurate charge fitting.*

        8.2 **Assign Atom Types with `atomtype`**
        - **Command:** `atomtype -i atomtype_input.ac0 -o atomtype_output.inf -p gaff`
        - *Assigns atom types using the GAFF force field, necessary for consistent parameterization.*

        8.3 **Generate RESP Input Files with `espgen`**
        - **Command:** `espgen -o antechamber_output_file.ESP -i input_log_file.gcrt.log`
        - *Creates RESP (Restrained ElectroStatic Potential) input files from Gaussian log data for charge fitting.*

        8.4 **Execute `respgen` for RESP Fitting**
        - **Command:** `respgen -i resp_input.ac -o resp_input1.IN -f resp1 -e 1`
        - **Command:** `respgen -i resp_input.ac -o resp_input2.IN -f resp2 -e 1`
        - *Generates RESP fitting input files for charge assignment.*

        8.5 **Run RESP Charge Fitting - Stage 1**
        - **Command:** `resp -O -i resp_input1.IN -o antechamber_resp_output1.OUT -e antechamber_output_file.ESP -t qout`
        - *Performs the first stage of RESP charge fitting to assign partial atomic charges.*

        8.6 **Run RESP Charge Fitting - Stage 2**
        - **Command:** `resp -O -i resp_input2.IN -o antechamber_resp_output2.OUT -e antechamber_output_file.ESP -q qout -t QOUT`
        - *Completes the second stage of RESP charge fitting for refined charge assignments.*

        8.7 **Generate PrePI Files with `prepgen`**
        - **Command:** `prepgen -i prep_input.ac -f int -o output_file.prepi -rn "MOL" -rf molecule.res`
        - *Creates PrePI (Pre-parameter Input) files required for parameter generation based on RESP fitting results.*
     - ---

9.  **Write Symmetry Information**
    - **Command:** `python write_symmetry.py resp_input1.IN symmetry_output.txt`
    - *Extracts and writes symmetry information from RESP input files and outputs them in a human-readable fashion to assist in accurate charge assignments, that considers symmetry correctly.*

10. **Edit RESP Input Files Based on Equivalence Groups**
    - **Command:** `python edit_resp.py resp_input1.IN equiv_groups.json resp_input1_MOD.IN`
    - *Modifies RESP input files using equivalence groups to ensure consistent and accurate charge assignment based on predefined symmetry.*
    - *The script enforce symmetry constraints based on user defined equivalence groups to improve the accuracy of the RESP fitting, thereby accounting for unidentified or falsely assigned symmetry through antechamber.*

11. **Write Modified Symmetry Information**
    - **Command:** `python write_symmetry.py resp_input1_MOD.IN symmetry_resp_MOD.txt`
    - *Updates symmetry information based on modified RESP input files.*

12. **Run RESP Charge Fitting with Modified Inputs - Stage 1**
    - **Command:** `resp -O -i resp_input1_MOD.IN -o resp_output1.OUT -p output_file.pch -t output_file1.chg -e output_file.esp`
    - *Performs RESP charge fitting using modified input files for initial charge assignment.*

13. **Run RESP Charge Fitting with Modified Inputs - Stage 2**
    - **Command:** `resp -O -i resp_input2.IN -o resp_output2.OUT -p output_file2.pch -t output_file2.chg -e output_file.esp -q output_file1.chg`
    - *Completes RESP charge fitting with modified inputs for refined charge assignments.*

14. **Update MOL2 File with New Charges**
    - **Command:** `python update_mol2.py input_file.mol2 charges_file.chg output_file_resp.mol2`
    - *Incorporates newly fitted RESP charges into the MOL2 molecular structure file.*

15. **Execute `acpype` for Topology and Parameter Generation**
    - **Command:** `acpype -i output_file_resp.mol2 -n <charge> -a <atom_type> -c user`
    - *Generates comprehensive topology and parameter files compatible with the <atom_type> force field, specifying net charge and user-defined charge assignments.*
    - *Where "-a <atom_type>" specifies the atom type (e.g., gaff, gaff2) and "-c user" tells acpype to use user-provided charges instead of calculating new ones.*

16. **Report `acpype` Output**
    - *Logs the creation of new topology and parameter files and confirms the completion of `acpype` processing.*

17. **Identify Output Files for Transfer**
    - *Lists all output files generated in the scratch directory for copying back to the original directory.*

18. **Copy Output Files Back to Original Directory**
    - *Transfers generated output files from the scratch directory back to the original working directory, ensuring that existing files are not overwritten by renaming originals if necessary.*

19. **Remove Scratch Directory**
    - *Cleans up by deleting the temporary scratch directory and all its contents after processing is complete.*

---

## **Notes**

- **Temporary Workspace Management**: Establishes and cleans up a dedicated environment to ensure that intermediate processing does not interfere with updates on main directories.
  
- **Input and Output Handling**: Ensures secure transfer of necessary files for processing and safely returns generated outputs without overwriting existing files.
  
- **Computational Chemistry Tools Execution**: Utilizes specialized tools (`espgen`, `antechamber`, `bondtype`, `atomtype`, `respgen`, `resp`, `prepgen`, `acpype`) for charge parameterization and fitting essential for accurate simulations.
  
- **Custom Scripting for Data Processing**: Employs user-defined Python scripts to handle specific tasks like symmetry information extraction, RESP input modification, and MOL2 file updates, ensuring tailored processing steps.
  
- **Logging and Verification**: Maintains detailed logs throughout the process for transparency, debugging, and verification of each step's successful execution.

---
