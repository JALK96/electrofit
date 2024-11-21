# Electrofit

Electrofit is a Python package that automates the partial charge fitting process for small molecules, streamlining parameterization for molecular simulations.
By performing multiple RESP fitting iterations across different conformational states, Electrofit generates a charge distribution from which average charges can be easily obtained.
Users can validate their charge fittings through integrated molecular dynamics simulations.

Electrofit simplifies the entire workflow, including:

- Quantum mechanical energy optimization of molecular structures
- Execution of molecular simulations
- Extraction and analysis of diverse conformations
- Calculation of atomic charges with symmetry handling and averaging
- Preparation of final simulations with optimized charges

## Table of Contents

- [Electrofit](#electrofit)
	- [Table of Contents](#table-of-contents)
	- [Features](#features)
	- [Installation](#installation)
	- [Usage](#usage)
		- [Workflow Overview](#workflow-overview)
		- [Configuration](#configuration)
		- [Step-by-Step Guide](#step-by-step-guide)
			- [0. Setup Process Directory](#0-setup-process-directory)
			- [1. Process Initial Structure](#1-process-initial-structure)
			- [2. Setup Simulation Directories](#2-setup-simulation-directories)
			- [3. Start Initial Simulation](#3-start-initial-simulation)
			- [4. Extract Conformations](#4-extract-conformations)
			- [5. Process Conformations](#5-process-conformations)
			- [6. Extract and Average Charges](#6-extract-and-average-charges)
			- [7. Setup Final Simulation Directory](#7-setup-final-simulation-directory)
			- [8. Start Final Simulation](#8-start-final-simulation)
	- [Module Descriptions](#module-descriptions)
		- [Main Executable Scripts](#main-executable-scripts)
		- [Helper Modules](#helper-modules)
			- [config\_parser.py](#config_parserpy)
			- [file\_manipulation.py](#file_manipulationpy)
			- [plotting.py](#plottingpy)
			- [set\_logging.py](#set_loggingpy)
			- [setup\_finalize\_scratch.py](#setup_finalize_scratchpy)
		- [Execution Scripts](#execution-scripts)
			- [edit\_resp.py](#edit_resppy)
			- [update\_mol2.py](#update_mol2py)
			- [plot\_results.py](#plot_resultspy)
		- [Bash Scripts](#bash-scripts)
			- [pis.sh](#pissh)
			- [gmx.sh](#gmxsh)
			- [pc.sh](#pcsh)
	- [Folder Structure](#folder-structure)
	- [Logging](#logging)
	- [License](#license)
	- [Contact](#contact)
	- [Examples](#examples)
		- [Scenario 1: Using AM1-BCC Charges Without Symmetry Adjustments](#scenario-1-using-am1-bcc-charges-without-symmetry-adjustments)
		- [Scenario 2: Full Optimization with Symmetry Adjustments](#scenario-2-full-optimization-with-symmetry-adjustments)
		- [Scenario 3: AM1-BCC Charges with Ignored Symmetry and Group Averaging](#scenario-3-am1-bcc-charges-with-ignored-symmetry-and-group-averaging)
	- [Dependencies](#dependencies)
	- [Acknowledgments](#acknowledgments)
	- [Additional Notes](#additional-notes)

---

## Features

- **Automated Partial Charge Fitting** : Streamlines the process of fitting partial charges using RESP, considering multiple conformations and symmetry constraints.

- **Simulation Workflow Automation** : Automates the setup and execution of molecular dynamics (MD) simulations with GROMACS.

- **Symmetry Handling** : Allows for user-defined symmetry adjustments during RESP fitting and charge averaging.

- **Flexible Protocols** : Supports different protocols (`bcc` and `opt`) for charge assignment and optimization.
  - bcc: Uses AM1-BCC charges assigned by Antechamber, suitable for quick parameterization.
  - opt: Performs full geometry optimization with Gaussian, followed by RESP fitting for the initial structure.

- **Visualization** : Generates plots to visualize charge distributions and analyze results, allowing for interpretation of charge fitting outcomes.

- **Logging** : Integrated logging system for monitoring processes and debugging, providing detailed information about each step of the workflow.

- **Scratch Directory Management** : Efficient handling of temporary files and directories during computations to optimize storage usage and avoid conflicts with regular backups, thereby preventing clogging the central backup
infrastructure.

---

## Installation

To install Electrofit, clone the repository and install the package using `pip`:

```bash
git clone https://github.com/JALK96/electrofit.git
cd electrofit
pip install .
```

---

## Usage

### Workflow Overview

The Electrofit workflow automates the process of parameterizing a molecule by fitting partial charges and testing the results through MD simulations. The workflow consists of the following main steps:

1. **Initial Setup and Structure Processing**

2. **Initial Simulation**

3. **Conformation Extraction and Processing**

4. **Charge Extraction and Averaging**

5. **Final Simulation Setup and Execution**

### Configuration

Before running the workflow, you need to set up the input directory (`data/input/your_molecule`), which must be located in the project root directory (`/electrofit`). 
This directory should contain:

- The molecular structure file of your molecule (`your_molecule.mol2`).
- A configuration file (`input.ef`) with parameters that control the behavior of the workflow.
- Optional: An equivalence groups JSON file (`equiv_groups.json`) if you are using symmetry adjustments.

**Example `input.ef` File:**

```ini
# Global parameters
Protocol: bcc
MoleculeName: IP_010111
ResidueName: I23
Charge: -8
AdjustSymmetry: TRUE
IgnoreSymmetry: TRUE
CalculateGroupAverage: TRUE
BaseScratchDir: /scratch/user_name/tmp

# Simulation parameters
AtomType: gaff2
BoxType: dodecahedron
BoxEdgeDistance: 1.2
Cation: NA
Anion: CL
IonConcentration: 0.15
```

**Parameter Descriptions:**  

- **Protocol** : Specifies the simulation protocol (`bcc` or `opt`).

- **MoleculeName** : Name of the molecule (must match the MOL2 structure file name).

- **ResidueName** : Short, three-character residue name.

- **Charge** : Total molecular charge.

- **AdjustSymmetry** : Whether to apply user-defined symmetry adjustments (`TRUE` or `FALSE`).

- **IgnoreSymmetry** : If `TRUE`, ignores symmetry during RESP fitting.

- **CalculateGroupAverage** : If `TRUE`, calculates group averages for symmetric atoms during charge averaging.

- **BaseScratchDir** : Base directory for scratch space.

- **AtomType** : Atom type used in simulations (`gaff`, `gaff2`, or `amber`).

- **BoxType** : Shape of the simulation box.

- **BoxEdgeDistance** : Distance between the molecule and the box edge (in nm).

- **Cation** : Name of the cation to add to the simulation.

- **Anion** : Name of the anion to add to the simulation.

- **IonConcentration** : Ion concentration for solvent replacement.

### Step-by-Step Guide

#### 0. Setup Process Directory

Prepare the initial processing directory by copying necessary input files.

```bash
python main_executable/0_setup_process_dir.py
```

- Copies input.ef, your_molecule.mol2, and optional symmetry files to the process/ directory.
- Ensures that the initial structure is ready for processing.


#### 1. Process Initial Structure

Run initial preprocessing of the molecular structure.

```bash
python main_executable/1_pis.py
```

This step involves:

- **Protocol `bcc`** : Assigns AM1-BCC charges using Antechamber and generates GROMACS input files using ACPYPE.

  - Assigns AM1-BCC charges using Antechamber.
  - Generates GROMACS input files using ACPYPE.

- **Protocol `opt`** : Performs Gaussian optimization and RESP fitting with optional symmetry adjustments.
  
  - Performs Gaussian optimization.
  - Generates electrostatic potential (ESP) files.
  - Performs RESP fitting with optional symmetry adjustments.


#### 2. Setup Simulation Directories

Prepare the simulation directories with necessary input files.

```bash
python main_executable/2_setup_sim_dir.py
```

- Organizes files for the initial MD simulation.
- Copies necessary topology and parameter files.

#### 3. Start Initial Simulation

Initiate the initial MD simulation using GROMACS.

```bash
python main_executable/3_start_sim.py
```

- Runs energy minimization, equilibration, and production MD.
- Uses the initial charges assigned in the previous steps.

#### 4. Extract Conformations

Extract conformations from the MD trajectory for further processing.

```bash
python main_executable/4_extract_conforms.py
```

- Extracts snapshots at specified intervals.
- Prepares individual directories for each conformation.

#### 5. Process Conformations

Process each extracted conformation to perform charge fitting.

```bash
python main_executable/5_process_conforms.py
```

- Performs Gaussian single point calculations on each conformation.
- Generates ESP files and performs RESP fitting.
- Applies symmetry adjustments if specified.

#### 6. Extract and Average Charges

Extract charges from processed conformations and calculate averages.

```bash
python main_executable/6_extract_charges.py
```

- Collects charges from all conformations.
- Averages charges, applying symmetry constraints if specified.

#### 7. Setup Final Simulation Directory

Set up the final simulation directories with updated charges.

```bash
python main_executable/7_setup_final_sim_dir.py
```

- Incorporates averaged charges into simulation input files.
- Prepares for the final MD simulation.


#### 8. Start Final Simulation

Run the final MD simulation using the updated charges.

```bash
python main_executable/8_start_final_sim.py
```

- Performs the final production MD simulation.
- Uses the refined charges for improved accuracy.

---

## Module Descriptions

### Main Executable Scripts

Located in the `main_executable/` directory, these scripts orchestrate the entire workflow.

- **0_setup_process_dir.py** : Sets up the processing directory.

- **1_pis.py** : Processes the initial molecular structure.

- **2_setup_sim_dir.py** : Sets up simulation directories.

- **3_start_sim.py** : Starts the initial MD simulation.

- **4_extract_conforms.py** : Extracts conformations from trajectories.

- **5_process_conforms.py** : Processes each conformation.

- **6_extract_charges.py** : Extracts and averages charges.

- **7_setup_final_sim_dir.py** : Sets up the final simulation directories.

- **8_start_final_sim.py** : Starts the final MD simulation.

### Helper Modules

Located in the `helper/` directory, these modules provide utility functions.

#### config_parser.py

- **Purpose** : Parses the `input.ef` configuration file.

- **Usage** :

```python
from electrofit.helper.config_parser import ConfigParser
config = ConfigParser('input.ef')
```

#### file_manipulation.py

- **Purpose** : Manages file operations and symmetry adjustments.

- **Key Functions** :
  - `update_mol2_charges()`: Updates charges in MOL2 files.

  - `edit_resp_input()`: Edits RESP input files to apply symmetry constraints.

  - `extract_charges_from_subdirectories()`: Collects charges from multiple conformations.

#### plotting.py

- **Purpose** : Generates plots for charge distributions.

- **Key Functions** :
  - `plot_charges_by_atom()`: Plots charge distribution colored by atom type.

  - `plot_charges_by_symmetry()`: Plots charge distribution colored by symmetry groups.

#### set_logging.py

- **Purpose** : Configures logging for scripts.

- **Usage** :

```python
from electrofit.helper.set_logging import setup_logging
setup_logging('process.log')
```

#### setup_finalize_scratch.py

- **Purpose** : Manages scratch directories.

- **Key Functions** :
  - `setup_scratch_directory()`: Sets up a scratch directory for temporary files and moves input files to scratch.

  - `finalize_scratch_directory()`: Cleans up and moves output files back to the original directory.

### Execution Scripts

Located in the `execution/` directory, these scripts perform specific tasks.

#### edit_resp.py

- **Purpose** : Edits RESP input files to incorporate user defined symmetry constraints, that needs to be specified in a json file (equiv_groups.json). 

- **Usage** :

```bash
python execution/edit_resp.py input_RESP_file.resp equivalence_groups.json output_RESP_file.resp --ignore_sym
```

#### update_mol2.py

- **Purpose** : Updates MOL2 files with new charges.

- **Usage** :

```bash
python execution/update_mol2.py input.mol2 charges.chg output.mol2
```

#### plot_results.py

- **Purpose** : Plots charge distributions.

- **Usage** :

```bash
python execution/plot_results.py -ic initial_charges.json -c charges_dict.json -s symmetry_groups.json -d output_directory
```

### Bash Scripts

Located in the `bash/` directory, these scripts facilitate execution on remote servers.
To make this work on your machine you must enable ssh-keygen to automatically switch to the remote server without the need of a password.

#### pis.sh

- **Purpose** : Runs initial structure processing on a remote server.

- **Usage** :

```bash
./bash/pis.sh
```

#### gmx.sh

- **Purpose** : Runs GROMACS simulations on a remote server.

- **Usage** :

```bash
./bash/gmx.sh
```

#### pc.sh

- **Purpose** : Processes conformations on a remote server.

- **Usage** :

```bash
./bash/pc.sh
```

---

## Folder Structure

```plaintext
electrofit/
├── bash/
│   ├── gmx.sh
│   ├── pc.sh
│   └── pis.sh
├── commands/
│   └── run_commands.py
├── execution/
│   ├── edit_resp.py
│   ├── plot_results.py
│   └── update_mol2.py
├── helper/
│   ├── config_parser.py
│   ├── file_manipulation.py
│   ├── plotting.py
│   ├── set_logging.py
│   └── setup_finalize_scratch.py
├── main/
│   ├── gmx_simulation.py
│   ├── process_conform.py
│   └── process_initial_structure.py
├── main_executable/
│   ├── 0_setup_process_dir.py
│   ├── 1_pis.py
│   ├── 2_setup_sim_dir.py
│   ├── 3_start_sim.py
│   ├── 4_extract_conforms.py
│   ├── 5_process_conforms.py
│   ├── 6_extract_charges.py
│   ├── 7_setup_final_sim_dir.py
│   └── 8_start_final_sim.py
├── data/input/your_molecule/
│   └── input.ef
│   └── your_molecule.mol2
│   └── equiv_groups.json
├── process/
│   └── (processing files)
├── README.md
├── setup.py
└── LICENSE
```

---

## Logging

Electrofit uses Python’s built-in `logging` module to record important events and errors. Logs are typically written to `process.log` in the working directory.
**Logging Setup Example:**

```python
from electrofit.helper.set_logging import setup_logging
setup_logging('process.log')
```

- **Features:**

  - Logs messages with different severity levels (DEBUG, INFO, WARNING, ERROR).
  - Outputs logs to both the console and log files.
  - Helps in tracking the progress and diagnosing issues.
---

## License

This project is licensed under the MIT License - see the [[LICENSE](#(https://mit-license.org))] file for details.

---

## Contact

For questions or suggestions, please open an issue on the GitHub repository or contact the maintainer at [johann.arthurlaux@gmail.com]() .

---

## Examples

### Scenario 1: Using AM1-BCC Charges Without Symmetry Adjustments

**Configuration:**

```ini
AdjustSymmetry: False
Protocol: bcc
```

**Workflow Highlights:**  

- **Charges Assigned** : AM1-BCC charges from Antechamber.

- **Symmetry Handling** : Default symmetry from Antechamber is used during all stages of the fitting process.

- **Charge Averaging** : No user defined symmetry constraints applied during averaging. Symmetry constraints enforced by Antechamber during RESP fitting result in equal final charges for by antechamber determined symmetric atoms. **CAUTION:** The symmetry assigned by Antechamber does not necessarily reflect the true symmetry of your molecule. 

### Scenario 2: Full Optimization with Symmetry Adjustments

**Configuration:**

```ini
AdjustSymmetry: True
Protocol: opt
```

**Workflow Highlights:**  

- **Charges Assigned** : RESP charges after full Gaussian optimization.

- **Symmetry Handling** : User-defined symmetry adjustments applied during initial and conformational RESP fitting, leading to equal charges for symmetric atoms.

- **Charge Averaging** : Symmetry constraints enforced during RESP fitting result in equal final charges for symmetric atoms.

### Scenario 3: AM1-BCC Charges with Ignored Symmetry and Group Averaging

**Configuration:**

```ini
AdjustSymmetry: True
Protocol: bcc
IgnoreSymmetry: True
CalculateGroupAverage: True
```

**Workflow Highlights:**  

- **Charges Assigned** : AM1-BCC charges; symmetry considered by antechamber initially.

- **Symmetry Handling** : All atoms treated independently during RESP fitting of single conformations.

- **Charge Averaging** : User-defined symmetry applied during charge averaging. The individually sampled charges for symmetric atoms are averaged and subsequently assigned to all symmetric atoms.

---

## Dependencies

- **Python 3.6 or higher**

- **Antechamber**  (part of AmberTools)

- **Gaussian** (Commercial software that requires a license. Ensure that Gaussian is installed and accessible in your environment.)

- **GROMACS** (Open-source MD simulation software)

- **ACPYPE** (to generate GROMACS input)

- **Open Babel**  (for file conversions)

- **Matplotlib**  (for plotting)

- **NumPy**

- **SciPy**

- **pandas**

**Note:** Install AmberTools23 via conda as a dedicated environment. Into this, install all other packages used.

---

## Acknowledgments

This package was developed to automate and streamline the process of molecular parameterization and simulation. Special thanks to all contributors and the developers of the underlying tools integrated into Electrofit.

- **AmberTools:** For Antechamber and RESP.
- **Gaussian:** For quantum chemical calculations.
- **GROMACS:** For molecular dynamics simulations.
- **ACPYPE:** For converting Amber parameters to GROMACS format.
- **Open Babel:** For file format conversions.

---

## Additional Notes

- **Remote Execution:** The bash scripts (`pis.sh`, `gmx.sh`, `pc.sh`) are designed to run on remote servers. Make sure to configure the scripts with the correct remote hostnames. You need to enable a passwordless connection to the remote server via ssh-keygen.
- **Scratch Directory:** The `BaseScratchDir` in input.ef should point to a directory with sufficient space and write permissions, as well as no backup protocols.
- **MDP Files:** Ensure that your simulation parameter files (.mdp files) are correctly set up in the `data/MDP/` directory.
- **Equivalence Groups:** When using symmetry adjustments, the equiv_groups.json file should be properly formatted and located in the input directory alongside the initial structure file (MOL2).

---

**Note:** Replace *your_molecule* with the actual name of your molecule in directory paths and filenames.
