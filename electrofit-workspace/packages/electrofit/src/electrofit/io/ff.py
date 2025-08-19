import os
import logging





def replace_posres_in_file(file_path):
    """Opens a topology file and replaces"POSRES_LIG" with "POSRES"

    Args:
        file_path (file): Topology file to eddit.
    """
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if POSRES_LIG is present, then replace it
        if "POSRES_LIG" in content:
            updated_content = content.replace("POSRES_LIG", "POSRES")
        else:
            logging.error("No 'POSRES_LIG' found in the file.")
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(f"Successfully replaced 'POSRES_LIG' with 'POSRES' in {file_path}")

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def include_tip3p(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Avoid duplicate insertion
            include_line = f'#include "{forcefield}/tip3p.itp"'
            if include_line in content:
                logging.info(f"tip3p.itp already included in {file_path}")
                return

            # Insert the tip3p topology block before [ system ]
            updated_content = content.replace(
                "[ system ]",
                f"""
; Include water topology
#include "{forcefield}/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

[ system ]
""",
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include tip3p topology (tip3p.itp)."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(
            f"Successfully included tip3p.itp and position restraints in {file_path}"
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def include_ions(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Avoid duplicate insertion
            include_line = f'#include "{forcefield}/ions.itp"'
            if include_line in content:
                logging.info(f"ions.itp already included in {file_path}")
                return

            # Insert the ions topology block before [ system ]
            updated_content = content.replace(
                "[ system ]",
                f"""
; Include ion topology
#include "{forcefield}/ions.itp"

[ system ]
""",
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include ions.itp."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(f"Successfully included ions.itp in {file_path}")

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def remove_defaults_section_lines(file_path):
    """
    Removes the '[ defaults ]' section and the next three lines from a GROMACS topology file.

    Parameters:
    - file_path (str): Path to the GROMACS topology (.top) file.

    Returns:
    - None

    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - Exception: For any other unforeseen errors during file processing.
    """
    try:
        # Check if the file exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        # Read the original content of the file
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Initialize variables to track the removal process
        new_lines = []
        skip_next = 0  # Counter to skip lines after [ defaults ]

        for index, line in enumerate(lines):
            stripped_line = line.strip()

            # Detect the start of the '[ defaults ]' section
            if stripped_line == "[ defaults ]" and skip_next == 0:
                logging.info(
                    "Found '[ defaults ]' section. Initiating removal of this section and the next three lines."
                )
                skip_next = 3  # [ defaults ] + 2 subsequent lines
                continue  # Skip the '[ defaults ]' line

            if skip_next > 0:
                logging.debug(f"Skipping line {index + 1}: {line.strip()}")
                skip_next -= 1
                continue  # Skip the lines within the [ defaults ] section

            # Retain all other lines
            new_lines.append(line)

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.writelines(new_lines)

        logging.info(
            f"Successfully removed '[ defaults ]' section and the following three lines from '{file_path}'."
        )

    except FileNotFoundError as fnf_error:
        logging.error(fnf_error)
        raise
    except Exception as e:
        logging.error(f"An error occurred while removing '[ defaults ]' section: {e}")
        raise

def include_ff(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read all lines
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Define the lines to insert
        include_comment = "; Include forcefield\n"
        include_line = f'#include "{forcefield}/forcefield.itp"\n'

        # Check if the include line already exists to prevent duplication
        if any(forcefield in line and "forcefield.itp" in line for line in lines):
            logging.info(
                f'"{forcefield}/forcefield.itp" is already included in {file_path}.'
            )
            return

        # Insert the include lines as the second and third lines
        if lines:
            lines.insert(1, include_comment)
            lines.insert(2, include_line)
        else:
            # If the file is empty, add the include lines at the beginning
            lines.append(include_comment)
            lines.append(include_line)

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.writelines(lines)

        logging.info(
            f'Successfully included "{forcefield}/forcefield.itp" in {file_path}.'
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")