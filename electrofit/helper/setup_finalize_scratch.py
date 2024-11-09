import os
import shutil
import logging
from electrofit.helper.set_logging import setup_logging


def setup_scratch_directory(input_files, base_scratch_dir="/scratch/johannal96/tmp"):
    """
    Sets up a scratch directory and copies essential input files.

    Parameters:
    - input_files (list): List of input file and directory names to copy into the scratch directory.
    - base_scratch_dir (str): Base directory for scratch space.

    Returns:
    - scratch_dir (str): Path to the created scratch directory.
    - original_dir (str): Path to the original working directory.
    """
    fullpath = os.getcwd()
    calcdir = os.path.basename(fullpath)
    parent_dir = os.path.dirname(fullpath)
    parent_folder_name = os.path.basename(parent_dir)
    scratch_dir = os.path.join(base_scratch_dir, parent_folder_name, calcdir)

    # Initialize logging with log file in original_dir
    log_file_path = os.path.join(fullpath, "process.log")
    setup_logging(log_file_path)
    # Print only if logging is not already initilized
    if not logging.getLogger().hasHandlers():
        logging.info(f"Logging initialized. Log file: {log_file_path}")

    # Create scratch directory
    os.makedirs(scratch_dir, exist_ok=True)
    logging.info(f"Created scratch directory: {scratch_dir}")

    # Copy input files
    for file in input_files:
        src = os.path.join(fullpath, file)
        dst = os.path.join(scratch_dir, file)
        if os.path.exists(src):
            if os.path.isfile(src):
                shutil.copy2(src, dst)
                logging.info(f"Copied file '{file}' to scratch directory.")
            elif os.path.isdir(src):
                shutil.copytree(src, dst)
                logging.info(f"Copied directory '{file}' to scratch directory.")
            else:
                logging.warning(
                    f"Input item '{file}' is neither a file nor a directory. Skipping copy."
                )
        else:
            logging.warning(
                f"Input file or directory '{file}' not found. Skipping copy."
            )

    # Log the contents of the scratch directory
    scratch_contents = os.listdir(scratch_dir)
    logging.info(f"Scratch directory contents after setup: {scratch_contents}")

    return scratch_dir, fullpath

def finalize_scratch_directory(
    original_dir, scratch_dir, input_files, output_files=None
):
    """
    Copies output files and directories back to the original directory and cleans up the scratch directory.

    Parameters:
    - original_dir (str): Path to the original directory.
    - scratch_dir (str): Path to the scratch directory.
    - input_files (list): List of input file and directory names to exclude from copying back.
    - output_files (list, optional): Specific list of output files/directories to copy back.
                                     If None, all items in scratch_dir excluding input_files are copied.
    """
    if output_files is None:
        # List all items in scratch_dir excluding input_files
        output_files = [f for f in os.listdir(scratch_dir) if f not in input_files]

    # Log the output files that will be copied back
    logging.info(f"Output files to be copied back: {output_files}")

    for item in output_files:
        src = os.path.join(scratch_dir, item)
        dst = os.path.join(original_dir, item)

        if not os.path.exists(src):
            logging.warning(
                f"Output item '{item}' does not exist in scratch directory."
            )
            continue

        # Determine if the item is a file or directory
        if os.path.isfile(src):
            # Handle file overwrites
            if os.path.exists(dst):
                base, ext = os.path.splitext(dst)
                counter = 1
                new_dst = f"{base}_copy{counter}{ext}"
                while os.path.exists(new_dst):
                    counter += 1
                    new_dst = f"{base}_copy{counter}{ext}"
                shutil.copy2(src, new_dst)
                logging.info(
                    f"File '{item}' already exists. Renamed to '{os.path.basename(new_dst)}'."
                )
            else:
                shutil.copy2(src, dst)
                logging.info(f"Copied file '{item}' back to original directory.")
        elif os.path.isdir(src):
            # Handle directory overwrites
            if os.path.exists(dst):
                base = dst
                counter = 1
                new_dst = f"{base}_copy{counter}"
                while os.path.exists(new_dst):
                    counter += 1
                    new_dst = f"{base}_copy{counter}"
                shutil.copytree(src, new_dst)
                logging.info(
                    f"Directory '{item}' already exists. Renamed to '{os.path.basename(new_dst)}'."
                )
            else:
                shutil.copytree(src, dst)
                logging.info(f"Copied directory '{item}' back to original directory.")
        else:
            logging.warning(
                f"Item '{item}' is neither a file nor a directory. Skipping."
            )

    # Clean up scratch directory
    try:
        shutil.rmtree(scratch_dir)
        logging.info(f"Removed scratch directory: {scratch_dir}")
    except Exception as e:
        logging.error(f"Failed to remove scratch directory '{scratch_dir}': {e}")