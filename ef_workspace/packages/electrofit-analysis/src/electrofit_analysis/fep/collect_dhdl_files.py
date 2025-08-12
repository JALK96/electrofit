#!/usr/bin/env python3

import os
import shutil
import sys


def main():
    """
    Collect dhdl.xvg files from lambda_* subdirectories and rename them
    to dhdl.X.xvg, placing them all in a single folder.

    Usage:
        python collect_dhdl_files.py <root_dir> [destination_dir]

    Example:
        python collect_dhdl_files.py /path/to/transitions all_dhdl
    """
    # --- Parse command-line arguments ---
    if len(sys.argv) < 2:
        print("Usage: python collect_dhdl_files.py <root_dir> [destination_dir]")
        sys.exit(1)

    root_dir = sys.argv[1]

    # You can optionally specify a custom destination folder:
    if len(sys.argv) > 2:
        destination_dir = sys.argv[2]
    else:
        destination_dir = "all_dhdl"

    # --- Create the destination folder if it doesn't exist ---
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # --- Iterate over entries in the root directory ---
    for entry in os.listdir(root_dir):
        entry_path = os.path.join(root_dir, entry)

        # Check if the entry is a directory and starts with 'lambda_'
        if os.path.isdir(entry_path) and entry.startswith("lambda_"):
            # Extract the index number from the directory name
            parts = entry.split("_")  # e.g. lambda_0 -> ["lambda", "0"]
            if len(parts) == 2 and parts[1].isdigit():
                lambda_index = parts[1]

                # Construct the path to dhdl.xvg inside this lambda folder
                dhdl_file = os.path.join(entry_path, "dhdl.xvg")

                # If the file exists, copy and rename
                if os.path.exists(dhdl_file):
                    new_name = f"dhdl.{lambda_index}.xvg"
                    destination_path = os.path.join(destination_dir, new_name)
                    shutil.copy2(dhdl_file, destination_path)
                    print(f"Copied {dhdl_file} -> {destination_path}")
                else:
                    print(f"Warning: No dhdl.xvg found in {entry_path}")

    print("Done collecting dhdl.xvg files.")


if __name__ == "__main__":
    main()
