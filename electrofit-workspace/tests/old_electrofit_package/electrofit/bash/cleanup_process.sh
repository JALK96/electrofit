#!/bin/bash

# ======================================================================
# Script Name: cleanup_process.sh
# Description: Deletes unwanted files and directories within the 'process/'
#              directory tree while preserving specified files and directories.
# Author: Arthur Laux
# Date: 06.09.2024
# ======================================================================

# Base directory
BASE_DIR="/home/johannal96/MA/electrofit/process"

# Verify that the base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory '$BASE_DIR' does not exist."
    exit 1
fi

# Iterate over each IP_XXXXXX directory
for ip_dir in "$BASE_DIR"/IP_*; do
    # Check if it's a directory
    if [ -d "$ip_dir" ]; then
        echo "Processing directory: $ip_dir"

        # Define the run_gau_create_gmx_in directory
        run_dir="$ip_dir/run_gau_create_gmx_in"

        if [ -d "$run_dir" ]; then
            echo "  Cleaning run_gau_create_gmx_in/..."

            # 1. Delete unwanted files in run_gau_create_gmx_in/
            echo "    Deleting unwanted files..."

            # Preserve only IPXX.mol2, delete IPXX_resp.mol2 and other unwanted files
            find "$run_dir" -type f ! \( \
                -name 'pis.sh' \
                -o -name '*.mol2' \
                -o -name '*.gcrt' \
                -o -name '*.gcrt.log' \
                -o -name '*.gesp' \
                -o -name '*.ef' \
                -o -name '*.json' \
                -o -name 'molecule.chk' \
            \) -print -delete

            # Explicitly delete IPXX_resp.mol2 if it exists
            find "$run_dir" -type f -name 'IP*_resp*.mol2' -print -delete

            echo "    Unwanted files deleted in run_gau_create_gmx_in/."

            # 2. Delete the IPXX_resp.acpype/ directory
            echo "    Deleting IPXX_resp.acpype/ directory if it exists..."
            resp_dir_pattern="$run_dir"/IP*_resp.acpype

            # Check if any directories match the pattern
            resp_dirs=( $resp_dir_pattern )
            if [ -d "${resp_dirs[0]}" ]; then
                for resp_dir in "${resp_dirs[@]}"; do
                    if [ -d "$resp_dir" ]; then
                        echo "      Deleting directory: $resp_dir"
                        rm -rf "$resp_dir"
                        echo "        Deleted $resp_dir"
                    fi
                done
            else
                echo "      No IP*_resp.acpype/ directories found."
            fi
        else
            echo "  Warning: Directory '$run_dir' does not exist."
        fi

        # Define the simulation directory
        simulation_dir="$ip_dir/run_gmx_simulation"

        if [ -d "$simulation_dir" ]; then
            echo "  Deleting run_gmx_simulation/ directory..."
            rm -rf "$simulation_dir"
            echo "    run_gmx_simulation/ directory deleted."
        else
            echo "  Warning: Directory '$simulation_dir' does not exist."
        fi

        echo "----------------------------------------"
    else
        echo "Warning: '$ip_dir' is not a directory."
    fi
done

echo "Cleanup process completed."