#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script Directory: $SCRIPT_DIR"

# Extract the path up to and including 'electrofit' from SCRIPT_DIR
PROJECT_ROOT="${SCRIPT_DIR%%electrofit*}electrofit"
echo "Project Root: $PROJECT_ROOT"

# Name of the Python script to execute
PYTHON_SCRIPT="$PROJECT_ROOT/electrofit/execution/run_gmx.py"  # Adjusted to avoid duplication
echo "Python Script: $PYTHON_SCRIPT"

# Remote machine details
REMOTE_MACHINE="qcm04"
echo "Remote Machine: $REMOTE_MACHINE"

# Generate a unique Screen session name to prevent conflicts
SCREEN_SESSION="$(basename "$(dirname "$SCRIPT_DIR")")_$(date +%Y%m%d%H%M%S)"
echo "Screen Session Name: $SCREEN_SESSION"

# Execute the Python script within a detached Screen session on the remote machine
ssh "$REMOTE_MACHINE" "cd \"$SCRIPT_DIR\" && source ~/.bashrc && screen -dmS \"$SCREEN_SESSION\" bash -c \"conda activate AmberTools23; python '$PYTHON_SCRIPT'; exec bash\""

# Check if SSH and Screen commands were successful
if [ $? -eq 0 ]; then
    echo "Successfully started run_gmx.py in Screen session '$SCREEN_SESSION' on $REMOTE_MACHINE in conda environment AmberTools23."
else
    echo "Error: Failed to execute run_gmx.py in Screen session on $REMOTE_MACHINE."
    exit 1
fi