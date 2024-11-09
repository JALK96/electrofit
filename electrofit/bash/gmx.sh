#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Get the parent folder name for the screen session
SCREEN_SESSION="$(basename "$(dirname "$SCRIPT_DIR")")"

# Extract the path up to and including 'electrofit' from SCRIPT_DIR
PROJECT_ROOT="${SCRIPT_DIR%%electrofit*}electrofit"
# Name of the Python script to execute
PYTHON_SCRIPT="$PROJECT_ROOT/electrofit/execution/run_gmx.py"

# Build and execute the command on the remote server
ssh "qcm04" "cd \"$SCRIPT_DIR\" && source ~/.bashrc && screen -dmS \"$SCREEN_SESSION\" bash -c \"conda activate AmberTools23; python '$PYTHON_SCRIPT'; exec bash\""