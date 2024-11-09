#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Name of the screen session: parent of parent directory basename and current directory basename
SCREEN_SESSION="$(basename "$(dirname "$(dirname "$SCRIPT_DIR")")")_$(basename "$SCRIPT_DIR")"

# Extract the path up to and including 'electrofit' from SCRIPT_DIR
PROJECT_ROOT="${SCRIPT_DIR%%electrofit*}electrofit"

# Name of the Python script to execute
PYTHON_SCRIPT="$PROJECT_ROOT/electrofit/execution/run_process_conform.py"

# Remote server details
REMOTE_HOST="qcm08"

# Execute the command on the remote server via SSH
ssh "$REMOTE_HOST" "cd \"$SCRIPT_DIR\" && source ~/.bashrc && screen -dmS \"$SCREEN_SESSION\" bash -c \"conda activate AmberTools23; python '$PYTHON_SCRIPT'; exec bash\""