#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# -----------------------------
# Configuration and Setup
# -----------------------------

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define the log file
LOG_FILE="$SCRIPT_DIR/process.log"

# Function to log messages with timestamps
log() {
    local MESSAGE="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $MESSAGE" | tee -a "$LOG_FILE"
}

# -----------------------------
# Start of the Script
# -----------------------------

log "=== Script Started ==="


log "Script Directory: $SCRIPT_DIR"

# Name of the screen session: parent of parent directory basename and current directory basename
PARENT_DIR="$(basename "$(dirname "$SCRIPT_DIR")")"
CURRENT_DIR="$(basename "$SCRIPT_DIR")"
SCREEN_SESSION="${PARENT_DIR}_${CURRENT_DIR}"
log "Screen Session Name: $SCREEN_SESSION"

# Extract the path up to and including 'electrofit' from SCRIPT_DIR
PROJECT_ROOT="${SCRIPT_DIR%%electrofit*}electrofit"
log "Project Root: $PROJECT_ROOT"

# Name of the Python script to execute
PYTHON_SCRIPT="$PROJECT_ROOT/electrofit/execution/run_process_conform.py"
log "Python Script Path: $PYTHON_SCRIPT"

# Remote server details
REMOTE_HOST="qcm05"
log "Remote Host: $REMOTE_HOST"

# -----------------------------
# Function Definitions
# -----------------------------

# Function to execute SSH command
execute_ssh_command() {
    local HOST="$1"
    local CMD="$2"

    log "Connecting to $HOST to execute the command."
    ssh "$HOST" "$CMD"
    SSH_STATUS=$?
    
    if [ $SSH_STATUS -ne 0 ]; then
        log "Error: SSH command failed with status $SSH_STATUS."
        exit 1
    else
        log "SSH command executed successfully on $HOST."
    fi
}

# -----------------------------
# Build the SSH Command
# -----------------------------

SSH_COMMAND="cd \"$SCRIPT_DIR\" && \
source ~/.bashrc && \
screen -dmS \"$SCREEN_SESSION\" bash -c \"conda activate AmberTools23; python '$PYTHON_SCRIPT'; exec bash\""

log "Built SSH Command: $SSH_COMMAND"

# -----------------------------
# Execute the SSH Command
# -----------------------------

execute_ssh_command "$REMOTE_HOST" "$SSH_COMMAND"

log "Successfully started run_process_conform.py in Screen session '$SCREEN_SESSION' on $REMOTE_HOST."

# -----------------------------
# End of the Script
# -----------------------------

log "=== Script Completed Successfully ==="