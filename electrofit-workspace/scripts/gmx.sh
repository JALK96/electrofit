#!/bin/bash
set -e

# -----------------------------
# Setup & logging
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$SCRIPT_DIR/process.log"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"; }
touch "$LOG_FILE"
log "=== Script Started ==="
log "Script Directory: $SCRIPT_DIR"

# -----------------------------
# Load env (written by Step 2)
# -----------------------------
ENV_FILE="$SCRIPT_DIR/.electrofit_env"
if [[ -f "$ENV_FILE" ]]; then
  # shellcheck disable=SC1090
  . "$ENV_FILE"
  log "Loaded env from $ENV_FILE"
else
  log "No .electrofit_env found; using defaults."
fi

# Defaults if not set
: "${REMOTE_HOST:=qcm04}"
: "${REMOTE_USER:=$USER}"
: "${REMOTE_SHELL_INIT:=~/.bashrc}"
: "${REMOTE_CONDA_ENV:=AmberTools23}"
: "${USE_SCREEN:=1}"
: "${SCREEN_NAME_PREFIX:=ef}"
: "${PYTHON_ENTRYPOINT:=python -m electrofit.external.gromacs}"

log "REMOTE_HOST=$REMOTE_HOST"
log "REMOTE_USER=$REMOTE_USER"
log "REMOTE_CONDA_ENV=$REMOTE_CONDA_ENV"
log "USE_SCREEN=$USE_SCREEN"
log "PYTHON_ENTRYPOINT=$PYTHON_ENTRYPOINT"

# Screen session name
PARENT_DIR="$(basename "$(dirname "$SCRIPT_DIR")")"
SCREEN_SESSION="${SCREEN_NAME_PREFIX}_${PARENT_DIR}_$(date +%Y%m%d%H%M%S)"
log "Screen Session Name: $SCREEN_SESSION"

# -----------------------------
# Build remote command
# -----------------------------
# Run from the directory where this script lives so relative paths resolve.
REMOTE_CMD="cd \"$SCRIPT_DIR\" && \
source \"$REMOTE_SHELL_INIT\" && \
conda activate \"$REMOTE_CONDA_ENV\" && \
$PYTHON_ENTRYPOINT"

if [[ "$USE_SCREEN" == "1" ]]; then
  REMOTE_CMD="cd \"$SCRIPT_DIR\" && source \"$REMOTE_SHELL_INIT\" && conda activate \"$REMOTE_CONDA_ENV\" && screen -dmS \"$SCREEN_SESSION\" bash -lc '$PYTHON_ENTRYPOINT; exec bash'"
fi

log "Remote Command: $REMOTE_CMD"

# -----------------------------
# Execute over SSH
# -----------------------------
SSH_TARGET="${REMOTE_USER}@${REMOTE_HOST}"
log "Connecting to $SSH_TARGET"
set +e
ssh "$SSH_TARGET" "$REMOTE_CMD"
SSH_STATUS=$?
set -e

if [[ $SSH_STATUS -ne 0 ]]; then
  log "Error: SSH command failed with status $SSH_STATUS."
  exit $SSH_STATUS
fi

if [[ "$USE_SCREEN" == "1" ]]; then
  log "Started in screen session '$SCREEN_SESSION' on $SSH_TARGET."
else
  log "Command completed on $SSH_TARGET."
fi

log "=== Script Completed Successfully ==="