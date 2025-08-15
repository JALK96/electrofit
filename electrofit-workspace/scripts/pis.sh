#!/usr/bin/env bash
set -euo pipefail
[[ "${EF_DEBUG:-0}" == "1" ]] && set -x

# --- tiny logger ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$SCRIPT_DIR/process.log"
: >"$LOG_FILE"
log(){ echo "$(date '+%F %T') - $*" | tee -a "$LOG_FILE"; }
fail(){ log "FATAL: $*"; exit 2; }
trap 's=$?; [[ $s -ne 0 ]] && log "ERROR: pis.sh exited with status $s"; exit $s' EXIT

log "=== pis.sh started ==="
log "SCRIPT_DIR=$SCRIPT_DIR"

# --- project + config paths ---
# NOTE: when called directly, three levels up is the project root (…/tests/integration)
PROJECT_DIR="${ELECTROFIT_PROJECT_PATH:-"$(cd "$SCRIPT_DIR/../../.." && pwd)"}"
CONFIG_FILE="${ELECTROFIT_CONFIG_PATH:-"$SCRIPT_DIR/electrofit.toml"}"
[[ -f "$CONFIG_FILE" ]] || fail "config not found: $CONFIG_FILE"
log "PROJECT_DIR=$PROJECT_DIR"
log "CONFIG_FILE=$CONFIG_FILE"

# --- defaults (overridden by TOML) ---
REMOTE_HOST=""
CONDA_ENV="AmberTools23"
USE_SCREEN=1

# --- read [compute] from TOML ---
PY_IN="${PYTHON_EXE:-python3}"
set +e
readarray -t CFG < <("$PY_IN" - "$CONFIG_FILE" <<'PY'
import sys
try:
    try: import tomllib
    except ModuleNotFoundError: import tomli as tomllib
    with open(sys.argv[1],'rb') as f: d=tomllib.load(f)
    c=d.get('compute') or {}
    print(c.get('remote_host',''))
    print(c.get('conda_env','AmberTools23'))
    v=c.get('use_screen',True)
    print('1' if (v is True or str(v).lower() in ('1','true','yes','y','on')) else '0')
except Exception as e:
    print('__ERR__:'+repr(e))
PY
)
st=$?
set -e
if [[ $st -ne 0 || "${CFG[0]}" =~ ^__ERR__ ]]; then
  log "TOML parse warning: ${CFG[0]#__ERR__:}; using defaults."
else
  REMOTE_HOST="${CFG[0]}"
  CONDA_ENV="${CFG[1]}"
  USE_SCREEN="${CFG[2]}"
fi

# Normalize markers for local
case "${REMOTE_HOST,,}" in ""|"local"|"<local>") REMOTE_HOST="";; esac

log "compute.remote_host=${REMOTE_HOST:-<local>}"
log "compute.conda_env=$CONDA_ENV"
log "compute.use_screen=$USE_SCREEN"

# --- helpers ---
conda_init() {
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1090
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1090
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/.bashrc" ]]; then
    # shellcheck disable=SC1090
    source "$HOME/.bashrc"
  fi
}

ensure_conda_env() {
  conda_init
  if ! command -v conda >/dev/null 2>&1; then
    fail "conda not found; cannot activate env '$CONDA_ENV'"
  fi
  log "Activating conda env: $CONDA_ENV"
  set +e
  conda activate "$CONDA_ENV"
  rc=$?
  set -e
  [[ $rc -ne 0 ]] && fail "failed to activate conda env '$CONDA_ENV'"
}

# --- local launch ---
run_local() {
  log "Mode: LOCAL"
  cd "$SCRIPT_DIR"
  ensure_conda_env

  if [[ "$USE_SCREEN" == "1" && ! $(command -v screen) ]]; then
    log "screen not found locally; running in foreground"
    USE_SCREEN=0
  fi

  if [[ "$USE_SCREEN" == "1" ]]; then
    SESSION="ef_pis_local_$(date +%Y%m%d%H%M%S)"
    CMD="cd '$SCRIPT_DIR' && electrofit run-process-initial-structure --project '$PROJECT_DIR' --config '$CONFIG_FILE'"
    log "Launching in screen session: $SESSION"
    screen -dmS "$SESSION" bash -lc "$CMD; exec bash"
    log "Started. Attach with: screen -r $SESSION"
  else
    log "Running in foreground…"
    electrofit run-process-initial-structure --project "$PROJECT_DIR" --config "$CONFIG_FILE"
  fi
}

# --- remote launch (shared FS assumed for now) ---
run_remote() {
  log "Mode: REMOTE (shared FS)"
  [[ -n "$REMOTE_HOST" ]] || fail "remote_host empty"

  SESSION="ef_pis_$(basename "$(dirname "$SCRIPT_DIR")")_$(date +%Y%m%d%H%M%S)"
  log "SSH target: $REMOTE_HOST"
  log "Screen session: $SESSION"

  # Build a robust remote command without fragile inlined quotes:
  # 1) create a small launcher on the remote; 2) start it in screen -dmS
  ssh "$REMOTE_HOST" bash -lc "$(printf 'SDIR=%q; PROJ=%q; CFILE=%q; CENV=%q; SESSION=%q; ' \
    "$SCRIPT_DIR" "$PROJECT_DIR" "$CONFIG_FILE" "$CONDA_ENV" "$SESSION"; cat <<'BASH'
set -e
# Init conda in this shell so 'conda activate' works
if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/.bashrc" ]]; then
  source "$HOME/.bashrc"
fi

# Create the launcher script in the run directory
mkdir -p "$SDIR"
cat > "$SDIR/.ef_launch.sh" <<'LAUNCH'
#!/usr/bin/env bash
set -e
if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/.bashrc" ]]; then
  source "$HOME/.bashrc"
fi
conda activate "$CENV"
cd "$SDIR"
exec electrofit run-process-initial-structure --project "$PROJ" --config "$CFILE"
LAUNCH
chmod +x "$SDIR/.ef_launch.sh"

# Prefer screen; fall back to nohup if needed
if command -v screen >/dev/null 2>&1; then
  screen -dmS "$SESSION" bash -lc "$SDIR/.ef_launch.sh; exec bash"
  echo "SCREEN_OK $SESSION"
else
  nohup bash -lc "$SDIR/.ef_launch.sh" > "$SDIR/nohup.out" 2>&1 < /dev/null &
  echo "NO_SCREEN"
fi
BASH
)"

  log "Remote launch issued. You can attach with: ssh $REMOTE_HOST 'screen -r $SESSION' (if screen is available)."
}

# --- dispatch ---
if [[ -z "${REMOTE_HOST}" ]]; then
  log "Dispatch: run_local()"
  run_local
else
  log "Dispatch: run_remote()"
  run_remote
fi

log "=== pis.sh done ==="