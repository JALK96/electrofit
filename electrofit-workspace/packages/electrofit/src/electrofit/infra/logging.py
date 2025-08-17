"""Logging Infrastruktur (aus Top-Level logging.py verschoben)."""
from __future__ import annotations

import logging
import os
import subprocess
from logging.handlers import WatchedFileHandler
from pathlib import Path

try:
    from electrofit import __version__ as _electrofit_version
except Exception:  # pragma: no cover
    _electrofit_version = "unknown"


def setup_logging(log_path, also_console: bool = True, suppress_initial_message: bool = False) -> None:
    path = Path(log_path).resolve()
    root = logging.getLogger()
    env_level = os.getenv("ELECTROFIT_LOG_LEVEL", "INFO").upper()
    desired_level = getattr(logging, env_level, logging.INFO)
    if root.level > desired_level:
        root.setLevel(desired_level)
    effective_level = logging.getLevelName(root.level)
    to_remove = []
    existing_same = False
    for h in root.handlers:
        if isinstance(h, logging.FileHandler):
            try:
                existing_path = Path(getattr(h, "baseFilename", ""))
                if existing_path.resolve() == path:
                    existing_same = True
                else:
                    to_remove.append(h)
            except Exception:
                to_remove.append(h)
    for h in to_remove:
        root.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    if also_console:
        if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in root.handlers):
            fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            ch = logging.StreamHandler()
            ch.setLevel(root.level)
            ch.setFormatter(fmt)
            root.addHandler(ch)
    else:
        for h in [h for h in root.handlers if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)]:
            root.removeHandler(h)
            try:
                h.close()
            except Exception:
                pass
    if not existing_same:
        fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        fh = WatchedFileHandler(path, mode="a", encoding="utf-8", delay=False)
        fh.setLevel(root.level)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        if not suppress_initial_message:
            root.info(f"Logging initialized. Log file: {path} (level={effective_level})")


def _git_commit_short() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL, timeout=1)
        return out.decode().strip()
    except Exception:  # pragma: no cover
        return ""


def log_run_header(step_name: str):
    commit = _git_commit_short()
    parts = [f"electrofit { _electrofit_version }", f"step={step_name}"]
    if commit:
        parts.append(f"git={commit}")
    header = " | ".join(parts)
    logger = logging.getLogger()
    logger.info(header)
    try:
        for h in logger.handlers:
            if hasattr(h, 'baseFilename'):
                fp = getattr(h, 'baseFilename')
                if not fp:
                    continue
                try:
                    with open(fp, 'r+') as f:
                        lines = f.read().splitlines()
                        if header not in lines:
                            f.write(header + "\n")
                except Exception:
                    pass
    except Exception:  # pragma: no cover
        pass


def reset_logging():
    logging.shutdown()
    root_logger = logging.getLogger()
    handlers = root_logger.handlers[:]
    for handler in handlers:
        root_logger.removeHandler(handler)
        handler.close()
    logger_dict = logging.Logger.manager.loggerDict
    for logger_name in logger_dict:
        logger = logging.getLogger(logger_name)
        if hasattr(logger, "handlers"):
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
                handler.close()
            logger.filters = []
