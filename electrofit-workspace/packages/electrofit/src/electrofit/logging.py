import logging
from logging.handlers import WatchedFileHandler
from pathlib import Path


def setup_logging(log_path, also_console: bool = True) -> None:
    """
    Ensure the root logger writes to the given log file.
    Idempotent: if a FileHandler to this log_path already exists, do nothing.
    Optionally add a console handler once.
    """
    path = Path(log_path).resolve()
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # If we already have a FileHandler for this exact path, bail out
    for h in root.handlers:
        if isinstance(h, logging.FileHandler):
            try:
                if Path(getattr(h, "baseFilename", "")).resolve() == path:
                    return
            except Exception:
                pass

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # File handler (monitors file replacement/rotation safely)
    fh = WatchedFileHandler(path, mode="a", encoding="utf-8", delay=False)
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    # Add a single console handler if none exists yet
    if also_console and not any(isinstance(h, logging.StreamHandler) for h in root.handlers):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(fmt)
        root.addHandler(ch)

    # Log once that we’re wired to this file
    root.info(f"Logging initialized. Log file: {path}")


def reset_logging():
    """
    Reset the logging configuration by removing all handlers.
    This allows reinitialization of logging as if starting fresh.
    """
    # Shutdown logging and flush any pending messages
    logging.shutdown()

    # Remove handlers from the root logger
    root_logger = logging.getLogger()
    handlers = root_logger.handlers[:]
    for handler in handlers:
        root_logger.removeHandler(handler)
        handler.close()

    # Optionally, reset configurations of all existing loggers
    # This is useful if you have child loggers that also need resetting
    logger_dict = logging.Logger.manager.loggerDict
    for logger_name in logger_dict:
        logger = logging.getLogger(logger_name)
        if hasattr(logger, "handlers"):
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
                handler.close()
            # Remove filters if any
            logger.filters = []
