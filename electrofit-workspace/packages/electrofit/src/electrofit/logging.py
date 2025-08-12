import logging


def setup_logging(log_file_path):
    """
    Sets up logging to output to both console and a log file.

    Parameters:
    - log_file_path (str): Path to the log file.
    """
    # Prevent adding multiple handlers if this function is called multiple times
    if not logging.getLogger().hasHandlers():
        # Create a custom logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        # Define the log message format
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

        # Console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)

        # File handler
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)

        # Add handlers to the logger
        logger.addHandler(ch)
        logger.addHandler(fh)


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
