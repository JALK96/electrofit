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
