import io
import os
import logging
import subprocess
from contextlib import redirect_stderr, redirect_stdout, contextmanager
from pathlib import Path



def _snapshot(path: str):
    p = Path(path)
    # record names, sizes, mtimes (you can extend as needed)
    return { (f.name, f.stat().st_size, int(f.stat().st_mtime)) for f in p.iterdir() }

@contextmanager
def _pushd(path: str | None):
    prev = os.getcwd()
    if path:
        os.chdir(path)
    try:
        yield
    finally:
        if path:
            os.chdir(prev)

def run_python(func, *args, cwd=None, log_level=logging.INFO, **kwargs):
    """
    Run a Python callable with the same ergonomics as `run_command`:
    - capture stdout/stderr into the logs,
    - log newly created files in `cwd`.
    """
    logging.info(f"Executing python: {func.__module__}.{func.__name__}(*args, **kwargs)")
    before = _snapshot(cwd or ".")
    buf = io.StringIO()
    with _pushd(cwd), redirect_stdout(buf), redirect_stderr(buf):
        result = func(*args, **kwargs)
    output = buf.getvalue()
    if output:
        logging.log(log_level, output.rstrip())

    after = _snapshot(cwd or ".")
    new = after - before
    if new:
        # List only file names; you can format with sizes/mtimes if useful
        logging.info("New files/folders created: " + ", ".join(sorted(n for (n, *_ ) in new)))
    else:
        logging.info("No new files/folders created.")
    return result


def check_gaussian_convergence(gaussian_log_path):
    """
    Checks if Gaussian SCF calculation converged.

    Parameters:
    - gaussian_log_path (str): Path to the Gaussian log file.

    Returns:
    - bool: True if converged, False otherwise.
    """
    try:
        with open(gaussian_log_path, "r") as log_file:
            for line in log_file:
                if "SCF Done" in line and "Converged" in line:
                    return True
        return False
    except FileNotFoundError:
        logging.error(f"Gaussian log file '{gaussian_log_path}' not found.")
        return False


def run_command(command, cwd=None):
    """
    Runs a shell command and handles errors, logging any new files or directories created.

    Parameters:
    - command (str): The shell command to execute.
    - cwd (str, optional): The working directory in which to execute the command.

    Returns:
    - output (str): Combined standard output from the command.

    Raises:
    - subprocess.CalledProcessError: If the command exits with a non-zero status.
    """
    try:
        logging.info(f"Executing command: {command}")

        # Before running the command, list files and directories in cwd
        if cwd:
            before_items = set(os.listdir(cwd))
        else:
            before_items = set(os.listdir())

        # Initialize Popen with stdout and stderr as pipes for real-time logging
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
        )

        # Collecting output and logging in real-time
        output_lines = []
        for stdout_line in iter(process.stdout.readline, ""):
            logging.info(stdout_line.strip())
            output_lines.append(stdout_line)
        for stderr_line in iter(process.stderr.readline, ""):
            logging.debug(stderr_line.strip())
            output_lines.append(stderr_line)

        # Ensure process has finished and capture the return code
        process.stdout.close()
        process.stderr.close()
        return_code = process.wait()

        # Handle exit status
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)

        # After running the command, list files and directories in cwd
        if cwd:
            after_items = set(os.listdir(cwd))
        else:
            after_items = set(os.listdir())

        # Determine new files/folders created
        new_items = after_items - before_items
        if new_items:
            logging.info(
                f"New files/folders created by '{command}': {', '.join(new_items)}"
            )
        else:
            logging.info(f"No new files/folders created by '{command}'.")

        return "".join(output_lines)

    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing command '{command}': {e.stderr}")
        raise


def run_and_log(func, *args, log_level=logging.INFO, **kwargs):
    """
    Runs the given function while capturing its stdout and stderr output,
    then logs that output at the specified log level.

    Parameters
    ----------
    func : callable
        The function to execute.
    *args, **kwargs :
        Arguments and keyword arguments to pass to the function.
    log_level : int, optional
        The logging level to use for the captured output (default is logging.INFO).

    Returns
    -------
    result
        The return value of the function.

    If the function raises an exception, its output will be logged before the exception is re-raised.
    """
    f = io.StringIO()
    try:
        with redirect_stdout(f), redirect_stderr(f):
            result = func(*args, **kwargs)
    except Exception:
        output = f.getvalue()
        if output:
            logging.log(log_level, f"Output before exception:\n{output}")
        logging.error("Function raised an exception", exc_info=True)
        raise
    else:
        output = f.getvalue()
        if output:
            logging.log(log_level, f"Captured output:\n{output}")
        return result


def create_gaussian_input(
    mol2_file, molecule_name, net_charge, scratch_dir, atom_type="gaff2"
):
    """
    Runs Antechamber to generate Gaussian input file from mol2-file.

    Parameters:
    - mol2_file (str): Path to the input mol2 file.
    - molecule_name (str): Name of the molecule to be used in output files.
    - net_charge (int): Net charge of the molecule.
    - scratch_dir (str): Directory where the command is executed.
    """

    command = (
        f"antechamber -i {mol2_file} -fi mol2 -nc {net_charge} -at {atom_type} "
        f"-o {molecule_name}.gcrt -fo gcrt -gv 1 -ge {molecule_name}.gesp"
    )
    run_command(command, cwd=scratch_dir)
    logging.info("Antechamber processing completed.")

    # Check if expected files were created
    expected_files = [f"{molecule_name}.gcrt"]
    created_files = []
    for file in expected_files:
        file_path = os.path.join(scratch_dir, file)
        if os.path.isfile(file_path):
            created_files.append(file)
        else:
            logging.warning(f"Expected file '{file}' was not created.")

    if created_files:
        logging.debug(
            f"Expected files created by antechamber: {', '.join(created_files)}"
        )
    else:
        logging.error("No expected files were created by antechamber.")


def gaussian_out_to_prepi(g_out, scratch_dir, prepi_file=None):
    """
    Converts Gaussian output to prepi file using antechamber.

    Parameters:
    - g_out (str): Path to Gaussian output file.
    - scratch_dir (str): Directory where the command is executed.
    - prepi_file (str, optional): Desired name for the prepi file. If None, derived from g_out.
    """
    # Convert g_out to a Path object for easier manipulation
    g_out_path = Path(g_out)

    if prepi_file is None:
        # Derive prepi_file name by replacing the extension with .prepi
        prepi_file = g_out_path.with_suffix(".prepi").name
    else:
        # Ensure prepi_file is a string (in case a Path is passed)
        prepi_file = str(prepi_file)

    command = f"antechamber -i {g_out} -fi gout -o {prepi_file} -fo prepi -c resp -s 2"
    run_command(command, cwd=scratch_dir)
    logging.info(f"Prepi file '{prepi_file}' generated successfully.")


def generate_mol2_with_resp_charges(g_out, mol2_output, scratch_dir, atom_type="gaff2"):
    """
    Generates mol2 file with RESP charges.

    Parameters:
    - g_out (str): Path to Gaussian output file.
    - mol2_output (str): Desired name for the output mol2 file with RESP charges.
    - scratch_dir (str): Directory where the command is executed.
    """
    command = f"antechamber -i {g_out} -fi gout -o {mol2_output} -fo mol2 -c resp -s 2 -at {atom_type}"
    run_command(command, cwd=scratch_dir)
    logging.info(f"Mol2 file '{mol2_output}' with RESP charges generated successfully.")


def run_gaussian_calculation(input_file, molecule_name, scratch_dir):
    """
    Runs Gaussian calculation on the input file.

    Parameters:
    - input_file (str): Gaussian input file.
    - scratch_dir (str): Directory where the command is executed.
    """
    # Run Gaussian
    command = f"rung16 {input_file}"
    run_command(command, cwd=scratch_dir)
    logging.info("Gaussian calculation completed.")

    # Check for Gaussian output files
    log_file = f"{input_file}.log"
    chk_file = "molecule.chk"
    gesp_file = f"{molecule_name}.gesp"
    expected_files = [log_file, chk_file, gesp_file]

    log_file_path = os.path.join(scratch_dir, log_file)
    check_gaussian_convergence(log_file_path)

    created_files = []
    for file in expected_files:
        file_path = os.path.join(scratch_dir, file)
        if os.path.isfile(file_path):
            created_files.append(file)
        else:
            logging.warning(f"Expected file '{file}' was not created.")

    if created_files:
        logging.debug(f"Expected files created by Gaussian: {', '.join(created_files)}")
    else:
        logging.error("No expected files were created by Gaussian.")


def run_espgen(gesp_file, esp_file, scratch_dir):
    """
    Runs espgen to convert .gesp to .esp.

    Parameters:
    - gesp_file (str): Input .gesp file.
    - esp_file (str): Output .esp file.
    - scratch_dir (str): Directory where the command is executed.
    """
    command = f"espgen -i {gesp_file} -o {esp_file}"
    run_command(command, cwd=scratch_dir)
    logging.info("espgen processing completed.")


def run_acpype(mol2_file, net_charge, scratch_dir, atom_type="gaff2", charges="bcc"):
    """
    Runs acpype to create GROMACS input.

    Parameters:
    - mol2_file (str): Input mol2 file with RESP charges.
    - net_charge (int): Net charge of the molecule.
    - scratch_dir (str): Directory where the command is executed.
    - atom_type (str): Atom Type for simulation input (gaff2, amber, gaff).
    - charges (str): Charge Method: gas, bcc (default), user (user's charges in mol2 file).
    """

    command = (
        f"acpype -i {mol2_file} -n {net_charge} -a {atom_type} -c {charges} -o gmx"
    )
    run_command(command, cwd=scratch_dir)
    logging.info("acpype processing completed.")


def run_resp(
    respin_file,
    esp_file,
    resp_output,
    resp_pch,
    resp_chg,
    scratch_dir,
    previous_charges=None,
):
    """
    Runs RESP fitting using the given respin file.

    Parameters:
    - respin_file (str): RESP input file.
    - esp_file (str): ESP file.
    - resp_output (str): RESP output file.
    - resp_pch (str): RESP point charge file.
    - resp_chg (str): RESP charge file.
    - previous_charges (str): RESP output of first stage resp fitting.
    - scratch_dir (str): Directory where the command is executed.
    """
    # Always use basenames in command to avoid RESP issues with long/absolute paths
    rf_abs = respin_file
    ef_abs = esp_file
    pq_abs = previous_charges
    rf = os.path.basename(respin_file)
    ef = os.path.basename(esp_file)
    pq = os.path.basename(previous_charges) if previous_charges else None
    if rf_abs != rf or ef_abs != ef:
        logging.debug(
            "[resp-path] Using basenames for RESP command: -i %s -e %s (orig -i %s -e %s)",
            rf, ef, rf_abs, ef_abs,
        )
    if pq_abs and pq_abs != pq:
        logging.debug(
            "[resp-path] Using basename for previous charges: -q %s (orig %s)", pq, pq_abs
        )
    # Sanity: ensure files exist in cwd
    for needed in [rf, ef] + ([pq] if pq else []):
        if not os.path.isfile(os.path.join(scratch_dir, needed)):
            logging.error("[resp-path] Required file '%s' not found in scratch; aborting RESP", needed)
            raise FileNotFoundError(needed)
    command = (
        f"resp -O -i {rf} -o {resp_output} "
        f"-p {resp_pch} -t {resp_chg} -e {ef}"
    )
    if pq:
        command += f" -q {pq}"
    run_command(command, cwd=scratch_dir)
    # Validate expected output artifacts
    missing = [f for f in (resp_output, resp_pch, resp_chg) if not os.path.isfile(os.path.join(scratch_dir, f))]
    if missing:
        # One retry (only if first run used absolute original path forms, but we already switched to basename) -> list directory
        listing = sorted(os.listdir(scratch_dir))
        logging.error(
            "[resp-validate] Missing expected RESP output file(s): %s (cwd listing: %s)",
            ", ".join(missing), ", ".join(listing)
        )
        raise RuntimeError(f"RESP failed to produce: {', '.join(missing)}")
    logging.info(f"RESP fitting using '{respin_file}' completed.")
