"""Common logging functions."""

import sys
import logging
import tempfile
import subprocess
import pylib.util as util

LOGGER = None  # Global logger so we can switch between processes
FORMATTER = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
NAME = 'pdc_logger'


def setup(log_file):
    """Logger setup."""
    global LOGGER  # pylint: disable=global-statement

    if not LOGGER:
        handler = logging.FileHandler(log_file)
        handler.setFormatter(FORMATTER)
        handler.setLevel(logging.DEBUG)

        stream = logging.StreamHandler()
        stream.setFormatter(FORMATTER)
        stream.setLevel(logging.INFO)

        LOGGER = logging.getLogger(log_file)
        LOGGER.setLevel(logging.DEBUG)
        LOGGER.addHandler(handler)
        LOGGER.addHandler(stream)

        info('#' * 80)
        info('{} version: {}'.format(util.__TITLE__, util.__VERSION__))
        info('Python version: {}'.format(' '.join(sys.version.split())))
        info(' '.join(sys.argv[:]))


def subcommand(cmd, temp_dir=None, timeout=None):
    """
    Call a subprocess and log the output.

    Note: stdout=PIPE is blocking and large logs cause a hang.
    So we don't use it.
    """
    LOGGER.info(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:
        try:
            subprocess.check_call(
                cmd,
                shell=True,
                timeout=timeout,
                stdout=log_output,
                stderr=log_output)
        finally:
            with open(log_output.name) as log_input:
                for line in log_input:
                    line = line.strip()
                    if line:
                        LOGGER.debug(line)


def capture(cmd, out_path, temp_dir=None, timeout=None):
    """
    Call a subprocess and capture the output.

    Note: stdout=PIPE is blocking and large logs cause a hang.
    So we don't use it.
    """
    LOGGER.info(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:
        with open(out_path, 'w') as out_file:
            try:
                subprocess.check_call(
                    cmd,
                    shell=True,
                    timeout=timeout,
                    stdout=out_file,
                    stderr=log_output)
            finally:
                with open(log_output.name) as log_input:
                    for line in log_input:
                        line = line.strip()
                        if line:
                            LOGGER.debug(line)


def info(msg):
    """Log an info message."""
    LOGGER.info(msg)


def warn(msg):
    """Log a warning message."""
    LOGGER.warning(msg)


def error(msg):
    """Log an error message."""
    LOGGER.error(msg)


def fatal(msg):
    """Log an error message and exit."""
    error(msg)
    sys.exit(1)
