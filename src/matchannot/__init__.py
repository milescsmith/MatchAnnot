import typer

try:
    from matchannot.__about__ import __author__, __email__, __version__
except ImportError:
    __author__ = "Tom Skelly"
    __email__ = "thomas.skelly@nih.gov"
    __version__ = "unknown"

import logging

from matchannot.logger import setup_logging

matchannot_logger = setup_logging("matchannot")


def version_callback(value: bool) -> None:
    """Prints the version of the package."""
    if value:
        print(f"{__name__} version: {__version__}")
        raise typer.Exit()


def set_verbosity(v: int = 3) -> None:
    if v == 0:
        matchannot_logger.setLevel(logging.CRITICAL)
    elif v == 1:
        matchannot_logger.setLevel(logging.ERROR)
    elif v == 2:
        matchannot_logger.setLevel(logging.WARNING)
    elif v == 3:
        matchannot_logger.setLevel(logging.INFO)
    elif v == 4:
        matchannot_logger.setLevel(logging.DEBUG)
