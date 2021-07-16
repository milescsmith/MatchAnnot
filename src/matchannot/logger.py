import logging
from typing import Optional

import coloredlogs


def setup_logging(name: Optional[str] = None, level="DEBUG") -> logging.Logger:
    coloredlogs.DEFAULT_FIELD_STYLES = {
        "asctime": {"color": "green"},
        "levelname": {"bold": True, "color": "red"},
        "module": {"color": 73},
        "funcName": {"color": 74},
        "lineno": {"bold": True, "color": "green"},
        "message": {"color": "yellow"},
    }

    if name is None:
        name = __name__
    logger = logging.getLogger(name)
    logger.propagate = True

    if level == "DEBUG":
        logger_level = logging.DEBUG
    elif level == "INFO":
        logger_level = logging.INFO
    elif level == "WARNING":
        logger_level = logging.WARNING
    elif level == "ERROR":
        logger_level = logging.ERROR
    elif level == "CRITICAL":
        logger_level = logging.CRITICAL

    logger.setLevel(logger_level)
    coloredlogs.install(
        level=level,
        fmt="[%(funcName)s(%(lineno)d)]: %(message)s",
        logger=logger,
    )

    if name:
        fh = logging.FileHandler(f"{name}.log")
    else:
        fh = logging.FileHandler(f"{__name__}")
    formatter = logging.Formatter(
        "[%(asctime)s] {%(module)s:%(funcName)s():%(lineno)d} %(levelname)s - %(message)s"
    )
    fh.setLevel(logger_level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger
