r"""*Single-source version number for* ``matchannot``.
"""
try:
    from importlib.metadata import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    from importlib_metadata import metadata

__author__ = metadata("matchannot")["Author"]
__email__ = metadata("matchannot")["Author-email"]
__version__ = metadata("matchannot")["Version"]
