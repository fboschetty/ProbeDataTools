"""Tools for processing electron-microprobe analyses and mineral compositions."""

from importlib.metadata import version

__version__ = version("probedatatools")

from .probe import ProbeData

__all__ = [
    "ProbeData",
]