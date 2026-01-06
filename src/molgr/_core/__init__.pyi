"""

MolGR Core C++ Implementation
-----------------------------
Exposes optimized C++ algorithms for MolOP Graph Reconstruction.

"""

from __future__ import annotations

import typing

from . import consts, metal, reconstruct, scoring, utils

__all__: list[str] = [
    "DEBUG",
    "ERROR",
    "INFO",
    "LogLevel",
    "OFF",
    "WARN",
    "consts",
    "free_obmol_ptr",
    "metal",
    "reconstruct",
    "scoring",
    "set_log_level",
    "utils",
]

class LogLevel:
    """
    Members:

      DEBUG

      INFO

      WARN

      ERROR

      OFF
    """

    DEBUG: typing.ClassVar[LogLevel]  # value = <LogLevel.DEBUG: 0>
    ERROR: typing.ClassVar[LogLevel]  # value = <LogLevel.ERROR: 3>
    INFO: typing.ClassVar[LogLevel]  # value = <LogLevel.INFO: 1>
    OFF: typing.ClassVar[LogLevel]  # value = <LogLevel.OFF: 4>
    WARN: typing.ClassVar[LogLevel]  # value = <LogLevel.WARN: 2>
    __members__: typing.ClassVar[
        dict[str, LogLevel]
    ]  # value = {'DEBUG': <LogLevel.DEBUG: 0>, 'INFO': <LogLevel.INFO: 1>, 'WARN': <LogLevel.WARN: 2>, 'ERROR': <LogLevel.ERROR: 3>, 'OFF': <LogLevel.OFF: 4>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

def free_obmol_ptr(arg0: typing.SupportsInt) -> None:
    """
    Manually delete the OBMol pointer
    """

def set_log_level(level: LogLevel) -> None:
    """
    Set the logging level for the C++ core (DEBUG=0, INFO=1, WARN=2, ERROR=3, OFF=4)
    """

DEBUG: LogLevel  # value = <LogLevel.DEBUG: 0>
ERROR: LogLevel  # value = <LogLevel.ERROR: 3>
INFO: LogLevel  # value = <LogLevel.INFO: 1>
OFF: LogLevel  # value = <LogLevel.OFF: 4>
WARN: LogLevel  # value = <LogLevel.WARN: 2>
