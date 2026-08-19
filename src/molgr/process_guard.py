"""Reject fork children that inherited an initialized MolGR/Open Babel state."""

from __future__ import annotations

import os


_INITIAL_PROCESS_ID = os.getpid()


def ensure_current_process(api: str = "MolGR") -> None:
    """Ensure this call is running in the process that imported MolGR.

    ``spawn`` imports MolGR again and receives a new baseline PID.  A POSIX
    ``fork`` after import keeps the old baseline and is rejected before any
    Python Open Babel object is touched.
    """

    if os.getpid() != _INITIAL_PROCESS_ID:
        raise RuntimeError(
            f"{api} cannot run in a forked child after MolGR/Open Babel was "
            "initialized; use multiprocessing start method 'spawn' (or fork "
            "before importing MolGR)."
        )


__all__ = ["ensure_current_process"]
