"""
Author: TMJ
Date: 2026-06-28 22:54:16
LastEditors: TMJ
LastEditTime: 2026-06-28 22:56:58
Description: 请填写简介
"""
from __future__ import annotations

import signal
from collections.abc import Iterator
from contextlib import contextmanager


class CaseTimeoutError(TimeoutError):
    pass


@contextmanager
def case_timeout(seconds: float | None, label: str) -> Iterator[None]:
    if seconds is None or seconds <= 0 or not hasattr(signal, "setitimer"):
        yield
        return

    def _handle_timeout(signum, frame):  # noqa: ARG001
        raise CaseTimeoutError(f"{label} timed out after {seconds:.3f}s")

    previous_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _handle_timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous_handler)
