"""
Author: TMJ
Date: 2026-04-19 15:46:50
LastEditors: TMJ
LastEditTime: 2026-04-19 22:15:02
Description: 请填写简介
"""

from __future__ import annotations

from functools import lru_cache
from typing import Callable, NamedTuple, Optional, cast

from typing_extensions import ParamSpec, Protocol, TypeVar


P = ParamSpec("P")
R_co = TypeVar("R_co", covariant=True)


class CacheInfo(NamedTuple):
    hits: int
    misses: int
    maxsize: Optional[int]
    currsize: int


class CachedFunc(Protocol[P, R_co]):
    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R_co: ...
    def cache_info(self) -> CacheInfo: ...
    def cache_clear(self) -> None: ...


def typed_lru_cache(
    maxsize: int = 128, typed: bool = False
) -> Callable[[Callable[P, R_co]], CachedFunc[P, R_co]]:
    def decorator(func: Callable[P, R_co]) -> CachedFunc[P, R_co]:
        cached = lru_cache(maxsize=maxsize, typed=typed)(func)
        return cast(CachedFunc[P, R_co], cached)

    return decorator


__all__ = ["CacheInfo", "typed_lru_cache"]
