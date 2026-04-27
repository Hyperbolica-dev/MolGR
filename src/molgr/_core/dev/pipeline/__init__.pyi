"""
Development-only pipeline helpers
"""

from __future__ import annotations

from . import reconstruct_with_metals, reconstruct_without_metals, resonance

__all__: list[str] = [
    "clear_force_field_evaluation_cache",
    "clear_resonance_move_score_cache",
    "clear_uff_atom_typing_cache",
    "get_uff_atom_typing_cache_info",
    "reconstruct_with_metals",
    "reconstruct_without_metals",
    "resonance",
]

def clear_force_field_evaluation_cache() -> None:
    """
    Clear the C++ force-field evaluation cache.
    """

def clear_resonance_move_score_cache() -> None:
    """
    Clear the C++ resonance move score cache.
    """

def clear_uff_atom_typing_cache() -> None:
    """
    Clear the C++ UFF atom-typing assignment cache.
    """

def get_uff_atom_typing_cache_info() -> dict:
    """
    Return C++ UFF atom-typing cache hit/miss/size counters.
    """
