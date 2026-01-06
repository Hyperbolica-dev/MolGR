"""
Chemical constants
"""
from __future__ import annotations
import typing
__all__: list[str] = ['get_possible_metal_radicals']
def get_possible_metal_radicals(metal: str, valence: typing.SupportsInt) -> ...:
    """
                Get possible radical electron counts for a metal given its valence.
                
                Args:
                    metal (str): The chemical symbol (e.g., "Fe").
                    valence (int): The oxidation state.
                
                Returns:
                    set[int]: A set of possible unpaired electron counts.
    """
