"""
Pipeline-level helpers
"""

from __future__ import annotations

import collections.abc
import typing

import molgr
import molgr._core.utils

__all__: list[str] = [
    "get_last_run_timing_breakdown_ms",
    "get_radical_resonances_ptr",
    "process_resonance_ptr",
    "reconstruct_with_metals",
    "reconstruct_without_metals",
    "resonance",
    "smiles_token_ptr",
]

class reconstruct_with_metals:
    """
    Fallback-aligned reconstruction helpers with metals
    """
    class MetalAtomPosition:
        symbol: str
        def __repr__(self) -> str: ...
        @property
        def element_idx(self) -> int: ...
        @element_idx.setter
        def element_idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
        @property
        def idx(self) -> int: ...
        @idx.setter
        def idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
        @property
        def position_x(self) -> float: ...
        @position_x.setter
        def position_x(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
        @property
        def position_y(self) -> float: ...
        @position_y.setter
        def position_y(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
        @property
        def position_z(self) -> float: ...
        @position_z.setter
        def position_z(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
        @property
        def radical_num(self) -> int: ...
        @radical_num.setter
        def radical_num(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
        @property
        def valence(self) -> int: ...
        @valence.setter
        def valence(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...

    class MetalHandler:
        @staticmethod
        def combine_metal_with_mol(
            mol_ptr: typing.SupportsInt | typing.SupportsIndex,
            metals: collections.abc.Sequence[reconstruct_with_metals.MetalAtomPosition],
        ) -> None: ...
        def __init__(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None: ...
        def generate_combinations(
            self, total_radical_electrons: typing.SupportsInt | typing.SupportsIndex
        ) -> list[list[reconstruct_with_metals.MetalAtomPosition]]: ...
        def strip_metals(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str: ...

    @staticmethod
    def build_metal_states_ptr(
        mol_ptr: typing.SupportsInt | typing.SupportsIndex,
        atom_idx: typing.SupportsInt | typing.SupportsIndex,
    ) -> list[reconstruct_with_metals.MetalAtomPosition]: ...
    @staticmethod
    def combine_metal_with_omol_ptr(
        mol_ptr: typing.SupportsInt | typing.SupportsIndex,
        metals: collections.abc.Sequence[reconstruct_with_metals.MetalAtomPosition],
    ) -> None: ...
    @staticmethod
    def get_possible_metal_radicals(
        metal: str, valence: typing.SupportsInt | typing.SupportsIndex
    ) -> set[int]: ...
    @staticmethod
    def xyz2omol(
        xyz_block: str,
        total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
        total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> molgr._core.utils.MoleculeData: ...

class reconstruct_without_metals:
    """
    Fallback-aligned no-metal reconstruction helpers
    """
    @staticmethod
    def xyz_to_omol_no_metal(
        xyz_block: str,
        total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
        total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> molgr._core.utils.MoleculeData: ...

class resonance:
    """
    Fallback-aligned resonance helpers
    """
    @staticmethod
    def get_radical_resonances_smi(smiles: str) -> list[str]: ...
    @staticmethod
    def process_resonance_smi(
        smiles: str, charge: typing.SupportsInt | typing.SupportsIndex
    ) -> tuple[str, int]: ...

def get_last_run_timing_breakdown_ms() -> dict: ...
def get_radical_resonances_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> list[int]: ...
def process_resonance_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    charge: typing.SupportsInt | typing.SupportsIndex,
) -> tuple: ...
def smiles_token_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str: ...
