"""
Pipeline-level helpers
"""

from __future__ import annotations

import collections.abc
import typing

from .utils import MoleculeData

__all__: list[str] = [
    "reconstruct_without_metals",
    "reconstruct_with_metals",
    "resonance",
    "get_last_run_timing_breakdown_ms",
    "get_radical_resonances_ptr",
    "process_resonance_ptr",
    "smiles_token_ptr",
]

class _ReconstructWithoutMetalsModule:
    def xyz_to_omol_no_metal(
        self,
        xyz_block: str,
        total_charge: typing.SupportsInt,
        total_radical_electrons: typing.SupportsInt,
    ) -> MoleculeData | None: ...

class _ReconstructWithMetalsModule:
    class MetalAtomPosition:
        symbol: str
        def __repr__(self) -> str: ...
        @property
        def element_idx(self) -> int: ...
        @element_idx.setter
        def element_idx(self, arg0: typing.SupportsInt) -> None: ...
        @property
        def idx(self) -> int: ...
        @idx.setter
        def idx(self, arg0: typing.SupportsInt) -> None: ...
        @property
        def radical_num(self) -> int: ...
        @radical_num.setter
        def radical_num(self, arg0: typing.SupportsInt) -> None: ...
        @property
        def valence(self) -> int: ...
        @valence.setter
        def valence(self, arg0: typing.SupportsInt) -> None: ...
        @property
        def position_x(self) -> float: ...
        @position_x.setter
        def position_x(self, arg0: typing.SupportsFloat) -> None: ...
        @property
        def position_y(self) -> float: ...
        @position_y.setter
        def position_y(self, arg0: typing.SupportsFloat) -> None: ...
        @property
        def position_z(self) -> float: ...
        @position_z.setter
        def position_z(self, arg0: typing.SupportsFloat) -> None: ...

    class MetalHandler:
        @staticmethod
        def combine_metal_with_mol(
            mol_ptr: typing.SupportsInt,
            metals: collections.abc.Sequence[_ReconstructWithMetalsModule.MetalAtomPosition],
        ) -> None: ...
        def __init__(self, mol_ptr: typing.SupportsInt) -> None: ...
        def generate_combinations(
            self, total_radical_electrons: typing.SupportsInt
        ) -> list[list[_ReconstructWithMetalsModule.MetalAtomPosition]]: ...
        def strip_metals(self, mol_ptr: typing.SupportsInt) -> str: ...

    def get_possible_metal_radicals(
        self,
        metal: str,
        valence: typing.SupportsInt,
    ) -> set[int]: ...
    def build_metal_states_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        atom_idx: typing.SupportsInt,
    ) -> list[_ReconstructWithMetalsModule.MetalAtomPosition]: ...
    def combine_metal_with_omol_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        metals: collections.abc.Sequence[_ReconstructWithMetalsModule.MetalAtomPosition],
    ) -> None: ...
    def xyz2omol(
        self,
        xyz_block: str,
        total_charge: typing.SupportsInt,
        total_radical_electrons: typing.SupportsInt,
    ) -> MoleculeData | None: ...

class _ResonanceModule:
    def get_radical_resonances_smi(self, smiles: str) -> list[str]: ...
    def process_resonance_smi(self, smiles: str, charge: typing.SupportsInt) -> tuple[str, int]: ...

reconstruct_without_metals: _ReconstructWithoutMetalsModule
reconstruct_with_metals: _ReconstructWithMetalsModule
resonance: _ResonanceModule

def get_last_run_timing_breakdown_ms() -> dict[str, float]: ...
def get_radical_resonances_ptr(mol_ptr: typing.SupportsInt) -> list[int]: ...
def process_resonance_ptr(
    mol_ptr: typing.SupportsInt, charge: typing.SupportsInt
) -> tuple[int, int]: ...
def smiles_token_ptr(mol_ptr: typing.SupportsInt) -> str: ...
