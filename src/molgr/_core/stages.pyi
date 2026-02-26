"""
Stage-level parity helpers.
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "preprocess",
    "fresh",
    "eliminate",
    "clean",
    "break_bond",
]

class _PreprocessModule:
    def make_connections_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        factor: float = 1.4,
    ) -> None: ...
    def pre_clean_ptr(self, mol_ptr: typing.SupportsInt) -> None: ...
    def validate_omol_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        total_charge: typing.SupportsInt,
        total_radical: typing.SupportsInt,
    ) -> bool: ...

class _FreshModule:
    def assign_radical_dots_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        atom_idx: typing.SupportsInt,
    ) -> int: ...
    def assign_charge_radical_for_atom_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        atom_idx: typing.SupportsInt,
    ) -> None: ...
    def fresh_omol_charge_radical_ptr(self, mol_ptr: typing.SupportsInt) -> None: ...

class _EliminateModule:
    def eliminate_1_3_dipole_ptr(
        self, mol_ptr: typing.SupportsInt, given_charge: typing.SupportsInt
    ) -> int: ...
    def eliminate_positive_charges_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_negative_charges_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_nnn_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_high_positive_charge_atoms_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_cn_in_doubt_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_carboxyl_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_carbene_neighbor_heteroatom_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...
    def eliminate_charge_spliting_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt,
    ) -> int: ...

class _CleanModule:
    def clean_neighbor_radicals_ptr(self, mol_ptr: typing.SupportsInt) -> None: ...
    def clean_carbene_neighbor_unsaturated_ptr(self, mol_ptr: typing.SupportsInt) -> None: ...
    def clean_resonances_ptr(self, mol_ptr: typing.SupportsInt) -> None: ...

class _BreakBondModule:
    def break_deformed_ene_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt = 0,
        given_radical: typing.SupportsInt = 0,
        tolerance: float = 5.0,
    ) -> None: ...
    def break_one_bond_ptr(
        self,
        mol_ptr: typing.SupportsInt,
        given_charge: typing.SupportsInt = 0,
        given_radical: typing.SupportsInt = 0,
    ) -> int: ...

preprocess: _PreprocessModule
fresh: _FreshModule
eliminate: _EliminateModule
clean: _CleanModule
break_bond: _BreakBondModule
