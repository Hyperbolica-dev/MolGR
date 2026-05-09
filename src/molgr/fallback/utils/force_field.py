from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, Mapping, Optional, Sequence, Tuple, Union, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import ForceFieldConfig, MolGRConfig, resolve_config
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.tools import typed_lru_cache


if TYPE_CHECKING:
    from molgr.fallback.state import ReconstructionState


ForceFieldChoice = Literal["auto", "uff"]
ForceFieldAtomKey = Tuple[int, int, int, int, int, int, bool]
ForceFieldBondKey = Tuple[int, int, int, bool]
ForceFieldScoreKey = Tuple[Tuple[ForceFieldAtomKey, ...], Tuple[ForceFieldBondKey, ...]]
MetalStateKey = Tuple[Tuple[int, int, int, int, int, int, int], ...]
BondOrderOverrides = Mapping[Tuple[int, int], int]
ForceFieldInput = Union[pybel.Molecule, "ReconstructionState", "OmolForceFieldContext"]
_DEFAULT_FORCE_FIELD_CACHE_MAXSIZE = 4096
_COORDINATE_SCALE = 1_000_000
_FORCE_FIELD_UNIT_TO_KJ_MOL = {
    "kj/mol": 1.0,
    "kcal/mol": 4.184,
}


def _quantized_coordinate(value: float) -> int:
    return int(round(float(value) * _COORDINATE_SCALE))


def _build_score_key(omol: pybel.Molecule) -> ForceFieldScoreKey:
    obmol = cast(ob.OBMol, omol.OBMol)
    atom_keys: list[ForceFieldAtomKey] = []
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        atom_keys.append(
            (
                int(obatom.GetAtomicNum()),
                int(obatom.GetFormalCharge()),
                int(obatom.GetSpinMultiplicity()),
                _quantized_coordinate(obatom.GetX()),
                _quantized_coordinate(obatom.GetY()),
                _quantized_coordinate(obatom.GetZ()),
                bool(obatom.IsAromatic()),
            )
        )

    bond_keys: list[ForceFieldBondKey] = []
    for bond in ob.OBMolBondIter(obmol):
        obbond = cast(ob.OBBond, bond)
        begin_idx = int(cast(ob.OBAtom, obbond.GetBeginAtom()).GetIdx())
        end_idx = int(cast(ob.OBAtom, obbond.GetEndAtom()).GetIdx())
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bond_keys.append(
            (begin_idx, end_idx, int(obbond.GetBondOrder()), bool(obbond.IsAromatic()))
        )
    return tuple(atom_keys), tuple(bond_keys)


def build_force_field_score_key(omol_or_state: ForceFieldInput) -> ForceFieldScoreKey:
    return OmolForceFieldContext.from_input(omol_or_state).score_key


def build_metal_state_key(metal_states: Sequence[MetalAtomPosition]) -> MetalStateKey:
    return tuple(
        (
            metal_state.idx,
            metal_state.element_idx,
            metal_state.valence,
            metal_state.radical_num,
            _quantized_coordinate(metal_state.position_x),
            _quantized_coordinate(metal_state.position_y),
            _quantized_coordinate(metal_state.position_z),
        )
        for metal_state in metal_states
    )


class OmolForceFieldContext:
    """Hashable context wrapper around one omol instance."""

    __slots__ = ("omol", "_score_key", "_score_key_hash", "_contains_metals")

    def __init__(
        self,
        omol: pybel.Molecule,
        *,
        score_key: ForceFieldScoreKey | None = None,
        contains_metals: bool | None = None,
    ) -> None:
        self.omol = omol
        self._score_key = score_key
        self._score_key_hash: int | None = None
        self._contains_metals = contains_metals

    @classmethod
    def from_input(cls, omol_or_state: ForceFieldInput) -> OmolForceFieldContext:
        if isinstance(omol_or_state, OmolForceFieldContext):
            return omol_or_state
        if isinstance(omol_or_state, pybel.Molecule):
            return cls(omol_or_state)

        from molgr.fallback.state import ReconstructionState

        if isinstance(omol_or_state, ReconstructionState):
            cached_key = omol_or_state.metadata.get("force_field_score_key")
            cached_score_key = cached_key if isinstance(cached_key, tuple) else None
            return cls(omol_or_state.omol, score_key=cached_score_key)

        raise TypeError(f"Unsupported force-field input: {type(omol_or_state)!r}")

    @property
    def score_key(self) -> ForceFieldScoreKey:
        if self._score_key is None:
            self._score_key = _build_score_key(self.omol)
        return self._score_key

    @property
    def contains_metals(self) -> bool:
        if self._contains_metals is None:
            self._contains_metals = _contains_metal_atoms(self.omol)
        return self._contains_metals

    def __hash__(self) -> int:
        if self._score_key_hash is None:
            self._score_key_hash = hash(self.score_key)
        return self._score_key_hash

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OmolForceFieldContext):
            return NotImplemented
        return self.score_key == other.score_key


@dataclass(frozen=True)
class ForceFieldEvaluation:
    requested_force_field: str
    resolved_force_field: str
    selection_reason: str
    raw_energy: float
    raw_unit: str
    energy_kj_mol: float
    atom_count: int
    heavy_atom_count: int
    contains_metals: bool


def _contains_metal_atoms(omol: pybel.Molecule) -> bool:
    obmol = cast(ob.OBMol, omol.OBMol)
    return any(cast(ob.OBAtom, atom).IsMetal() for atom in ob.OBMolAtomIter(obmol))


def _strip_metal_atoms(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = ob.OBMol(cast(ob.OBMol, omol.OBMol))
    metal_indices = [
        cast(int, cast(ob.OBAtom, atom).GetIdx())
        for atom in ob.OBMolAtomIter(obmol)
        if cast(ob.OBAtom, atom).IsMetal()
    ]
    if not metal_indices:
        return pybel.Molecule(obmol)
    obmol.BeginModify()
    try:
        for atom_idx in reversed(metal_indices):
            atom = cast(Optional[ob.OBAtom], obmol.GetAtom(atom_idx))
            if atom is not None:
                obmol.DeleteAtom(atom)
    finally:
        obmol.EndModify()
    return pybel.Molecule(obmol)


def _count_heavy_atoms(obmol: ob.OBMol) -> int:
    return sum(1 for atom in ob.OBMolAtomIter(obmol) if cast(ob.OBAtom, atom).GetAtomicNum() != 1)


def _resolve_force_field_config(config: MolGRConfig | None = None) -> ForceFieldConfig:
    return resolve_config(config).force_field


def _normalized_force_field_order(force_fields: Sequence[str]) -> Tuple[ForceFieldChoice, ...]:
    normalized_order: list[ForceFieldChoice] = []
    for force_field in force_fields:
        normalized_force_field = _normalize_force_field_choice(force_field)
        if normalized_force_field in normalized_order:
            continue
        normalized_order.append(normalized_force_field)
    if not normalized_order:
        raise ValueError("Auto force-field candidate order cannot be empty.")
    return tuple(normalized_order)


def _normalize_force_field_choice(force_field: str) -> ForceFieldChoice:
    normalized = force_field.strip().lower()
    if normalized not in {"auto", "uff"}:
        raise ValueError("Unsupported force field choice. Expected one of 'auto' or 'uff'.")
    return cast(ForceFieldChoice, normalized)


def _force_field_energy_to_kj_mol(raw_energy: float, raw_unit: str) -> float:
    normalized_unit = raw_unit.strip().lower()
    factor = _FORCE_FIELD_UNIT_TO_KJ_MOL.get(normalized_unit)
    if factor is None:
        raise ValueError(f"Unsupported force-field energy unit: {raw_unit!r}")
    return raw_energy * factor


def _build_force_field_evaluation(
    context: OmolForceFieldContext,
    requested_force_field: str,
    *,
    auto_force_fields_metal_free: Sequence[str],
    auto_force_fields_with_metals: Sequence[str],
) -> ForceFieldEvaluation:
    normalized_force_field = _normalize_force_field_choice(requested_force_field)
    contains_metals = context.contains_metals
    working_obmol = ob.OBMol(cast(ob.OBMol, context.omol.OBMol))
    working_obmol.SetAromaticPerceived(False)
    atom_count = int(working_obmol.NumAtoms())
    heavy_atom_count = _count_heavy_atoms(working_obmol)

    candidates: list[Tuple[str, str]] = []
    if normalized_force_field == "auto":
        auto_force_fields = (
            _normalized_force_field_order(auto_force_fields_with_metals)
            if contains_metals
            else _normalized_force_field_order(auto_force_fields_metal_free)
        )
        for idx, candidate_name in enumerate(auto_force_fields):
            if idx == 0:
                selection_reason = f"auto_prefer_{candidate_name}"
            else:
                selection_reason = f"auto_fallback_to_{candidate_name}"
            candidates.append((candidate_name, selection_reason))
    else:
        candidates.append((normalized_force_field, "explicit_request"))

    failures: list[str] = []
    for candidate_force_field_name, selection_reason in candidates:
        force_field = ob.OBForceField.FindForceField(candidate_force_field_name)
        if not force_field:
            failures.append(f"{candidate_force_field_name}: unavailable")
            continue
        candidate_obmol = ob.OBMol(working_obmol)
        if not bool(force_field.Setup(candidate_obmol)):
            failures.append(f"{candidate_force_field_name}: setup_failed")
            continue
        raw_energy = float(force_field.Energy())
        raw_unit = str(force_field.GetUnit())
        energy_kj_mol = _force_field_energy_to_kj_mol(raw_energy, raw_unit)
        return ForceFieldEvaluation(
            requested_force_field=normalized_force_field,
            resolved_force_field=candidate_force_field_name,
            selection_reason=selection_reason,
            raw_energy=raw_energy,
            raw_unit=raw_unit,
            energy_kj_mol=energy_kj_mol,
            atom_count=atom_count,
            heavy_atom_count=heavy_atom_count,
            contains_metals=contains_metals,
        )

    raise ValueError(
        f"Could not evaluate force-field energy with {normalized_force_field!r}: {failures!r}"
    )


def _build_force_field_evaluation_cached_impl(
    context: OmolForceFieldContext,
    requested_force_field: ForceFieldChoice,
    config: MolGRConfig,
) -> ForceFieldEvaluation:
    force_field_config = _resolve_force_field_config(config)
    return _build_force_field_evaluation(
        context,
        requested_force_field,
        auto_force_fields_metal_free=_normalized_force_field_order(
            force_field_config.auto_force_fields_metal_free
        ),
        auto_force_fields_with_metals=_normalized_force_field_order(
            force_field_config.auto_force_fields_with_metals
        ),
    )


_force_field_evaluation_cached = typed_lru_cache(
    maxsize=_DEFAULT_FORCE_FIELD_CACHE_MAXSIZE,
    typed=True,
)(_build_force_field_evaluation_cached_impl)


def force_field_evaluation_cache_info() -> Tuple[int, int, int]:
    info = _force_field_evaluation_cached.cache_info()
    return info.hits, info.misses, info.currsize


def force_field_evaluation_cache_clear() -> None:
    _force_field_evaluation_cached.cache_clear()


def _force_field_evaluation_from_context(
    context: OmolForceFieldContext,
    requested_force_field: ForceFieldChoice,
    *,
    config: MolGRConfig | None = None,
) -> ForceFieldEvaluation:
    resolved_config = resolve_config(config)
    return _force_field_evaluation_cached(
        context,
        requested_force_field,
        resolved_config,
    )


def force_field_evaluation(
    omol_or_state: ForceFieldInput,
    *,
    force_field: ForceFieldChoice = "auto",
    config: MolGRConfig | None = None,
) -> ForceFieldEvaluation:
    context = OmolForceFieldContext.from_input(omol_or_state)
    requested_force_field = _normalize_force_field_choice(force_field)
    return _force_field_evaluation_from_context(
        context,
        requested_force_field,
        config=config,
    )


def force_field_energy(
    omol_or_state: ForceFieldInput,
    *,
    force_field: ForceFieldChoice = "auto",
    config: MolGRConfig | None = None,
) -> float:
    return force_field_evaluation(
        omol_or_state,
        force_field=force_field,
        config=config,
    ).energy_kj_mol


def organic_force_field_evaluation(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> ForceFieldEvaluation:
    context = OmolForceFieldContext.from_input(omol_or_state)
    if context.contains_metals:
        raise ValueError(
            "Organic force-field evaluation only supports metal-free molecules; "
            "use selection_force_field_evaluation(...) for the default metal-aware policy or "
            "combined_force_field_evaluation(...) for raw full-molecule UFF diagnostics."
        )
    requested_force_field = _normalize_force_field_choice(
        _resolve_force_field_config(config).organic_force_field
    )
    return _force_field_evaluation_from_context(
        context,
        requested_force_field,
        config=config,
    )


def organic_force_field_energy(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> float:
    return organic_force_field_evaluation(omol_or_state, config=config).energy_kj_mol


def combined_force_field_evaluation(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> ForceFieldEvaluation:
    requested_force_field = _normalize_force_field_choice(
        _resolve_force_field_config(config).combined_force_field
    )
    return force_field_evaluation(
        omol_or_state,
        force_field=requested_force_field,
        config=config,
    )


def combined_force_field_energy(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> float:
    return combined_force_field_evaluation(omol_or_state, config=config).energy_kj_mol


def selection_force_field_evaluation(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> ForceFieldEvaluation:
    context = OmolForceFieldContext.from_input(omol_or_state)
    requested_force_field = _normalize_force_field_choice(
        _resolve_force_field_config(config).selection_force_field
    )
    if context.contains_metals:
        organic_context = OmolForceFieldContext(
            _strip_metal_atoms(context.omol), contains_metals=False
        )
        return _force_field_evaluation_from_context(
            organic_context,
            requested_force_field,
            config=config,
        )
    return _force_field_evaluation_from_context(
        context,
        requested_force_field,
        config=config,
    )


def selection_force_field_energy(
    omol_or_state: ForceFieldInput,
    *,
    config: MolGRConfig | None = None,
) -> float:
    return selection_force_field_evaluation(omol_or_state, config=config).energy_kj_mol


__all__ = [
    "BondOrderOverrides",
    "ForceFieldChoice",
    "ForceFieldEvaluation",
    "ForceFieldInput",
    "ForceFieldScoreKey",
    "MetalStateKey",
    "OmolForceFieldContext",
    "build_force_field_score_key",
    "build_metal_state_key",
    "combined_force_field_energy",
    "combined_force_field_evaluation",
    "force_field_energy",
    "force_field_evaluation",
    "force_field_evaluation_cache_clear",
    "force_field_evaluation_cache_info",
    "organic_force_field_energy",
    "organic_force_field_evaluation",
    "selection_force_field_energy",
    "selection_force_field_evaluation",
]
