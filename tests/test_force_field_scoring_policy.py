from __future__ import annotations

import pytest
from openbabel import pybel

from molgr.fallback.utils import force_field as force_field_module
from molgr.fallback.utils.force_field import ForceFieldEvaluation


def _dummy_evaluation() -> ForceFieldEvaluation:
    return ForceFieldEvaluation(
        raw_energy=1.0,
        raw_unit="kj/mol",
        energy_kj_mol=1.0,
        atom_count=2,
        heavy_atom_count=2,
        contains_metals=False,
    )


def test_organic_force_field_evaluation_uses_fixed_uff_request(monkeypatch) -> None:
    calls = 0

    def fake_cached(context):
        nonlocal calls
        calls += 1
        return _dummy_evaluation()

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.organic_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.energy_kj_mol == pytest.approx(1.0)
    assert calls == 1


def test_force_field_evaluation_uses_fixed_uff_for_organic_input(monkeypatch) -> None:
    force_field_module.force_field_evaluation_cache_clear()
    setup_attempts: list[str] = []
    setup_atom_counts: list[int] = []

    class FakeForceField:
        def __init__(self, name: str) -> None:
            self.name = name

        def Setup(self, obmol) -> bool:
            setup_attempts.append(self.name)
            setup_atom_counts.append(int(obmol.NumAtoms()))
            return self.name == "uff"

        def Energy(self) -> float:
            return 2.0

        def GetUnit(self) -> str:
            return "kj/mol"

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(
        force_field_module.ob.OBForceField,
        "FindForceField",
        lambda name: FakeForceField(name),
    )

    evaluation = force_field_module.force_field_evaluation(pybel.readstring("smi", "CO"))

    assert setup_attempts == ["uff", "uff"]
    assert setup_atom_counts[0] == 0
    assert setup_atom_counts[1] > 0
    assert evaluation.energy_kj_mol == pytest.approx(2.0)


def test_organic_force_field_evaluation_rejects_metals(monkeypatch) -> None:
    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: True)

    with pytest.raises(ValueError, match="metal-free molecules"):
        force_field_module.organic_force_field_evaluation(pybel.readstring("smi", "CO"))


def test_combined_force_field_evaluation_uses_uff_request(monkeypatch) -> None:
    calls = 0

    def fake_force_field_evaluation(omol_or_state):
        nonlocal calls
        calls += 1
        return _dummy_evaluation()

    monkeypatch.setattr(force_field_module, "force_field_evaluation", fake_force_field_evaluation)

    evaluation = force_field_module.combined_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.energy_kj_mol == pytest.approx(1.0)
    assert calls == 1


def test_force_field_evaluation_uses_raw_uff_for_metal_input(
    monkeypatch,
) -> None:
    omol = pybel.readstring(
        "xyz",
        """2
LiO
Li 0.0 0.0 0.0
O 2.0 0.0 0.0
""",
    )
    li_atom = omol.OBMol.GetAtom(1)
    o_atom = omol.OBMol.GetAtom(2)
    li_atom.SetFormalCharge(1)
    o_atom.SetFormalCharge(-1)

    monkeypatch.setattr(
        force_field_module,
        "_force_field_evaluation_cached",
        lambda context: ForceFieldEvaluation(
            raw_energy=10.0,
            raw_unit="kj/mol",
            energy_kj_mol=10.0,
            atom_count=2,
            heavy_atom_count=2,
            contains_metals=True,
        ),
    )

    evaluation = force_field_module.force_field_evaluation(omol)

    assert evaluation.energy_kj_mol == pytest.approx(10.0)


def test_selection_force_field_evaluation_uses_fixed_uff_for_metal_free_input(monkeypatch) -> None:
    calls = 0

    def fake_cached(context):
        nonlocal calls
        calls += 1
        return _dummy_evaluation()

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.selection_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.energy_kj_mol == pytest.approx(1.0)
    assert calls == 1


def test_selection_force_field_evaluation_uses_fixed_uff_on_stripped_organic_part_for_metal_input(
    monkeypatch,
) -> None:
    calls = 0

    def fake_cached(context):
        nonlocal calls
        calls += 1
        return _dummy_evaluation()

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: True)
    monkeypatch.setattr(
        force_field_module,
        "_strip_metal_atoms",
        lambda omol: pybel.readstring("smi", "CO"),
    )
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.selection_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.energy_kj_mol == pytest.approx(1.0)
    assert calls == 1


def test_force_field_evaluation_cache_is_independent_of_config(monkeypatch) -> None:
    force_field_module.force_field_evaluation_cache_clear()
    setup_attempts = 0
    setup_atom_counts: list[int] = []

    class FakeForceField:
        def Setup(self, obmol) -> bool:
            nonlocal setup_attempts
            setup_attempts += 1
            setup_atom_counts.append(int(obmol.NumAtoms()))
            return True

        def Energy(self) -> float:
            return 2.0

        def GetUnit(self) -> str:
            return "kj/mol"

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(
        force_field_module.ob.OBForceField,
        "FindForceField",
        lambda name: FakeForceField(),
    )

    omol = pybel.readstring("smi", "CO")
    first = force_field_module.selection_force_field_energy(omol)
    second = force_field_module.selection_force_field_energy(omol)

    assert first == second == pytest.approx(2.0)
    assert setup_attempts == 2
    assert setup_atom_counts[0] == 0
    assert setup_atom_counts[1] > 0


def test_force_field_resets_shared_uff_when_exact_setup_key_changes() -> None:
    force_field_module.force_field_evaluation_cache_clear()
    neutral = pybel.readstring("smi", "CS(C)=O")
    neutral.addh()
    force_field_module.ob.OBBuilder().Build(neutral.OBMol)

    charged = pybel.Molecule(force_field_module.ob.OBMol(neutral.OBMol))
    charged.OBMol.GetAtom(2).SetFormalCharge(-2)

    assert force_field_module._build_openbabel_setup_key(
        neutral.OBMol
    ) == force_field_module._build_openbabel_setup_key(charged.OBMol)
    assert force_field_module._build_force_field_setup_key(
        neutral.OBMol
    ) != force_field_module._build_force_field_setup_key(charged.OBMol)

    neutral_energy = force_field_module.selection_force_field_energy(neutral)
    charged_energy = force_field_module.selection_force_field_energy(charged)

    assert neutral_energy != pytest.approx(charged_energy)
    assert charged_energy == pytest.approx(380.44249055506697)
