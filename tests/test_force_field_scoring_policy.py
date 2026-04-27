from __future__ import annotations

from dataclasses import replace

import pytest
from openbabel import pybel

from molgr.config import get_config, make_default_config, set_config
from molgr.fallback.utils import force_field as force_field_module
from molgr.fallback.utils.force_field import ForceFieldEvaluation


def _dummy_evaluation(requested: str, resolved: str) -> ForceFieldEvaluation:
    return ForceFieldEvaluation(
        requested_force_field=requested,
        resolved_force_field=resolved,
        selection_reason="test",
        raw_energy=1.0,
        raw_unit="kj/mol",
        energy_kj_mol=1.0,
        atom_count=2,
        heavy_atom_count=2,
        contains_metals=False,
    )


def test_organic_force_field_evaluation_uses_auto_request(monkeypatch) -> None:
    calls: list[str] = []

    def fake_cached(context, requested_force_field, config):
        calls.append(requested_force_field)
        return _dummy_evaluation(requested_force_field, "uff")

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.organic_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.requested_force_field == "auto"
    assert evaluation.resolved_force_field == "uff"
    assert calls == ["auto"]


def test_auto_force_field_uses_uff_for_organic_input(monkeypatch) -> None:
    setup_attempts: list[str] = []

    class FakeForceField:
        def __init__(self, name: str) -> None:
            self.name = name

        def Setup(self, _obmol) -> bool:
            setup_attempts.append(self.name)
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

    assert setup_attempts == ["uff"]
    assert evaluation.requested_force_field == "auto"
    assert evaluation.resolved_force_field == "uff"
    assert evaluation.selection_reason == "auto_prefer_uff"


def test_organic_force_field_evaluation_rejects_metals(monkeypatch) -> None:
    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: True)

    with pytest.raises(ValueError, match="metal-free molecules"):
        force_field_module.organic_force_field_evaluation(pybel.readstring("smi", "CO"))


def test_combined_force_field_evaluation_uses_uff_request(monkeypatch) -> None:
    calls: list[str] = []

    def fake_force_field_evaluation(omol_or_state, *, force_field="auto", config=None):
        calls.append(force_field)
        return _dummy_evaluation(force_field, "uff")

    monkeypatch.setattr(force_field_module, "force_field_evaluation", fake_force_field_evaluation)

    evaluation = force_field_module.combined_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.requested_force_field == "uff"
    assert evaluation.resolved_force_field == "uff"
    assert calls == ["uff"]


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
        lambda context, requested_force_field, config: ForceFieldEvaluation(
            requested_force_field=requested_force_field,
            resolved_force_field="uff",
            selection_reason="test",
            raw_energy=10.0,
            raw_unit="kj/mol",
            energy_kj_mol=10.0,
            atom_count=2,
            heavy_atom_count=2,
            contains_metals=True,
        ),
    )

    evaluation = force_field_module.force_field_evaluation(omol, force_field="uff")

    assert evaluation.requested_force_field == "uff"
    assert evaluation.resolved_force_field == "uff"
    assert evaluation.energy_kj_mol == pytest.approx(10.0)


def test_selection_force_field_evaluation_uses_auto_for_metal_free_input(monkeypatch) -> None:
    calls: list[str] = []

    def fake_cached(context, requested_force_field, config):
        calls.append(requested_force_field)
        return _dummy_evaluation(requested_force_field, "uff")

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: False)
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.selection_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.requested_force_field == "auto"
    assert calls == ["auto"]


def test_selection_force_field_evaluation_uses_auto_on_stripped_organic_part_for_metal_input(
    monkeypatch,
) -> None:
    calls: list[str] = []

    def fake_cached(context, requested_force_field, config):
        calls.append(requested_force_field)
        return _dummy_evaluation(requested_force_field, "uff")

    monkeypatch.setattr(force_field_module, "_contains_metal_atoms", lambda omol: True)
    monkeypatch.setattr(
        force_field_module,
        "_strip_metal_atoms",
        lambda omol: pybel.readstring("smi", "CO"),
    )
    monkeypatch.setattr(force_field_module, "_force_field_evaluation_cached", fake_cached)

    evaluation = force_field_module.selection_force_field_evaluation(pybel.readstring("smi", "CO"))

    assert evaluation.requested_force_field == "auto"
    assert evaluation.resolved_force_field == "uff"
    assert calls == ["auto"]


def test_force_field_config_rejects_mmff94_candidate_order() -> None:
    original_config = get_config()
    base_config = make_default_config()
    invalid_config = replace(
        base_config,
        force_field=replace(
            base_config.force_field,
            auto_force_fields_metal_free=("mmff94",),
        ),
    )

    try:
        set_config(invalid_config)
        with pytest.raises(ValueError, match="Expected one of 'auto' or 'uff'"):
            force_field_module.force_field_evaluation(pybel.readstring("smi", "CO"))
    finally:
        set_config(original_config)


def test_force_field_cache_uses_config_in_key_and_execution(monkeypatch) -> None:
    base_config = make_default_config()
    config_a = replace(
        base_config,
        force_field=replace(
            base_config.force_field,
            auto_force_fields_metal_free=("uff",),
        ),
    )
    config_b = replace(
        base_config,
        force_field=replace(
            base_config.force_field,
            auto_force_fields_metal_free=("uff", "uff"),
        ),
    )
    seen_configs = []

    def fake_resolve_force_field_config(config=None):
        seen_configs.append(config)
        assert config is not None
        return config.force_field

    force_field_module.force_field_evaluation_cache_clear()
    monkeypatch.setattr(
        force_field_module,
        "_resolve_force_field_config",
        fake_resolve_force_field_config,
    )

    omol = pybel.readstring("smi", "CO")
    first = force_field_module.force_field_evaluation(omol, config=config_a)
    second = force_field_module.force_field_evaluation(omol, config=config_a)
    third = force_field_module.force_field_evaluation(omol, config=config_b)

    assert first.energy_kj_mol == second.energy_kj_mol
    assert first.resolved_force_field == second.resolved_force_field
    assert third.resolved_force_field == "uff"
    assert seen_configs == [config_a, config_b]
