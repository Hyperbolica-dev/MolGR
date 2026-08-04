from __future__ import annotations

from typing import get_type_hints

from molgr import config as root_config
from molgr.fallback import utils as fallback_utils


def test_fallback_utils_reexports_config_from_single_module() -> None:
    assert fallback_utils.CONFIG is root_config.CONFIG
    assert fallback_utils.ResonanceTraversalScore is root_config.ResonanceTraversalScore


def test_config_string_options_use_literal_aliases() -> None:
    resonance_hints = get_type_hints(root_config.ResonanceConfig)
    interface_hints = get_type_hints(root_config.PythonInterfaceConfig)

    assert resonance_hints["traversal_score"] is root_config.ResonanceTraversalScore
    assert interface_hints["reconstruction_failure_policy"] is (
        root_config.ReconstructionFailurePolicy
    )


def test_config_no_longer_exposes_force_field_selection() -> None:
    config = root_config.MolGRConfig()

    assert all("force_field" not in name for name in root_config.__all__)
    assert not hasattr(config, "force_field")
