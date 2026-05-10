from __future__ import annotations

from typing import get_type_hints

from molgr import config as root_config
from molgr.fallback import utils as fallback_utils
from molgr.fallback.utils import config as fallback_config


def test_fallback_utils_reexports_config_from_single_module() -> None:
    assert fallback_utils.config is fallback_config
    assert fallback_utils.resolve_config is root_config.resolve_config
    assert fallback_utils.DEFAULT_MOLGR_CONFIG is fallback_config.DEFAULT_MOLGR_CONFIG
    assert fallback_utils.ResonanceTraversalScore is root_config.ResonanceTraversalScore


def test_config_string_options_use_literal_aliases() -> None:
    resonance_hints = get_type_hints(root_config.ResonanceConfig)
    interface_hints = get_type_hints(root_config.PythonInterfaceConfig)

    assert resonance_hints["traversal_score"] is root_config.ResonanceTraversalScore
    assert interface_hints["reconstruction_failure_policy"] is (
        root_config.ReconstructionFailurePolicy
    )


def test_config_no_longer_exposes_force_field_selection() -> None:
    config = root_config.make_default_config()

    assert all("force_field" not in name for name in root_config.__all__)
    assert not hasattr(config, "force_field")
