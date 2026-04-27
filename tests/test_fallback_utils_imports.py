from __future__ import annotations

from molgr import config as root_config
from molgr.fallback import utils as fallback_utils
from molgr.fallback.utils import config as fallback_config


def test_fallback_utils_reexports_config_from_single_module() -> None:
    assert fallback_utils.config is fallback_config
    assert fallback_utils.resolve_config is root_config.resolve_config
    assert fallback_utils.DEFAULT_MOLGR_CONFIG is fallback_config.DEFAULT_MOLGR_CONFIG
