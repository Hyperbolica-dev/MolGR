from __future__ import annotations

import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = APP_DIR.parents[1]
sys.path.insert(0, str(APP_DIR))

from project_runtime import validate_project_runtime  # noqa: E402


def test_review_runtime_uses_current_checkout_and_built_cpp_extension() -> None:
    runtime = validate_project_runtime(REPO_ROOT)

    assert Path(runtime["molgr_source"]).resolve().relative_to(REPO_ROOT / "src" / "molgr")
    assert Path(runtime["cpp_extension"]).is_file()
    assert runtime["python_prefix"] == sys.prefix
