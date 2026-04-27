from __future__ import annotations

import shutil
from pathlib import Path


CORE_ROOT = Path("src/molgr/_core")
PIPELINE_ROOT = CORE_ROOT / "pipeline"

MODULE_TO_PACKAGE_INIT = ((CORE_ROOT / "pipeline.pyi", PIPELINE_ROOT / "__init__.pyi"),)

REQUIRED_PATHS = (
    CORE_ROOT / "__init__.pyi",
    CORE_ROOT / "utils.pyi",
    PIPELINE_ROOT / "__init__.pyi",
    PIPELINE_ROOT / "reconstruct_with_metals.pyi",
    PIPELINE_ROOT / "reconstruct_without_metals.pyi",
)

STALE_PATHS = (
    CORE_ROOT / "dev.pyi",
    CORE_ROOT / "dev",
    CORE_ROOT / "stages.pyi",
    CORE_ROOT / "stages",
    PIPELINE_ROOT / "resonance.pyi",
)


def _promote_module_stub(module_stub: Path, package_init: Path) -> None:
    if not module_stub.exists():
        return

    if package_init.exists():
        raise RuntimeError(
            "found both module-form and package-form stubs; clean up stubgen output first: "
            f"{module_stub} and {package_init}"
        )

    package_init.parent.mkdir(parents=True, exist_ok=True)
    module_stub.replace(package_init)


def _remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
        return
    if path.exists():
        path.unlink()


def _require_path(path: Path) -> None:
    if not path.exists():
        raise RuntimeError(f"missing expected stub: {path}")


def main() -> None:
    for module_stub, package_init in MODULE_TO_PACKAGE_INIT:
        _promote_module_stub(module_stub, package_init)

    for path in REQUIRED_PATHS:
        _require_path(path)

    for path in STALE_PATHS:
        _remove_path(path)


if __name__ == "__main__":
    main()
