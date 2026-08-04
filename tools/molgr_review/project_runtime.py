"""Validate that review reconstruction uses the current project checkout."""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import molgr


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _cpp_source_files(repo_root: Path) -> Iterable[Path]:
    for path in (repo_root / "CMakeLists.txt", repo_root / "pyproject.toml"):
        if path.is_file():
            yield path
    for directory in (repo_root / "src" / "cpp", repo_root / "src" / "bindings"):
        if directory.is_dir():
            yield from (path for path in directory.rglob("*") if path.is_file())


def _git_runtime(repo_root: Path) -> tuple[str, bool]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=repo_root,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        return "", False
    return revision, dirty


def validate_project_runtime(repo_root: Path) -> dict[str, Any]:
    """Return runtime provenance or fail when the editable/C++ build is stale."""

    repo_root = repo_root.resolve()
    package_file = Path(str(molgr.__file__)).resolve()
    expected_package_dir = (repo_root / "src" / "molgr").resolve()
    if not _is_within(package_file, expected_package_dir):
        raise RuntimeError(
            "molgr is not imported from the current checkout: "
            f"{package_file}. Run `uv pip install -e . -v --no-build-isolation`."
        )

    core = importlib.import_module("molgr._core")
    core_file = Path(str(core.__file__)).resolve()
    if not core_file.is_file():
        raise RuntimeError(f"MolGR C++ extension is unavailable: {core_file}")

    sources = list(_cpp_source_files(repo_root))
    latest_source = max(sources, key=lambda path: path.stat().st_mtime) if sources else None
    if latest_source is not None and latest_source.stat().st_mtime > core_file.stat().st_mtime:
        raise RuntimeError(
            "MolGR C++ extension is older than the current source checkout: "
            f"{latest_source}. Rebuild with "
            "`uv pip install -e . -v --no-build-isolation` before reviewing."
        )

    revision, dirty = _git_runtime(repo_root)
    return {
        "python": sys.executable,
        "python_prefix": sys.prefix,
        "molgr_source": str(package_file),
        "cpp_extension": str(core_file),
        "latest_cpp_source": str(latest_source or ""),
        "git_revision": revision,
        "git_dirty": dirty,
    }


__all__ = ["validate_project_runtime"]
