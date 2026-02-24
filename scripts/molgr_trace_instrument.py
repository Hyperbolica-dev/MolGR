from __future__ import annotations

import functools
import importlib
import inspect
import os
import pkgutil
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from types import ModuleType
from typing import Any, Optional

from molgr import fallback
from molgr.fallback.pipeline import reconstruct_with_metals, reconstruct_without_metals
from molgr.fallback.utils import scoring


try:
    from scripts.molgr_trace_schema import TraceWriter
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from molgr_trace_schema import TraceWriter


_TRACE_WRAPPER_FLAG = "__molgr_trace_wrapper__"
_STAGE_FUNCTION_SKIPLIST = {
    ("molgr.fallback.stages.fresh", "assign_radical_dots"),
    ("molgr.fallback.stages.fresh", "assign_charge_radical_for_atom"),
}
_CALLER_MODULES_TO_HEAL = (
    "molgr.fallback.pipeline.reconstruct_without_metals",
    "molgr.fallback.pipeline.resonance",
    "molgr.fallback.pipeline.reconstruct_with_metals",
)


def _safe_smiles(value: Any) -> Optional[str]:
    if value is None:
        return None
    writer = getattr(value, "write", None)
    if not callable(writer):
        return None
    try:
        smiles_raw = writer("smi")
    except Exception:
        return None
    if not isinstance(smiles_raw, str):
        return None
    smiles_text = smiles_raw.strip()
    if not smiles_text:
        return None
    return smiles_text.split()[0]


def _extract_molecule(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, "OBMol") and callable(getattr(value, "write", None)):
        return value
    if isinstance(value, tuple) and value:
        return _extract_molecule(value[0])
    if isinstance(value, list) and value:
        return _extract_molecule(value[0])
    return None


def _extract_metrics(value: Any) -> Optional[dict[str, Any]]:
    mol = _extract_molecule(value)
    if mol is None:
        return None
    atoms = getattr(mol, "atoms", None)
    if atoms is None:
        return None
    try:
        atom_list = list(atoms)
    except Exception:
        return None

    atom_count = len(atom_list)
    charge = 0
    radical = 0
    bond_count: Optional[int] = None

    for atom in atom_list:
        obatom = getattr(atom, "OBAtom", None)
        if obatom is None:
            continue
        get_charge = getattr(obatom, "GetFormalCharge", None)
        if callable(get_charge):
            try:  # noqa: SIM105
                charge += int(get_charge())  # pyright: ignore[reportArgumentType]
            except Exception:
                pass
        get_spin = getattr(obatom, "GetSpinMultiplicity", None)
        if callable(get_spin):
            try:
                spin = int(get_spin())  # pyright: ignore[reportArgumentType]
                if spin > 1:
                    radical += spin - 1
            except Exception:
                pass

    obmol = getattr(mol, "OBMol", None)
    if obmol is not None:
        num_bonds = getattr(obmol, "NumBonds", None)
        if callable(num_bonds):
            try:
                bond_count = int(num_bonds())  # pyright: ignore[reportArgumentType]
            except Exception:
                bond_count = None

    return {
        "charge": charge,
        "radical": radical,
        "atom_count": atom_count,
        "bond_count": bond_count,
    }


def _error_text(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def _source_meta(fn: Callable[..., Any]) -> dict[str, Any]:
    try:
        target = inspect.unwrap(fn)
    except Exception:
        target = fn

    src_qualname = getattr(target, "__qualname__", None) or getattr(target, "__name__", None)
    src_module = getattr(target, "__module__", None)

    src_path: Optional[str] = None
    src_line: Optional[int] = None
    try:
        raw_path = inspect.getsourcefile(target) or inspect.getfile(target)
        if isinstance(raw_path, str) and raw_path.strip() != "":
            src_path = os.path.abspath(raw_path)
    except Exception:
        src_path = None

    if src_path is not None:
        try:
            _lines, start_line = inspect.getsourcelines(target)
            if isinstance(start_line, int) and start_line > 0:
                src_line = start_line
        except Exception:
            src_line = None

    meta: dict[str, Any] = {}
    if isinstance(src_path, str):
        meta["src_path"] = src_path
    if isinstance(src_line, int):
        meta["src_line"] = src_line
    if isinstance(src_qualname, str) and src_qualname.strip() != "":
        meta["src_qualname"] = src_qualname
    if isinstance(src_module, str) and src_module.strip() != "":
        meta["src_module"] = src_module
    return meta


def _patch_targets() -> list[tuple[Any, str, str]]:
    targets: list[tuple[Any, str, str]] = []
    exact = {
        "make_connections",
        "pre_clean",
        "fresh_omol_charge_radical",
        "omol_score",
        "break_deformed_ene",
        "break_one_bond",
        "get_radical_resonances",
        "process_resonance",
        "validate_omol",
        "xyz_to_omol_no_metal",
        "xyz2omol",
    }

    for name in dir(reconstruct_without_metals):
        if not (name in exact or name.startswith("eliminate_") or name.startswith("clean_")):
            continue
        value = getattr(reconstruct_without_metals, name, None)
        module_name = getattr(value, "__module__", "")
        if isinstance(module_name, str) and module_name.startswith("molgr.fallback.stages."):
            continue
        if callable(value):
            targets.append((reconstruct_without_metals, name, "reconstruct_without_metals"))

    for name in ("xyz2omol", "xyz_to_omol_no_metal"):
        value = getattr(reconstruct_with_metals, name, None)
        if callable(value):
            targets.append((reconstruct_with_metals, name, "reconstruct_with_metals"))

    if callable(getattr(scoring, "omol_score", None)):
        targets.append((scoring, "omol_score", "scoring"))

    if callable(getattr(fallback, "xyz2omol", None)):
        targets.append((fallback, "xyz2omol", "fallback"))

    return targets


def _discover_stage_modules() -> list[ModuleType]:
    stages_pkg = importlib.import_module("molgr.fallback.stages")
    stage_modules: list[ModuleType] = []
    package_paths = getattr(stages_pkg, "__path__", None)
    if package_paths is None:
        return stage_modules
    for module_info in pkgutil.walk_packages(package_paths, f"{stages_pkg.__name__}."):
        stage_modules.append(importlib.import_module(module_info.name))
    return stage_modules


def _iter_stage_defined_functions(
    stage_module: ModuleType,
) -> Iterator[tuple[str, Callable[..., Any]]]:
    for name, value in vars(stage_module).items():
        if (stage_module.__name__, name) in _STAGE_FUNCTION_SKIPLIST:
            continue
        if not callable(value):
            continue
        if getattr(value, "__module__", None) != stage_module.__name__:
            continue
        yield name, value


def _heal_stage_callsites(
    *,
    originals: list[tuple[Any, str, Any]],
    wrapped_by_original: dict[Callable[..., Any], Callable[..., Any]],
) -> None:
    for module_name in _CALLER_MODULES_TO_HEAL:
        try:
            caller_module = importlib.import_module(module_name)
        except Exception:
            continue
        for attr_name, attr_value in vars(caller_module).items():
            try:
                replacement = wrapped_by_original.get(attr_value)
            except TypeError:
                continue
            if replacement is None:
                continue
            if attr_value is replacement:
                continue
            originals.append((caller_module, attr_name, attr_value))
            setattr(caller_module, attr_name, replacement)


@contextmanager
def trace_monkeypatch(tracer: TraceWriter) -> Iterator[None]:
    originals: list[tuple[Any, str, Any]] = []
    wrapped_by_original: dict[Callable[..., Any], Callable[..., Any]] = {}

    def wrap(fn: Callable[..., Any], *, op: str, phase: str) -> Callable[..., Any]:
        existing = wrapped_by_original.get(fn)
        if existing is not None:
            return existing
        if getattr(fn, _TRACE_WRAPPER_FLAG, False):
            return fn

        source_meta = _source_meta(fn)

        @functools.wraps(fn)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            input_mol = _extract_molecule(args[0]) if args else None
            smiles_in = _safe_smiles(input_mol)
            metrics_in = _extract_metrics(input_mol)
            start_id = tracer.span_start(
                op=op,
                phase=phase,
                smiles_in=smiles_in,
                meta=dict(source_meta),
                metrics=metrics_in,
            )
            status = "ok"
            error_text: Optional[str] = None
            result: Any = None
            try:
                result = fn(*args, **kwargs)
                return result
            except Exception as exc:
                status = "error"
                error_text = _error_text(exc)
                raise
            finally:
                output_mol = _extract_molecule(result)
                if output_mol is None:
                    output_mol = input_mol
                smiles_out = _safe_smiles(output_mol)
                metrics_out = _extract_metrics(output_mol)
                changed = (smiles_in != smiles_out) and bool((smiles_in or "") + (smiles_out or ""))
                meta: dict[str, Any] = {"changed": changed}
                if "omol_score" in op.lower() and isinstance(result, (int, float)):
                    meta["omol_score"] = float(result)
                tracer.span_end(
                    start_event_id=start_id,
                    op=op,
                    status=status,
                    error=error_text,
                    smiles_out=smiles_out,
                    meta=meta,
                    metrics=metrics_out,
                )

        setattr(wrapped, _TRACE_WRAPPER_FLAG, True)
        cache_clear = getattr(fn, "cache_clear", None)
        if callable(cache_clear):
            wrapped.cache_clear = cache_clear  # type: ignore[assignment]
        cache_info = getattr(fn, "cache_info", None)
        if callable(cache_info):
            wrapped.cache_info = cache_info  # type: ignore[assignment]
        wrapped_by_original[fn] = wrapped
        return wrapped

    try:
        for module_obj, name, phase in _patch_targets():
            original = getattr(module_obj, name)
            wrapped = wrap(original, op=name, phase=phase)
            if original is wrapped:
                continue
            originals.append((module_obj, name, original))
            setattr(module_obj, name, wrapped)

        for stage_module in _discover_stage_modules():
            for name, original in _iter_stage_defined_functions(stage_module):
                op = f"{stage_module.__name__}.{original.__qualname__}"
                wrapped = wrap(original, op=op, phase="stage")
                if original is wrapped:
                    continue
                originals.append((stage_module, name, original))
                setattr(stage_module, name, wrapped)

        _heal_stage_callsites(originals=originals, wrapped_by_original=wrapped_by_original)
        yield
    finally:
        for module_obj, name, original in reversed(originals):
            setattr(module_obj, name, original)
