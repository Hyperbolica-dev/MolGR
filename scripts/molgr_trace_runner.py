from __future__ import annotations

from pathlib import Path
from typing import Any

from rdkit import Chem

from molgr import fallback
from molgr.fallback.pipeline import reconstruct_without_metals, resonance
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


try:
    from scripts.molgr_debug_html import render_trace_report
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from molgr_debug_html import render_trace_report

try:
    from scripts.molgr_trace_instrument import trace_monkeypatch
    from scripts.molgr_trace_schema import TraceWriter
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from molgr_trace_instrument import trace_monkeypatch
    from molgr_trace_schema import TraceWriter


def _error_text(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def _clear_caches() -> None:
    clear_no_metal = getattr(reconstruct_without_metals.xyz_to_omol_no_metal, "cache_clear", None)
    if callable(clear_no_metal):
        clear_no_metal()
    clear_resonance = getattr(resonance.process_resonance, "cache_clear", None)
    if callable(clear_resonance):
        clear_resonance()


def _pybel_smiles(omol: Any) -> str | None:
    if omol is None:
        return None
    writer = getattr(omol, "write", None)
    if not callable(writer):
        return None
    try:
        smiles_raw = writer("smi")
    except Exception:
        return None
    if not isinstance(smiles_raw, str):
        return None
    smiles_text = smiles_raw.strip()
    if smiles_text == "":
        return None
    return smiles_text.split()[0]


def _run_equivalence_check(
    *,
    tracer: TraceWriter,
    case_idx: int,
    ground_truth_rdmol: Chem.Mol,
    ground_truth_smiles: str,
    predicted_omol: Any,
    use_chirality: bool,
    max_resonance: int,
) -> tuple[bool, str | None]:
    start_id = tracer.span_start(
        op="equivalence.check_equivalence",
        phase="equivalence",
        smiles_in=ground_truth_smiles,
        meta={
            "case_idx": case_idx,
            "equivalent": False,
            "method": None,
            "use_chirality": use_chirality,
            "max_resonance": max_resonance,
        },
        metrics=None,
    )

    status = "ok"
    error_text = None
    smiles_out = _pybel_smiles(predicted_omol)
    meta = {
        "case_idx": case_idx,
        "equivalent": False,
        "method": None,
        "use_chirality": use_chirality,
        "max_resonance": max_resonance,
    }

    try:
        if predicted_omol is None:
            raise ValueError("fallback.xyz2omol returned None")

        predicted_rdmol = pybel_to_rdmol(predicted_omol)
        smiles_out = Chem.MolToSmiles(
            predicted_rdmol,
            canonical=True,
            isomericSmiles=use_chirality,
        )

        equivalent, info = check_equivalence(
            ground_truth_rdmol,
            predicted_rdmol,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
        )
        meta["equivalent"] = bool(equivalent)
        method = info.method
        meta["method"] = method.value if method is not None else None

        if not equivalent:
            status = "error"
            error_text = info.reason
    except Exception as exc:
        status = "error"
        error_text = _error_text(exc)
    finally:
        tracer.span_end(
            start_event_id=start_id,
            op="equivalence.check_equivalence",
            status=status,
            error=error_text,
            smiles_out=smiles_out,
            meta=meta,
            metrics=None,
        )

    return status == "ok", error_text


def _run_case(
    case: dict[str, Any],
    case_dir: Path,
    *,
    use_chirality: bool,
    max_resonance: int,
) -> bool:
    trace_path = case_dir / "trace.jsonl"
    tracer = TraceWriter(trace_path)
    case_idx = int(case.get("case_idx", 0))
    input_smiles = case.get("input_smiles")
    smiles_text = input_smiles if isinstance(input_smiles, str) else None

    root_start_id = tracer.span_start(
        op="case",
        phase="backtest",
        smiles_in=smiles_text,
        meta={"case_idx": case_idx, "input_smiles": smiles_text},
        metrics=None,
    )
    status = "ok"
    error_text = None
    try:
        with trace_monkeypatch(tracer):
            _clear_caches()
            provider_error = case.get("provider_error")
            if isinstance(provider_error, str) and provider_error.strip() != "":
                raise ValueError(provider_error)

            xyz_block = case.get("xyz_block")
            total_charge = case.get("total_charge")
            total_radical_electrons = case.get("total_radical_electrons")

            if not isinstance(xyz_block, str) or xyz_block.strip() == "":
                raise ValueError("missing xyz_block in case")
            if not isinstance(total_charge, int):
                raise ValueError("missing total_charge in case")
            if not isinstance(total_radical_electrons, int):
                raise ValueError("missing total_radical_electrons in case")

            predicted_omol = fallback.xyz2omol(
                xyz_block,
                total_charge=total_charge,
                total_radical_electrons=total_radical_electrons,
            )

            ground_truth_rdmol = case.get("ground_truth_rdmol")
            if isinstance(ground_truth_rdmol, Chem.Mol):
                ground_truth_smiles = case.get("ground_truth_smiles")
                smiles_in = (
                    ground_truth_smiles
                    if isinstance(ground_truth_smiles, str) and ground_truth_smiles.strip() != ""
                    else Chem.MolToSmiles(
                        ground_truth_rdmol,
                        canonical=True,
                        isomericSmiles=use_chirality,
                    )
                )
                eq_ok, eq_error = _run_equivalence_check(
                    tracer=tracer,
                    case_idx=case_idx,
                    ground_truth_rdmol=ground_truth_rdmol,
                    ground_truth_smiles=smiles_in,
                    predicted_omol=predicted_omol,
                    use_chirality=use_chirality,
                    max_resonance=max_resonance,
                )
                if not eq_ok:
                    raise ValueError(eq_error or "equivalence check failed")
    except Exception as exc:
        status = "error"
        error_text = _error_text(exc)
    finally:
        tracer.span_end(
            start_event_id=root_start_id,
            op="case",
            status=status,
            error=error_text,
            smiles_out=None,
            meta={"case_idx": case_idx, "input_smiles": smiles_text},
            metrics=None,
        )
        tracer.close()
        render_trace_report(trace_dir=case_dir)

    return status == "ok"


def run_trace_cases(
    cases: list[dict[str, Any]],
    out_root: Path,
    *,
    use_chirality: bool = False,
    max_resonance: int = 50,
) -> int:
    out_root.mkdir(parents=True, exist_ok=True)

    all_ok = True
    for idx, case in enumerate(cases, start=1):
        case_dir = out_root / f"run-case{idx:03d}"
        case_dir.mkdir(parents=True, exist_ok=True)
        case_ok = _run_case(
            case=case,
            case_dir=case_dir,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
        )
        all_ok = all_ok and case_ok

    return 0 if all_ok else 1
