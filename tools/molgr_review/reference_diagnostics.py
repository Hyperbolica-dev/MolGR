from __future__ import annotations

import json
import re
from collections import Counter
from typing import Any, Iterable, Mapping


_EQUIVALENCE_TIMEOUT = re.compile(r"\bequivalence\s+\d+\s+timed out after\s+[0-9.]+s$")
_CANDIDATE_REPARSE_FAILURE = "predicted_smiles could not be reparsed"


def comparison_skip_reasons(error: Any) -> list[str]:
    """Extract retained benchmark comparison-skip reasons without altering source rows."""

    if not error:
        return []
    payload: Any = error
    if isinstance(payload, str):
        try:
            payload = json.loads(payload)
        except json.JSONDecodeError:
            return []
    if not isinstance(payload, Mapping):
        return []
    rows = payload.get("rows")
    if not isinstance(rows, Mapping):
        return []
    reasons: list[str] = []
    for row in rows.values():
        if not isinstance(row, Mapping):
            continue
        skipped = str(row.get("comparison_skipped") or "").strip().lower()
        if skipped not in {"1", "true", "yes", "y"}:
            continue
        reason = str(row.get("comparison_skip_reason") or "").strip()
        if reason:
            reasons.append(reason)
    return reasons


def classify_reference_problem(
    *,
    reference_smiles: Any,
    skip_reasons: Iterable[str] = (),
    formula_status: Any = "",
) -> tuple[str, str]:
    """Return a diagnostic category separate from formula validity and review verdicts."""

    if not str(reference_smiles or "").strip():
        return "missing_reference", "Reference SMILES is unavailable."

    normalized_reasons = [str(reason).strip() for reason in skip_reasons if str(reason).strip()]
    if normalized_reasons and all(
        _EQUIVALENCE_TIMEOUT.search(reason) for reason in normalized_reasons
    ):
        return (
            "equivalence_timeout",
            "Candidate and Reference parsed, but their equivalence check timed out.",
        )
    if normalized_reasons and all(
        _CANDIDATE_REPARSE_FAILURE in reason for reason in normalized_reasons
    ):
        return (
            "candidate_reparse_failure",
            "Reconstruction succeeded, but its generated Candidate SMILES could not be reparsed.",
        )
    if str(formula_status or "").strip().lower() == "formula_mismatch":
        return "formula_mismatch", "The Reference and XYZ molecular formulae differ."
    if normalized_reasons:
        return "comparison_skipped", normalized_reasons[0]
    return "ok", ""


def reference_metal_charges(reference_smiles: str | None) -> list[tuple[str, int]] | None:
    """Return [(element, formal_charge)] for metal atoms in the reference SMILES.

    Returns None when the reference cannot be parsed, and [] when it parses but
    contains no metal atoms.
    """
    if not reference_smiles or not str(reference_smiles).strip():
        return None
    try:
        from rdkit import Chem

        from molgr.utils.equivalence import _NON_METAL_ATOMIC_NUMBERS
    except ImportError:
        return None
    try:
        mol = Chem.MolFromSmiles(str(reference_smiles))
    except Exception:
        return None
    if mol is None:
        return None
    out: list[tuple[str, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _NON_METAL_ATOMIC_NUMBERS:
            out.append((atom.GetSymbol(), int(atom.GetFormalCharge())))
    return out


def site_metal_valences(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
) -> dict[str, list[int]] | None:
    """Return {element: sorted candidate valences} for the XYZ metal site.

    Uses MolGR's real per-site enumeration path
    (``prepare_metal_state`` -> ``available_valence_radical_states``) rather than a
    hardcoded historical valence set. Returns None if the code path cannot run.
    """
    try:
        from molgr.fallback.utils.metals import preparation
    except ImportError:
        return None
    try:
        prepared = preparation.prepare_metal_state(
            xyz_block,
            int(total_charge or 0),
            int(total_radical_electrons or 0),
        )
    except Exception:
        return None
    valences: dict[str, set[int]] = {}
    for position in prepared.available_valence_radical_states:
        for state in position:
            symbol = str(state.symbol)
            valences.setdefault(symbol, set()).add(int(state.valence))
    if not valences:
        return None
    return {element: sorted(values) for element, values in valences.items()}


def _xyz_formula(xyz_block: str) -> Counter[str]:
    lines = xyz_block.splitlines()
    try:
        atom_count = int(str(lines[0]).strip())
    except (ValueError, IndexError):
        return Counter()
    atom_lines = lines[2 : 2 + atom_count]
    return Counter(str(line.split()[0]) for line in atom_lines if line.split())


def reference_formula_conserved(
    reference_smiles: str,
    xyz_block: str,
) -> tuple[bool, str]:
    """Return (conserved, mismatch_detail) comparing reference vs XYZ atom counts."""
    from rdkit import Chem

    xyz_counts = _xyz_formula(xyz_block)
    try:
        reference_mol = Chem.MolFromSmiles(str(reference_smiles))
    except Exception:
        reference_mol = None
    if reference_mol is None:
        return False, "reference SMILES could not be parsed"
    reference_with_h = Chem.AddHs(reference_mol)
    reference_counts: Counter[str] = Counter(
        atom.GetSymbol()
        for atom in reference_with_h.GetAtoms()  # pyright: ignore[reportCallIssue]
    )
    if xyz_counts == reference_counts:
        return True, ""
    detail = "; ".join(
        f"{symbol}:xyz={xyz_counts.get(symbol, 0)},ref={reference_counts.get(symbol, 0)}"
        for symbol in sorted(set(xyz_counts) | set(reference_counts))
        if xyz_counts.get(symbol, 0) != reference_counts.get(symbol, 0)
    )
    return False, detail


__all__ = [
    "classify_reference_problem",
    "comparison_skip_reasons",
    "reference_metal_charges",
    "reference_formula_conserved",
    "site_metal_valences",
]
