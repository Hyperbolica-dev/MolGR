from __future__ import annotations

from typing import Any


STRONG_MAPPING_CONFIDENCES = frozenset({"exact", "unique_graph_mapping"})


def mapping_allows_strong(mapping: dict[str, Any]) -> bool:
    return (
        mapping.get("confidence") in STRONG_MAPPING_CONFIDENCES
        and not mapping.get("enumeration_truncated", False)
        and not mapping.get("timeout", False)
        and mapping.get("signature_count") == 1
    )


def mapping_fallback_bucket(diff_labels: list[str], confidence: str) -> str:
    if confidence == "failed":
        return "unknown"
    if "metal-coordination-edge difference" in diff_labels:
        return "metal_coordination_ambiguous"
    if "hydrogen-assignment difference" in diff_labels or "multiple differences" in diff_labels:
        return "complex_multi_difference"
    return "unknown"


def mapping_csv_fields(mapping: dict[str, Any]) -> dict[str, Any]:
    return {
        "mapping_confidence": mapping.get("confidence", "failed"),
        "mapping_count_examined": mapping.get("mapping_count_examined", 0),
        "mapping_truncated": bool(mapping.get("enumeration_truncated", False)),
        "mapping_signature_count": mapping.get("signature_count", 0),
    }
