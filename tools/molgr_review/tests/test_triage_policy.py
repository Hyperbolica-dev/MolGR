from tools.molgr_review.triage_policy import (
    mapping_allows_strong,
    mapping_csv_fields,
    mapping_fallback_bucket,
)


def test_only_proven_mapping_can_enter_strong_bucket() -> None:
    assert mapping_allows_strong(
        {
            "confidence": "unique_graph_mapping",
            "enumeration_truncated": False,
            "timeout": False,
            "signature_count": 1,
        }
    )
    assert not mapping_allows_strong(
        {
            "confidence": "ambiguous",
            "enumeration_truncated": False,
            "timeout": False,
            "signature_count": 1,
        }
    )
    assert not mapping_allows_strong(
        {
            "confidence": "unique_graph_mapping",
            "enumeration_truncated": True,
            "timeout": False,
            "signature_count": 1,
        }
    )


def test_ambiguous_mapping_routes_by_relevant_diff() -> None:
    assert (
        mapping_fallback_bucket(["metal-coordination-edge difference"], "ambiguous")
        == "metal_coordination_ambiguous"
    )
    assert (
        mapping_fallback_bucket(["hydrogen-assignment difference"], "ambiguous")
        == "complex_multi_difference"
    )
    assert mapping_fallback_bucket(["unknown"], "failed") == "unknown"


def test_mapping_csv_fields_are_stable() -> None:
    assert mapping_csv_fields(
        {
            "confidence": "ambiguous",
            "mapping_count_examined": 128,
            "enumeration_truncated": True,
            "signature_count": 2,
        }
    ) == {
        "mapping_confidence": "ambiguous",
        "mapping_count_examined": 128,
        "mapping_truncated": True,
        "mapping_signature_count": 2,
    }
