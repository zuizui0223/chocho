"""Information must come from strata with nonzero conditional variance."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_butterfly_239_GBIF_occurrence_host_opportunity_20261019 import (
    stratum_information_audit,
)


def test_zero_or_all_positive_strata_do_not_count_as_information():
    strata = [
        {"N": 8, "A": 3, "K": 0},
        {"N": 8, "A": 3, "K": 8},
        {"N": 8, "A": 3, "K": 2},
    ]
    result = stratum_information_audit(strata)
    assert result["mixed_resource_strata"] == 3
    assert result["positive_strata_legacy_definition"] == 2
    assert result["nonzero_conditional_variance_strata"] == 1
    assert result["zero_positive_strata"] == 1
    assert result["all_positive_strata"] == 1
    assert result["status"] == "BELOW_PROTOCOL_INFORMATION_THRESHOLD"


def test_fifty_informative_strata_are_required_even_with_many_empty_strata():
    varied = {"N": 8, "A": 3, "K": 2}
    empty = {"N": 8, "A": 3, "K": 0}
    assert stratum_information_audit([varied] * 49 + [empty] * 200)["status"] == (
        "BELOW_PROTOCOL_INFORMATION_THRESHOLD"
    )
    result = stratum_information_audit([varied] * 50)
    assert result["protocol_minimum_informative_strata"] == 50
    assert result["status"] == "MEETS_PROTOCOL_INFORMATION_THRESHOLD"
