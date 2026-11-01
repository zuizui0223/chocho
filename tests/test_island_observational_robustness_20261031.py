"""Source-locked, independent audit of 51-island compositional association."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031"
PROGRAM = DIR / "audit_pair_observation.py"
RECEIPT = DIR / "island_pair_sensitivity_audit.json"


def run_program():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island/port source files are not in independent GEB review pack")
    spec = importlib.util.spec_from_file_location("island_obs_20261031", PROGRAM)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def close(actual, expected, tolerance=1e-10):
    assert math.isclose(actual, expected, abs_tol=tolerance, rel_tol=0)


def test_observation_is_not_just_double_zero_pairs():
    m = run_program()
    result = m.analyse(repetitions=99)
    observed = result["observed"]
    assert result["eligible_islands"] == 51
    assert result["island_pairs_not_independent"] == 1275
    close(observed["host_new_rho"], 0.29980576483377025)
    close(observed["host_prior_rho"], 0.4178255833828514)
    close(observed["partial_distance_prior"], 0.2153535842700987)
    close(observed["partial_distance_prior_richness"], 0.22080814781853753)
    assert observed["pairs_host_positive"] == 711
    assert observed["pairs_new_positive"] == 600
    assert observed["pairs_both_zero"] == 358
    assert observed["within_positive_subsets"]["both_positive"]["n"] == 394
    close(observed["within_positive_subsets"]["both_positive"]["rho"],
          0.3261260753902157)
    assert [d["n"] for d in observed["within_distance_quartiles"]] == [319, 319, 319, 318]
    assert all(d["rho"] > 0 for d in observed["within_distance_quartiles"])
    assert observed["host_prior_rho"] > observed["host_new_rho"]


def test_permutation_receipt_reproduced_from_source_without_independent_pairs():
    m = run_program()
    fresh = m.analyse(repetitions=4999)
    receipt = json.loads(RECEIPT.read_text(encoding="utf8"))
    assert fresh["source_git_blobs"] == receipt["source_git_blobs"]
    assert fresh["interpretation_limits"] == receipt["interpretation_limits"]
    assert fresh["status"] == receipt["status"]
    for key, number in receipt["observed"].items():
        actual = fresh["observed"][key]
        if isinstance(number, (float, int)) and not isinstance(number, bool):
            close(actual, number, tolerance=1e-9)
    assert set(fresh["permutation_sensitivity"]) == set(receipt["permutation_sensitivity"])
    for name, expected in receipt["permutation_sensitivity"].items():
        actual = fresh["permutation_sensitivity"][name]
        for key in ("permutations", "seed", "groups", "largest_group", "movable_islands"):
            assert actual[key] == expected[key], (name, key)
        for key in ("observed", "one_sided_positive_p", "two_sided_abs_p",
                    "null_mean", "null_q025", "null_q975"):
            close(actual[key], expected[key], tolerance=2e-4)
    # Do NOT describe the residual-label p-values as confirmatory/causal.
    assert "not preregistered" in " ".join(fresh["interpretation_limits"])
