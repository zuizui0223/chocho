"""Offline tests for non-submission exploratory port nutrient pilot."""
from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

FILE = Path(__file__).resolve().parents[1] / "scripts" / "exploratory_port_nutrient_ellenberg.py"
spec = importlib.util.spec_from_file_location("port_nutrient_ellenberg_pilot", FILE)
assert spec and spec.loader
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def test_nutrient_parser_is_strict():
    assert pilot.parse_value("5") == 5
    assert pilot.parse_value("(5)") == 5
    assert pilot.parse_value("5~") == 5
    assert pilot.parse_value("  (7)  ") == 7
    assert pilot.parse_value("0") is None
    assert pilot.parse_value("0", allow_zero=True) == 0
    assert pilot.parse_value("x") is None
    assert pilot.parse_value("5)") is None
    assert pilot.parse_value("(5") is None
    assert pilot.parse_value("10") is None


def test_exact_matching_does_not_impute_missing_traits(tmp_path):
    host_file = tmp_path / "hosts.csv"
    trait_file = tmp_path / "traits.tsv"
    with host_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["host_name", "ports"])
        writer.writerows([["Plantago major", 30], ["Rumex crispus", 5], ["NotInTraits species", 1]])
    with trait_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Name", "L", "F", "N", "S"])
        writer.writerows([["Plantago major", "8", "6", "(8)", "0"],
                          ["Rumex crispus", "7", "5", "6", "x"]])
    rows, meta = pilot.load_data(host_file, trait_file)
    assert len(rows) == 2
    assert meta["n_port_candidate_host_species"] == 3
    assert meta["n_exact_species_matches_with_N"] == 2
    assert meta["n_unmatched_host_species"] == 1
    assert meta["median_ports_matched"] == 17.5
    assert meta["median_ports_unmatched"] == 1
    assert rows[0]["S"] == 0
    assert rows[1]["S"] is None


def test_tied_ranks_and_association():
    assert pilot.ranks([2, 2, 4]) == [1.5, 1.5, 3]
    assert math.isclose(pilot.spearman([1, 2, 3], [10, 20, 30]), 1)
    assert math.isclose(pilot.spearman([1, 2, 3], [30, 20, 10]), -1)


def test_within_genus_permutation_preserves_genus_and_does_not_claim_effect():
    sample = [
        {"genus": "A", "N": 1, "ports": 1},
        {"genus": "A", "N": 6, "ports": 22},
        {"genus": "B", "N": 3, "ports": 12},
        {"genus": "B", "N": 9, "ports": 15},
        {"genus": "C", "N": 5, "ports": 99},
    ]
    out1 = pilot.within_genus_test(sample, permutations=99)
    out2 = pilot.within_genus_test(sample, permutations=99)
    assert out1 == out2
    assert out1["eligible_genera"] == 2
    assert out1["eligible_species"] == 4
    assert 0 < out1["permutation_two_sided_p"] <= 1
    assert len(out1["null_95_percentiles"]) == 2
