"""Hand-derived fixtures protecting overlap from misleading marginal positives."""
import importlib.util
from pathlib import Path


def module():
    path = Path(__file__).parents[1] / "scripts/analyze_island_host_identifiability_20261029.py"
    assert path.exists(), "The support diagnostic must exist"
    spec = importlib.util.spec_from_file_location("island_support", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_marginal_island_positives_do_not_supply_within_island_host_contrast():
    # A bug pooling different species would incorrectly turn these into support.
    rows = [
        {"species": "A a", "genus": "A", "stratum": "1", "cells": [
            {"country": "CU", "added": 1, "positive": 1},
            {"country": "GT", "added": 0, "positive": 0}]},
        {"species": "B b", "genus": "B", "stratum": "1", "cells": [
            {"country": "CU", "added": 0, "positive": 1},
            {"country": "GT", "added": 1, "positive": 0}]}]
    result, _, _ = module().analyze_support(rows, "fixture")
    island = result["within_class_host_support"]["wholly_island_country"]
    assert result["descriptive_margins_confounded"]["wholly_island_country"]["positive_cells"] == 2
    assert island["host_contrast_strata"] == 0
    assert island["informative_strata"] == 0


def test_four_group_overlap_without_any_detection_is_not_informative():
    # A bug treating exposure overlap alone as information would fail here.
    rows = [{"species": "A a", "genus": "A", "stratum": "1", "cells": [
        {"country": "CU", "added": 1, "positive": 0},
        {"country": "DO", "added": 0, "positive": 0},
        {"country": "GT", "added": 1, "positive": 0},
        {"country": "HN", "added": 0, "positive": 0}]}]
    result, _, _ = module().analyze_support(rows, "fixture")
    assert result["all_four_groups_support"]["strata"] == 1
    assert result["all_four_groups_support"]["cells"] == 4
    assert result["all_four_groups_support"]["informative_strata"] == 0
    rows[0]["cells"][0]["positive"] = 1
    result, _, _ = module().analyze_support(rows, "fixture")
    assert result["all_four_groups_support"]["informative_strata"] == 1
    assert result["within_class_host_support"]["wholly_island_country"]["informative_strata"] == 1


def test_late_only_risk_set_excludes_prior_detections_instead_of_recoding_them():
    # A bug keeping early-positive cells would invent eligible late-only trials.
    rows = [{"species": "A a", "genus": "A", "stratum": "1", "cells": [
        {"country": "CU", "added": 1, "early": 1, "late": 1},
        {"country": "DO", "added": 0, "early": 0, "late": 1}]}]
    result = module().directional_rows(rows, "late_only")
    assert len(result[0]["cells"]) == 1
    assert result[0]["cells"][0]["country"] == "DO"
    assert result[0]["cells"][0]["positive"] == 1
