"""Figure-source integrity: show observed versus structural expectations honestly."""
from __future__ import annotations

import importlib.util
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/render_butterfly_resource_novelty_figure.py"
spec = importlib.util.spec_from_file_location("resource_novelty_figure", SCRIPT)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BASE = ROOT / "provenance/reviewer_defenses/results"
NO_CLIMATE = BASE / "butterfly_plant_range_structural_null_two_chains_v0.2.json"
CLIMATE = BASE / "butterfly_climate_matched_absorption_v0.1.json"
LOCAL = BASE / "butterfly_spatial_null_longchain_qualified_v0.2.json"


def test_plot_requires_all_original_source_identifiers() -> None:
    data = mod.extract_results(
        mod.load(NO_CLIMATE), mod.load(CLIMATE), mod.load(LOCAL)
    )
    assert data["total_candidate_incidences"] == 75346
    assert [p["n"] for p in data["accounting"]] == [46387, 14406, 14553]
    assert sum(p["n"] for p in data["accounting"]) == 75346
    assert data["focal_butterflies"] == 239
    assert data["resource_regions"] == 365
    assert data["controls"][0]["a_excess"] == [1086, 1089]
    assert data["controls"][0]["u_deficit"] == [831, 832]
    assert data["controls"][1]["a_excess"] == [552, 558]
    assert data["controls"][1]["u_deficit"] == [426, 431]
    assert data["controls"][2]["a_excess"] == [271, 276]
    assert data["controls"][2]["u_deficit"] == [264, 267]
    assert [c["draws_per_chain"] for c in data["controls"]] == [499, 299, 99]


def test_renderer_cannot_promote_unqualified_spatial_null() -> None:
    b, c, s = mod.load(NO_CLIMATE), mod.load(CLIMATE), mod.load(LOCAL)
    local = deepcopy(s)
    local["qualification"]["checks"][0]["pass"] = False
    with pytest.raises(ValueError, match="not qualified"):
        mod.extract_results(b, c, local)


def test_renderer_refuses_incidence_unit_mismatch() -> None:
    b, c, s = mod.load(NO_CLIMATE), mod.load(CLIMATE), mod.load(LOCAL)
    wrong = deepcopy(b)
    wrong["accounting"]["observed_unique_added_butterfly_region_cells"] += 1
    with pytest.raises(ValueError, match="identity has changed"):
        mod.extract_results(wrong, c, s)


def test_chart_describes_what_it_does_not_measure() -> None:
    b, c, s = mod.load(NO_CLIMATE), mod.load(CLIMATE), mod.load(LOCAL)
    data = mod.extract_results(b, c, s)
    assert "not local larval use" in data["definition"]
    assert "colonization" in data["definition"]
