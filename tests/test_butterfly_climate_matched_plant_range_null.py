"""Synthetic unit tests for climate-conditioned range-shuffle design only."""
from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/analyze_butterfly_climate_matched_plant_range_null.py"
spec = importlib.util.spec_from_file_location("climate_plant_null", PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def _climate(values, code: str, n: int = 12) -> str:
    return (
        f"test,{code},{n},"
        + ",".join(str(value) for value in values)
        + "\n"
    )


def test_climate_duplicates_are_consistent_across_butterflies(
    tmp_path: Path,
) -> None:
    f = tmp_path / "units.csv"
    f.write_text(
        "species,wgsrpd3_code,climate_sample_points,bio1,bio7,bio12,bio15\n"
        + _climate([12, 6, 700, 80], "AAA")
        + _climate([12, 6, 700, 80], "AAA")
        + _climate([8, 10, 400, 60], "BBB"),
        encoding="utf-8",
    )
    units, metadata = m.read_unit_climate(f)
    assert len(units) == 2
    assert units["AAA"] == (12.0, 6.0, 700.0, 80.0)
    assert metadata["duplicate_rows_validated"] == 1
    assert len(metadata["source_file_sha256"]) == 64

    f.write_text(
        "species,wgsrpd3_code,climate_sample_points,bio1,bio7,bio12,bio15\n"
        + _climate([12, 6, 700, 80], "AAA")
        + _climate([12.2, 6, 700, 80], "AAA"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="nonidentical"):
        m.read_unit_climate(f)


def test_missing_climate_is_not_imputed(tmp_path: Path) -> None:
    f = tmp_path / "units.csv"
    f.write_text(
        "species,wgsrpd3_code,climate_sample_points,bio1,bio7,bio12,bio15\n"
        + _climate([12, 6, 700, 80], "AAA")
        + "B,BBB,0,,,,\n",
        encoding="utf-8",
    )
    out, meta = m.read_unit_climate(f)
    assert set(out) == {"AAA"}
    assert meta["regions_only_missing_climate"] == 1


def test_categorical_blocks_do_not_exchange_cross_climate() -> None:
    codes = ["R1", "R2", "R3", "R4"]
    level1 = {i: "X" for i in range(4)}
    sig = {
        "R1": (0., 10., 100., 0.),
        "R2": (1., 10., 110., 0.),
        "R3": (9., 10., 900., 0.),
        "R4": (10., 10., 910., 0.),
    }
    blocks, meta = m.climate_strata(
        sig, codes, level1, "temperature_precipitation_2x2"
    )
    assert blocks[0] == blocks[1]
    assert blocks[2] == blocks[3]
    assert blocks[0] != blocks[2]
    assert meta["variables"] == ["bio1", "bio12"]
    source = {
        "h1": {0, 2},
        "h2": {1, 3},
    }
    groups, cumulative, diag = m.grouped_introduced_edges(
        source, {"h1": "Fabaceae", "h2": "Fabaceae"}, level1, blocks
    )
    assert len(groups) == 2
    assert diag["total_introduced_plant_region_edges"] == 4
    assert diag["viable_edges"] == 4
    assert diag["climate_covered_fraction"] == 1
    shuffled = {host: set(v) for host, v in source.items()}
    accepted = m.attempt_range_swaps(
        groups, cumulative, {"h1": set(), "h2": set()},
        shuffled, 250, random.Random(12)
    )
    assert accepted > 0
    assert {h: len(units) for h, units in shuffled.items()} == {
        h: 2 for h in shuffled
    }
    assert all(len(units & {0, 1}) == 1 for units in shuffled.values())
    assert all(len(units & {2, 3}) == 1 for units in shuffled.values())
    assert all(sum(r in x for x in shuffled.values()) == 1 for r in range(4))


def test_plant_family_and_continent_are_always_preserved() -> None:
    intro = {"F1": {0}, "G1": {1}, "F2": {2}, "F3": {3}}
    climate_classes = {i: (0, 0) for i in range(4)}
    level1 = {0: "X", 1: "X", 2: "Y", 3: "Y"}
    families = {"F1": "F", "F2": "F", "F3": "F", "G1": "G"}
    groups, cumulative, diag = m.grouped_introduced_edges(
        intro, families, level1, climate_classes
    )
    assert diag["viable_edges"] == 2
    assert diag["fixed_edges"] == 2
    assert len(groups) == 1
    assert set(groups[0]) == {("F2", 2), ("F3", 3)}
    assert cumulative == [2]


def test_insufficient_climate_coverage_is_not_evaluable() -> None:
    intro = {"a": {0, 1, 2}, "b": {0, 1, 2}}
    groups, edges, result = m.grouped_introduced_edges(
        intro,
        {"a": "F", "b": "F"},
        {0: "X", 1: "X", 2: "X"},
        {0: (0, 0)},
    )
    assert result["climate_covered_fraction"] == pytest.approx(1/3)
    assert result["viable_edges"] == 0
    assert edges == []


def test_quantile_cuts_repeat_when_values_are_equal() -> None:
    assert m.quantile_cutpoints([3.,3.,3.,3.], 2) == (3.,)
    assert m.quantile_cutpoints([0.,10.,20.,30.], 2) == (15.,)


def test_source_baselines_are_fixed_by_exact_incidences() -> None:
    assert m.EXPECTED["introduced_host_butterfly_region_incidences"] == 75346
    assert m.EXPECTED["introduced_incidences_in_native_supported_cells"] == 46387
    assert m.EXPECTED["added_units"] == 14553
    assert m.EXPECTED["introduced_incidences_in_native_supported_cells"] + (
        m.EXPECTED["introduced_incidences_in_newly_supported_cells"]
    ) == m.EXPECTED["introduced_host_butterfly_region_incidences"]
    assert m.EXPECTED["introduced_incidences_in_newly_supported_cells"] - (
        m.EXPECTED["redundant_introduced_host_incidences_in_new_cells"]
    ) == m.EXPECTED["added_units"]
