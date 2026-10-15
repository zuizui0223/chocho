"""Structural checks and frozen-input replication of exploratory redundancy audit."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analyze_butterfly_species_redundancy.py"
spec = importlib.util.spec_from_file_location("redundancy_audit", SCRIPT)
assert spec and spec.loader
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_tied_average_ranks() -> None:
    assert audit.average_ranks([2, 2, 1, 9]).tolist() == [2.5, 2.5, 1, 4]
    with pytest.raises(ValueError):
        audit.average_ranks([2, float("nan")])


def test_no_zero_fill_when_additions_absent(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text(
        "species,introduced_added_units,resolved_host_species,host_family_count,"
        "fraction_added_units_with_multiple_contributing_hosts,"
        "native_resource_units,host_taxonomy_lower_bound_adequate\n"
        "A,0,3,1,0.0,15,True\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="zero-added"):
        audit.read_frozen_rows(path)

    path.write_text(
        "species,introduced_added_units,resolved_host_species,host_family_count,"
        "fraction_added_units_with_multiple_contributing_hosts,"
        "native_resource_units,host_taxonomy_lower_bound_adequate\n"
        "A,5,3,1,,15,True\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="lacks"):
        audit.read_frozen_rows(path)


def test_partial_rank_is_not_a_marginal_rank_correlation() -> None:
    n = np.array([1., 2., 3., 4., 5., 6., 7., 8.])
    x = np.array([1., 4., 2., 5., 3., 6., 8., 7.])
    y = np.array([2., 1., 4., 3., 6., 5., 7., 8.])
    raw = audit.spearman(x, y)
    conditional = audit.partial_rank_correlation(x, y, (n,))
    assert np.isfinite(raw) and np.isfinite(conditional)
    assert abs(raw - conditional) > 1e-3


def test_frozen_species_results_are_exactly_grounded() -> None:
    rows, digest = audit.read_frozen_rows(
        ROOT / "data/frozen/figure_sources/host_contribution_metrics.csv"
    )
    assert digest == "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa"
    result = audit.describe(rows)
    assert result["n_focal"] == 239
    assert result["n_expanded"] == 206
    assert result["added_butterfly_region_units"] == 14553
    assert result["weighted_single_host_supported_units"] == 8574
    assert result["weighted_single_host_supported_fraction"] == pytest.approx(0.589157, abs=1e-6)
    assert result["species_median_single_host_supported_fraction"] == pytest.approx(0.688259, abs=1e-6)
    assert [result["strata"][k]["butterflies"] for k in ("1", "2-4", "5-9", "10+")] == [14, 67, 45, 80]
    assert result["strata"]["10+"]["majority_single_supported_butterflies"] == 37
    assert result["rich_portfolio_probe"]["butterflies"] == 78
    assert result["rich_portfolio_probe"]["majority_single_supported"] == 35
    assert result["rich_portfolio_probe"]["adequate_butterflies"] == 78
    effects = result["rank_associations"]
    assert effects["spearman_families_multi_support"] == pytest.approx(0.290752, abs=1e-6)
    assert effects["spearman_host_species_multi_support"] == pytest.approx(0.668770, abs=1e-6)
    assert effects["partial_families_given_host_species"] == pytest.approx(0.015444, abs=1e-6)
    assert effects["partial_families_given_host_species_native_and_added_area"] == pytest.approx(0.035279, abs=1e-6)
    assert effects["partial_host_species_given_families"] == pytest.approx(0.629567, abs=1e-6)
