"""Regression checks for the isolated, post-hoc geographic redundancy analysis."""
from __future__ import annotations

import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = runpy.run_path(str(ROOT / "scripts/analyze_butterfly_resource_redundancy_hostcount.py"))


def test_average_tie_ranks_and_seed():
    assert ANALYSIS["ranks"]([1.0, 1.0, 3.0, 2.0]) == [1.5, 1.5, 4.0, 3.0]
    rng = ANALYSIS["Mulberry32"](20261011)
    assert rng.random() == pytest.approx(0.34649015543982387)
    assert rng.random() == pytest.approx(0.8141545492690057)


def test_frozen_hostcount_diagnostic_receipt():
    expansion = ANALYSIS["csv_index"](
        ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    )
    contributions = ANALYSIS["csv_index"](
        ROOT / "data/frozen/figure_sources/host_contribution_metrics.csv"
    )
    rows = ANALYSIS["prepare"](expansion, contributions)
    assert len(rows) == 191
    got = ANALYSIS["analyze"](rows, permutations=9999, seed=20261011, bootstraps=10000)
    assert got["marginal_spearman_family_breadth_vs_multiple_host_fraction"] == pytest.approx(
        0.3730510474391473
    )
    assert got["exact_host_species_count_adjusted_rank_correlation"] == pytest.approx(
        -0.00041154985111599223
    )
    assert got["permutation_p_two_sided"] == pytest.approx(0.9954)
    assert got["n_informative_species"] == 158
    assert [v["species"] for v in got["host_count_bins"]] == [32, 37, 42, 80]
    assert [v["single_host_added_units"] for v in got["host_count_bins"]] == [
        1653, 1426, 1831, 2823
    ]
    assert got["host_count_bins"][-1]["single_host_added_share"] == pytest.approx(
        0.4535668380462725
    )
    assert got["rich_host_ge10_species_bootstrap_q025"] == pytest.approx(
        0.40752465658398124, abs=1e-10
    )
    assert got["rich_host_ge10_species_bootstrap_q975"] == pytest.approx(
        0.5009126742730914, abs=1e-10
    )
