"""Audit global 2010 vs five-year index and country-node sensitivity."""
import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031"
SCRIPT = BASE / "audit_full_lsbci2010_sensitivity.py"


def load():
    spec = importlib.util.spec_from_file_location("audit_global_lsbci", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_2010_country_index_is_not_port_network_effect():
    if "_anonymous_stage" in ROOT.parts or not (BASE/"full_lsbci2010_island_pairs.csv").exists():
        pytest.skip("Exploratory full UNCTAD source results not part of GEB anonymous review")
    d = load().evaluate()
    assert (d["island_pairs"], d["unique_islands"], d["countries"], d["country_pairs"]) == (1130, 49, 26, 325)
    c = d["main_2010"]
    assert math.isclose(c["ship_host_raw_rho"], 0.16676837789933527, abs_tol=1e-8)
    assert math.isclose(c["ship_butterfly_raw_rho"], 0.04053820312186473, abs_tol=1e-8)
    assert math.isclose(c["ship_butterfly_adjusted_host"], -0.02537308792904556, abs_tol=1e-8)
    assert math.isclose(c["host_butterfly_adjusted_shipping"], 0.17541155606848266, abs_tol=1e-8)
    c5 = d["mean_2006_2010"]
    assert math.isclose(c5["ship_butterfly_adjusted_host"], -0.0033090715355358537, abs_tol=1e-8)
    assert math.isclose(c5["host_butterfly_adjusted_shipping"], 0.1743978226034486, abs_tol=1e-8)
    loo = d["leave_one_country_out"]
    assert loo["n_countries"] == 26
    assert loo["ship_adjusted_negative_folds"] == 22
    assert loo["host_adjusted_positive_folds"] == 26
    assert loo["ship_adjusted_min_max"][0] < 0 < loo["ship_adjusted_min_max"][1]
    assert all(x > 0 for x in loo["host_adjusted_min_max"])
    assert "NOT_PORT_VOYAGES" in d["status"]
