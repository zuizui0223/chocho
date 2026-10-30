"""Consequential checks for finite-record richness and supported comparisons."""
import importlib.util
import itertools
import math
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
PATH = Path(__file__).resolve().parents[1] / "scripts/analyze_expanded_island_recordings_20261029.py"


def module():
    assert PATH.exists(), "Prespecified recording analysis is not implemented"
    spec = importlib.util.spec_from_file_location("expanded_analysis", PATH)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_rarefaction_matches_all_finite_subsets_not_a_random_draw():
    m = module()
    pool = ["a", "a", "a", "b", "b", "c"]
    expected = np.mean([len(set(s)) for s in itertools.combinations(pool, 3)])
    assert m.rarefied_richness([3, 2, 1], 3) == pytest.approx(expected, abs=1e-12)
    assert m.rarefied_richness([3, 2, 1], 6) == 3
    assert m.rarefied_richness([1] * 100, 20) == pytest.approx(20)
    assert m.rarefied_richness([1000000], 20) == 1


def test_rarefaction_rejects_insufficient_or_noninteger_record_pools():
    m = module()
    for counts, size in [([2, 1], 4), ([1.5, 2], 2), ([-1, 3], 1), ([3], 0)]:
        with pytest.raises(ValueError):
            m.rarefied_richness(counts, size)


def test_roster_clusters_retain_bridge_and_handle_dateline():
    m = module()
    rows = [{"island_key": key, "representative_lon": lon, "representative_lat": 0}
            for key, lon in [("a", 179), ("bridge", -179), ("c", -177), ("far", 0)]]
    groups, details = m.fixed_clusters(rows, 250)
    assert groups["a"] == groups["c"] == groups["bridge"]
    assert groups["far"] != groups["a"]
    assert max(x["diameter_km"] for x in details) > 250
    assert m.fixed_clusters(list(reversed(rows)), 250)[0] == groups


def specimen(**changes):
    return {"query_status": "COMPLETE", "species_assigned_records": 20,
            "species_assignment_fraction": .9, "positive_years": 3,
            "human_observation_fraction": .95, "dominant_dataset_share": .5, **changes}


def test_zero_missing_and_incomplete_assignment_are_not_zero_richness():
    m = module()
    assert m.eligibility_reasons(specimen(), 20) == []
    assert "query_not_complete" in m.eligibility_reasons(specimen(query_status="MISSING_QUERY_FAILED"), 20)
    assert "too_few_species_assigned_records" in m.eligibility_reasons(specimen(species_assigned_records=0), 20)
    assert "species_assignment_below_90pct" in m.eligibility_reasons(specimen(species_assignment_fraction=.899), 20)
    assert "fewer_than_required_recorded_years" in m.eligibility_reasons(specimen(positive_years=2), 20)


def regression_rows():
    rng = np.random.default_rng(42)
    result = []
    for i in range(48):
        gapclass = "under100" if i % 2 else "atleast500"
        gap = 10 + i if gapclass == "under100" else 500 + 37*i
        area = math.exp(float(rng.uniform(3, 9)))
        lat = float(rng.uniform(1, 50))
        region = f"longitude_{i % 4}"
        logy = .2 + .08*math.log10(area) - .11*math.log10(1+gap) + .04*(lat/10) + .03*(i % 4)
        result.append({"island_key": str(i), "area_km2": area, "representative_lat": lat,
                       "mainland_coast_gap_upper_km": gap, "gap_class": gapclass,
                       "port_class": "within50" if (i//2) % 2 else "beyond100",
                       "longitude_block": region, "cluster_500": f"g{i//2}",
                       "rarefied_species_20": math.exp(logy)})
    # Each of two paired longitude blocks must actually span both gap classes.
    for i, row in enumerate(result):
        row["longitude_block"] = f"longitude_{(i//2) % 4}"
        row["rarefied_species_20"] *= math.exp(.03*((i//2) % 4) - .03*(i % 4))
    return result


def test_fit_returns_known_gap_effect_only_after_support_gates():
    m = module()
    rows = regression_rows()
    result = m.gated_model(rows, 20)
    assert result["status"] == "FIT_EXPLORATORY_ASSOCIATION"
    assert result["effects"]["log10_one_plus_gap"]["ratio"] == pytest.approx(math.exp(-.11), abs=1e-9)
    assert result["cluster_df"] == 23
    assert result["geographic_clusters"] == 24
    restricted = [dict(r, gap_class="under100") for r in rows]
    denied = m.gated_model(restricted, 20)
    assert denied["status"] == "NOT_FIT_SUPPORT_GATE"
    assert "fewer_than_10_islands_per_gap_class" in denied["gate_failures"]
    assert "effects" not in denied


def test_port_fit_refuses_empty_gap_by_port_cell():
    m = module()
    rows = regression_rows()
    for r in rows:
        r["port_class"] = "within50" if r["gap_class"] == "under100" else "beyond100"
    result = m.gated_model(rows, 20, with_port=True)
    assert result["status"] == "NOT_FIT_SUPPORT_GATE"
    assert "fewer_than_3_islands_per_gap_port_cell" in result["gate_failures"]


def test_unrelated_taxonomic_namespace_is_rejected_before_facet_analysis():
    m = module()
    wrong_source = {"taxonomy": {"checklistKey": "unrelated-checklist", "taxonKey": "5G9"}}
    with pytest.raises(ValueError, match="taxonomy"):
        m.build_metrics({}, wrong_source)
