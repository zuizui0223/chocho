"""Do not turn harbor-node indices, missing links or GIS rasters into routes."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031/port_pair_network_gate.py"


def module():
    spec = importlib.util.spec_from_file_location("port_pair_network_gate", SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_link_coverage_unknown_zero_and_measured_modes():
    m = module()
    mapping = [
        dict(port_id="WPI10", island_id="a", source_port_id="PORT10",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="True"),
        dict(port_id="WPI11", island_id="b", source_port_id="PORT11",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="True"),
        dict(port_id="WPI12", island_id="c", source_port_id="PORT12",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="False"),
        dict(port_id="WPI13", island_id="d", source_port_id="PORT13",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="True"),
    ]
    edges = [
        dict(from_port_id="PORT10", to_port_id="PORT11", year="2019",
             weight="3", measurement_type=m.OBSERVED),
        dict(from_port_id="PORT11", to_port_id="PORT10", year="2020",
             weight="2", measurement_type=m.OBSERVED),
        dict(from_port_id="PORT10", to_port_id="PORT11", year="2019",
             weight="3500", measurement_type=m.MODELED),
    ]
    p = lambda a, b: dict(island_a=a, island_b=b, host_jaccard=".2",
                         new_butterfly_jaccard=".3", earlier_butterfly_jaccard=".1",
                         port_distance_km="123")
    result = m.scan(mapping, edges, 2019, 2020,
                    [p("a", "b"), p("a", "c"), p("a", "d")])
    rows = result["pair_rows"]
    assert result["covered_pairs"] == 2
    assert result["uncovered_pairs"] == 1
    assert rows[0]["observed_direct_link_weight"] == 5.0
    assert rows[0]["modeled_freight_link_weight"] == 3500.0
    assert rows[1]["island_pair_covered"] is False
    assert rows[1]["observed_direct_link_weight"] is None
    assert rows[2]["island_pair_covered"] is True
    assert rows[2]["observed_direct_link_weight"] == 0.0
    assert result["matched_link_rows"] == 3


def test_node_score_is_not_an_od_edge():
    m = module()
    mapping = [
        dict(port_id="X", island_id="a", source_port_id="X",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="True"),
        dict(port_id="Y", island_id="b", source_port_id="Y",
             coverage_start_year="2019", coverage_end_year="2020",
             complete_observation="True"),
    ]
    pair = [dict(island_a="a", island_b="b", host_jaccard=".1",
                 new_butterfly_jaccard=".2",
                 earlier_butterfly_jaccard=".3",
                 port_distance_km="100")]
    bad = [dict(from_port_id="X", to_port_id="Y", year="2019",
                weight="10", measurement_type="PLSCI_PORT_SCORE")]
    with pytest.raises(ValueError, match="cannot be a route"):
        m.scan(mapping, bad, 2019, 2020, pair)
    good = [{**bad[0], "measurement_type": m.OBSERVED}]
    with pytest.raises(ValueError, match="multiple islands"):
        m.scan(mapping + [{**mapping[1], "island_id": "c"}],
               good, 2019, 2020, pair)
    with pytest.raises(ValueError, match="reversed"):
        m.scan(mapping, good, 2020, 2019, pair)
