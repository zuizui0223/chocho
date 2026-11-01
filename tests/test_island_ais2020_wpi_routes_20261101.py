"""Public 2020 AIS source port endpoint coverage, positive-only island voyage tests."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/derive_ais_2020_island_port_links.py"


def load():
    pytest.importorskip("numpy")
    spec=importlib.util.spec_from_file_location("ais2020_WPI",SRC)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_no_ambiguous_wpi_port_and_never_impute_unseen_routes():
    m=load()
    # WPI lat,lng. B and C ports are 2.2 km apart; an AIS site at
    # their midpoint is ambiguous even though both are plausible harbors.
    original=np.asarray([[0.,0.],[0.,1.],[0.,1.02]])
    site={"A":(0.,0.),"B":(1.,0.),"ambiguous":(1.01,0.),
          "too_far":(6.,0.)}
    match,states,dist=m.match_source_sites(site,original,5)
    assert match["A"]==0
    assert "B" not in match
    assert "ambiguous" not in match
    assert "too_far" not in match
    assert states["unique_WPI_coordinate_match"]==1
    assert states["ambiguous_multiple_WPI_ports_in_radius"]==2
    assert states["not_locatable_to_frozen_194_WPI_ports"]==1
    with pytest.raises(ValueError,match="predeclared"):
        m.match_source_sites(site,original,20)


def test_observed_direct_AIS_link_retained_missing_links_unknown():
    m=load()
    ports=[
        {"port_id":"1","landmass_id":"10"},
        {"port_id":"2","landmass_id":"20"},
        {"port_id":"3","landmass_id":"30"},
    ]
    # The real source accepts exactly the frozen 51 choose 2 cohort;
    # synthetic tests preserve that scientific grain.
    islands=["10","20","30"]+[str(i) for i in range(100,148)]
    pairs=[
        {"island_a":a,"island_b":b}
        for i,a in enumerate(islands)
        for b in islands[i+1:]
    ]
    routes={
        ("SITE_A","SITE_B"):4,
        ("SITE_B","SITE_A"):3,
    }
    matched={"SITE_A":0,"SITE_B":1,"SITE_C":2}
    observed,edges,a=m.to_derived_island_network(
        {},routes,{"raw_source_route_rows":7},matched,ports,pairs,5,32500111)
    assert len(observed)==1275
    assert observed[0]["ais_2020_direct_voyages_observed_in_sample"]==7
    assert observed[0]["has_positive_2020_AIS_observed_from_to_pair"]==1
    assert observed[1]["ais_2020_direct_voyages_observed_in_sample"]==""
    assert observed[1]["island_2020_AIS_route_without_link_is_UNKNOWN"]==1
    assert a["source_undirected_island_pairs_observed"]==1
    assert a["matched_voyages_between_different_islands"]==7
    assert edges[0]["ais_2020_sampled_voyage_rows"]==7


def test_actual_archival_source_columns_stay_dated_from_to():
    m=load()
    assert "lon_from" in m.SOURCE_COLUMNS
    assert "lon_to" in m.SOURCE_COLUMNS
    assert "date-out" in m.SOURCE_COLUMNS
    assert m.EXPECTED_MD5=="af27e1bef1f706e2e29ed53aee47ff10"
    assert m.EXPECTED_CSV_BYTES==113764323
