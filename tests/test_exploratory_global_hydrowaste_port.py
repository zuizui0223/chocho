"""Global hydro-sewage and indexed host port geography sanity tests (offline)."""
import importlib.util,math,sys
from pathlib import Path
import numpy as np
import pytest
p=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_global_hydrowaste_port.py"
spec=importlib.util.spec_from_file_location("hydrowaste_port",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_spherical_wwtp_distance_and_neighbor_radius():
    pts=m.xyz(np.array([0.,0.,0.]),np.array([0.,1.,10.]))
    assert np.isclose(np.linalg.norm(pts[0]-pts[1]),m.ball_radius(111.195),atol=.0001)
    from scipy.spatial import cKDTree
    tree=cKDTree(pts)
    assert len(tree.query_ball_point(m.xyz([0.],[0.])[0],r=m.ball_radius(25)))==1
    assert len(tree.query_ball_point(m.xyz([0.],[0.])[0],r=m.ball_radius(200)))==2

def test_2019_port_not_equivalent_to_historical_outfall():
    assert m.KNOWN_EXPOSURES[0]=="wwtp_25km"
    assert "outfall_25km" in m.KNOWN_EXPOSURES
    assert 0<m.ball_radius(10)<m.ball_radius(50)<.02

def test_quality_wwtp_status_filters_are_not_inferred_from_missing():
    import inspect
    src=inspect.getsource(m.exposure)
    assert '"Operational","Not Reported"' in src
    assert "outfall" in src
    assert "qual" in src.lower()

def test_original_host_inventory_has_country_positive_and_incomplete_population():
    hosts=m.read_host_inventory(Path(__file__).resolve().parents[1]/
        "provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
    assert len(hosts)==1596
    assert (hosts.contemporary_not_native_in_catalogue>0).sum()==1150
