"""Exact WPI original source inverse solver: unique solution or fail closed."""
import pytest
for x in ("numpy","pandas"):pytest.importorskip(x)
import importlib.util,math
import numpy as np
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"recover_island51_exact_original_port_ids.py"
s=importlib.util.spec_from_file_location("solve51",f)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def test_unique_port_source_completion_uses_archive_constraints_not_jaccard():
    loc={"1":[0,0],"2":[0,2],"3":[22,2],"4":[15,0],"5":[-40,33]}
    anchor={"A":m.mean_xyz(["3"],loc),"B":m.mean_xyz(["4"],loc),
            "C":m.mean_xyz(["5"],loc)}
    true=m.mean_xyz(["1","2"],loc)
    distances={k:m.angular_distance(true,v) for k,v in anchor.items()}
    r=m.candidates_one_missing(["1"],["1","2","4"],2,loc,distances,anchor)
    assert r["n_candidates_uniquely_satisfying_all_known_distances"]==1
    assert r["winner"]["extra_port_ids"]==["2"]
    assert r["winner"]["max_source_distance_error_km"]<1e-6

def test_absent_original_port_id_never_inferred_by_nearest_guess():
    loc={"1":[0,0],"2":[3,0],"3":[21,0]}
    anchor={"A":m.mean_xyz(["3"],loc)}
    true=m.mean_xyz(["1","2"],loc)
    r=m.candidates_one_missing(["1"],["1","3"],2,loc,
          {"A":m.angular_distance(true,anchor["A"])},anchor)
    assert r["winner"] is None
    assert r["n_candidates_uniquely_satisfying_all_known_distances"]==0
