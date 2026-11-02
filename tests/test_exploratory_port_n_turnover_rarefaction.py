"""Offline negative controls on GBIF recording-effort thinning."""
import importlib.util,sys,math
from pathlib import Path
F=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_n_turnover_rarefaction.py"
spec=importlib.util.spec_from_file_location("port_n_rarefaction",F)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_rarefy_without_replacement_by_species_key():
    import random
    x={"1":10,"2":1,"3":3}
    assert m.sampled_keys_from_records(x,random.Random(1),14)=={1,2,3}
    assert m.sampled_keys_from_records(x,random.Random(1),5) <= {1,2,3}

def test_exact_focal_vs_other_mean_requires_scored_species():
    traits={i:{"N":v} for i,v in enumerate([1,3,5,7,9,8,7,6],start=1)}
    focal={1,2}
    assert m.mean_group_N(set(traits),traits,focal)==(1+3)/2-(5+7+9+8+7+6)/6
    assert m.mean_group_N({1,3,4,5,6,7,8},traits,focal) is None

def test_exact_early_late_turnover_no_fake_colonization():
    x={}
    for side in m.SIDES:
        for per,keys in [(m.PERIODS[0],{"1":2,"2":7,"3":1}),
                         (m.PERIODS[1],{"1":3,"3":5,"4":9})]:
            x[("10",1,side,per)]={"faceted_specieskey_counts":keys}
    traits={j:{"N":j} for j in range(1,5)}
    rows,summary=m.exact_set_turnover(x,traits,{1,2},[("10",1)])
    h=next(r for r in rows if r["group"]=="focal_host_N" and r["side"]=="port")
    assert h["shared_species"]==1 and h["newly_indexed_species"]==0
    assert h["lost_species"]==1
    bg=next(r for r in rows if r["group"]=="other_plant_N" and r["side"]=="port")
    assert bg["newly_indexed_species"]==1 and bg["lost_species"]==0

def test_full_cohort_rarefaction_is_deterministic_given_rng():
    import random
    vals=m.sampled_keys_from_records({"1":2,"2":3},random.Random(9),2)
    assert vals == m.sampled_keys_from_records({"1":2,"2":3},random.Random(9),2)
