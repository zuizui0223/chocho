import importlib.util,sys,math
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_allflora_N_offline_audit.py"
spec=importlib.util.spec_from_file_location("allflora_robust",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def observation(side,period,host,other,plant=100,coverage=.5):
    return {
        "island_id":"10","port_id":5,"side":side,"period":period,
        "plant_records":plant,"N_scored_candidate_host_species":3,
        "N_scored_other_plant_species":10,
        "N_scored_fraction_of_all_indexed_plant_species":coverage,
        "mean_N_candidate_hosts":host,
        "mean_N_other_scored_plants":other
    }
def test_paired_and_temporal_exact_key():
    a="early_2001_2010";b="later_2011_2024"
    obs=[observation("port",a,7,5),observation("nonport",a,5,6),
         observation("port",b,9,6),observation("nonport",b,5,6)]
    records=m.paired({"observations":obs})
    assert len(records)==2
    assert m.temporal_same_pairs(records)[0]["target_minus_background"]==1
    assert m.stats(m.temporal_same_pairs(records),"target_minus_background")["island_mean"]==1
def test_coverage_and_effort_exclusions():
    a="early_2001_2010"
    obs=[observation("port",a,7,6,plant=100,coverage=.19),
         observation("nonport",a,6,5,plant=30,coverage=.25)]
    data={"observations":obs}
    assert len(m.paired(data,min_coverage=.2))==0
    assert len(m.paired(data,min_effort_ratio=.3))==1
    assert len(m.paired(data,min_effort_ratio=.4))==0
def test_exact_sign_flip_six_islands():
    assert math.isclose(m.signflip_islands([1,2,3,4,5,6]),2/64)
