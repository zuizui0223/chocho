"""Offline tests: port versus nonport N of hosts and other plants."""
import importlib.util, math, sys
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_allplant_N_background.py"
sys.path.insert(0,str(f.parent))
spec=importlib.util.spec_from_file_location("allflorapilot",f)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_strict_binomial_and_value_parser(tmp_path):
    p=tmp_path/"N.tsv"
    p.write_text("Name\tN\nPlantago major\t(8)\nRumex crispus\t6\nRumex crispus subsp. tinctorius\t7\nSomething sp.\tx\nTrifolium repens\t7~\n")
    assert {x["name"]:x["N"] for x in m.read_trait_names(p)}=={
        "Plantago major":8,"Rumex crispus":6,"Trifolium repens":7}

def test_source_facet_is_species_weighted_not_occurrence_weighted():
    catalog={"accepted_unambiguous_by_speciesKey":{
        "1":{"name":"A a","N":9},"2":{"name":"B b","N":3},
        "3":{"name":"C c","N":7},"4":{"name":"D d","N":5}},
        "focal_host_by_speciesKey":{"1":{},"2":{}}}
    facet={"faceted_specieskey_counts":{"1":18,"2":2,"3":7,"4":3,"5":10},
           "indexed_plant_records":40,"indexed_species_keys":5,"records_with_speciesKey":40}
    result=m.summarize_facet(facet,catalog)
    assert result["N_scored_candidate_host_species"]==2
    assert result["N_scored_other_plant_species"]==2
    assert result["mean_N_hosts_minus_other"]==0
    assert result["N_scored_fraction_of_all_indexed_plant_species"]==.8

def site(side,period,host,bg,nb=10):
    return {"island_id":"7","port_id":41,"side":side,"period":period,
        "plant_records":60,"N_scored_candidate_host_species":5,
        "N_scored_other_plant_species":nb,
        "N_scored_fraction_of_all_indexed_plant_species":.6,
        "mean_N_candidate_hosts":host,"mean_N_other_scored_plants":bg}

def test_same_sites_and_times_background_differencing():
    a="early_2001_2010";b="later_2011_2024"
    obs=[site("port",a,7,5),site("nonport",a,5,6),
         site("port",b,7.5,6),site("nonport",b,5.5,6.5)]
    early=next(x for x in m.pair_rows(obs) if x["period"]==a)
    assert early["host_N_and_other_N_supported"]
    assert early["ports_minus_nonport_background_N"]==-1
    assert early["ports_minus_nonport_target_host_N"]==2
    assert early["difference_in_differences_hosts_minus_other_N"]==3
    later=m.analyze(obs)["identical_pair_period_difference"]
    assert later["shared_pairs"]==1
    assert math.isclose(later["host_minus_background_N"]["island_equal_weight_mean"],-.5)
    missing=m.pair_rows([site("port",a,7,5,nb=1),site("nonport",a,5,6)])
    assert missing[0]["difference_in_differences_hosts_minus_other_N"] is None
