"""Offline validation for the port versus non-port botanical GBIF pilot."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_nonport_botany.py"
spec=importlib.util.spec_from_file_location("coast_botany",SCRIPT)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fake_pair(island,port,records):
    return {
        "island_id":str(island),"port_id":port,"status":"PAIRED",
        "prior_2001_2010_5km_indexed_botanical_records":records,
        "port_sample_lat":0,"port_sample_lon":0,
        "nonport_lat":0.1,"nonport_lon":0.1,
    }


def test_pilot_selection_stratifies_but_does_not_duplicate_islands():
    x=[fake_pair(i,i,10*i) for i in range(1,37)]
    chosen=mod.choose_pilot(x,12)
    assert len(chosen)==12
    assert len({p["island_id"] for p in chosen})==12
    assert len({p["port_id"] for p in chosen})==12
    assert mod.choose_pilot(x,12)==chosen
    assert len(mod.choose_pilot(x,0))==36
    assert len(mod.choose_pilot(x,-1))==36


def test_plant_summary_is_species_presence_not_record_count():
    data={"indexed_plant_records":30,"indexed_species_keys":3,
          "records_with_speciesKey":28,"faceted_specieskey_counts":{1:18,2:1,3:9}}
    traits={1:{"N":9},2:{"N":1}}
    o=mod.plant_summary(data,traits)
    assert o["candidate_hosts_with_N_indexed_species"]==2
    assert o["candidate_hosts_N_high_7_9"]==1
    assert o["candidate_hosts_N_low_1_3"]==1
    assert o["candidate_hosts_mean_N_unweighted_by_occurrences"]==5
    assert o["candidate_hosts_N_high_fraction"]==0.5
    assert o["not_observed_does_not_mean_absent"]


def test_paired_data_sufficiency_is_strict_and_auditable():
    def make(side,period,n,records,mean):
        return {"island_id":"100","port_id":24,"period":period,
                "site_type":side,"all_plant_records":records,
                "candidate_hosts_with_N_indexed_species":n,
                "candidate_hosts_mean_N_unweighted_by_occurrences":mean}
    valid="early_2001_2010"
    later="later_2011_2024"
    receipt={"observations":[make("port",valid,5,100,7),make("nonport",valid,4,90,6),
                             make("port",later,5,120,8),make("nonport",later,1,15,4)]}
    o=mod.summarise(receipt,min_plant_records=30,min_trait_species=2)
    assert o["periods"][valid]["pairs_trait_comparable"]==1
    assert o["periods"][valid]["island_weighted_mean_N_difference"]==1
    assert o["periods"][valid]["independent_islands_trait_comparable"]==1
    assert o["periods"][later]["pairs_trait_comparable"]==0
    assert o["periods"][later]["island_weighted_mean_N_difference"] is None
    assert not o["periods"][valid]["data_gate_passed"]
