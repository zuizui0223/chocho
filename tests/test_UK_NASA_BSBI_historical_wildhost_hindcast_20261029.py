"""No-network causal firewall and pure covariate tests."""
import math,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_UK_NASA_BSBI_historical_wildhost_hindcast_20261029 import (
    exposure,unrelated_plant_portfolios,covariates,site_species_weights)

def toy():
    return {
      "pre1970_observed_bsbi_plant_names_by_site":{"S1":frozenset(("Urtica dioica","Lotus corniculatus")),"S2":frozenset(("Lotus corniculatus",))},
      "exact_source_host_plant_names_by_butterfly":{"Aglais urticae":frozenset(("Urtica dioica",)),"Polyommatus icarus":frozenset(("Lotus corniculatus",))},
      "all_source_plant_name_pool":frozenset(("Urtica dioica","Lotus corniculatus","Rhamnus cathartica","Trifolium pratense")),
      "original_pre1970_site_botanical_effort":{"S1":430,"S2":390},
    }

def test_original_host_positive_and_unknown_absence_not_confused():
    d=toy()
    assert exposure(d,"S1","Aglais urticae")==1
    assert exposure(d,"S2","Aglais urticae")==0
    assert d["original_pre1970_site_botanical_effort"]["S2"]>0

def test_loose_source_nonhost_placebo_never_contains_original_known_larval_hosts():
    d=toy()
    p=unrelated_plant_portfolios(d,1)
    assert not p["Aglais urticae"] & d["exact_source_host_plant_names_by_butterfly"]["Aglais urticae"]
    assert len(p["Aglais urticae"])==len(d["exact_source_host_plant_names_by_butterfly"]["Aglais urticae"])
    assert p==unrelated_plant_portfolios(d,1)

def test_site_species_weight_and_training_effort_covariate_no_outcome_used():
    d=toy()
    raw=[
      ("S1","Aglais urticae","early",1,[1.,2.],[4.,5.]),
      ("S1","Aglais urticae","early",0,[1.,2.],[4.,5.]),
      ("S2","Aglais urticae","early",0,[1.,2.],[4.,5.])
    ]
    b=covariates(raw,[0,1,2],d,None)
    m=covariates(raw,[0,1,2],d,"own")
    assert b.shape==(3,5) and m.shape==(3,6)
    assert m[0,-1]==pytest.approx(math.log(2))
    assert m[2,-1]==0.
    w=site_species_weights(raw,[0,1,2])
    assert sum(w)==pytest.approx(3.)
    assert w[0]==w[1] and w[2]>w[1]
