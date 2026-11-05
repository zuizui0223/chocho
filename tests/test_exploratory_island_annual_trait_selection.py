"""Source-only plant trait definition and island-level population checks."""
from __future__ import annotations
import pytest
for dep in ("numpy","pandas","scipy","pyreadr"):
    pytest.importorskip(dep)
import importlib.util,sys
from pathlib import Path
import numpy as np
import pandas as pd

file=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_island_annual_trait_selection.py"
spec=importlib.util.spec_from_file_location("island_trait",file)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_exact_accepted_species_only_and_annual_biennial_separation():
    assert m.as_bins("Poa annua")=="poa annua"
    assert m.as_bins("Poa annua subsp. maritima") is None
    assert m.as_bins("Rosa × hybrid") is None
    assert m.as_bins("Poa annua L.") is None
    species=[{"work_ID":1,"work_species":"Poa annua"},
             {"work_ID":2,"work_species":"Plantago major"},
             {"work_ID":3,"work_species":"Secale cereale"}]
    raw={"lifecycle":[{"work_ID":1,"trait_value":"annual","agreement":.99},
                      {"work_ID":2,"trait_value":"perennial","agreement":.99},
                      {"work_ID":3,"trait_value":"annual","agreement":.50}],
         "growth_form":[{"work_ID":1,"trait_value":"herb","agreement":1},
                        {"work_ID":2,"trait_value":"herb","agreement":1},
                        {"work_ID":3,"trait_value":"herb","agreement":1}]}
    mapping,counts,coverage=m.traits_from_gift(species,raw)
    assert mapping["poa annua"]=={"lifecycle":"annual","growth_form":"herb"}
    assert mapping["plantago major"]["lifecycle"]=="perennial"
    assert "lifecycle" not in mapping["secale cereale"]
    assert counts["lifecycle_low_agreement"]==1

def test_wcvp_source_rejects_conflicting_native_alien_region_status():
    names=pd.DataFrame({"plant_name_id":["1","2"],"taxon_key":["poa annua","plantago major"]})
    d=pd.DataFrame({"plant_name_id":["1","1","2"],
        "area_code_l3":["ALD","ALD","ALD"],
        "introduced":[0,1,0],"extinct":[0,0,0],"location_doubtful":[0,0,0]})
    with pytest.raises(ValueError,match="Conflicting"):
        m.source_table(names,d,{"poa annua":{"lifecycle":"annual"}},{"ALD"})

def test_trait_unknown_is_not_negative_and_port_permutation_requires_support():
    recs=[]
    for k in range(10):
        recs.append({"code":str(k),"level1":"1","log2ports":float(k),
             "log_area":float(k),"GBIF_allplant_indexed":100.0,
             "introduced_N":40,"native_N":100,
             "introduced_annual_fraction_known":k/10,
             "native_annual_fraction_known":(10-k)/10,
             "introduced_lifecycle_coverage":.5,
             "native_lifecycle_coverage":.5})
    df=pd.DataFrame(recs)
    r=m.port_perm(df,"introduced_annual_fraction_known",
          ["introduced_lifecycle_coverage","native_lifecycle_coverage"],B=29)
    assert r["status"]=="INSUFFICIENT_EXCHANGEABLE_GLOBAL_TRAIT_SUPPORT"
