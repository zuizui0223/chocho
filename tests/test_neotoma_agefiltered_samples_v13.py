"""Neotoma sample taxon, fossil proxy and age selection strict source gates."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_neotoma_agefiltered_plant_samples_v13 import (
    check_sample,record_source_metadata,source_identity,CANDIDATE)
ROOT=Path(__file__).resolve().parents[1]

def row(taxid=1299,name="Alnus glutinosa",datasettype="plant macrofossil",age=21400):
    return {"occid":1234,"sample":{"taxonid":taxid,"taxonname":name},
            "site":{"datasettype":datasettype,"siteid":193,"datasetid":52,
                    "location":"POINT(-60.0 40.0)"},
            "age":{"age":age,"ageyounger":20000,"ageolder":22200}}

def test_candidate_fossil_source_never_becomes_ancient_occupancy():
    q=row()
    assert check_sample(q,1299,"Alnus glutinosa","plant macrofossil")==CANDIDATE
    s=record_source_metadata(q,1299,"Alnus glutinosa","plant macrofossil")
    assert s["source_grade"]==CANDIDATE
    assert s["age_model_calibration_independently_verified"] is False
    assert s["pollen_species_level_identification_independently_verified"] is False
    assert s["source_datasetid"]==52 and s["source_siteid"]==193

def test_species_and_datasettype_negative_controls_fail_closed():
    assert check_sample(row(taxid=3461),1299,"Alnus glutinosa","plant macrofossil")=="WRONG_TAXON_ID_FILTER_FAILURE"
    assert check_sample(row(name="Alnus"),1299,"Alnus glutinosa","plant macrofossil")=="NOT_EXACT_SPECIES_IN_SOURCE"
    assert check_sample(row(datasettype="pollen"),1299,"Alnus glutinosa","plant macrofossil")=="WRONG_DATASETTYPE_FILTER_FAILURE"
    assert check_sample(row(age=1000),1299,"Alnus glutinosa","plant macrofossil")=="SAMPLE_CENTER_AGE_NOT_18_24KA"
    assert check_sample(row(age=None),1299,"Alnus glutinosa","plant macrofossil")=="UNKNOWN_SAMPLE_AGE_DO_NOT_INFER_ABSENCE"

def test_original_source_verified_numeric_taxonids_and_age_contract():
    c=ROOT/"docs/supporting/lgm_neotoma_idage_sample_candidates_v13_contract.json"
    v11=ROOT/"results/lgm_neotoma_taxon_filter_integrity_v11/source_filter_integrity_summary.json"
    v12=ROOT/"results/lgm_neotoma_taxonid_sample_age_v12/summary.json"
    spec=source_identity(c,v11,v12)
    assert spec["status"]=="BOUNDED_SOURCE_SAMPLE_VERIFICATION_NOT_LGM_OCCUPANCY"
    assert len(spec["cases"])==2
    assert spec["queries"]["limit"]==30
    assert spec["queries"]["datasettypes"]==["plant macrofossil","pollen"]
    assert spec["queries"]["ageold"]==24000
