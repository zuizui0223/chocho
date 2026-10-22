"""Neotoma v2 taxon-ID path sample/chronology response field diagnostics."""
import json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_neotoma_taxonid_sample_age_v12 import (
    frozen_contract,source_http,fossil_age_site_preview,value_shape,run_probe
)
ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"docs/supporting/lgm_neotoma_taxonid_age_source_v12_contract.json"
PREVIOUS=ROOT/"results/lgm_neotoma_taxon_filter_integrity_v11/source_filter_integrity_summary.json"

def test_v11_source_has_definitively_invalid_named_occurrence_filter():
    q=frozen_contract(CONTRACT,PREVIOUS)
    assert q["status"]=="SOURCE_SAMPLE_LEVEL_CHRONOLOGY_PRECHECK_NOT_LGM_PRESENCE"
    assert q["exact_taxon_records"]==[
        {"accepted":"Alnus glutinosa","neotoma_taxonid":1299},
        {"accepted":"Betula pendula","neotoma_taxonid":3461}]
    assert len(q["generic_endpoint_negative_controls"])==4

def test_age_record_shapes_are_not_calibrated_fossil_species_proofs():
    z=fossil_age_site_preview({"occid":324,
         "age":{"ageold":24200,"ageyoung":18450,"model":None},
         "sample":{"sampleid":19},"site":{"siteid":120}})
    assert z["occurrence_id"]=="324"
    assert z["sample_id"]==19 and z["site_id"]==120
    assert z["age_shape"]["kind"]=="object"
    assert "ageold" in z["age_shape"]["keys"]
    assert "calibrated_21ka_species_occupancy" not in z

class FakeResponse:
    status_code=200
    headers={"Content-Type":"application/json"}
    url="https://api.neotomadb.org/mock"
    def __init__(self,v): self.v=v;self.content=json.dumps(v).encode()
    def raise_for_status(self): pass
    def json(self):return self.v
class FakeSession:
    def get(self,url,params,**kwargs):
        if "/taxa/1299/" in url:
            return FakeResponse({"data":[{"occid":123,"age":{"age":1100},"site":{"siteid":1}}]})
        if "/taxa/3461/" in url:
            return FakeResponse({"data":[{"occid":456,"age":{"age":19000},"site":{"siteid":2}}]})
        return FakeResponse({"data":[{"occid":99999,"age":{"age":19000}}]})

def test_taxonid_paths_distinguish_occurrence_ids_while_generic_filters_fail():
    q=json.loads(CONTRACT.read_text())
    outputs,preview,report=run_probe(FakeSession(),q)
    assert len(outputs)==8 and len(preview)==4
    assert report["status"]=="TAXON_ID_PATH_DIFFERENT_OCCURRENCE_RECORDS_CHRONOLOGY_UNVERIFIED"
    assert len(report["generic_negative_controls_identical_payload_pairs"])==6
    assert report["source_any_original_species_18_24ka_fossil_occupancy_verified"] is False
    assert report["source_direct_occurrence_ids_covered_by_bounded_pages"]=={
        "Alnus glutinosa":1,"Betula pendula":1}

def test_malformed_json_without_data_list_is_source_failure_not_absence():
    x,y=source_http(FakeSession(),"https://api.neotomadb.org/v2.0/data/taxa/1299/occurrences",{})
    assert y["status"]=="SOURCE_HTTP200_NOT_FOSSIL_VALIDATION"
    assert x[0]["occid"]==123
