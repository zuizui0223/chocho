"""Neotoma v2 source integrity tests: no fossil claim from HTTP 200 or fake filters."""
import csv,hashlib,json,sys
from pathlib import Path
import pytest
import requests
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_taxon_filter_integrity_v11 import (
    original_replay,strict_registry_taxon_id,audit_live,
    clean_json_response,live_request,NAMES,BAD
)

ROOT=Path(__file__).resolve().parents[1]

def test_original_15_taxa_30_queries_have_only_two_unique_success_payloads():
    base=ROOT/"results/lgm_neotoma_15_plant_fossil_source_pilot_v01"
    x=original_replay(base/"all_15_taxa_source_query_status.csv",
                      base/"source_pilot_summary.json")
    assert x["prior_proxy_queries"]==30
    assert x["prior_http200_json"]==12
    assert x["unique_http200_payload_sha256"]==2
    assert sorted(a["query_count"] for a in x["repeated_payload_clusters"])==[5,7]
    assert all(a["distinct_taxon_names"]==a["query_count"]
               for a in x["repeated_payload_clusters"])
    assert x["not_a_confirmed_fossil_occurrence"] is True

def test_neotoma_taxon_registry_rejects_ambiguous_or_mismatched_name():
    correct={"data":[{"taxonname":"Alnus glutinosa","taxonid":123}]}
    wrong={"data":[{"taxonname":"Alnus sp.","taxonid":123}]}
    assert strict_registry_taxon_id(correct,"Alnus glutinosa")==123
    assert strict_registry_taxon_id(wrong,"Alnus glutinosa") is None
    assert strict_registry_taxon_id({"data":[]},"Alnus glutinosa") is None
    assert strict_registry_taxon_id({"data":[{"taxonname":"Betula pendula","taxonid":888}]},"Alnus glutinosa") is None

class FakeResp:
    status_code=200
    url="https://api.neotomadb.org/fake"
    headers={"Content-Type":"application/json"}
    def __init__(self,obj):
        self.payload=obj
        self.content=json.dumps(obj,sort_keys=True).encode()
    def json(self):return self.payload

class FakeSession:
    def get(self,url,params,**kwargs):
        name=params.get("taxonname")
        if "/data/taxa/" in url and "/occurrences" in url:
            return FakeResp({"data":[{"occid":1001 if "/101/" in url else 1002}]})
        if "/data/taxa" in url:
            if name==BAD:return FakeResp({"data":[]})
            return FakeResp({"data":[{"taxonname":name,"taxonid":101 if name==NAMES[0] else 102}]})
        # Simulate Neotoma taxonname query silently ignored by server:
        return FakeResp({"data":[{"occid":1},{"occid":2}]})

def test_negative_control_detects_fossil_api_name_filter_failure_even_on_http200():
    from pathlib import Path
    spec=json.loads((ROOT/"docs/supporting/lgm_neotoma_taxon_filter_integrity_v11_contract.json").read_text())
    rows,q=audit_live(FakeSession(),spec)
    assert q["status"]=="FOSSIL_TAXONNAME_FILTER_INVALID_SENTINEL_MATCHES_REAL_TAXON"
    assert q["success_occurrence_query_n"]==3
    assert len(q["identical_distinct_name_response_pairs"])==3
    assert q["fossil_species_21ka_sample_occurrence_verified"] is False
    assert q["exact_taxon_registry_ids"]=={"Alnus glutinosa":101,"Betula pendula":102}
    assert len(rows)==8

def test_no_data_or_invalid_json_does_not_verify_fossil_absence():
    p,meta=live_request(FakeSession(),"https://api.neotomadb.org/v2.0/data/taxa",{"taxonname":BAD})
    assert meta["status"]=="HTTP200_SOURCE_CONTENT_UNVERIFIED"
    assert meta["record_count"]==0
    assert meta["has_nonempty_data"] is False
    assert "fossil" not in meta

def test_api_record_identity_must_not_be_substituted_for_sample_age_or_species():
    info=clean_json_response(FakeResp({"data":[{"occid":1234,"age":{"age":21000}}]}))[1]
    assert info["record_count"]==1
    assert info["record_ids"]==["1234"]
    assert "sample_age_model_verified" not in info
