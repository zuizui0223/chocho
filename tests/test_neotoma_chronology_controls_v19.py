"""Regression tests for chronology source only; labelled ages never certify LGM trees."""
import json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from fetch_neotoma_Alnus_chronology_controls_v19 import (
    response_scan,source_gate,fetch,CHRONID_KEYS
)

def test_chronology_object_reads_original_identity_and_age_control_metadata_only():
    r={"data":[{"chronologyid":13034,"chronologyname":"historic author fit",
                 "ageboundyounger":0,"ageboundolder":27000,
                 "chroncontrols":[{"chroncontrolid":5,"depth":300,"age":20000,
                                  "agelimityounger":19000,"agelimitolder":21000}]}]}
    x=response_scan(r)
    assert x['explicit_chronologyids']==['13034']
    assert x['chroncontrol_container_count']==1
    assert x['chroncontrol_record_count']==1
    assert x['reliability_bound_key_count']==2
    assert x['scanned_nodes']>=2

def test_neotoma_source_without_id_not_verified_by_an_echoed_url():
    fake={"data":[{"message":"Chronology not found","request":"/chronologies/987654321"}]}
    q=response_scan(fake)
    assert q['explicit_chronologyids']==[]
    assert q['chroncontrol_record_count']==0

def test_second_chronology_same_sample_not_independent_sample():
    from collections import Counter
    expected={"parent_source":{
      "source_dataset_to_all_chronology_ids":{"4396":[30953],"14556":[7725,7726],
      "22662":[13010,13011],"22683":[13033,13034],
      "24197":[14551,14552],"41351":[25461,25462]},
      "source_dataset_to_labelled_candidate_chronology_id":{
        "14556":7725,"22662":13010,"22683":13034,"24197":14552,"41351":25462}},
      "status":"SOURCE_CHRONOLOGY_PREFLIGHT_FROZEN_BEFORE_CONTROL_EVIDENCE"}
    age=[]
    per={"14556":9,"22662":1,"22683":1,"24197":3,"41351":7}
    for did,n in per.items():
        for index in range(n):
            age.append({"datasetid":did,"sampleid":f"{did}0{index}",
                 "source_chronologyid":str(expected["parent_source"]["source_dataset_to_labelled_candidate_chronology_id"][did]),
                 "source_agetype_label":"Calibrated radiocarbon years BP",
                 "datum_positive":"True","source_labelled_18_24ka_cal_BP_candidate":"True"})
    while len(age)<788:
        age.append({"datasetid":"4396","sampleid":str(len(age)),
                    "datum_positive":"True","source_labelled_18_24ka_cal_BP_candidate":"False"})
    source={"source_original_datasets":6,"source_labelled_candidate_age_rows_n":21,
            "accepted_21ka_Alnus_host_occupied_population_sites":0,
            "source_sha256":{k:"fixed" for k in expected["parent_source"]["source_dataset_to_all_chronology_ids"]}}
    got,ids=source_gate(expected,age,source)
    assert len(got)==21 and len({(x['datasetid'],x['sampleid']) for x in got})==21
    assert len({i for v in ids.values() for i in v})==11
    age[0]["source_chronologyid"]="7726"
    with pytest.raises(ValueError):source_gate(expected,age,source)

class FakeResponse:
    def __init__(self,body,ctype="application/json",code=200):
        self._body=body
        self.headers={"Content-Type":ctype}
        self.status_code=code
        self.url="https://api.neotomadb.org/v2.0/data/chronologies/13034"
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def iter_content(self,chunk_size=65536):yield self._body

class FakeSession:
    def __init__(self,body,ctype="application/json",code=200):
        self.body=body;self.ctype=ctype;self.code=code
    def get(self,url,**kwargs):return FakeResponse(self.body,self.ctype,self.code)

def test_only_json_source_bodies_and_no_silent_error_to_absence():
    data=json.dumps({"data":[{"chronologyid":13034,"chroncontrols":[]}]}).encode()
    meta,raw=fetch(FakeSession(data),"https://api.neotomadb.org/v2.0/data/chronologies/13034",10000,(2,5),1)
    assert meta["status"]=="SOURCE_RECEIVED_UNCLASSIFIED"
    assert meta["source_structure"]["explicit_chronologyids"]==["13034"]
    assert raw==data
    bad,raw=fetch(FakeSession(b"<html>error</html>","text/html"),"https://test.invalid/",10000,(2,5),1)
    assert bad["status"]=="SOURCE_UNAVAILABLE_NOT_CHRONOLOGY_ABSENCE"
    assert raw is None

def test_false_source_response_with_id_but_no_real_controls_not_independent_calibration():
    x=response_scan({"data":[{"chronologyid":13034}]})
    assert x["explicit_chronologyids"]==["13034"]
    assert x["chroncontrol_record_count"]==0
    assert x["reliability_bound_key_count"]==0
