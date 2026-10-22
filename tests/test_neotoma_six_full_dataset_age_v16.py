"""Source-only Neotoma full-dataset chronology proof gates.

Synthetic tests prohibit treating HTTP-200, dataset-level age ranges or any
co-located but unassigned chronology IDs as confirmed 21 ka plant occupancy.
"""
from pathlib import Path
import sys,json,hashlib
from contextlib import AbstractContextManager
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from neotoma_six_full_dataset_download_provenance_v16 import (
    tree_field_summary,query_source,source_gate,EXPECTED
)

def test_explicit_source_taxon_and_sample_chronology_structure_distinct():
    source={"data":[{"datasetid":14556,
       "chronologies":[{"chronologyid":200,"agetype":"Radiocarbon years BP"},
                       {"chronologyid":201,"agetype":"Calibrated radiocarbon years BP"}],
       "samples":[{"sampleid":12,"taxonid":1299,"age":20456},
                  {"sampleid":13,"chronologyid":201,"occid":2164809,
                   "age":{"age":20456,"ageolder":20900,"ageyounger":20000}}]}]}
    q=tree_field_summary(source["data"])
    assert q["all_explicit_datasetids"]==["14556"]
    assert set(q["all_explicit_chronologyids"])=={"200","201"}
    assert q["all_explicit_occurrenceids"]==["2164809"]
    assert set(q["all_explicit_sampleids"])=={"12","13"}
    assert q["colocated_sample_chronology_count"]==1
    # No result returns confirmed chronology model/ancient occupancy.
    assert "sample_age_calibration_verified" not in q

def test_pollen_taxon_string_without_exact_dataset_id_is_not_validated():
    ambiguous={"data":[{"datasettype":"pollen","taxonname":"Alnus glutinosa",
         "site":{"siteid":9834},"sample":{"age":21000}}]}
    q=tree_field_summary(ambiguous)
    assert "14556" not in q["all_explicit_datasetids"]
    assert q["all_explicit_occurrenceids"]==[]
    assert q["all_explicit_chronologyids"]==[]

def test_no_artificial_chronology_from_parentage_or_arbitrary_age_number():
    q=tree_field_summary([{"chronology":{"chronologyid":110,"agetype":"cal BP"},
                           "sample":{"sampleid":10,"age":20500}}])
    assert q["colocated_sample_chronology_count"]==0

def test_original_six_dataset_ids_and_26_source_occurrences_are_exactly_pinned():
    contract={"status":"SOURCE_ONLY_DETERMINISTIC_FROZEN_BEFORE_FULL_DATASET_DOWNLOADS",
              "source":{"frozen_datasetids":list(map(int,EXPECTED))}}
    candidates=[]
    for i,(did,sid) in enumerate(EXPECTED.items()):
        for j in range(4 if i<3 else 5):
            candidates.append({"source_datasetid":did,"source_siteid":sid,
                "occid":str(100+i*10+j),"source_sample_taxonid":"1299",
                "source_sample_taxonname":"Alnus glutinosa","source_datasettype":"pollen"})
    # 12 + 15 =27, with exactly 26 IDs after duplicating one within a dataset.
    candidates[-1]["occid"]=candidates[-2]["occid"]
    flags=[{"occid":x["occid"]} for x in candidates]
    lookup=source_gate(contract,candidates,flags)
    assert len(lookup)==6 and sum(map(len,lookup.values()))==26
    candidates[0]["source_sample_taxonid"]="3461"
    with pytest.raises(ValueError):source_gate(contract,candidates,flags)

class Response:
    def __init__(self,payload,mime="application/json",status=200):
        self.content=json.dumps(payload).encode()
        self.headers={"Content-Type":mime,"Content-Length":str(len(self.content))}
        self.status_code=status
        self.url="https://api.neotomadb.org/v2.0/data/downloads/4396"
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def raise_for_status(self):
        if self.status_code!=200:
            raise ValueError("HTTP fake failure")
    def iter_content(self,chunk_size):
        yield self.content

class Session:
    def __init__(self,result):self.result=result
    def get(self,url,**kwargs):return self.result

def test_neotoma_http200_response_requires_parseable_data():
    payload={"data":[{"datasetid":4396,"samples":[{"sampleid":10}]}]}
    response,body,data=query_source(Session(Response(payload)),
        "https://api.neotomadb.org/v2.0/data/downloads/4396",100000,(2,3))
    assert response["status"]=="HTTP_JSON_SOURCE_RECEIVED"
    assert response["source_sha256"]==hashlib.sha256(body).hexdigest()
    assert "4396" in tree_field_summary(data)["all_explicit_datasetids"]

def test_generic_http200_metadata_without_data_never_proves_fossil():
    response,body,data=query_source(Session(Response({"success":True,"results":[{"taxonname":"Alnus glutinosa"}]})),
        "https://api.neotomadb.org/v2.0/data/downloads/4396",100000,(2,3))
    assert response["status"]=="SOURCE_UNAVAILABLE_OR_INVALID_NOT_FOSSIL_ABSENCE"
    assert body is None and data is None

def test_response_size_cap_does_not_infer_absence():
    response,body,data=query_source(Session(Response({"data":[{"datasetid":4396}]})),
        "https://api.neotomadb.org/v2.0/data/downloads/4396",4,(2,3))
    assert response["status"]=="SOURCE_UNAVAILABLE_OR_INVALID_NOT_FOSSIL_ABSENCE"
    assert body is None

def test_missing_or_unresolved_source_never_claims_21ka_occupied_host():
    text=(Path(__file__).resolve().parents[1]/"scripts"/
          "neotoma_six_full_dataset_download_provenance_v16.py").read_text()
    assert '"Alnus_21ka_fossil_site_presence_accepted":0' in text
    assert '"sample_assigned_age_depth_models_verified":0' in text
    assert '"historical_host_occurrence_or_genetic_response_used":False' in text
