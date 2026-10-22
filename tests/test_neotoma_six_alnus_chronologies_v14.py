"""Original 6-source Alnus pollen/Neotoma chronology gate and metadata source tests."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_six_alnus_chronologies_v14 import gate_source,run,fetch_metadata
ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"docs/supporting/lgm_neotoma_Alnus_six_dataset_chronologies_v14_contract.json"
PREVIOUS=ROOT/"results/lgm_neotoma_idage_18_24ka_source_candidates_v13"

def test_six_datasets_fixed_from_original_27_taxon_and_sample_age_filtered_source_rows():
    q,rows=gate_source(CONTRACT,PREVIOUS/"bounded_18_24ka_sample_taxon_proxy_checks.csv",
                       PREVIOUS/"summary.json")
    assert len(rows)==27
    assert q["source"]["unique_source_occids"]==26
    assert {x["datasetid"] for x in q["source"]["candidate_dataset_sites"]}=={
        4396,14556,22662,22683,24197,41351}
    assert {x["siteid"] for x in q["source"]["candidate_dataset_sites"]}=={
        3392,9834,14474,14485,15750,23855}

class Fake:
    status_code=200
    headers={"Content-Type":"application/json"}
    url="https://api.neotomadb.org/source"
    def __init__(self,payload):self.obj=payload;self.content=json.dumps(payload).encode()
    def json(self):return self.obj
    def raise_for_status(self):pass
class FakeSession:
    def get(self,url,**kwargs):
        if "9999999" in url:return Fake({"data":[]})
        if "chronologies" in url:
            return Fake({"data":[{"chronid":234,"chronologyname":"Age model unknown"}]})
        if "publications" in url:
            return Fake({"data":[{"publicationid":222,"doi":None}]})
        return Fake({"data":[{"datasetid":1,"datasettype":"pollen"}]})

def test_dataset_chronology_publications_are_metadata_not_ancient_occurrence():
    spec=json.loads(CONTRACT.read_text())
    rows,q=run(FakeSession(),spec)
    assert len(rows)==19
    assert q["candidate_source_datasets"]==6 and q["candidate_source_sites"]==6
    assert q["metadata_coverage"]["chronologies"]["http_available"]==6
    assert q["metadata_coverage"]["publications"]["http_available"]==6
    assert q["metadata_coverage"]["dataset"]["http_available"]==6
    assert q["sample_ages_independently_chronology_validated"] is False
    assert q["pollen_morphology_species_level_validated"] is False
    assert q["original_21ka_host_plant_populations_inferred"] is False

def test_invalid_dataset_chronology_is_not_fossil_absence():
    _,q=run(FakeSession(),json.loads(CONTRACT.read_text()))
    assert q["status"]=="CHRONOLOGY_AND_PUBLICATION_METADATA_AVAILABILITY_AUDITED_NO_FOSSIL_OCCUPANCY"
    assert q["genetic_response_used"] is False
