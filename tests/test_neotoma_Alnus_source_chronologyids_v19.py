"""Tests for Neotoma frozen exact 21 candidate source chronology IDs and age-type metadata."""
import sys,json
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_Alnus_21_frozen_chronologyids_v19 import (
    source_candidates,parse_chronologies,source_chron_meta,fetch_chronologies,norm_label
)

DATA={"14556":["7725"],"22662":["13010"],"22683":["13034"],
      "24197":["14552"],"41351":["25462"]}

def contract():
    return {"sources":{"required_source_21_candidate_chron_ids":DATA}}

def synth_rows():
    src=[]
    for did,cs in DATA.items():
        num={"14556":9,"22662":1,"22683":1,"24197":3,"41351":7}[did]
        for j in range(num):
            src.append({"datasetid":did,"sampleid":str(j+1),
             "source_labelled_18_24ka_cal_BP_candidate":"True",
             "datum_positive":"True","source_agetype_label":"Calibrated radiocarbon years BP",
             "source_age_within_18_24ka":"True","source_chronologyid":cs[0]})
    src.extend([{"datasetid":"4396","sampleid":str(i),
                "source_labelled_18_24ka_cal_BP_candidate":"False"} for i in range(767)])
    assert len(src)==788
    return src

def test_original_21_ages_and_five_chronologyids_frozen_without_age_selection():
    a=source_candidates(synth_rows(),contract())
    assert len(a)==21
    assert set(x["datasetid"] for x in a)==set(DATA)
    a[0]["source_chronologyid"]="9999"
    src=synth_rows()
    src[0]["source_chronologyid"]="9999"
    with pytest.raises(ValueError):source_candidates(src,contract())
    src=synth_rows()
    src[0]["datum_positive"]="False"
    with pytest.raises(ValueError):source_candidates(src,contract())

def test_official_chronology_parser_accepts_matching_id_and_rejects_cross_datasets():
    rec=[{"chronology":{"chronologyid":7725,"datasetid":14556,
         "agetype":"Calibrated radiocarbon years BP","controls":[{"controltype":"Radiocarbon"}]}}]
    out=parse_chronologies(rec,"14556")
    assert set(out)=={"7725"}
    assert source_chron_meta(out["7725"])["source_chronology_controls_total"]==1
    bad=json.loads(json.dumps(rec));bad[0]["chronology"]["datasetid"]=9999
    with pytest.raises(ValueError):parse_chronologies(bad,"14556")
    with pytest.raises(ValueError):parse_chronologies(rec+rec,"14556")
    with pytest.raises(ValueError):parse_chronologies([{"chronology":{}}],"14556")

def test_chronology_agetype_conflict_cannot_be_silently_calibrated():
    q=source_chron_meta({"chronologyid":30953,
       "agetype":"Radiocarbon years BP",
       "controls":[{"controltype":"Radiocarbon",
           "geochronagetype":"Radiocarbon years BP",
           "controllimits":{"ageolder":24000,"ageyounger":22000}}]})
    assert q["chronology_source_agetype"]=="Radiocarbon years BP"
    assert q["source_chronology_controls_with_numerical_bounds"]==1
    assert norm_label(q["chronology_source_agetype"])!="calibrated radiocarbon years bp"
    assert q["original_model_independent_calibration_not_verified"] is True

class Response:
    def __init__(self,obj,status=200,kind="application/json"):
        self.data=json.dumps(obj).encode()
        self.status_code=status
        self.headers={"content-type":kind,"content-length":str(len(self.data))}
        self.url="https://api.neotomadb.org/v2.0/data/datasets/14556/chronologies"
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def raise_for_status(self):
        if self.status_code!=200:raise ValueError("fake upstream error")
    def iter_content(self,chunk_size):yield self.data
class Session:
    def __init__(self,r):self.r=r
    def get(self,*args,**kwargs):return self.r

def test_source_fetch_returns_data_without_promoting_ancient_plant_occupancy():
    src={"data":[{"chronology":{"chronologyid":7725,
          "agetype":"Calibrated radiocarbon years BP"}}]}
    result,body,data=fetch_chronologies(Session(Response(src)),"http://test",100000,(2,3))
    assert result["status"]=="SOURCE_JSON_AVAILABLE"
    assert parse_chronologies(data,"14556")["7725"]["agetype"].startswith("Calibrated")
    assert body is not None
    failure,body,data=fetch_chronologies(Session(Response(src,kind="text/html")),
        "http://test",100000,(2,3))
    assert failure["status"]=="SOURCE_UNAVAILABLE_NOT_FOSSIL_ABSENCE" and body is None

def test_v19_never_promotes_sample_chronology_link_to_occupied_host():
    s=(Path(__file__).resolve().parents[1]/"scripts"/
       "audit_neotoma_Alnus_21_frozen_chronologyids_v19.py").read_text()
    assert '"source_glacial_occupied_Alnus_populations_accepted":0' in s
    assert '"sample_age_models_independently_calibrated_n":0' in s
    assert '"insect_genetic_responses_read":False' in s
