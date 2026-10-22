"""Tests: raw Neotoma sample taxon datum, multiple original age models and calibrated label boundaries."""
import sys,gzip,json,hashlib
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_direct_Alnus_source_sampleages_v18 import (
    source_samples,source_pollen_datum,classify_source_age,validate_source,STATUS)

DATASETS={"4396":"3392","14556":"9834","22662":"14474",
          "22683":"14485","24197":"15750","41351":"23855"}

def contract(sourcehash):
    return {"status":"PREDECLARED_SAMPLE_LEVEL_SOURCE_CHECK_NO_ACCEPTED_LGM_OCCUPANCY",
        "source":{"dataset_sha256":sourcehash,"dataset_site":DATASETS,
                  "taxonid":1299,"source_taxon_label":"Alnus glutinosa",
                  "age_min_calibrated_BP":18000,"age_max_calibrated_BP":24000}}

def sample(sid,age=21234,count=3):
    return {"sampleid":sid,"analysisunitid":10000+sid,"depth":float(sid),
        "datum":[{"taxonid":1299,"value":count},{"taxonid":7777,"value":45}],
        "ages":[{"age":age,"chronologyid":900,"agetype":"Calibrated radiocarbon years BP",
                 "ageolder":22000,"ageyounger":20000},
                {"age":18500,"chronologyid":901,"agetype":"Radiocarbon years BP",
                 "ageolder":None,"ageyounger":None}]}

def full_source(did,sit,samples):
    return {"data":[{"site":{"siteid":int(sit),"collectionunit":{
        "dataset":{"datasetid":int(did),"datasettype":"pollen",
                   "samples":samples}}}}]}

def test_sample_is_exact_source_dataset_and_not_neotoma_occurrenceid():
    q=full_source("4396","3392",[sample(1)])
    out=source_samples(json.dumps(q).encode(),"4396","3392")
    assert out[0]["sampleid"]==1
    with pytest.raises(ValueError):
        source_samples(json.dumps(q).encode(),"14556","3392")
    with pytest.raises(ValueError):
        source_samples(json.dumps(q).encode(),"4396","9834")

def test_positive_datum_must_be_numeric_and_gt_zero():
    assert source_pollen_datum(sample(1,count=12))["total_positive_numeric_datum_value"]==12
    assert source_pollen_datum(sample(1,count=0))["source_positive_numeric_datums"]==0
    assert source_pollen_datum(sample(1,count=None))["source_nonnumeric_or_negative_datums"]==1
    assert source_pollen_datum(sample(1,count="trace"))["source_nonnumeric_or_negative_datums"]==1
    assert source_pollen_datum(sample(1,count=-1))["source_nonnumeric_or_negative_datums"]==1
    assert source_pollen_datum({"datum":[{"taxonid":111,"value":40}]})["source_taxon_1299_datum_entries"]==0
    with pytest.raises(ValueError):
        source_pollen_datum({"datum":None})

def test_calibrated_label_and_explicit_chronology_are_both_required():
    c=contract({})
    old=classify_source_age(sample(1)["ages"][0],c)
    assert old["source_labelled_18_24ka_cal_BP_candidate"] is True
    assert old["has_both_source_age_bounds"]
    raw=classify_source_age(sample(1)["ages"][1],c)
    assert raw["source_labelled_18_24ka_cal_BP_candidate"] is False
    without_chron=sample(1)["ages"][0].copy()
    without_chron["chronologyid"]=None
    assert classify_source_age(without_chron,c)["source_labelled_18_24ka_cal_BP_candidate"] is False
    without_bounds=without_chron.copy();without_bounds["chronologyid"]=900
    without_bounds["ageyounger"]=None
    result=classify_source_age(without_bounds,c)
    assert result["source_labelled_18_24ka_cal_BP_candidate"] is True
    assert not result["has_both_source_age_bounds"]

def test_uncalibrated_18_24ka_numeric_age_is_not_a_candidate():
    q=classify_source_age({"age":21000,"agetype":"Radiocarbon years BP",
                           "chronologyid":222},contract({}))
    assert not q["source_labelled_18_24ka_cal_BP_candidate"]
    q=classify_source_age({"age":25000,"agetype":"Calibrated radiocarbon years BP",
                           "chronologyid":222},contract({}))
    assert not q["source_labelled_18_24ka_cal_BP_candidate"]
    with pytest.raises(ValueError):
        classify_source_age({"age":21000,"ageolder":18000,"ageyounger":24000,
            "agetype":"Calibrated radiocarbon years BP","chronologyid":222},contract({}))

def test_complete_six_source_dataset_inventory_cannot_choose_best_chronology(tmp_path):
    digests={}
    for did,sit in DATASETS.items():
        one=sample(1)
        two=sample(2,age=20125,count=0)
        raw=json.dumps(full_source(did,sit,[one,two])).encode()
        digests[did]=hashlib.sha256(raw).hexdigest()
        with gzip.open(tmp_path/f"dataset_{did}.json.gz","wb") as f:f.write(raw)
    plants,ages,datasets,q=validate_source(contract(digests),tmp_path)
    assert len(datasets)==6 and len(plants)==12 and len(ages)==24
    assert q["source_positive_datum_samples_total"]==6
    assert q["source_positive_datum_samples_with_labelled_cal_BP_18_24ka_n"]==6
    assert q["source_labelled_candidate_age_rows_n"]==6
    assert q["accepted_21ka_Alnus_host_occupied_population_sites"]==0
    assert q["status"]==STATUS
    # A duplicate source sampleid must not be counted twice under a single site.
    raw=json.dumps(full_source("4396","3392",[sample(1),sample(1)])).encode()
    digests["4396"]=hashlib.sha256(raw).hexdigest()
    with gzip.open(tmp_path/"dataset_4396.json.gz","wb") as f:f.write(raw)
    with pytest.raises(ValueError):
        validate_source(contract(digests),tmp_path)

def test_no_source_species_morphology_or_genetic_promotion():
    text=(Path(__file__).resolve().parents[1]/"scripts"/
          "audit_neotoma_direct_Alnus_source_sampleages_v18.py").read_text()
    assert '"accepted_21ka_Alnus_host_occupied_population_sites":0' in text
    assert '"source_Alnus_pollen_morphology_independently_validated":0' in text
    assert 'No join to original v13 occurrence ID' in text
