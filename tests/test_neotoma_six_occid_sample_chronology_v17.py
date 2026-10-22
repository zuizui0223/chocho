"""Tests for explicit Neotoma 26-occurrence / 6 full dataset source sample chronology crosswalk."""
import sys, json
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from neotoma_six_occid_exact_sample_chronology_v17 import (
    source_datasets,inspect_download,match_explicit_sample,gate_original,
    ids,SITES
)

def fixture_dataset():
    return {"data":[{"site":{"siteid":9834,"collectionunit":{
      "dataset":{"datasetid":14556,"samples":[
        {"sampleid":12,"analysisunitid":91,
         "ages":[{"age":20456,"agetype":"Radiocarbon years BP","chronologyid":100},
                 {"age":22311,"agetype":"Calibrated radiocarbon years BP","chronologyid":101}],
         "datum":[{"taxonid":1299,"value":14}]},
        {"sampleid":13,"analysisunitid":92,
         "ages":[{"age":20456,"agetype":"Calibrated radiocarbon years BP","chronologyid":101}],
         "datum":[{"taxonid":1299,"value":7}]}],
       "chronologies":[{"chronologyid":100},{"chronologyid":101}]}}}]}

def test_original_six_full_dataset_ids_must_align_exactly():
    d=source_datasets(fixture_dataset(),"14556","9834")
    assert len(d["samples"])==2
    with pytest.raises(ValueError):
        source_datasets(fixture_dataset(),"14557","9834")
    with pytest.raises(ValueError):
        source_datasets(fixture_dataset(),"14556","9999")

def test_same_numeric_age_not_a_valid_occurrence_join():
    d=source_datasets(fixture_dataset(),"14556","9834")
    occurrence={"occid":2164809,"sample":{"age":20456},"siteid":9834,"datasetid":14556}
    matches,fields=match_explicit_sample(occurrence,d)
    assert len(matches)==0 and fields["occurrence_explicit_sample_ids"]==[]
    assert fields["occurrence_explicit_analysisunit_ids"]==[]

def test_explicit_id_join_does_not_silently_calibrate_or_select_an_age_model():
    d=source_datasets(fixture_dataset(),"14556","9834")
    occurrence={"occid":2164809,"sample":{"sampleid":12},"site":{"siteid":9834}}
    matches,fields=match_explicit_sample(occurrence,d)
    assert len(matches)==1
    assert matches[0]["explicit_id_join_modes"]==["EXACT_SAMPLEID"]
    assert matches[0]["sample_chronologyids"]==["100","101"]
    assert len(matches[0]["sample_age_metadata"])==2
    assert matches[0]["candidate_pollen_taxon_1299_in_download_datum"]
    assert "age_model_independently_calibrated" not in matches[0]

def test_analysisunitid_alone_is_distinct_from_sampleid():
    d=source_datasets(fixture_dataset(),"14556","9834")
    occurrence={"analysisunitid":92,"ages":[{"age":20456}],"occid":21}
    matches,_=match_explicit_sample(occurrence,d)
    assert len(matches)==1 and matches[0]["analysisunitid"]=="92"
    assert matches[0]["explicit_id_join_modes"]==["EXACT_ANALYSISUNITID"]

def test_ambiguous_original_occurrence_cannot_be_auto_assigned():
    d=source_datasets(fixture_dataset(),"14556","9834")
    occurrence={"sampleids":[{"sampleid":12},{"sampleid":13}]}
    matches,_=match_explicit_sample(occurrence,d)
    assert len(matches)==2

def test_all_original_27_taxonid_1299_pollen_rows_and_6_datasets_fixed():
    rows=[]
    for j,(did,site) in enumerate(SITES.items()):
        for i in range(4 if j<3 else 5):
            rows.append({"source_datasetid":did,"source_siteid":site,
                "source_sample_taxonid":"1299","source_sample_taxonname":"Alnus glutinosa",
                "source_datasettype":"pollen","occid":str(10000+j*100+i)})
    rows[-1]["occid"]=rows[-2]["occid"]
    c={"status":"SOURCE_ONLY_FROZEN_SIX_FIRST_OCCURRENCE_ID_JOIN_TEST_BEFORE_NEW_HTTP"}
    select=gate_original(c,rows)
    assert len(select)==6
    for k in SITES:assert select[k]==str(10000+list(SITES).index(k)*100)
    rows[0]["source_sample_taxonname"]="Quercus robur"
    with pytest.raises(ValueError):gate_original(c,rows)

def test_raw_fossil_source_sample_object_reports_both_ages_not_one_chosen_model():
    d=source_datasets(fixture_dataset(),"14556","9834")
    r=inspect_download(d)
    assert r["source_full_samples_n"]==2
    assert r["source_sample_age_records_n"]==3
    assert r["source_unique_sample_ids_n"]==2
    assert r["source_distinct_age_chronology_id_n"]==2
    assert r["source_sample_taxon1299_with_datum_n"]==2

def test_zero_glacial_occupancy_and_genetic_response_claims_in_probe_code():
    s=(Path(__file__).resolve().parents[1]/"scripts"/
       "neotoma_six_occid_exact_sample_chronology_v17.py").read_text()
    assert '"LGM_alnus_population_occupancy_accepted":0' in s
    assert '"source_only_no_genetic_response_or_21ka_occupancy":True' in s
