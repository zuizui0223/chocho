"""Source-only chronological default-model comparison, never ancient occupancy inference."""
import sys,json
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_Alnus_model_choice_v21 import (
    choose_calibrated_default,evaluate,readcsv,yes,numeric
)

def models():
    return {
      "7725":{"source_datasetid":"14556","isdefault":"False",
         "source_reliable_younger_BP":"-40","source_reliable_older_BP":"34335",
         "source_control_depth_min_cm":"0","source_control_depth_max_cm":"450",
         "source_age_controls_n":"7"},
      "7726":{"source_datasetid":"14556","isdefault":"True",
         "source_reliable_younger_BP":"-28","source_reliable_older_BP":"11172",
         "source_control_depth_min_cm":"0","source_control_depth_max_cm":"450",
         "source_age_controls_n":"7"},
      "13034":{"source_datasetid":"22683","isdefault":"True",
         "source_reliable_younger_BP":"-50","source_reliable_older_BP":"12540",
         "source_control_depth_min_cm":"0","source_control_depth_max_cm":"406",
         "source_age_controls_n":"7"}
    }

def original_sample(chron,default,age,depth,did,alternates):
    return {"datasetid":did,"sampleid":"sample1","siteid":"100",
       "source_chronologyid":chron,"original_chronology_isdefault":default,
       "sample_center_age_BP":str(age),"sample_depth_cm":str(depth),
       "source_age_type":"Calibrated radiocarbon years BP",
       "source_sample_has_both_age_bounds":"False",
       "original_Alnus_pollen_positive_datum_value":"4.0",
       "alternative_original_source_age_models":json.dumps(alternates)}

def test_same_original_source_sample_nondefault_21ka_becomes_9ka_with_true_default():
    row=original_sample("7725","False",21087,528,"14556",[
       {"chronologyid":"7726","original_age_numeric":9003,
        "original_age_type":"Calibrated radiocarbon years BP",
        "default_model_flag_from_chronology":True}])
    q=choose_calibrated_default(row,models())
    assert q["original_selected_source_age_BP"]==21087
    assert q["default_calibrated_source_age_BP"]==9003
    assert not q["default_calibrated_center_inside_18_24ka"]
    assert not q["sample_depth_in_default_returned_control_extent"]

def test_pollen_age_inside_window_but_older_than_own_source_reliable_limit():
    row=original_sample("13034","True",20368,392,"22683",[
       {"chronologyid":"13033","original_age_numeric":17830,
        "original_age_type":"Calibrated radiocarbon years BP",
        "default_model_flag_from_chronology":False}])
    q=choose_calibrated_default(row,models())
    assert q["default_calibrated_center_inside_18_24ka"]
    assert q["sample_depth_in_default_returned_control_extent"]
    assert q["default_center_within_reported_reliable_span"] is False

def test_dont_convert_raw_radiocarbon_years_to_calibrated_BP():
    row=original_sample("13034","True",20368,392,"22683",[
       {"chronologyid":"13033","original_age_numeric":20000,
        "original_age_type":"Radiocarbon years BP",
        "default_model_flag_from_chronology":True}])
    q=choose_calibrated_default(row,models())
    assert q["default_calibrated_source_age_BP"]==20368

def test_missing_or_inconsistent_default_model_fails_closed():
    r=original_sample("7725","False",21087,528,"14556",[])
    with pytest.raises(ValueError):choose_calibrated_default(r,models())
    r=original_sample("7725","False",21087,528,"14556",[
       {"chronologyid":"7726","original_age_numeric":9003,
        "original_age_type":"Calibrated radiocarbon years BP",
        "default_model_flag_from_chronology":True}])
    m=models();m["7726"]["isdefault"]="False"
    with pytest.raises(ValueError):choose_calibrated_default(r,m)

def test_fixed_real_21_neotoma_samples_have_nine_metadata_only_survivors():
    root=Path(__file__).resolve().parents[1]/"results"/"lgm_neotoma_21_Alnus_sample_chronology_control_qc_v20"
    samples=readcsv(root/"all_21_original_Alnus_age_depth_caveats_v20.csv")
    models_source=readcsv(root/"all_11_original_chronology_metadata_v20.csv")
    summary=json.loads((root/"summary.json").read_text())
    result,q=evaluate(samples,models_source,summary)
    assert len(result)==21
    counts=q["sample_counts"]
    assert counts["default_calibrated_center_inside_18_24ka"]==11
    assert counts["default_calibrated_center_outside_18_24ka"]==10
    assert counts["default_calibrated_center_in_18_24ka_and_reported_reliable_span"]==10
    assert counts["default_calibrated_center_18_24ka_both_reliable_and_controls_extent"]==9
    assert counts["source_missing_individual_uncertainty_bounds"]==20
    assert q["accepted_independently_calibrated_21ka_samples"]==0
