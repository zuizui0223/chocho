import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from qualify_ttf_native_host_reference_v01 import qualifications

def test_both_side_reference_threshold_and_native_host_completeness():
    native={"10":[{"gbif_id":"1"}],"20":[{"gbif_id":"2"}]}
    points=[{"wcvp_id":"10"}]*30+[{"wcvp_id":"20"}]*20
    a=qualifications({"10","20"},40,points,native)
    assert a["paired"] and a["all_native_hosts_observed"]
    assert a["all_native_hosts_selected"]
    assert a["source_native_host_ids"]==2

def test_missing_current_native_record_is_not_inferred_absent():
    native={"10":[{"gbif_id":"1"}]}
    points=[{"wcvp_id":"10"}]*100
    a=qualifications({"10","20"},100,points,native)
    assert a["paired"]
    assert not a["all_native_hosts_observed"]
    assert not a["all_native_hosts_selected"]

def test_insect_count_still_required_after_native_host_screen():
    native={"10":[{"gbif_id":"1"}]}
    p=[{"wcvp_id":"10"}]*200
    assert not qualifications({"10"},29,p,native)["paired"]
    assert qualifications({"10"},30,p,native)["paired"]

def test_too_few_native_host_positions_cannot_be_backfilled():
    native={"10":[{"gbif_id":"1"}]}
    p=[{"wcvp_id":"10"}]*29
    assert not qualifications({"10"},200,p,native)["paired"]
