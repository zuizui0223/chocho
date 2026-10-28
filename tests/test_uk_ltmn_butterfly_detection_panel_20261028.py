"""Identity and absence-interpetation tests for exact visited transect data."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK,site_effort,count_at_visits,tally

def visit(s,y,m=6,d=2,t="720"):
    return {"UKBMS_SITENO":s,"YEAR":y,"MONTH":m,"DAY":d,"START_TIME":t}

def counted(s,y,sp,number=2,m=6,d=2,t="720"):
    return {**visit(s,y,m,d,t),"WEEKNO":9,"COMMON_NAME":sp,"VISIT_TOTAL":number}

def test_immutable_uk_common_names_are_unique_and_not_unresolved_wall():
    assert len(CROSSWALK)==27 and len(set(CROSSWALK.values()))==27
    assert "Lasiommata megera" not in CROSSWALK
    assert CROSSWALK["Maniola jurtina"]=="Meadow Brown"

def test_other_species_denominator_never_conditions_on_focal_record():
    focal="Meadow Brown"
    ev=("A",1988,6,2,"720")
    panel={"A"}
    ve={"A":{"early":{ev},"late":{("A",2015,6,2,"720")}}}
    observations={
        ev:{focal,"Small White","Comma","Peacock"},
        ("A",2015,6,2,"720"):{"Small White","Comma","Peacock"},
    }
    out=tally(ve,panel,observations,3)
    row=next(x for x in out if x["original_UKBMS_common_name"]==focal)
    assert row["target_specific_eligible_source_visits_all_panel_by_era"]=={"early":1,"late":1}
    assert row["periods"]["early"]["focal_detected_visits"]==0
    assert row["information_gate"]=="INSUFFICIENT_FOCAL_DETECTION_SUPPORT"

def test_target_group_names_gated_and_site_not_excluded_when_zero_count():
    events=[visit("B",1986),visit("B",2016)]
    original,panel,period=site_effort(iter(events),{"B":(2,4,("grassland",""))})
    assert len(original)==2 and panel==set()

def test_exact_visit_join_disallows_unvisited_and_moths():
    v=visit("B",1986)
    ev=("B",1986,6,2,"720")
    data=[counted("B",1986,"Meadow Brown"),counted("B",1986,"Silver Y"),counted("B",1986,"Peacock",1),
          counted("B",1986,"No butterflies seen",2)]
    out,qc=count_at_visits(iter(data),{ev},{"B"})
    assert out[ev]=={"Meadow Brown","Peacock"}
    assert qc["n_positive_original_named_visit_species_pairs"]==2
