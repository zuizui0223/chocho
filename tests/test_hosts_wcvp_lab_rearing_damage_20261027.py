"""Frozen source-independent QA for original HOSTS Lab Rearing/Damage flags by WCVP origin."""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from collections import Counter
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_hosts_wcvp_lab_rearing_damage_20261027 as m


def triples():
    return {
        ("Butterfly test","P1","Hawaii"): {
            "lab":{"Y","<EMPTY>"},
            "damage":{"on leaves","<EMPTY>"},"rows":2},
        ("Butterfly test","P2","USA"): {
            "lab":{"?"}, "damage":{"adult association"},"rows":1},
        ("Butterfly test","P3","Nearctic"): {
            "lab":{"<EMPTY>"},"damage":{"<EMPTY>"},"rows":1},
        ("Butterfly test","P4","West Malaysia"): {
            "lab":{"Y"},"damage":{"oviposition only"},"rows":1},
        ("Butterfly test","P5","Japan"): {
            "lab":{"<EMPTY>"},"damage":{"on leaves"},"rows":1},
    }


def test_raw_source_flag_to_accepted_triple_collapsing():
    values=triples()
    statuses={
        ("Butterfly test","P1","Hawaii"):"INTRODUCED_CONFIRMED",
        ("Butterfly test","P2","USA"):"NATIVE_CONFIRMED",
        ("Butterfly test","P3","Nearctic"):"UNMAPPED_OR_UNREVIEWED_LOCATION",
        ("Butterfly test","P4","West Malaysia"):"INTRODUCED_CONFIRMED",
        ("Butterfly test","P5","Japan"):"WCVP_UNLISTED",
    }
    result=m.summarize(values,statuses)
    introduced=result["counts_by_source_origin"]["INTRODUCED_CONFIRMED"]
    assert introduced["source_distinct_butterfly_host_Location_triples"]==2
    assert introduced["underlying_original_HOSTS_rows"]==3
    assert introduced["triple_contains_any_Lab_Rearing_Y"]==2
    assert introduced["triple_contains_any_nonempty_Damage"]==2
    assert introduced["triple_contains_multiple_rearing_literal_values"]==1
    assert introduced.get("triple_all_Lab_Rearing_blank",0)==0
    assert result["introduced_area_top_damage_literal_values"]==[
        {"original_literal_damage":"on leaves","distinct_triples":1},
        {"original_literal_damage":"oviposition only","distinct_triples":1}
    ]
    native=result["counts_by_source_origin"]["NATIVE_CONFIRMED"]
    assert native["triple_contains_any_Lab_Rearing_question_mark"]==1
    assert native["triple_contains_any_Lab_Rearing_Y"]==0
    assert result["counts_by_source_origin"]["UNMAPPED_OR_UNREVIEWED_LOCATION"][
        "triple_all_Lab_Rearing_blank"]==1


def test_false_nonempty_damage_as_juvenile_rearing_not_inferred():
    vals={("B b","P1","USA"):{"lab":{"<EMPTY>"},
                                 "damage":{"adult association","oviposition only"},
                                 "rows":2}}
    output=m.summarize(vals,{("B b","P1","USA"):"INTRODUCED_CONFIRMED"})
    counts=output["counts_by_source_origin"]["INTRODUCED_CONFIRMED"]
    assert counts["triple_contains_any_nonempty_Damage"]==1
    assert counts.get("triple_contains_any_Lab_Rearing_Y",0)==0
    assert "confirmed_larval_survival" not in str(output)


def test_missing_source_status_fails_closed():
    with pytest.raises(ValueError,match="Missing expected source"):
        m.summarize(triples(),{})


def test_accept_taxon_host_raw_location_and_no_row_duplication(tmp_path,monkeypatch):
    fields=["Insect Genus","Insect Species","Hostplant Genus","Hostplant Species",
            "Location","Lab Rearing","Damage"]
    cases=[
      ["Butterfly","test","Plant","one","Hawaii","Y","on leaves"],
      ["Butterfly","test","Plant","one","Hawaii","",""],
      ["Butterfly","test","Plant","two","USA","?","adult association"],
      ["Wrong","taxon","Plant","one","Hawaii","Y","on leaves"],
      ["Butterfly","test","Plant","","USA","Y","on leaves"],
    ]
    p=tmp_path/"original_source.csv"
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(fields);w.writerows(cases)
    # Deliberately synthetic source margins; verify actual full sizes in CI
    # via source SHA checks and the fixed original counters.
    monkeypatch.setattr(m,"ORIGINAL_SOURCE_ROW_COUNTERS",(5,4,3,3,2))
    t,flow=m.collect(p,{
        ("Butterfly test","Plant one"):"accepted1",
        ("Butterfly test","Plant two"):"accepted2"
    },{"Butterfly test"})
    assert len(t)==2
    assert t["Butterfly test","accepted1","Hawaii"]["lab"]=={"Y","<EMPTY>"}
    assert t["Butterfly test","accepted1","Hawaii"]["rows"]==2
    assert flow["accepted_focal_source_rows"]==3


def test_unknown_rearing_symbol_cannot_be_silently_positive(tmp_path,monkeypatch):
    fields=["Insect Genus","Insect Species","Hostplant Genus","Hostplant Species",
            "Location","Lab Rearing","Damage"]
    p=tmp_path/"original_source.csv"
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(fields)
        w.writerow(["Butterfly","test","Plant","one","Hawaii","yes",""])
    monkeypatch.setattr(m,"ORIGINAL_SOURCE_ROW_COUNTERS",(1,1,1,1,1))
    with pytest.raises(ValueError,match="Unknown original HOSTS laboratory"):
        m.collect(p,{("Butterfly test","Plant one"):"P1"},{"Butterfly test"})
