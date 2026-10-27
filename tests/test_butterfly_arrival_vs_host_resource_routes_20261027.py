"""No-network negative controls for the two independently transportable channels."""
import csv
import json
import math
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_butterfly_arrival_vs_host_resource_routes_20261027 import (
    TOP10,MODES,TRADE,added_resource_credits_for_one_species,
    load_route_pilot,audit115,join_top10
)

LEDGER=ROOT/"data/exploratory/arrival_resource_top10_propagation_route_pilot_20261027.csv"


def test_added_resource_unit_is_fractionally_split_not_double_counted():
    native={"A":{"oldA"},"B":{"oldB"}}
    current={"A":{"oldA","newB"},"B":{"oldB","newB"}}
    added,credit=added_resource_credits_for_one_species(
        {"A","B"},native,current)
    assert added=={"newB"}
    assert credit["A"]==.5
    assert credit["B"]==.5
    assert math.isclose(sum(credit.values()),len(added))
    current["C"]={"newB"}
    added,c=added_resource_credits_for_one_species(
        {"A","B","C"},native,current)
    assert len(added)==1 and math.isclose(c["A"],1/3)
    assert math.isclose(sum(c.values()),1.0)


def test_original_top10_plant_tags_allow_overlapping_modes_but_no_trade_observation():
    rows=load_route_pilot(LEDGER)
    assert len(rows)==10
    assert tuple(x["accepted_botanical_name"] for x in rows)==TOP10
    assert all(row["documented_seed_planting_or_seed_spread"]=="1" for row in rows)
    assert sum(row["documented_live_sod_sprigs_or_rooted_material"]=="1"
               for row in rows)==2
    assert sum(row["documented_container_nursery_weed_association"]=="1"
               for row in rows)==2
    assert sum(row["documented_local_seedling_transplantation"]=="1"
               for row in rows)==1
    for row in rows:
        assert all(row[key]=="NOT_VERIFIED" for key in TRADE)


def test_top10_propagation_is_not_an_exclusive_trade_partition():
    ledger=load_route_pilot(LEDGER)
    names=[r["accepted_botanical_name"] for r in ledger]
    fake=[{"accepted_botanical_name":name,
           "fraction_of_all_14553_original_added_units":.02,
           "source_fractional_added_resource_credit":2.,
           "original_rank":i+1,"family":"Synthetic"}
          for i,name in enumerate(names)]
    result=join_top10(fake,ledger,{"top_10_host_fractional_credit_fraction":.2})
    categories=result["mode_coverage_of_14553_original_added_resource_units_not_mutually_exclusive"]
    assert math.isclose(categories[MODES[0]],.2)
    assert math.isclose(categories[MODES[1]],.04)
    assert math.isclose(categories[MODES[2]],.04)
    assert math.isclose(categories[MODES[3]],.02)
    assert sum(categories.values())>.2
    assert result["not_a_multinomial_mode_partition"] is True
    assert all(p["actual_species_level_international_plant_shipment_observed"] is None
               for p in result["original_top10"])


def test_original_115_are_all_observed_and_do_not_contain_49_absences(tmp_path):
    receipt=tmp_path/"receipt.json"
    receipt.write_text(json.dumps({
        "outside_native_species_x_units":115,
        "outside_native_species_x_units_recovered_by_introduced_hosts":66
    }))
    species=tmp_path/"species.csv"
    with species.open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["species","observed_outside_native",
                         "outside_native_units_explained_by_introduced_host_ranges"])
        for i in range(32):
            writer.writerow([f"Species{i}",115 if i==0 else 0,
                             66 if i==0 else 0])
    result=audit115(receipt,species)
    assert result["outside_native_butterfly_observed_units"]==115
    assert result["observed_units_with_introduced_resource_overlap"]==66
    assert result["observed_units_not_recovered_by_any_reconstructed_current_host_range"]==49
    assert result["butterfly_occurrence_outcome_by_original_115_units"]=="ALL_115_OBSERVED"
    assert result["no_unobserved_units_supplied_for_arrival_presence_absence_comparison"]


def test_source_115_report_not_a_presence_absence_outcome(tmp_path):
    receipt=tmp_path/"receipt.json"
    receipt.write_text(json.dumps({
        "outside_native_species_x_units":115,
        "outside_native_species_x_units_recovered_by_introduced_hosts":66
    }))
    csvfile=tmp_path/"panel.csv"
    with csvfile.open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["species","observed_outside_native",
                         "outside_native_units_explained_by_introduced_host_ranges"])
        for i in range(32):
            writer.writerow([f"Species{i}",1,1])
    with pytest.raises(ValueError,match="source margins changed"):
        audit115(receipt,csvfile)


def test_unsourced_trade_must_never_become_zero_transport(tmp_path):
    rows=load_route_pilot(LEDGER)
    rows[0]["crossborder_commodity_transit_by_this_form"]="0"
    p=tmp_path/"misleading.csv"
    with p.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError,match="NO species"):
        load_route_pilot(p)
