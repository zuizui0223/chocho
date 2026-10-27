"""Non-network independent tests for fixed-source 24 recipient host pruning."""
import csv
import json
import sys
from pathlib import Path
from collections import Counter

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_original24_host_pruning_independent_route_cases_20261027 import (
    SOURCE_PROFILE_SHA,SOURCE_PROFILE_SCHEMA,EXPECTED_COUNTS,
    SCHEMA,load_profile,count_local_host_source,load_cases,
)


def strict_row(i, *, intro=3, native=0, H=10, sp=None):
    return {
        "species":sp or f"Butterfly species{i}",
        "independent_2024_recipient_zone":"Hawaii",
        "geographic_status":"MAIN_COMPARABLE",
        "accepted_original_known_global_larval_hosts":H,
        "introduced_only_wild_source_hosts":intro,
        "native_only_wild_source_hosts":native,
        "source_plant_native_and_introduced_in_country_atoms":0,
        "unlisted_or_below_area_gate_known_hosts":H-intro-native,
    }


def coarse_row(i, *, grade="UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY"):
    return {
        "species":f"Unknown species{i}",
        "recipient_zone":"Europe",
        "geographic_status":grade,
        "source_mapping_status":"UNSCORABLE_SUPRA_COUNTRY_REGION",
    }


def test_original_24_species_local_host_count_conservation_and_no_failures():
    rows=[strict_row(i,intro=1+(i%6),native=int(i%2==0),H=11)
          for i in range(24)]+[coarse_row(i) for i in range(11)]
    summary, strict, nonmain, repeats=count_local_host_source(rows)
    assert summary["n_original_strict_geographic_source_positive_cases"]==24
    assert summary["n_distinct_original_butterfly_taxa_in_24_cases"]==24
    assert summary["n_with_at_most_3_introduced_wild_hosts"]==12
    assert summary["source_local_intro_species_median"]==3.5
    assert summary["n_with_any_native_wild_host_alternative"]==12
    assert all(r["source_host_range_not_identical_to_historical_plant_distribution"]
               for r in strict)
    assert all(r["not_scored_as_absent"] for r in nonmain)
    assert repeats==[]


def test_24_same_species_in_multiple_regions_varies_with_constant_original_hosts():
    rows=[strict_row(i,H=40) for i in range(24)]
    rows[0]["species"]="Butterfly alpha"
    rows[0]["independent_2024_recipient_zone"]="Hawaii"
    rows[0]["introduced_only_wild_source_hosts"]=3
    rows[0]["unlisted_or_below_area_gate_known_hosts"]=37
    rows[1]["species"]="Butterfly alpha"
    rows[1]["independent_2024_recipient_zone"]="New Zealand"
    rows[1]["introduced_only_wild_source_hosts"]=20
    rows[1]["unlisted_or_below_area_gate_known_hosts"]=20
    output,_,_,r=count_local_host_source(rows+[coarse_row(i) for i in range(11)])
    a=next(x for x in r if x["species"]=="Butterfly alpha")
    assert a["source_min_current_introduced_host_species"]==3
    assert a["source_max_current_introduced_host_species"]==20
    assert a["number_of_distinct_source_receiving_areas"]==2
    assert a["constant_original_species_HOSTS_list_across_receiving_areas"]


def test_host_species_count_drift_and_geography_coarseness_stop():
    a=[strict_row(i) for i in range(24)]
    a[0]["unlisted_or_below_area_gate_known_hosts"]=999
    with pytest.raises(ValueError,match="decomposition"):
        count_local_host_source(a+[coarse_row(i) for i in range(11)])


def test_wrong_profile_SHA_and_original35_population_fail_closed(tmp_path):
    path=tmp_path/"profile.json"
    path.write_text(json.dumps({
      "schema":SOURCE_PROFILE_SCHEMA,
      "source":{"original_2024_positive_species_recipient_records":34,
                "original_2024_matched_aliens":18,
                "original_239_butterflies":239,
                "original_accepted_butterfly_host_pairs":2596}
    }))
    with pytest.raises(ValueError,match="SHA changed"):
        load_profile(path)
    with pytest.raises(ValueError,match="source panel changed"):
        load_profile(path,enforce_source=False)


def test_source_pathway_first_report_does_not_identify_ship_or_actual_host(tmp_path):
    source=ROOT/"data/exploratory/recipient_island_pathway_resource_case_evidence_20261027.csv"
    assert source.exists()
    rows=[
        {"species":"Calpodes ethlius","recipient_zone":"Galapagos",
         "geographic_status":"MAIN_COMPARABLE",
         "introduced_only_wild_source_hosts":3,
         "native_only_wild_source_hosts":0,
         "accepted_original_known_global_larval_hosts":5},
        {"species":"Leptosia nina","independent_2024_recipient_zone":"Okinawa",
         "geographic_status":"SENSITIVITY_BOTANICAL_PROXY",
         "introduced_only_wild_source_hosts":0,
         "native_only_wild_source_hosts":2,
         "accepted_original_known_global_larval_hosts":10},
    ]
    case=load_cases(source,rows)
    assert len(case)==2
    gal=next(x for x in case if x["independently_established_recipient"]=="Galapagos")
    assert gal["source_first_record_of_consumer_year"]==1987
    assert gal["source_first_record_of_named_plant_year"]==1991
    assert gal["first_record_temporal_order_does_not_prove_absence"][
        "butterfly_report_precedes_this_botanical_first_report"]
    assert gal["carrier_host_identification_status"]=="UNVERIFIED_NOT_ZERO"
    ok=next(x for x in case if x["independently_established_recipient"]=="Okinawa")
    assert ok["source_evidence_grade"]=="BOTANICAL_NANSEI_PROXY_NOT_EXACT_OKINAWA"
    assert ok["source_first_record_of_named_plant_year"] is None


def test_case_source_must_not_call_carrier_host_species_verification_true(tmp_path):
    source=ROOT/"data/exploratory/recipient_island_pathway_resource_case_evidence_20261027.csv"
    with source.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        fieldnames=reader.fieldnames
        rows=list(reader)
    rows[0]["same_carrier_and_actual_larval_host_species_verified"]="YES"
    out=tmp_path/"false_confirmation.csv"
    with out.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader();writer.writerows(rows)
    with pytest.raises(ValueError,match="No original source evidence links"):
        load_cases(out,[
            {"species":"Calpodes ethlius","recipient_zone":"Galapagos",
             "geographic_status":"MAIN_COMPARABLE",
             "introduced_only_wild_source_hosts":3,
             "native_only_wild_source_hosts":0,
             "accepted_original_known_global_larval_hosts":5},
            {"species":"Leptosia nina","independent_2024_recipient_zone":"Okinawa",
             "geographic_status":"SENSITIVITY_BOTANICAL_PROXY",
             "introduced_only_wild_source_hosts":0,
             "native_only_wild_source_hosts":2,
             "accepted_original_known_global_larval_hosts":10},
        ])
