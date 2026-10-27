"""Fail-closed tests for original HOSTS 257-label scope taxonomy."""
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from build_hosts_location_scope_crosswalk_20261027 import classify, from_counts


def test_country_can_be_isolated_but_plant_status_must_be_unknown():
    assert classify("USA") == ("SINGLE_COUNTRY_OR_TERRITORY_NAME","US")
    assert classify("Taiwan") == ("SINGLE_COUNTRY_OR_TERRITORY_NAME","TW")
    assert classify("Puerto Rico") == ("SINGLE_COUNTRY_OR_TERRITORY_NAME","PR")
    result=from_counts({"USA":5,"Taiwan":2},"a"*64)
    assert len(result["source_verified_label_entries"])==2
    assert all(x["plant_native_or_introduced_status"]=="NOT_CLASSIFIED"
               for x in result["source_verified_label_entries"])
    assert all(x["wgsrpd3_source_region_code"] is None
               for x in result["source_verified_label_entries"])


def test_island_and_multi_region_labels_never_promoted_to_whole_country():
    assert classify("Hawaii") == ("SUBNATIONAL_NAME_NEEDS_GEOMETRY","US")
    assert classify("West Malaysia") == ("SUBNATIONAL_NAME_NEEDS_GEOMETRY","MY")
    assert classify("Nearctic")[0]=="BIOGEOGRAPHIC_OR_MULTICOUNTRY_NOT_PRIMARY"
    assert classify("Oriental")[0]=="BIOGEOGRAPHIC_OR_MULTICOUNTRY_NOT_PRIMARY"
    assert classify("British Isles")[0]=="BIOGEOGRAPHIC_OR_MULTICOUNTRY_NOT_PRIMARY"


def test_no_field_inferences_from_interceptions_or_provisional_localities():
    assert classify("USA [interception]")[0]=="INTERCEPTION_NOT_FIELD_HOST_USE"
    assert classify("Ecuador (prov)")[0]=="PROVISIONAL_LOCATION_NOT_PRIMARY"
    assert classify("Korea")[0]=="HISTORICAL_OR_MULTINATIONAL_NOT_PRIMARY"
    assert classify("Czechoslovakia")[0]=="HISTORICAL_OR_MULTINATIONAL_NOT_PRIMARY"
    assert classify("Made-Up-Location-Unit")[0]=="UNREVIEWED_LABEL_NOT_PRIMARY"


def test_area_stats_preserve_all_host_source_denominators():
    n={"Nearctic":24,"USA":15,"West Malaysia":6,"USA [interception]":2,
       "Unknown archaeological region":3}
    x=from_counts(n,"a"*64)
    assert x["source_species_level_annotated_rows"]==sum(n.values())
    assert sum(x["row_counts_by_scope"].values())==sum(n.values())
    assert x["distinct_location_labels"]==len(n)
