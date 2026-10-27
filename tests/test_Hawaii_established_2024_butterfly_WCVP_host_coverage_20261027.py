"""Exact positive established species and WCVP Hawaii host source controls."""
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_Hawaii_established_2024_butterfly_WCVP_host_coverage_20261027 import (
    classify_plants, classify_source_positive_species, SCHEMA
)


def test_WCVP_native_precedes_introduced_and_no_wild_is_missing_not_zero():
    a=classify_plants({"x","y","z"},
        native={"x":{"HAW"}} ,
        current={"x":{"HAW"},"y":{"HAW"}})
    assert a["NATIVE_HAW_WILD"]==1
    assert a["INTRODUCED_HAW_WILD"]==1
    assert a["NOT_RECORDED_HAW_WILD"]==1
    assert sum(a.values())==3


def case(species,n,year="1977"):
    return {
      "species":species,
      "original_239_host_species_breadth":n,
      "original_239_source_added_WCVP_resource_units":20,
      "source_documented_non_native_recipient_areas":["Hawaii"],
      "original_source_region_and_year_records":[{
         "recipient_zone":"Hawaii",
         "year_source_first_report":year,
         "source_intentional_release_flag":"No",
         "source_eradication_flag":"NA",
         "source_established_indoor_outdoor_flag":"Out"
      }]
    }


def test_distinct_species_positive_only_host_coverage_not_establishment_failure():
    hawaii=[
       case("Butterfly first",2),
       case("Butterfly second",1),
       case("Butterfly third",1,"1758"),
       case("Butterfly fourth",1),
    ]
    host={
        "Butterfly first":{"x","y"},
        "Butterfly second":{"x"},
        "Butterfly third":{"u"},
        "Butterfly fourth":{"y"},
    }
    native={"x":{"HAW"}}
    current={"x":{"HAW"},"y":{"HAW"}}
    x=classify_source_positive_species(hawaii,host,native,current)
    d={r["species"]:r for r in x}
    assert d["Butterfly first"]["accepted_host_support_class"]=="BOTH_NATIVE_AND_INTRODUCED_WILD_HOSTS"
    assert d["Butterfly second"]["accepted_host_support_class"]=="ONLY_WCVP_NATIVE_WILD_HOSTS"
    assert d["Butterfly third"]["accepted_host_support_class"]=="NO_HOST_IN_WCVP_WILD_HAWAII_MAP"
    assert d["Butterfly fourth"]["accepted_host_support_class"]=="ONLY_WCVP_INTRODUCED_WILD_HOSTS"
    assert d["Butterfly third"]["chronology_flag"].startswith("YEAR_<=1800")
    assert all(r["does_NOT_identify_actual_Hawaii_larval_host_use"] for r in x)


def test_source_host_richness_drift_stops():
    with pytest.raises(ValueError,match="per species host link count"):
        classify_source_positive_species(
            [case("Pieris rapae",3)],
            {"Pieris rapae":{"host1","host2"}},
            native={}, current={},
        )
