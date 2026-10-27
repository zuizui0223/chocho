"""Source-independent tests: WCVP locality correctness and positive-only 35 recipient records."""
import math
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_zenodo35_recipient_WCVP_host_opportunities_20261027 import (
    zone_mapping,country_plant,classify_plant,species_source_record,
    hypergeom,artificial_reference,host_bin,EXACT_ZONES,
    ALL_ZONE_NAMES,SCHEMA,
)


def geo():
    regions={
       "HAW":"Hawaii","OGA":"Ogasawara-shoto",
       "GAL":"Galápagos","KOR":"Korea",
       "NNS":"Nansei-shoto",
    }
    atoms={"CL":[(("CLC",),95.),(("CLN",),5.)],
           "JP":[(("JAP",),95.),(("NNS",),5.)]}
    metadata={
       "CL":{"is_geo_safe":True,"guard_coverage":.9496638677671182,
             "country_area_km2":105.3},
       "JP":{"is_geo_safe":True,"guard_coverage":.9721632427516506,
             "country_area_km2":102.86},
    }
    return regions,atoms,metadata


def test_11_fixed_receiving_areas_and_exact_botanical_source_names():
    assert len(ALL_ZONE_NAMES)==11
    regions,atoms,meta=geo()
    for name in ("Hawaii","Ogasawara","Galapagos","Korea"):
        r=zone_mapping(name,regions,atoms,meta)
        assert r["mapping_status"]=="EXACT_WCVP_LEVEL3_SOURCE_REGION"
        assert r["primary_geographic_evidence"]
    o=zone_mapping("Okinawa",regions,atoms,meta)
    assert o["source_region_code"]=="NNS"
    assert not o["primary_geographic_evidence"]
    for broad in ("Europe","North America"):
        assert zone_mapping(broad,regions,atoms,meta)[
            "mapping_status"]=="UNSCORABLE_SUPRA_COUNTRY_REGION"
    with pytest.raises(ValueError,match="Unknown independent"):
        zone_mapping("Random Country",regions,atoms,meta)


def test_Chile_94966pct_original_map_coverage_must_fail95_but_pass90():
    regions,atoms,meta=geo()
    a=zone_mapping("Chile",regions,atoms,meta,coverage=.95)
    assert a["mapping_status"]=="ORIGINAL_COUNTRY_SOURCE_AREA_COVERAGE_BELOW_GATE"
    a=zone_mapping("Chile",regions,atoms,meta,coverage=.90)
    assert a["mapping_status"]=="ORIGINAL_ISO_COUNTRY_ATOM_AREA_RECONSTRUCTION"


def test_exact_island_host_source_not_US_country_and_missing_not_failure():
    regions,atoms,meta=geo()
    m=zone_mapping("Hawaii",regions,atoms,meta)
    native={"x":{"HAW"}}
    current={"x":{"HAW"},"y":{"HAW"}}
    assert classify_plant("x",m,atoms,meta,native,current)=="NATIVE_WILD"
    assert classify_plant("y",m,atoms,meta,native,current)=="INTRODUCED_WILD"
    assert classify_plant("unknown",m,atoms,meta,native,current)==\
        "NO_WCVP_WILD_RECORDED_OR_BELOW_AREA_GATE"
    r=species_source_record("Butterfly one","Hawaii",{"x","y","unknown"},
       m,atoms,meta,native,current,[])
    assert r["source_plant_host_support_class"]==\
        "BOTH_NATIVE_AND_INTRODUCED_POTENTIAL_HOSTS"
    assert r["accepted_original_known_global_larval_hosts"]==3
    assert r["unlisted_or_below_area_gate_known_hosts"]==1
    assert r["source_global_HOSTS_link_not_direct_local_larval_feeding"]


def test_country_union_splits_native_and_introduced_and_guards_small_fragments():
    _,atoms,meta=geo()
    n={"plant":{"JAP"},"native":{"JAP"}}
    current={"plant":{"JAP","NNS"},"native":{"JAP"},"intro":{"NNS"}}
    status=country_plant("plant","JP",atoms,meta,n,current)
    assert status[0] and status[1]
    jp=zone_mapping("Japan",{"HAW":"Hawaii"},atoms,meta,coverage=.95)
    assert classify_plant("plant",jp,atoms,meta,n,current)==\
        "BOTH_NATIVE_AND_INTRODUCED_WILD"
    assert classify_plant("intro",jp,atoms,meta,n,current)==\
        "INTRODUCED_WILD"
    assert classify_plant("native",jp,atoms,meta,n,current)==\
        "NATIVE_WILD"


def test_source_positive_in_host_count_bin_matched_reference_is_not_p_value():
    original=[]
    for i in range(24):
        original.append({
            "accepted_original_known_global_larval_hosts":2 if i<12 else 25,
            "source_any_introduced_WCVP_wild_host": i%2==0
        })
    positive=[original[0],original[2],original[14],original[16]]
    a=artificial_reference(original,positive)
    assert a["source_established_butterflies_with_intro_wild_hosts_observed"]==4
    assert math.isclose(a["source_accepted_breadth_matched_artificial_expected"],2.)
    assert a["source_reference_95pct_interval"][0]<=4
    assert a["NO_butterfly_arrival_or_nonestablishment_data_for_this_reference"]
    assert [host_bin(x) for x in (1,5,6,20,21,50,51)]==[
        "1-5","1-5","6-20","6-20","21-50","21-50","51+"
    ]


def test_source_exact_hypergeom_not_fabricated_for_impossible_bins():
    p=hypergeom(4,2,2)
    assert math.isclose(p[0],1/6)
    assert math.isclose(p[1],4/6)
    assert math.isclose(p[2],1/6)
    with pytest.raises(ValueError,match="Impossible"):
        hypergeom(5,10,2)
