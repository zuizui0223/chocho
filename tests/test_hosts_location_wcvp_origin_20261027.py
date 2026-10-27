"""Conservative original HOSTS geographic source and WCVP provenance tests."""
import csv
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_hosts_location_wcvp_origin_20261027 as source


def synthetic_geo():
    area={
        "US":[(("HAW",),10.),(("NCA",),90.)],
        "MY":[(("MLY",),40.),(("SBH",),40.),(("SWK",),20.)],
        "BR":[(("BRA",),70.),(("BRL",),30.)],
    }
    meta={
        k:{"country_area_km2":100.,"guard_coverage":1.0,
           "area_identity_checks":"PASS","is_geo_safe":True}
        for k in area
    }
    return area,meta,{"HAW","NCA","MLY","SBH","SWK","BRA","BRL"}


def test_country_area_is_not_the_count_of_equal_botanical_regions():
    atoms,meta,regions=synthetic_geo()
    native={"plant1":{"NCA"}}
    current={"plant1":{"NCA","HAW"}}
    assert source.dominant_country(
        "plant1","US",threshold=.95,atoms=atoms,metadata=meta,
        native=native,current=current)=="MIXED_OR_PARTIAL"
    assert source.dominant_country(
        "plant1","US",threshold=.80,atoms=atoms,metadata=meta,
        native=native,current=current)=="NATIVE_CONFIRMED"
    # WCVP origin in the true Hawaiian subregion does NOT inherit US status.
    assert source.classify_location(
        "plant1","Hawaii",threshold=.95,atoms=atoms,
        atom_metadata=meta,original_regions=regions,
        native=native,current=current)=="INTRODUCED_CONFIRMED"


def test_west_malaysia_is_not_automatically_east_malaysia():
    atoms,meta,regions=synthetic_geo()
    native={"x":{"SBH","SWK"}}
    current={"x":{"MLY","SBH","SWK"}}
    assert source.classify_location("x","West Malaysia",threshold=.95,
        atoms=atoms,atom_metadata=meta,original_regions=regions,
        native=native,current=current)=="INTRODUCED_CONFIRMED"
    assert source.classify_location("x","East Malaysia",threshold=.95,
        atoms=atoms,atom_metadata=meta,original_regions=regions,
        native=native,current=current)=="SUBREGION_MAPPING_UNSAFE"
    assert source.classify_location("x","West Malaysia",threshold=.95,
        atoms=atoms,atom_metadata=meta,
        original_regions={"SBH","SWK"},native=native,current=current)== \
        "SUBREGION_MAPPING_UNSAFE"


def test_unknown_country_area_not_asserted_unlisted_plant_geography():
    atoms,meta,regions=synthetic_geo()
    meta["US"]["country_area_km2"]=105.
    meta["US"]["guard_coverage"]=100/105
    result=source.dominant_country(
        "unknown","US",threshold=.95,atoms=atoms,metadata=meta,
        native={},current={})
    assert result=="WCVP_UNLISTED"
    meta["US"]["country_area_km2"]=112.
    meta["US"]["guard_coverage"]=100/112
    result=source.dominant_country(
        "unknown","US",threshold=.8,atoms=atoms,metadata=meta,
        native={},current={})
    assert result=="COUNTRY_MAPPING_UNSAFE"
    meta["US"]["is_geo_safe"]=False
    assert source.dominant_country(
        "unknown","US",threshold=.8,atoms=atoms,metadata=meta,
        native={},current={})=="COUNTRY_MAPPING_UNSAFE"


def test_overlapping_country_atoms_of_conflicting_WCVP_class_are_ambiguous():
    atoms,meta,regions=synthetic_geo()
    atoms["BR"]=[(("BRA","BRL"),100.)]
    native={"x":{"BRA"}}
    current={"x":{"BRA","BRL"}}
    assert source.dominant_country(
        "x","BR",threshold=.8,atoms=atoms,metadata=meta,
        native=native,current=current)=="MIXED_OR_PARTIAL"


def test_broad_source_and_interception_are_not_local_larval_reports():
    atoms,meta,regions=synthetic_geo()
    for name in ("Nearctic","Oriental","East Africa",
                 "USA [interception]","British Isles"):
        assert source.classify_location(
            "plant1",name,threshold=.95,atoms=atoms,
            atom_metadata=meta,original_regions=regions,
            native={},current={})=="UNMAPPED_OR_UNREVIEWED_LOCATION"


def test_invasive_host_documentation_is_not_an_adoption_event(tmp_path,monkeypatch):
    rows=[
        ["Butterfly","test","Plant","test","USA"],
        ["Butterfly","test","Plant","test","Hawaii"],
        ["Butterfly","test","Plant","other","Hawaii"],
        ["Butterfly","test","Plant","other","Hawaii"], # exact duplicate
        ["Butterfly","test","Plant","unknown","Nearctic"],
    ]
    p=tmp_path/"original.csv"
    with p.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["Insect Genus","Insect Species","Hostplant Genus",
                    "Hostplant Species","Location"])
        w.writerows(rows)
    monkeypatch.setattr(source,"EXPECTED_HOSTS_ROWS",5)
    monkeypatch.setattr(source,"EXPECTED_SPECIES_HOST_ROWS",5)
    geo,meta,region=synthetic_geo()
    native={"x":{"NCA"},"y":{"NCA"}}
    current={"x":{"NCA","HAW"},"y":{"NCA","HAW"}}
    result=source.observations(
        p,{
            ("Butterfly test","Plant test"):"x",
            ("Butterfly test","Plant other"):"y",
            ("Butterfly test","Plant unknown"):"z"
        },{"Butterfly test"},native=native,current=current,
        atoms=geo,atom_metadata=meta,source_regions=region,families={})
    main=result["primary_and_sensitivity"]["0.95"]
    assert main["unique_species_host_location_triples"]==4
    assert main["triples_by_plant_origin"]["INTRODUCED_CONFIRMED"]==2
    assert main["introduced_triples_with_same_pair_documented_in_native_area"]==0
    assert main["introduced_triples_without_native_area_pair_report"]==2
    assert main["original_raw_records_by_plant_origin"]["INTRODUCED_CONFIRMED"]==3


def test_extant_source_presence_cannot_hide_native_subset_error(tmp_path):
    native=tmp_path/"native.csv";current=tmp_path/"current.csv"
    for path, rows in ((native,[("X","AAA")]),
                       (current,[("X","BBB")])):
        with path.open("w",newline="") as f:
            w=csv.writer(f);w.writerow(["accepted_plant_name_id","area_code_l3"])
            w.writerows(rows)
    with pytest.raises(ValueError,match="not subset"):
        source.wcvp_range(native,current)
