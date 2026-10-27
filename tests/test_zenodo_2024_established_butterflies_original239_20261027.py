"""Synthetic positive-only actual established alien insects cross original butterfly taxa.

No named Zenodo establishment gap is a known arrival failure. No
species-region independence assumptions. Exact host resource source checks.
"""
import csv,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_zenodo_2024_established_butterflies_original239_20261027 as m


def test_original_zenodo_content_hashes_and_table_grain_are_frozen():
    assert m.ZENODO_RECORD==11519262
    assert m.RAW_ROWS=={
        "occurrence_table_20230317.csv":13047,
        "attribute_table_20230317.csv":8514,
        "source_table_20230317.csv":13210,
    }
    assert all(len(x)==64 for x in m.DATA.values())


def test_accepted_species_binomial_cannot_be_auto_rescued_from_genus():
    assert m.binomial("Pieris rapae")=="Pieris rapae"
    assert m.binomial("  Pieris    rapae ")=="Pieris rapae"
    assert m.binomial("Pieris sp.") is None
    assert m.binomial("Pieris brassicae nepalensis") is None
    assert m.binomial("cf. Pieris rapae") is None


def positive_rows():
    return [
        {"genus_species":"Pieris rapae","family":"Pieridae",
         "order":"Lepidoptera","rank":"species",
         "year":"1860","region":"North America",
         "intentional_release":"No","established_indoors_or_outdoors":"outdoors",
         "eradicated":"NA"},
        {"genus_species":"Pieris rapae","family":"Pieridae",
         "order":"Lepidoptera","rank":"species",
         "year":"1930","region":"New Zealand",
         "intentional_release":"No","established_indoors_or_outdoors":"outdoors",
         "eradicated":"NA"},
        {"genus_species":"Thymelicus lineola","family":"Hesperiidae",
         "order":"Lepidoptera","rank":"species",
         "year":"1910","region":"North America",
         "intentional_release":"NA","established_indoors_or_outdoors":"NA",
         "eradicated":"NA"},
        {"genus_species":"Drosophila simulans","family":"Drosophilidae",
         "order":"Diptera","rank":"species","year":"1970",
         "region":"Hawaii","intentional_release":"NA",
         "established_indoors_or_outdoors":"NA","eradicated":"NA"},
    ]


def data():
    return {
        "occurrence_table_20230317.csv":positive_rows(),
        "attribute_table_20230317.csv":[
            {"genus_species":"Pieris rapae","order":"Lepidoptera",
             "plant_feeding":"Yes","total":"2"},
            {"genus_species":"Thymelicus lineola","order":"Lepidoptera",
             "plant_feeding":"Yes","total":"1"},
        ],
        "source_table_20230317.csv":[
            {"genus_species":"Pieris rapae","source":"Museum original"}
        ],
    }


def panel():
    return {
        "Pieris rapae":{"species":"Pieris rapae","host_family_count":"2",
             "resolved_host_species":"8","native_resource_units":"50",
             "introduced_added_units":"42","introduced_share_of_contemporary":"0.4"},
        "Thymelicus lineola":{"species":"Thymelicus lineola",
             "host_family_count":"1","resolved_host_species":"3",
             "native_resource_units":"80","introduced_added_units":"50",
             "introduced_share_of_contemporary":"0.3"},
        "Danaus plexippus":{"species":"Danaus plexippus",
             "host_family_count":"1","resolved_host_species":"4",
             "native_resource_units":"100","introduced_added_units":"20",
             "introduced_share_of_contemporary":"0.2"},
    }


def test_positive_nonnative_establishment_species_x_region_not_negative_arrival():
    result=m.analyze(data(),panel())
    assert result[
        "original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region"]==2
    assert result["original_239_butterflies_not_in_source_positive_list"]==1
    assert result[
        "original_239_nonlisted_species_cannot_be_called_failed_or_untransported"]
    entries=result["source_positive_original239_species"]
    p=next(x for x in entries if x["species"]=="Pieris rapae")
    assert p["source_documented_non_native_recipient_area_count"]==2
    assert p["source_earliest_documented_establishment_year"]==1860
    assert p["original_239_source_added_WCVP_resource_units"]==42
    assert p["source_attribute_plant_feeding_label"]=="Yes"
    assert result["positive_established_species_only_descriptive"][
        "status"]=="INSUFFICIENT_MATCHED_ESTABLISHMENT_TAXA"


def test_order_and_rank_cannot_turn_diptera_or_subspecies_into_butterflies():
    d=data()
    d["occurrence_table_20230317.csv"][3]["genus_species"]="Danaus plexippus"
    r=m.analyze(d,panel())
    assert r["original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region"]==2
    d["occurrence_table_20230317.csv"][3]["order"]="Lepidoptera"
    d["occurrence_table_20230317.csv"][3]["rank"]="subspecies"
    r=m.analyze(d,panel())
    assert r["original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region"]==2


def test_species_region_year_missing_kept_uncertain():
    d=data()
    d["occurrence_table_20230317.csv"][0]["year"]="NA"
    x=m.analyze(d,panel())
    p=next(z for z in x["source_positive_original239_species"]
           if z["species"]=="Pieris rapae")
    assert p["source_earliest_documented_establishment_year"]==1930
    assert p["original_source_region_and_year_records"][0][
        "year_source_first_report"] is None


def test_metadata_group_does_not_infer_plant_origin_from_insect_origin():
    r=m.analyze(data(),panel())
    assert "source_plant_status_NATIVE" not in str(r)
    assert "failed_intro" not in str(r)
