"""No-network source/denominator tests for original SA import laboratory Lepidoptera."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_sa_import_lepidoptera_detection_by_commodity_20261027 as m


def test_NA_not_tested_is_distinct_from_tested_0_and_positive():
    assert m.lepidoptera_result("NA")=="NOT_TESTED"
    assert m.lepidoptera_result("na")=="NOT_TESTED"
    assert m.lepidoptera_result("0")=="TESTED_NEGATIVE"
    assert m.lepidoptera_result("Plutella xylostella (Plutellidae)")=="DETECTED_LEPIDOPTERA"
    assert m.lepidoptera_result("Lepidoptera larva")=="DETECTED_LEPIDOPTERA"
    assert m.lepidoptera_result("")=="MISSING_ASSESSMENT"
    q=m.detection_count({"NOT_TESTED":7,"TESTED_NEGATIVE":10,"DETECTED_LEPIDOPTERA":2,
                         "MISSING_ASSESSMENT":1})
    assert q["original_raw_samples"]==20
    assert q["lepidoptera_tested_samples"]==12
    assert q["detected_Lepidoptera_sample_rows"]==2
    assert q["conditional_detection_per_100_TESTED"]==100*2/12


def test_transport_classes_are_nonexclusive_of_insect_stage():
    assert m.commodity_class("Cuttings: budwood")=="LIVE_PROPAGATION"
    assert m.commodity_class("In vitro plantlets")=="LIVE_PROPAGATION"
    assert m.commodity_class("Bulbs")=="LIVE_PROPAGATION"
    assert m.commodity_class("Plants")=="LIVE_PROPAGATION"
    assert m.commodity_class("Seeds")=="SEED_GRAIN"
    assert m.commodity_class("Grain")=="SEED_GRAIN"
    assert m.commodity_class("Beans")=="SEED_GRAIN"
    assert m.commodity_class("Nuts")=="SEED_GRAIN"
    assert m.commodity_class("Fresh fruit")=="FRESH_OR_DRY_PRODUCE"
    assert m.commodity_class("Fresh foliage")=="FRESH_OR_DRY_PRODUCE"
    assert m.commodity_class("Foreign commodity")=="OTHER_OR_UNRESOLVED"
    assert m.commodity_class("")=="OTHER_OR_UNRESOLVED"


def test_historical_reporting_era_does_not_mix_source_sample_years():
    assert m.sample_period("2001")=="1994-2004"
    assert m.sample_period("2019")=="2005-2019"
    assert m.sample_period("2023")=="MISSING_OR_OTHER_YEAR"
    assert m.sample_period("")=="MISSING_OR_OTHER_YEAR"


def fixture_sample(field,regime,commodity,year="2008",idnum="1"):
    return {
        "catalogNumber":idnum,"year":year,
        "sample.source":regime,"crop.commodity":commodity,
        "Insecta: Lepidoptera":field
    }


def test_binary_taxon_metadata_has_sufficiently_qualified_butterfly_identity():
    # A family stage label is NOT a butterfly species identity.
    rows=[
        {"order":"Lepidoptera","scientificName":"Plutella xylostella",
         "SA.occurenceStatus.historical":"absent",
         "SA.occurenceStatus.current":"present (2015)","Taxon":"Insects"},
        {"order":"Lepidoptera","scientificName":"Pieris rapae",
         "SA.occurenceStatus.historical":"absent",
         "SA.occurenceStatus.current":"absent","Taxon":"Insects"},
        {"order":"Coleoptera","scientificName":"Naupactus",
         "SA.occurenceStatus.historical":"absent",
         "SA.occurenceStatus.current":"present","Taxon":"Insects"},
    ]
    result=m.taxon_metadata(
        rows,{"Plutella xylostella (Plutellidae)","Pieris rapae",
              "Pyralidae larva","Lepidoptera larva","cf. Papilio demoleus"},
        {"Pieris rapae","Papilio demoleus"}
    )
    assert result["Lepidoptera_metadata_taxa_rows"]==2
    assert result["source_later_established_Lepidoptera_rows"][
        "reported_historically_absent_and_later_present_count"]==1
    assert result["exact_positive_species_names_matching_frozen_239_butterflies"]==["Pieris rapae"]
    assert "Lepidoptera larva" not in result[
        "unambiguous_positive_names_matching_original_metadata"]


def test_life_stage_inference_is_forbidden_when_only_taxon_positive():
    x=m.detection_count({"DETECTED_LEPIDOPTERA":1,"TESTED_NEGATIVE":99,
                         "NOT_TESTED":900})
    assert x["lepidoptera_tested_samples"]==100
    assert x["conditional_detection_per_100_TESTED"]==1
    assert x["original_raw_samples"]==1000
    assert "wild_established" not in x


def test_source_taxa_and_consignment_records_not_interchangeable():
    assert m.SCHEMA.startswith("chocho_SA_original_plant_import_lepidoptera")
    assert m.SOURCE_SHA["sample"]!=m.SOURCE_SHA["metadata"]
    assert m.URLS["sample"].endswith("25166393")
    assert m.URLS["metadata"].endswith("26235020")
    assert m.EXPECTED_STATUS=={
        "TESTED_NEGATIVE":15673,"NOT_TESTED":10443,
        "DETECTED_LEPIDOPTERA":178,"MISSING_ASSESSMENT":1
    }
