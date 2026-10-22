"""Neotoma six-source source-age and chronology-scale integrity checks."""
import json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_six_dataset_age_scale_v15 import (
    source_previews,parse_dataset_source,align_records,run,csvrows
)

ROOT=Path(__file__).resolve().parents[1]
CAND=ROOT/"results/lgm_neotoma_idage_18_24ka_source_candidates_v13"
META=ROOT/"results/lgm_neotoma_Alnus_six_pollen_datasets_chronology_v14"

def input_table():
    return source_previews(csvrows(META/"all_six_dataset_chronology_publication_source_queries.csv"))

def test_six_source_dataset_and_chronology_citations_are_available_and_exact():
    groups=input_table()
    assert len(groups)==6 and all(len(q)==3 for q in groups.values())
    items={k:parse_dataset_source(k,v) for k,v in groups.items()}
    assert items["4396"]["site_name"]=="Plesheevo Lake"
    assert "provisional" in items["4396"]["first_chronology_notes"].lower()
    assert "SOURCE_SAYS_CHRONOLOGY_VERY_PROVISIONAL" in items["4396"]["metadata_flags"]
    assert "MULTIPLE_INCOMPATIBLE_PUBLISHED_AGE_SCALES" in items["24197"]["metadata_flags"]
    assert "MULTIPLE_INCOMPATIBLE_PUBLISHED_AGE_SCALES" in items["41351"]["metadata_flags"]
    assert "OLD_CALIBRATION_CURVE_SOURCE_REVIEW" in items["22662"]["metadata_flags"]
    assert "OLD_CALIBRATION_CURVE_SOURCE_REVIEW" in items["22683"]["metadata_flags"]
    assert all(q["publication_citation"] for q in items.values())

def test_exact_source_rows_27_expose_four_numeric_dataset_envelope_mismatches():
    items={k:parse_dataset_source(k,v) for k,v in input_table().items()}
    rows,counts=align_records(csvrows(CAND/"bounded_18_24ka_sample_taxon_proxy_checks.csv"),items)
    assert len(rows)==27
    assert sum(x["numeric_age_outside_any_published_range"] for x in counts.values())==4
    assert counts["14556"]["numeric_age_outside_any_published_range"]==3
    assert counts["22662"]["numeric_age_outside_any_published_range"]==1
    assert sum("SAMPLE_AGE_UNCERTAINTY_NOT_IN_OCCURRENCE_RESPONSE" in q["flags"] for q in rows)==26
    assert sum("DUPLICATE_OCCURRENCE_NOT_INDEPENDENT_OBSERVATION" in q["flags"] for q in rows)==1
    assert all(not r["sample_age_calibration_model_verified"] and not r["fossil_species_occurrence_validated"] for r in rows)

def test_dataset_range_mismatch_is_review_only_not_automatic_invalidity():
    r=run(type("A",(),{
        "candidate":CAND/"bounded_18_24ka_sample_taxon_proxy_checks.csv",
        "candidate_summary":CAND/"summary.json",
        "dataset_metadata":META/"all_six_dataset_chronology_publication_source_queries.csv",
        "dataset_summary":META/"summary.json"
    })())
    rows,sites,q=r
    assert q["status"]=="SOURCE_AGE_SCALE_CONFLICTS_REQUIRING_FOSSIL_CHRONOLOGY_AND_TAXONOMY_REVIEW"
    assert q["source_sample_calibrated_18_24ka_ages_verified_n"]==0
    assert q["fossil_occurrences_admitted_to_LGM_host_occupancy_model_n"]==0
    assert q["candidate_rows_outside_all_published_dataset_age_envelopes"]==4
    assert q["insect_genetic_response_opened"] is False

def test_wrong_record_identity_fails_closed():
    items={k:parse_dataset_source(k,v) for k,v in input_table().items()}
    source=csvrows(CAND/"bounded_18_24ka_sample_taxon_proxy_checks.csv")
    source[0]=dict(source[0],source_sample_taxonid="9999999")
    with pytest.raises(ValueError):align_records(source,items)
