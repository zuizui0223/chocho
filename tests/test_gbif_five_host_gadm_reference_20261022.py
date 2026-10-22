"""No-network exact hypergeometric and conservative candidate-plant checks."""
from __future__ import annotations
import importlib.util
from fractions import Fraction
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=ROOT/"scripts/analyze_gbif_five_host_gadm_reference_20261022.py"
spec=importlib.util.spec_from_file_location("local_gadm",s)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def row(N,H,j,O):
    return {"all_butterfly_gadm1_records":N,
            "all_butterfly_records_in_plant_gadm1":H,
            "immature_source_labels":j,
            "immature_labels_in_plant_detected_gadm1":O}

def test_exact_conditional_hypergeom_matches_toy():
    r=[row(4,2,1,1)]
    d=m.exact_reference(r)
    assert d==[Fraction(1,2),Fraction(1,2)]
    v=m.summarize_reference(r)
    assert v["null_expected_labels"]==.5
    assert v["observed_minus_reference_pp"]==50
    assert v["conditional_finite_label_reassignment_upper_tail_NOT_ecological_p"]==.5

def test_two_stage_finite_margins_are_independently_conditional():
    r=[row(4,2,1,1),row(4,2,1,0)]
    d=m.exact_reference(r)
    assert d==[Fraction(1,4), Fraction(1,2),Fraction(1,4)]
    v=m.summarize_reference(r)
    assert v["observed_minus_reference_pp"]==0
    assert sum(d)==1

def test_reject_ineligible_sources_and_invalid_hypergeom():
    with pytest.raises(ValueError):
        m.analyse({"schema":"changed"},{},{})
    with pytest.raises(ValueError):
        m.exact_reference([row(3,4,1,1)])
    with pytest.raises(ValueError):
        m.exact_reference([row(3,2,4,1)])


def test_archived_gadm_five_source_replays_exactly_without_network():
    import json
    base=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"
    plant=json.loads((base/"gbif_five_host_plant_gadm_source_20261021.json").read_text(encoding="utf8"))
    stage=json.loads((base/"gbif_juvenile_gadm1_conditional_20261021.json").read_text(encoding="utf8"))
    butterfly=json.loads((base/"gbif_positive48_gadm_lifestage_facets_20261020.json").read_text(encoding="utf8"))
    saved=json.loads((base/"gbif_five_host_gadm_conditioned_reference_20261022.json").read_text(encoding="utf8"))
    out=m.analyse(plant,stage,butterfly)
    for key in ("schema","status","source_candidate_host_country_links",
                "plant_country_links_with_at_least_one_indexed_record",
                "source_indexed_botanical_records_across_20_links",
                "source_juvenile_label_count_in_five_cases",
                "source_juvenile_labels_in_candidate_plant_positive_gadm1",
                "posthoc_label_allocation_reference","cases"):
        assert out[key]==saved[key]
    assert out["exact_gbif_accepted_botanical_species"]==8
    assert out["plant_country_links_with_at_least_one_indexed_record"]==9
    assert out["source_indexed_botanical_records_across_20_links"]==38
    assert len(out["cases"])==5
    assert sum(r["immature_labels_in_plant_detected_gadm1"] for r in out["cases"])==16
    assert abs(out["posthoc_label_allocation_reference"]["all_five"]["conditional_finite_label_reassignment_upper_tail_NOT_ecological_p"]-0.2291464884)<1e-10
    assert out["posthoc_label_allocation_reference"]["excluding_papilio_dominican_republic"]["observed_minus_reference_pp"]<0
