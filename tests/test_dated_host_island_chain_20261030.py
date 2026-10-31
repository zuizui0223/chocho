"""Fixed five prior-host island pairs: chronology and geography evidence screen."""
import importlib.util
import json
import pathlib
import pytest

ROOT=pathlib.Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/audit_dated_host_island_chain.py"
spec=importlib.util.spec_from_file_location("host_island_chain",SCRIPT)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_missing_followup_year_is_unknown_not_later_or_absent():
    end=module.year_observation(None,"42")
    assert end["first_year_in_2011_2017_followup_table"] is None
    assert end["first_year_not_in_followup_table_is_not_evidence_of_absence"] is True
    with pytest.raises(ValueError):
        module.canonical_pair("7","7")
    with pytest.raises(ValueError):
        module.year_observation({
            "first_record_year":"2015","followup_year_list":"2015",
            "any_followup_record":"True","record_in_year5":"False"
        },"27")


def test_real_prior_host_island_linkage_is_not_claimed_as_dispersal():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island source absent by design from original GEB anonymous bundle")
    r=module.audit(ROOT)
    assert len(r["links"])==5
    assert r["unique_island_pairs"]==4
    assert r["shorter_than_500km_port_representative_pairs_counting_butterfly"]==3
    assert r["both_record_years_2011_2017_available"]==1
    assert r["ordered_and_upstream_re_recorded_before_second"]==1
    by_species={x["butterfly"]:x for x in r["links"]}
    v=by_species["Vanessa cardui"]
    assert v["port_representative_great_circle_km"]>18000
    assert v["earlier_record_island_if_both_verified"]=="183"
    assert v["upstream_has_re_record_before_other_island_if_both_verified"] is True
    for sp in ["Colias croceus","Colias hyale"]:
        x=by_species[sp]
        assert 100<x["port_representative_great_circle_km"]<110
        assert not x["both_initial_record_years_available_in_followup_table"]
        assert x["earlier_record_island_if_both_verified"] is None
        assert x["unobserved_or_unknown_endpoint_year_not_inferred"] is True
    assert r["null_event_not_proof_of_date_2018_2022"] is True
    saved=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/dated_host_island_chain_screen.json"
    assert json.loads(saved.read_text(encoding="utf8"))==r
