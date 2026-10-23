"""No-network privacy, evidence and source-drift checks for original 39 GBIF juvenile labels."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs,urlparse
import pytest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("junhost",
    ROOT/"scripts/audit_gbif_juvenile_host_annotations_20261023.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
BASE=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"

def original():
    return m.source_cohort(
        BASE/"gbif_juvenile_gadm1_conditional_20261021.json",
        BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json",
        BASE/"gbif_source_introduced_hosts_five_cases_20261021.json")

def fake_response(q):
    data=[{"key":i+1,"countryCode":q["country"],"lifeStage":q["stage"],
           "associatedTaxa":"Larva feeding on Triphasia trifolia"} for i in range(q["expected_count"])]
    return json.dumps({"offset":0,"limit":100,"count":len(data),
                        "endOfRecords":True,"results":data}).encode()

def test_frozen_39_queries_and_safe_response_url():
    q,meta=original()
    assert len(q)==15 and sum(x["expected_count"] for x in q)==39
    assert len({(x["species"],x["country"]) for x in q})==11
    assert sum(bool(x["candidate_names"]) for x in q)>0
    for x in q:
        z=parse_qs(urlparse(m.query_url(x)).query)
        assert z["taxonKey"]==[str(x["taxon_key"])]
        assert z["country"]==[x["country"]]
        assert z["year"]==["2010,2025"]
        assert z["lifeStage"]==[x["stage"]]
        assert z["limit"]==["100"] and z["offset"]==["0"]
        assert set(z["basisOfRecord"])==set(m.BASIS)

def test_exact_source_plant_mentions_not_raw_comment_leak():
    q,_=original()
    source=next(x for x in q if x["species"]=="Papilio demoleus" and x["country"]=="DO" and x["stage"]=="Larva")
    r=m.inspect_response(source,fake_response(source))
    assert r["frozen_stage_count"]==source["expected_count"]
    assert r["records_with_exact_source_candidate_binomial_mention"]==source["expected_count"]
    assert r["records_with_candidate_mention_and_host_related_words_unverified"]==source["expected_count"]
    assert "associatedTaxa" in r["field_nonempty_record_counts"]
    assert "observer" not in json.dumps(r).lower()
    assert "Larva feeding" not in json.dumps(r)

def test_source_drift_and_incomplete_paging_are_not_negative_evidence():
    q,_=original()
    row=q[0]
    with pytest.raises(ValueError,match="drift"):
        m.inspect_response(row,json.dumps({"offset":0,"limit":100,"count":0,
                           "endOfRecords":True,"results":[]}).encode())
    d=json.loads(fake_response(row))
    d["endOfRecords"]=False
    with pytest.raises(ValueError):
        m.inspect_response(row,json.dumps(d).encode())

def test_no_associated_fields_does_not_imply_no_host():
    q,_=original()
    x=q[0];d=json.loads(fake_response(x))
    for r in d["results"]:r.pop("associatedTaxa")
    o=m.inspect_response(x,json.dumps(d).encode())
    assert o["records_with_any_selected_text_field"]==0
    assert o["records_with_exact_source_candidate_binomial_mention"]==0
