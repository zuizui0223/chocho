"""No-network source-locked juvenile metadata/GADM1 contract tests."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse
import pytest

ROOT=Path(__file__).resolve().parents[1]
sp=ROOT/"scripts/audit_gbif_juvenile_gadm1_20261021.py"
spec=importlib.util.spec_from_file_location("juvenile_gadm",sp)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
source=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_gadm_lifestage_facets_20261020.json"

def sample():
    return json.dumps({"count":2,"limit":0,"results":[],
      "facets":[{"field":"GADM_LEVEL_1_GID","counts":[
          {"name":"DOM.1_1","count":1},{"name":"DOM.2_1","count":1}]}]}).encode()

def test_exact_original_locked_cohort_and_filters():
    req,cells,meta=m.cohort(source)
    assert len(req)==15 and len(cells)==48
    assert len({(r["species"],r["country"]) for r in req})==11
    assert sum(r["expected_stage_count"] for r in req)==39
    for r in req:
        qs=parse_qs(urlparse(m.make_url(r)).query)
        assert qs["taxonKey"]==[str(r["taxon_key"])]
        assert qs["country"]==[r["country"]]
        assert qs["year"]==["2010,2025"]
        assert qs["lifeStage"]==[r["life_stage"]]
        assert qs["facet"]==["gadmLevel1Gid"]
        assert qs["limit"]==["0"]
        assert set(qs["basisOfRecord"])==set(m.BASIS)

def test_complete_gadm_and_unmapped_counts():
    v=m.parse_gadm(sample())
    assert v["record_count"]==2 and v["gadm1_mapped"]==2
    assert v["gadm1_unmapped"]==0
    v=json.loads(sample())
    v["facets"][0]["counts"].pop()
    x=m.parse_gadm(json.dumps(v).encode())
    assert x["gadm1_unmapped"]==1

def test_facet_error_and_subset_fail_closed():
    a=json.loads(sample())
    a["facets"][0]["counts"].append(dict(a["facets"][0]["counts"][0]))
    with pytest.raises(ValueError,match="duplicate"):
        m.parse_gadm(json.dumps(a).encode())
    a=json.loads(sample())
    a["facets"][0]["counts"][0]["count"]=8
    with pytest.raises(ValueError,match="exceed"):
        m.parse_gadm(json.dumps(a).encode())
    req,cells,meta=m.cohort(source)
    fake=[{**r,"status":"API_UNAVAILABLE_NOT_ZERO"} for r in req]
    assert m.summarize(fake,cells,meta)["status"]=="INCOMPLETE_STAGE_GADM_DATA_NOT_ZERO"

def test_replay_requires_same_species_country_stage_identities(tmp_path):
    req,pairs,meta=m.cohort(source)
    fake=[{**r,"status":"API_UNAVAILABLE_NOT_ZERO"} for r in req]
    z=m.summarize(fake,pairs,meta)
    assert z["status"]=="INCOMPLETE_STAGE_GADM_DATA_NOT_ZERO"
    assert len(z["attempts"])==15

def test_full_saved_stage_conditional_source_matches_frozen_48():
    from pathlib import Path
    import json
    base=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"
    expected=json.loads((base/"gbif_juvenile_gadm1_conditional_20261021.json").read_text(encoding="utf-8"))
    qs,full,meta=m.cohort(source)
    replay=m.summarize(expected["attempts"],full,meta)
    assert replay["status"]=="COMPLETE_CONDITIONAL_JUVENILE_GADM1"
    assert replay["matched_conditionals"]==15
    assert replay["resolved_juvenile_records"]==replay["gadm1_mapped_juvenile_records"]==39
    assert replay["cell_gadm1_juvenile_summary"]==expected["cell_gadm1_juvenile_summary"]
    assert all(x["cross_category_gadm_mismatch"]=={}
               for x in replay["cell_gadm1_juvenile_summary"])
