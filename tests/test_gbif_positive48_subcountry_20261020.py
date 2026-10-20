"""No-network tests for 48 locked GBIF positive butterfly-country spatial facets."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from urllib.parse import urlparse,parse_qs
import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/audit_gbif_positive48_subcountry_20261020.py"
spec=importlib.util.spec_from_file_location("positive48_qc",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

STRATA=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json"
MANIFEST=ROOT/"data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json"
COUNTS=ROOT/"data/exploratory/gbif_country_20261020/gbif_country_counts.csv"


def positive():
    cells,src=m.load_locked_positive_cells(STRATA,MANIFEST,COUNTS)
    return cells,src


def sample(n=9):
    return json.dumps({"count":n,"limit":0,"results":[],
       "facets":[
          {"field":"GADM_LEVEL_1_GID","counts":[
             {"name":"DOM.1_1","count":6},
             {"name":"DOM.2_1","count":2}]},
          {"field":"LIFE_STAGE","counts":[{"name":"adult","count":1}]},
          {"field":"ESTABLISHMENT_MEANS","counts":[]}
       ]}).encode()


def test_frozen_48_and_two_pilots_and_filters():
    cells,src=positive()
    assert len(cells)==48
    assert src["source_distinct_species"]==20
    assert src["source_distinct_countries"]==23
    assert sum(c["added_resource_country"] for c in cells)==26
    pilots={ (x["species"],x["country"]):x for x in cells if (x["species"],x["country"]) in m.PILOT}
    assert len(pilots)==2
    assert pilots["Aguna asander","DO"]["source_late_count"]==9
    assert pilots["Papilio demoleus","IL"]["source_late_count"]==86
    for x in cells:
        query=parse_qs(urlparse(m.url_for(x)).query)
        assert query["taxonKey"]==[str(x["taxon_key"])]
        assert query["country"]==[x["country"]]
        assert query["year"]==["2010,2025"]
        assert set(query["basisOfRecord"])==set(m.BASIS)
        assert set(query["facet"])==set(m.FIELDS)
        assert query["limit"]==["0"]


def test_gadm1_counts_are_not_renormalized_to_mapped_part():
    v=m.parse_facets(sample())
    assert v["record_count"]==9
    assert v["facets"]["gadm1"]["assigned_records"]==8
    assert v["facets"]["gadm1"]["unassigned_records"]==1
    assert v["facets"]["life_stage"]["assigned_records"]==1
    assert v["facets"]["establishment_means"]["unassigned_records"]==9


def test_reject_overflow_duplicate_and_missing_gadm():
    bad=json.loads(sample())
    bad["facets"][0]["counts"][1]["count"]=20
    with pytest.raises(ValueError,match="exceeds"):
        m.parse_facets(json.dumps(bad).encode())
    bad=json.loads(sample())
    bad["facets"][0]["counts"].append(dict(bad["facets"][0]["counts"][0]))
    with pytest.raises(ValueError,match="duplicate"):
        m.parse_facets(json.dumps(bad).encode())
    bad=json.loads(sample())
    bad["facets"].pop(0)
    with pytest.raises(ValueError,match="No GADM1"):
        m.parse_facets(json.dumps(bad).encode())


def test_missing_optional_facets_never_are_counted_as_larvae():
    data=json.loads(sample())
    data["facets"]=data["facets"][:1]
    x=m.parse_facets(json.dumps(data).encode())
    assert x["missing_optional_facets"]==["establishment_means","life_stage"]
    assert "life_stage" not in x["facets"]


def test_positive_and_control_summaries_keep_original_country_cohort():
    cells,_=positive()
    p=cells[:2]
    tested=[{**x,"status":"COMPLETE_GADM1","record_count":x["source_pooled_count"],
             "source_pooled_delta":0,"facets":{"gadm1":{"categories":{"AAA.1_1":x["source_pooled_count"]},"assigned_records":x["source_pooled_count"],"unassigned_records":0}}} for x in p]
    r=m.summarize(tested,pilot=False)
    assert r["status"]=="PARTIAL_OR_INDEX_DRIFT_NOT_ZERO"
    assert r["gadm1_facets_complete"]==2
    assert r["cells_with_one_mapped_gadm_unit"]==2
    assert r["life_stage_metadata_records"]==0
