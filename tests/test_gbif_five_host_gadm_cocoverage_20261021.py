"""No-network checks for GBIF plant GADM1 source resolution and co-detection."""
from __future__ import annotations
import importlib.util,json
from pathlib import Path
from urllib.parse import parse_qs,urlparse
import pytest
ROOT=Path(__file__).resolve().parents[1]
f=ROOT/"scripts/audit_gbif_five_host_gadm_cocoverage_20261021.py"
spec=importlib.util.spec_from_file_location("cocoverage",f)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
BASE=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"
def test_locked_five_and_8_wcvp_plants_no_outcome_filter():
    cases,names,pairs,j,a,meta=m.load_fixed(BASE/"gbif_source_introduced_hosts_five_cases_20261021.json",
        BASE/"gbif_juvenile_gadm1_conditional_20261021.json",
        BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json")
    assert len(cases)==5 and len(names)==8 and len(pairs)==20
    assert meta["source_focal_case_n"]==5
    assert {r["species"] for r in pairs}=={"Anteos clorinde","Colobura dirce","Papilio demoleus"}
def test_species_resolver_does_not_accept_synonyms_or_fuzzy_matches():
    base={"status":"SOURCE_HTTP_OK","url":"https://api.gbif.org","raw_sha256":"0",
        "retrieved_utc":"2026-10-21","json":{"usageKey":123,"matchType":"EXACT",
        "rank":"SPECIES","status":"ACCEPTED","kingdom":"Plantae","canonicalName":"Cassia fistula"}}
    assert m.strict_match("Cassia fistula",base)["accepted_species_key"]==123
    assert m.strict_match("Cassia fistula",{**base,"json":{**base["json"],"matchType":"FUZZY"}})["accepted_species_key"] is None
    assert m.strict_match("Cassia fistula",{**base,"json":{**base["json"],"canonicalName":"Cassia acutifolia"}})["accepted_species_key"] is None
def test_source_filters_gadm_facets_and_unmapped_propagation():
    pair={"country":"DO","dataset_key":None}
    q=parse_qs(urlparse(m.occurrence_url(pair,987)).query)
    assert q["taxonKey"]==["987"] and q["country"]==["DO"]
    assert q["year"]==["2010,2025"] and q["facet"]==["gadmLevel1Gid"]
    assert set(q["basisOfRecord"])==set(m.BASIS)
    raw={"status":"SOURCE_HTTP_OK","retrieved_utc":"2026-10-21","url":"x","raw_sha256":"foo",
        "json":{"count":3,"limit":0,"results":[],"facets":[{"field":"GADM_LEVEL_1_GID","counts":[{"name":"DOM.1_1","count":2}]}]}}
    a=m.parse_plant_gadm(raw)
    assert a["indexed_country_plant_records"]==3
    assert a["gadm1_mapped_records"]==2 and a["gadm1_unmapped_records"]==1
    assert m.parse_plant_gadm({**raw,"json":{**raw["json"],"count":0,"facets":[]}})["indexed_country_plant_records"]==0
