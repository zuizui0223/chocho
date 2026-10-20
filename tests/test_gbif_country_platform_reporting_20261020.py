"""Check 11 frozen country/platform denominators without provider network calls."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/audit_gbif_country_platform_reporting_20261020.py"
spec=importlib.util.spec_from_file_location("reporting",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
FIXED=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_dataset_year_matrix_20261020.json"


def response(years):
    return json.dumps({"count":sum(years.values()),"limit":0,"results":[],
      "facets":[{"field":"YEAR","counts":[{"name":y,"count":n} for y,n in years.items()]}]}).encode()


def test_source_group_selection_and_exact_query_nested_filter():
    rows,meta=m.read_locked_dataset_year(FIXED)
    assert len(rows)==meta["source_group_count"]==11
    assert all(2018<=int(y)<=2025 for r in rows for y in r["source_years"])
    assert len([r for r in rows if r["country"]=="IL" and r["species"]=="Papilio demoleus"])==3
    for r in rows:
        for scope in m.FACETS:
            q=parse_qs(urlparse(m.url_for(r,scope)).query)
            assert q["datasetKey"]==[r["dataset_key"]]
            assert q["country"]==[r["country"]]
            assert q["year"]==["2018,2025"]
            assert q["limit"]==["0"] and q["facet"]==["year"]
            assert set(q["basisOfRecord"])==set(m.BASIS)
            if scope=="lepidoptera_legacy_797":
                assert q["orderKey"]==["797"]
            else:
                assert "orderKey" not in q


def test_facet_contract_and_taxonomic_nesting_gate():
    assert m.parse_years(response({"2024":7,"2025":12}))=={
        "index_count":19,"year_counts":{"2024":7,"2025":12}}
    a={"status":"COMPLETE","index_count":30,"year_counts":{"2024":20,"2025":10}}
    assert m.assess_nested({"2024":9,"2025":33},a)["status"]=="DENOMINATOR_NOT_NESTED_DO_NOT_INTERPRET"
    a["year_counts"]["2025"]=100
    assert m.assess_nested({"2024":9,"2025":33},a)["status"]=="NESTED_VALID"
    assert m.assess_nested({"2024":9},{"status":"UNAVAILABLE_NOT_ZERO"})["status"]=="UNAVAILABLE_NOT_ZERO"
    with pytest.raises(ValueError):
        m.parse_years(b"{")
    bad=json.loads(response({"2024":7}))
    bad["facets"][0]["counts"][0]["count"]=8
    with pytest.raises(ValueError, match="do not match"):
        m.parse_years(json.dumps(bad).encode())
    bad=json.loads(response({"2024":7}))
    bad["facets"][0]["counts"].append({"name":"2024","count":1})
    with pytest.raises(ValueError):
        m.parse_years(json.dumps(bad).encode())


def test_comparison_requires_same_provider_with_two_positive_years():
    rows,_=m.read_locked_dataset_year(FIXED)
    source_by_key={(r["dataset_key"],r["country"]):r for r in rows}
    # Synthetic reporting denominators reproduce the positive source counts,
    # but imply constant focal share in both years, therefore ratio=1.
    fetched={}
    for r in rows:
        for kind in m.FACETS:
            # All source years -> multiply by 10 to preserve nesting.
            v={y:int(n*10) for y,n in r["source_years"].items()}
            fetched[(r["dataset_key"],r["country"],kind)]={
                "status":"COMPLETE","index_count":sum(v.values()),"year_counts":v}
    report=m.assemble(rows,fetched,{"test":True})
    assert report["status"]=="ALL_BROAD_DENOMINATORS_VALID"
    for provider in ("INPA","iNaturalist"):
        for scope in m.FACETS:
            x=report["focal_summary"]["israel_papilio"][provider][scope]
            assert x["relative_focal_rate_2025_to_2024"]==1
