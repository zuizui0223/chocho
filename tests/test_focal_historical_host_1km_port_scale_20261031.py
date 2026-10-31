"""Validate earlier indexed butterfly host plants at different WPI buffer scales.

Biological scope: 3 focal butterflies, 2 named hosts, 4 host×island
endpoints. An indexed plant speciesKey within 5km is NOT evidence of a
port-local 1km source, original shipping movement or butterfly breeding.
"""
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"
F=D/"focal_host_5km_to_1km_source_audit.json"
CODE=D/"audit_focal_host_spatial_scale.py"

def test_named_5km_historical_host_vs_harbor_1km_records():
    if not F.is_file() or not CODE.is_file():
        pytest.skip("Historic island/port exploratory branch omitted from frozen anonymous GEB bundle")
    source=json.loads(F.read_text(encoding="utf8"))
    spec=importlib.util.spec_from_file_location("chocho_focal_historical_host_1km",CODE)
    assert spec is not None and spec.loader is not None
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    replay=module.analyze()
    assert replay==source
    assert source["status"]=="TEMPORALLY_PRIOR_5KM_HOST_INDEX_DOES_NOT_ESTABLISH_PORT_LOCAL_1KM_HOST"
    assert source["named_focal_butterfly_species"]==3
    assert source["distinct_candidate_hosts"]==2
    assert source["unique_host_by_island_endpoints"]==4
    assert source["butterfly_species_island_endpoints_counting_shared_hosts"]==6
    assert all(q["5km_prior_host_source_indexed"] is True
               for q in source["site_index_summary"])
    assert sum(q["1km_port_site_indexed_count"] for q in source["site_index_summary"])==0
    assert sum(q["1km_offport_site_indexed_count"] for q in source["site_index_summary"])==2
    # Medicago at the two earlier indexed islands has no 1km source GBIF
    # record in the present sample: do NOT recode these as biological absences.
    alfalfa=[q for q in source["site_index_summary"] if q["host"]=="Medicago sativa"]
    assert len(alfalfa)==2
    assert {q["island"] for q in alfalfa}=={"169","276"}
    assert all(q["1km_port_site_indexed_count"]==0 and
               q["1km_offport_site_indexed_count"]==0 for q in alfalfa)
    alder=[q for q in source["site_index_summary"] if q["host"]=="Alnus incana"]
    assert len(alder)==2
    assert next(q for q in alder if q["island"]=="276")["1km_offport_site_indexed_count"]==2
    assert next(q for q in alder if q["island"]=="715")["1km_offport_site_indexed_count"]==0
    assert all(q["1km_port_site_indexed_count"]==0 for q in alder)
    assert all(r["one_km_nonindex_NOT_biological_absence"] is True
               for r in source["butterfly_island_endpoints"])
    assert all(r["host_index_lead_years"]>=14 and r["host_index_lead_years"]<=19
               for r in source["butterfly_island_endpoints"])
    assert any("not imply host plant absent" in x for x in source["limits"])
    assert any("NOT independent replicates" in x for x in source["limits"])
