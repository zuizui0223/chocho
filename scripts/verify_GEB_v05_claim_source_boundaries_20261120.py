#!/usr/bin/env python3
"""Fail-closed GEB v0.5 main-claim versus auxiliary GBIF facet-source scope gate.

Confirms manuscript main figure assertions do NOT accidentally borrow from
the incomplete 1596-port botanical facet source or new 2026-11-20 index,
while retaining distinct 155/224 original download citation and genuine
source-reproducible fixed239 island opportunity analysis.

This is a source-dependency guard, NOT an ecological replication, DOI mint,
journal-submission or anonymous data archive.
"""
from __future__ import annotations
import json
from pathlib import Path

R = Path(__file__).resolve().parents[1]
M = R / "manuscript/integrated_butterfly_resource_access_20261109"
ARTICLE = M / "ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
AUDIT = M / "GEB_V05_EVIDENCE_DEPENDENCY_AND_DATA_ACCESS_HANDBACK_20261120.md"
BASE = R / "results/geb_v05_gbif_source_attribution_20261120"
PORTS = BASE / "NEW_GBIF_ORIGINAL1596_PORTS_DATASETKEY_INDEX_20261120.json"
RECENT = BASE / "PORT1596_SOURCE_EFFORT_DIAGNOSTIC_20261120.json"
SHARED = BASE / "VERIFIED_GBIF_155_224_ORIGINAL_DOWNLOAD_PUBLISHER_COUNTS.json"

def section(text: str, start: str, end: str) -> str:
    begin=text.find(start)
    finish=text.find(end,begin+len(start))
    if begin<0 or finish<0 or finish<=begin:
        raise AssertionError(f"Expected manuscript section missing: {start}")
    return text[begin:finish]

def audit():
    article=ARTICLE.read_text(encoding="utf-8")
    memo=AUDIT.read_text(encoding="utf-8")
    port=json.loads(PORTS.read_text(encoding="utf-8"))
    recent=json.loads(RECENT.read_text(encoding="utf-8"))
    gbif=json.loads(SHARED.read_text(encoding="utf-8"))
    if "Ports, the sea and two layers of island floras" not in article or "Coauthor review draft; unsubmitted" not in article:
        raise AssertionError("Manuscript identity / still-draft status drift")
    required={
       "Figure 1":"Fig1_resource_geography_scales_v03.svg",
       "Figure 2":"Fig2_ports_native_alien_floras_isolation_v05.svg",
       "Figure 3":"Fig3_global_host_contributors_at_ports_restructured_v04.svg",
       "Figure 4":"Fig4_independent_butterfly_record_predictability_v04.svg"
    }
    for name,fig in required.items():
        if fig not in article:
            raise AssertionError(name+" source-backed illustration link drift")
    for fragment in ["239 source-eligible butterflies × 369 WGSRPD3","126 botanical islands",
                     "237 mainland","155 butterfly species","224 butterflies",
                     "1,596 port buffers","eight of these ten leading species"]:
        if fragment not in article:
            raise AssertionError("Original exact source population/status changed: "+fragment)
    methods=section(article,"### 2.3. Port associations","### 2.4. Plant identities")
    fig2=section(article,"### 3.2. Ports","### 3.3.")
    fig3=section(article,"### 3.3.","### 3.4.")
    fig4=section(article,"### 3.4.","### 4.")
    for core in ["WPI port counts","whole-flora richness","MEM10","potential food geography"]:
        if core not in methods:
            raise AssertionError("Main WCVP/WPI/239 host opportunities source independence changed: "+core)
    if "1,596" not in methods or "Separate surveys" not in methods:
        raise AssertionError("Original 1596 botanical survey not separated in methods")
    if "whole introduced flora richness" not in fig2 or "Mainland isolation" not in fig2:
        raise AssertionError("Main Figure2 WCVP/WPI botanical-port and fixed239 isolation no longer distinguished")
    if "separate historical port-buffer observation" not in fig3 or "1,596" not in fig3:
        raise AssertionError("Auxiliary historical GBIF port data not separated from Fig3 main HOSTS/WCVP")
    if "155 butterflies" not in fig4 or "224 GBIF-indexed butterflies" not in fig4:
        raise AssertionError("Independent consumer cohorts lost in paper")
    if "10.15468/dl.umqmsf" not in article:
        raise AssertionError("Original GBIF source citation missing")
    if not all(v in memo for v in ("Fig.2A","Fig.2D","Fig.3の","Fig.4",
                                  "1,596港", "142群", "未閉鎖", "未作成")):
        raise AssertionError("Source dependency memo lacks established cohort and release gates")
    if "155 and 224 modelled species" not in memo or "same original GBIF occurrence download" not in memo:
        raise AssertionError("Original shared 155/224 source DOI not accurately described")
    if "1,482" not in memo or "3,138,675" not in memo:
        raise AssertionError("Separate NEW GBIF 1596 provider source results missing from scope memo")
    if any(v in article for v in ("3,138,675","1,482 datasetKey", "PORT1596_SOURCE_EFFORT_DIAGNOSTIC")):
        raise AssertionError("New different-snapshot and exploratory provider result silently introduced into manuscript")
    if len(port["new_source_port_index_facets"])!=1596 or port["new_source_datasetKey_contributors"]!=1482:
        raise AssertionError("Original/source newer 1596 port source sample drift")
    if port["source_original_2026_GBIF_contributing_datasetKeys_recovered"]:
        raise AssertionError("New provider index falsely attributed to original source")
    if recent["n_orig_WPI_complete_ports"]!=1596 or recent["new_GBIF_current_zero_index_ports"]!=135:
        raise AssertionError("Auxiliary observation source sample drift")
    if recent["status"]!="EXPLORATORY_CROSS_SNAPSHOT_SOURCE_EFFORT_DIAGNOSTIC_NOT_HISTORIC_PROVENANCE_OR_MAIN_WCVP_RESULT":
        raise AssertionError("Original vs current source scope warning erased")
    if gbif["original_DOI"]!="10.15468/dl.umqmsf" or gbif["provider_count"]!=1716 or gbif["contributors_record_sum_if_available"]!=14637152:
        raise AssertionError("Independent consumer GBIF original shared download DOI/publisher record attribution drift")
    print("GEB_V05_MAIN_SOURCE_SCOPE_AND_AUX_1596_GBIF_CITATION_FIREWALL_PASS")
    print(json.dumps({"figures":len(required),"current_GBIF_port_datasetKeys":1482,
                      "original_GBIF155_224_publishers":1716,
                      "observational_1596_source_is_auxiliary_only":True,
                      "derived_GBIF_DOI_invented":False,
                      "journal_anonymous_source_archive_approved":False},ensure_ascii=False))

if __name__=="__main__":
    audit()
