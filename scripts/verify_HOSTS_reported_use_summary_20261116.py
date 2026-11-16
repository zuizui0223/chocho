#!/usr/bin/env python3
"""Recheck frozen HOSTS+WCVP+FAO location report estimates without original archives.

This is a verification of archived research outputs, NOT independent
reanalysis of the 140485-row HOSTS input; full rerun is a separate pinned CI.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"results/hosts_island_local_use_20261116/VERIFIED_SUMMARY_ALL24_20261116.json"
MANUSCRIPT=ROOT/"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"


def main():
    j=json.loads(RESULT.read_text(encoding="utf8"))
    assert j["schema"]=="chocho_HOSTS_source_locked_reported_use_summary_20261116_v1"
    assert j["status"]=="SOURCE_VERIFIED_DESCRIPTIVE_NOT_FIELD_CONSUMPTION"
    assert len(j["scopes"])==24
    src=j["source"]
    assert src["hosts_commit"]=="808e0b869f9ec1adf8efff87cf6a395adda103e0"
    assert src["wcvp_commit"]=="65bed76bae9d644ccb6ad200c05f9f5071d89e05"
    assert src["original_HOSTS_rows"]==140485
    assert src["resolved_origin_rows"]==21011
    assert src["model_source_workflow"]==43881198606
    def scope(tax, crop):
        return j["scopes"][f"{tax}__strict_exact_or_composite__{crop}"]
    def close(actual,expected,tol=1e-5):
        assert math.isclose(float(actual),float(expected),abs_tol=tol,rel_tol=0),(actual,expected)

    butterfly="HOSTS_Papilionoidea"
    base=scope(butterfly,"all_hosts")
    fao=scope(butterfly,"exclude_FAO_exact")
    wild=scope(butterfly,"exclude_FAO_exact_plus_genus_wildcards")
    allgenus=scope(butterfly,"exclude_FAO_all_crop_genera")
    assert (base["island"]["localities"],base["mainland"]["localities"])==(13,15)
    close(base["island"]["median_share"],.42424242424242425)
    close(base["mainland"]["median_share"],.09090909090909091)
    close(base["island_mainland_p"],.009879678872791467)
    assert (fao["island"]["localities"],fao["mainland"]["localities"])==(12,14)
    close(fao["island"]["median_share"],.20900537634408603)
    close(fao["mainland"]["median_share"],.08116883116883117)
    close(fao["island_mainland_p"],.03486742478298251)
    assert (wild["island"]["localities"],wild["mainland"]["localities"])==(12,14)
    close(wild["island"]["median_share"],.1973870056497175)
    assert (allgenus["island"]["localities"],allgenus["mainland"]["localities"])==(12,12)
    close(allgenus["island"]["median_share"],.17467581998474446)
    close(allgenus["mainland"]["median_share"],.07142857142857142)
    assert base["oceanic"]["localities"]==4 and base["continental"]["localities"]==5
    close(base["oceanic"]["median_share"],.6122093023255815)
    close(base["continental"]["median_share"],.10638297872340426)
    assert base["port_bootstrap_95"][0]<0<base["port_bootstrap_95"][1]
    assert fao["port_bootstrap_95"][0]<0<fao["port_bootstrap_95"][1]
    assert base["port_beta_logodds"]<0
    fixed239=scope("fixed239_butterflies","all_hosts")
    assert (fixed239["island"]["localities"],fixed239["mainland"]["localities"])==(9,6)
    assert fixed239["island_mainland_p"]>.05
    all_lep=scope("all_HOSTS_lepidoptera","all_hosts")
    assert all_lep["island"]["localities"]!=base["island"]["localities"]
    text=MANUSCRIPT.read_text(encoding="utf8")
    for marker in ("42.4%", "9.1%", "20.9%", "8.1%", "−0.128", "21,011",
                   "22.6%", "14.5%", "48.1%", "0.088",
                   "not an independent field study", "same HOSTS catalogue"):
        assert marker in text,marker
    assert "not uniformly field-verified feeding" in text
    print("PINNED_HOSTS_REPORT_RESULT_SUMMARY_AND_MANUSCRIPT_GATES_PASS")


if __name__=="__main__":
    main()
