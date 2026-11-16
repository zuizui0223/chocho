#!/usr/bin/env python3
"""Keep historical HOSTS positives and exact flora-normalization limits in GEB."""
from __future__ import annotations

import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ARCH=ROOT/"results/hosts_use_vs_flora_20261116/VERIFIED_FLORISTIC_AND_FAMILY_BASELINE_20261116.json"
GEB=ROOT/"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"

def close(x,y,t=.00001):
    assert abs(x-y)<t,(x,y)

def main():
    j=json.loads(ARCH.read_text(encoding="utf8"))
    assert j["schema"]=="chocho_HOSTS_vs_whole_family_matched_WCVP_baseline_20261116_v1"
    assert j["status"]=="SOURCE_VERIFIED_OBSERVATIONAL_TAXON_OPPORTUNITY_BASELINE_NOT_FIELD_PREFERENCE"
    assert (j["source_run"],j["regions"],j["accepted_WCVP_species_region_records"])==(43887252612,139,809162)
    assert j["accepted_WCVP_global_species"]==360064
    assert j["original_WCVP_island_flora_validation"]["matching_island_regions"]==44
    assert j["original_WCVP_island_flora_validation"]["mismatches"]==0
    assert len(j["scopes"])==12
    b=j["scopes"]["HOSTS_Papilionoidea__all_hosts__min10"]
    c=j["scopes"]["HOSTS_Papilionoidea__exclude_FAO_exact__min10"]
    assert (b["n"],c["n"])==(24,22)
    assert (b["island"]["localities"],b["mainland"]["localities"])==(12,12)
    assert (c["island"]["localities"],c["mainland"]["localities"])==(11,11)
    def contrast(s,k):
        return s["contrasts"][k]
    close(contrast(b,"raw_report_share")["median_difference"],.3465160075329567)
    close(contrast(c,"raw_report_share")["median_difference"],.1382488479262673)
    close(contrast(b,"whole_flora_adjusted_excess")["median_difference"],.13885707075087386)
    close(contrast(c,"whole_flora_adjusted_excess")["median_difference"],.005594400050484519)
    close(contrast(b,"family_flora_adjusted_excess")["median_difference"],.07142381346699858)
    close(contrast(c,"family_flora_adjusted_excess")["median_difference"],-.020796688108226935)
    for k in ("whole_flora_adjusted_excess","family_flora_adjusted_excess"):
        lo,hi=contrast(c,k)["descriptive_locality_bootstrap95"]
        assert lo<0<hi,(k,lo,hi)
    assert contrast(c,"family_flora_adjusted_excess")["exploratory_mannwhitney_p"]>.9
    f=j["scopes"]["fixed239_butterflies__exclude_FAO_exact__min10"]
    assert f["n"]==13 and f["island"]["localities"]==8
    t=GEB.read_text(encoding="utf8")
    for v in ("+118.9%","+54.9%","4.12%","42.4%","9.1%","20.9%","8.1%",
              "+7.14","−2.08","−10.04","+12.69","809,162","44 matched",
              "not measures of host plant abundance","do not establish a surplus"):
        assert v in t,("Manuscript text omitted source claim/boundary",v)
    ab=t.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    import re
    words=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",ab))
    assert words<=300,words
    print("HOSTS_FLORA_FAMILY_SOURCE_AND_GEB_NONCAUSAL_CLAIM_VERIFIED",words)

if __name__=="__main__":
    main()
