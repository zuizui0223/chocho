#!/usr/bin/env python3
"""Source-lock the species-paired HOSTS island/mainland section of one GEB paper."""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"results/hosts_use_vs_flora_20261116/PAIRED_SAME_BUTTERFLY_FULL_SOURCE_20261117.json"
PAPER=ROOT/"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"


def close(a,b,tol=1e-9):
    if abs(a-b)>tol: raise AssertionError((a,b))


def main():
    s=json.loads(RESULT.read_text(encoding="utf8"))
    assert s["schema"]=="chocho_same_butterfly_in_island_mainland_reported_host_20261117_v1"
    assert s["status"]=="SOURCE_PAIRED_LITERATURE_REPORTS_NOT_CONFIRMED_WILD_HOST_SWITCHING"
    assert s["source_original"]["fixed_focal_butterflies"]==239
    assert s["source_original"]["HOSTS_Git_commit"]=="808e0b869f9ec1adf8efff87cf6a395adda103e0"
    assert s["source_original"]["WCVP_Git_commit"]=="65bed76bae9d644ccb6ad200c05f9f5071d89e05"
    assert len(s["cohorts"])==8
    all_=s["cohorts"]["HOSTS_Papilionoidea__all_hosts"]
    non=s["cohorts"]["HOSTS_Papilionoidea__exclude_FAO_exact"]
    fixed=s["cohorts"]["fixed239_butterflies__exclude_FAO_exact"]
    assert all_["shared_species_island_and_mainland"]==67
    assert non["shared_species_island_and_mainland"]==57
    assert fixed["shared_species_island_and_mainland"]==20
    assert all_["unique_insect_host_original_location_records"]==1928
    assert non["unique_insect_host_original_location_records"]==1729
    for c in (all_,non,fixed):
        n=c["shared_species_island_and_mainland"]
        assert len(c["source_matched_species_complete"])==n
        p=c["paired_comparison"]
        for x in p["difference_raw"],p["difference_family_excess"],p["difference_family_expectation"]:
            assert x["n_shared_butterflies"]==n
        close(p["difference_raw"]["mean_paired_difference"],
              p["difference_family_excess"]["mean_paired_difference"]+
              p["difference_family_expectation"]["mean_paired_difference"])
        for k in ("difference_raw","difference_family_excess"):
            lo,hi=p[k]["bootstrap95_species_level_mean_descriptive"]
            assert lo<0<hi,(c["crop"],k,lo,hi)
    assert all_["same_insect_accepted_host_indexed_both"]==29
    assert non["same_insect_accepted_host_indexed_both"]==19
    orig={
        (z["butterfly"],z["accepted_host_id"]) for z in all_["source_shared_host_pairs"]
        if z["source_pattern"]=="island_introduced_mainland_native"
    }
    subset=[
        z for z in non["source_shared_host_pairs"]
        if z["source_pattern"]=="island_introduced_mainland_native"
    ]
    assert len(orig)==13 and len(subset)==11
    assert {(z["butterfly"],z["accepted_host_id"]) for z in subset} <= orig
    assert len({z["butterfly"] for z in subset})==10
    assert len({z["accepted_host_id"] for z in subset})==7
    assert sum(z["mainlands"]==["Thailand"] for z in subset)==10
    assert sum(z["mainlands"]==["Finland"] for z in subset)==1
    all_genus=s["cohorts"]["HOSTS_Papilionoidea__exclude_FAO_all_crop_genera"]
    strict=[z for z in all_genus["source_shared_host_pairs"]
            if z["source_pattern"]=="island_introduced_mainland_native"]
    assert all_genus["shared_species_island_and_mainland"]==47
    assert all_genus["same_insect_accepted_host_indexed_both"]==15
    assert len(strict)==7
    assert {(z["butterfly"],z["accepted_host_id"]) for z in strict} <= {
        (z["butterfly"],z["accepted_host_id"]) for z in subset
    }
    assert sum(z["mainlands"]==["Thailand"] for z in strict)==6
    assert sum(z["mainlands"]==["Finland"] for z in strict)==1
    table=(ROOT/"results/hosts_use_vs_flora_20261116/SHARED_BUTTERFLY_HOST_PAIRS_TABLE_20261117.md")
    t=table.read_text(encoding="utf8")
    assert t.count(" | Yes |")==7 and t.count(" | No |")==4

    assert any(z["butterfly"]=="Vanessa cardui" and z["host_name"]=="Cirsium vulgare"
               and z["islands"]==["Hawaii"] for z in subset)
    assert any(z["butterfly"]=="Danaus plexippus" and
               z["host_name"]=="Calotropis gigantea" for z in subset)
    p=PAPER.read_text(encoding="utf8")
    for marker in ("67 species","57","−6.01","−18.49","+6.15",
                   "Ten of those eleven pairs","Thailand","Cirsium vulgare",
                   "seven","six of the seven",
                   "NOT","not independent","same HOSTS catalogue"):
        if marker=="NOT":continue
        assert marker in p,marker
    assert "## Structured abstract" in p
    print("HOSTS_SAME_BUTTERFLY_BOTH_REGIONS_SOURCE_FIGURES_PASS")


if __name__=="__main__":main()
