#!/usr/bin/env python3
"""Focal host species and FAO crop contributions within Papilionoidea locality rows.

Runs on source-locked derivative CSVs from the completed main HOSTS workflow.
No WCVP reinstall or new origin reclassification is performed.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
from analyze_HOSTS_reported_use_island_fao_ports_20261116 import build_rows


def report(rows):
    out={}
    for group in ("island","mainland"):
        subset=[r for r in rows if r["comparison_group"]==group]
        introduced=[r for r in subset if r["origin"]=="introduced"]
        native=[r for r in subset if r["origin"]=="native"]
        sites={r["location"] for r in subset}
        def rates(rr):
            n=len(rr)
            return {
                "n":n,
                "FAO_exact_crop_n":sum(r["crop_exact"] for r in rr),
                "FAO_exact_crop_fraction":(sum(r["crop_exact"] for r in rr)/n if n else None),
                "FAO_any_crop_genus_n":sum(r["crop_any_FAO_genus"] for r in rr),
                "FAO_any_crop_genus_fraction":(sum(r["crop_any_FAO_genus"] for r in rr)/n if n else None),
            }
        by_taxon=Counter(r["plant"] for r in introduced)
        by_pair=Counter((r["insect_species"],r["plant"]) for r in introduced)
        out[group]={
            "locality_count":len(sites),
            "known_origin_rows":len(subset),
            "introduced":rates(introduced),
            "native":rates(native),
            "introduced_top_host_plants":[{"accepted_name":x,"report_rows":n,
                        "insects_using":len({r["insect_species"] for r in introduced if r["plant"]==x})}
                    for x,n in by_taxon.most_common(30)],
            "introduced_distinct_butterfly_host_pairs":len(by_pair),
            "introduced_unique_accepted_host_plants":len(by_taxon),
        }
    return out


def main():
    p=argparse.ArgumentParser()
    for field in ("hosts","host-resolution","wcvp","locations","fao","ports",
                  "areas","butterflies239","out-dir"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args()
    rows,_,coverage,source=build_rows(a)
    output={}
    for cohort in ("HOSTS_Papilionoidea","fixed239_butterflies","all_HOSTS_lepidoptera"):
        for geo in ("strict_exact_or_composite","expanded_approximate"):
            sub=[r for r in rows if r["rearing"]!="Y" and
                 (geo=="expanded_approximate" or r["mapping_tier"] in ("exact","composite"))
                 and r["mapping_tier"]!="provisional" and
                 (cohort=="all_HOSTS_lepidoptera" or
                  (cohort=="fixed239_butterflies" and r["focal_butterfly239"]) or
                  (cohort=="HOSTS_Papilionoidea" and r["papilionoidea"]))]
            output[cohort+"__"+geo]=report(sub)
    payload={
        "schema":"chocho_HOSTS_crop_dependency_by_report_region_v01",
        "status":"LITERATURE_HOST_ASSOCIATIONS_NOT_DIRECT_FIELD_FEEDING",
        "source":source,
        "raw_origin_matched_rows":len(rows),
        "cohorts":output,
        "limits":[
            "A raw HOSTS literature association row is not a replicate butterfly or verified field larval consumption.",
            "FAO species and broad genera are independent taxonomic filters, not proof of agricultural use.",
            "The ratio of crop taxa among WCVP origin-known records excludes unknown and cannot represent all island host records.",
            "The source paper or experiment underlying a report is not resolved from these columns.",
            "Same HOSTS reference may contribute multiple rows; site record-based crop percentages reflect publication bias.",
        ],
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"HOSTS_crops_contribution_by_origin.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("HOSTS_CROP_DEPENDENCE_AUDIT",json.dumps({
      "primary":output["HOSTS_Papilionoidea__strict_exact_or_composite"],
      "fixed239":output["fixed239_butterflies__strict_exact_or_composite"],
      "all_lepidoptera":output["all_HOSTS_lepidoptera__strict_exact_or_composite"],
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
