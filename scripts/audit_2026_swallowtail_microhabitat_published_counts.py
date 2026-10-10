#!/usr/bin/env python3
"""Source-transcribed 2026 paper arithmetic consistency audit (NO ORIGINAL RAW DATA).

This script checks mutually inconsistent quadrat denominators in a recent
published swallowtail study. It does not reconstruct the unshared original
survey rows, rerun logistic regression, or test competition causality.

Jang et al (2026), Journal of Ecology and Environment 50:15
https://doi.org/10.5141/jee.26.017
https://www.e-jecoenv.org/journal/view.html?uid=1281&vmd=Full
Transcribed from article HTML on 2026-10-10.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

PUBLICATION={
    "doi":"10.5141/jee.26.017",
    "online_publication_date":"2026-08-06",
    "date_source_checked":"2026-10-10",
    "article_url":"https://www.e-jecoenv.org/journal/view.html?uid=1281&vmd=Full",
    "native_host_species_in_south_korea":"Aristolochia contorta",
    "site_quadrat_counts_as_published":{"AY":50,"CJ":51,"JM":57,"PT":57},
    "declared_total_quadrats_in_methods":255,
    "table1_quadrat_labels":{"neither":17,"S_only":155,"A_only":72,"both":11},
    "table1_quadrat_by_site":{
        "AY":{"neither":3,"S_only":26,"A_only":20,"both":1},
        "CJ":{"neither":4,"S_only":28,"A_only":16,"both":3},
        "JM":{"neither":6,"S_only":32,"A_only":14,"both":5},
        "PT":{"neither":4,"S_only":29,"A_only":22,"both":2},
    },
    "table1_quadrat_by_month":{
        "May":{"neither":4,"S_only":12,"A_only":17,"both":3},
        "June":{"neither":13,"S_only":21,"A_only":23,"both":6},
        "July":{"neither":0,"S_only":21,"A_only":14,"both":2},
        "August":{"neither":0,"S_only":31,"A_only":11,"both":0},
        "September":{"neither":0,"S_only":30,"A_only":7,"both":0},
    },
    "table1_ramet_status":{"neither":542,"S_only":204,"A_only":124,"both":5},
    "adjusted_relative_light_intensity_or_per_1pct":{"S_montela":1.009,"A_alcinous":0.979},
    "article_original_data_availability":"corresponding author on reasonable request; no public row-level file listed",
}
CATEGORIES=("neither","S_only","A_only","both")

def sum_categories(groups):
    return {k:sum(group[k] for group in groups.values()) for k in CATEGORIES}

def audit(source=PUBLICATION):
    site_total=sum(source["site_quadrat_counts_as_published"].values())
    site_categories=sum_categories(source["table1_quadrat_by_site"])
    month_categories=sum_categories(source["table1_quadrat_by_month"])
    ramets=source["table1_ramet_status"]
    n_ramets=sum(ramets.values())
    published_header=source["table1_quadrat_labels"]
    publish_total=sum(published_header.values())
    return {
      "schema":"chocho_2026_published_swallowtail_quadrat_consistency_v01",
      "source":{k:source[k] for k in ("doi","online_publication_date","article_url","date_source_checked")},
      "source_kind":"PUBLICATION_PRINTED_TABLES_ONLY_NOT_ORIGINAL_RECORDS",
      "original_raw_field_data_accessible_in_this_audit":False,
      "published_methods_total_quadrats":source["declared_total_quadrats_in_methods"],
      "published_site_total_quadrats":site_total,
      "published_table_header_total_quadrats":publish_total,
      "site_breakdown_by_presence":site_categories,
      "month_breakdown_by_presence":month_categories,
      "printed_table_header_by_presence":published_header,
      "inferred_internal_115_not_155_S_only_count":{
        "site_breakdown_S_only":site_categories["S_only"],
        "month_breakdown_S_only":month_categories["S_only"],
        "printed_header_S_only":published_header["S_only"],
        "discrepancy":published_header["S_only"]-site_categories["S_only"],
      },
      "site_totals_agree_with_published_table_rows":all(
          sum(rows.values())==source["site_quadrat_counts_as_published"][site]
          for site,rows in source["table1_quadrat_by_site"].items()),
      "site_category_and_month_category_sums_agree":site_categories==month_categories,
      "article_quadrat_header_is_internally_consistent":(site_total==source["declared_total_quadrats_in_methods"]==publish_total==sum(site_categories.values())),
      "source_based_plausible_quadrat_denominator":sum(site_categories.values()),
      "denominator_requires_author_confirmation":True,
      "original_ramets_by_presence":ramets,
      "original_ramet_total":n_ramets,
      "ramets_with_both_species":ramets["both"],
      "both_species_ramet_fraction":ramets["both"]/n_ramets,
      "quadrat_both_count_from_all_source_breakdowns":site_categories["both"],
      "both_species_quadrat_fraction_if_internal_215_denominator":site_categories["both"]/site_total,
      "ramet_and_quadrat_statistics_not_independent_samples":True,
      "ramet_overlap_is_not_causal_competition":True,
      "geographic_transport_to_kyoto_A_debilis_not_identified":True,
      "new_ecological_effect_estimated":False,
      "GEB_PR38_untouched":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    results=audit()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(results,indent=2,ensure_ascii=False))

if __name__=="__main__":main()
