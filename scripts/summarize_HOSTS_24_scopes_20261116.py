#!/usr/bin/env python3
"""Review-ready overview from the original 24 frozen HOSTS locality analysis scopes.

Reads only the completed source-locked JSON, never re-queries or reweights records.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path

COHORTS={
 "HOSTS_Papilionoidea":"Papilionoidea (butterfly families)",
 "fixed239_butterflies":"Original fixed 239 chocho butterflies",
 "all_HOSTS_lepidoptera":"All HOSTS Lepidoptera (including moths)",
}
CROPS={
 "all_hosts":"All accepted known hosts",
 "exclude_FAO_exact":"FAO exact crop species excluded",
 "exclude_FAO_exact_plus_genus_wildcards":"FAO exact species and explicit spp. genus excluded",
 "exclude_FAO_all_crop_genera":"All genera represented in FAO crops excluded (over-conservative)",
}
GEOGRAPHY=("strict_exact_or_composite","expanded_approximate")


def pct(x):
    return "" if x is None else f"{100*float(x):.1f}%"


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input-json",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    j=json.loads(a.input_json.read_text(encoding="utf8"))
    if j["status"]!="SOURCE_ANALYSIS_NOT_FIELD_CONFIRMED_LARVAL_USE" or len(j["results"])!=24:
        raise RuntimeError("Incomplete source/24-scope model provenance")
    table=[]
    for cohort in COHORTS:
        for geo in GEOGRAPHY:
            for crop in CROPS:
                key="__".join((cohort,geo,crop))
                if key not in j["results"]:raise ValueError("Missing required source scope "+key)
                item=j["results"][key]
                for threshold in ("min10","min20"):
                    t=item[threshold]
                    med=t["tests"].get("island_vs_mainland",{})
                    ocean=t["tests"].get("oceanic_vs_continental",{})
                    lm=t["adjusted_model_island_vs_mainland"]
                    portcoef=lm.get("coef_log_odds",{}).get("log2_1plus_WPI_ports")
                    pci=lm.get("source_locality_bootstrap_95pct_descriptive",{}).get("log2_1plus_WPI_ports")
                    table.append({
                        "cohort":cohort,"geography_tier":geo,"crop_rule":crop,"min_records":threshold,
                        "eligible_original_rows":item["eligible_original_HOSTS_rows"],
                        "distinct_insects":item["distinct_insects"],
                        "distinct_host_species":item["distinct_accepted_hosts"],
                        "island_localities":t["island"]["localities"],
                        "mainland_localities":t["mainland"]["localities"],
                        "island_median":t["island"]["median_share"],
                        "mainland_median":t["mainland"]["median_share"],
                        "island_weighted":t["island"]["record_weighted_share"],
                        "mainland_weighted":t["mainland"]["record_weighted_share"],
                        "island_distinct_pair_median":t["island"]["unique_pair_median_share"],
                        "mainland_distinct_pair_median":t["mainland"]["unique_pair_median_share"],
                        "MW_p_exploratory":med.get("descriptive_Mann_Whitney_p"),
                        "median_diff_bootstrap95":med.get("locality_resample_median_difference_95pct"),
                        "oceanic_localities":t["oceanic"]["localities"],
                        "continental_localities":t["continental"]["localities"],
                        "oceanic_median":t["oceanic"]["median_share"],
                        "continental_median":t["continental"]["median_share"],
                        "oceanic_vs_continental_p":ocean.get("descriptive_Mann_Whitney_p"),
                        "port_beta_logodds":portcoef,
                        "port_bootstrap95":pci,
                    })
    a.out_dir.mkdir(parents=True,exist_ok=True)
    out=a.out_dir/"HOSTS_all_24_scopes_review_table.csv"
    with out.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)

    lines=["# Reproduced NHM HOSTS locality-use audit: all 24 source-locked scopes",
           "",
           "This analysis uses the original HOSTS Location field and WCVP v13",
           "regional introduced/native flags. **These are literature-attributed",
           "host-use reports, not uniformly validated field larval feeding.**",
           "Each location is a unit; crop and cohort variants are sensitivities",
           "rather than independent replications. Only complete WCVP region-origin",
           "matches enter denominators.",
           "",
           "| Butterfly/insect cohort | Locality mapping | Crop rule | Island median (n) | Mainland median (n) | Location-level p | Port coefficient (log odds) |",
           "|---|---|---|---:|---:|---:|---:|"]
    for x in table:
        if x["min_records"]!="min10":continue
        lines.append("| "+ " | ".join([
            COHORTS[x["cohort"]],
            x["geography_tier"],
            CROPS[x["crop_rule"]],
            f'{pct(x["island_median"])} ({x["island_localities"]})',
            f'{pct(x["mainland_median"])} ({x["mainland_localities"]})',
            f'{x["MW_p_exploratory"]:.3f}' if x["MW_p_exploratory"] is not None else "insufficient",
            f'{x["port_beta_logodds"]:+.3f}' if x["port_beta_logodds"] is not None else "unidentifiable",
        ])+" |")
    lines += ["",
              "## Explicit inference boundary",
              "",
              "- No common biological field-abundance denominator or butterfly observation-effort controls.",
              "- All HOSTS-Lepidoptera results include moths and must not be called butterfly-only.",
              "- Original fixed239 and Papilionoidea-family results have different taxon-selection and denominator sets.",
              "- WCVP introduced status can differ by WGSRPD3 unit and is not dated local species arrival.",
              "- Port coefficients reflect a partial locality-level association, not ship-mediated transport.",
              "- Original preliminary user figures are not treated as verified before this source-model result.",
              ""]
    (a.out_dir/"HOSTS_all_24_scopes_review_table.md").write_text("\n".join(lines),encoding="utf8")
    simplified=[{
      "cohort":v["cohort"],"geo":v["geography_tier"],"crop":v["crop_rule"],"nsite":[v["island_localities"],v["mainland_localities"]],
      "median_pct":[v["island_median"]*100 if v["island_median"] is not None else None,
                    v["mainland_median"]*100 if v["mainland_median"] is not None else None],
      "p":v["MW_p_exploratory"],"port_beta":v["port_beta_logodds"],
      "island_oceanic_continental":[v["oceanic_localities"],v["continental_localities"],v["oceanic_median"],v["continental_median"]]
    } for v in table if v["min_records"]=="min10"]
    print("HOSTS_24_SCOPE_REVIEW",json.dumps(simplified,ensure_ascii=False),flush=True)
    print("HOSTS_24_SCOPE_SUMMARY_VERIFIED",str(out),flush=True)


if __name__=="__main__":main()
