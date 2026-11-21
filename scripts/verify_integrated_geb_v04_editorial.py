#!/usr/bin/env python3
"""Fail-closed, source-backed audit of the four-figure GEB v0.4 editorial restructuring.

This verifies exactly four main figures, original v0.3 reference/supplement
survival, source values against immutable main evidence, SVG parseability,
and that supporting but non-independent positive results remain qualified.
Only renders PNG previews on request; does not fit biological models.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

REL = Path("manuscript/integrated_butterfly_resource_access_20261109")
MAIN = REL / "ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
CAPTIONS = REL / "FIGURE_CAPTIONS_V03.md"
SUP_G = REL / "SUPPLEMENTARY_INDEPENDENT_BUTTERFLY_RECORDS_20261119.md"
NUMBERS = Path("figures/integrated_island_port_butterfly_v03/FIGURE_V04_SOURCE_NUMBERS.json")
ORIGINAL_PREFIX = "## References (DOI identities checked; final source and style audit outstanding)"
LEGACY_ARCHIVE = REL / "SUPPLEMENTARY_LEGACY_REFERENCES_AND_APPENDICES_V03_UNEDITED.md"
MAIN_SOURCE = "3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"

def near(a,b,limit=1e-6):
    if not math.isfinite(float(a)) or abs(float(a)-float(b))>limit:
        raise AssertionError(f"Source figure value mismatch {a} versus {b}")

def load_j(p):
    return json.loads(p.read_text(encoding="utf8"))

def assert_source(receipt,source):
    """Source main is the immutable original numeric PR #180–#187 evidence snapshot."""
    if receipt["main_evidence_commit_at_read"] != MAIN_SOURCE:
        raise AssertionError("New figure numerical receipt pinned to other source")
    rows=[load_j(source / f) for f in receipt["source_paths"]]
    port,flora,forecast_time,forecast_space,forecast_distant,gbif,world,sa=rows
    fig2=receipt["figures"]["fig2"]
    if fig2["island_botanical_regions"]!=126 or fig2["whole_flora_jaccard"]["unorded_pairs"]!=7875:
        raise AssertionError("Mixed island geographic source denominator")
    selection=[r for r in port["results"]
               if r["cohort"]=="WCVP126" and not r["additional_flora_covariates"]]
    by_response={r["response"]:r for r in selection}
    for given,source_key in [("native","native"),("introduced","introduced"),
                              ("introduced_minus_native","alien_minus_native_log_ratio")]:
        q=fig2["port_coefficient"][given]
        orig=by_response[source_key]
        near(q["estimate"],orig["port_beta_per_log2_1plus_WPI_port"])
        for a,b in zip(q["descriptive_HC3_95"],orig["HC3_95pct_descriptive_interval"]):
            near(a,b)
    near(fig2["whole_flora_jaccard"]["native_median"],flora["median_native_Jaccard"])
    near(fig2["whole_flora_jaccard"]["introduced_median"],flora["median_alien_Jaccard"])
    if flora["total_unordered_strict_island_pairs"]!=7875:
        raise AssertionError("Whole flora dyad sample changed")
    f4=receipt["figures"]["fig4"]
    f155=f4["panel_A_155_independently_indexed_butterflies"]
    if f155["test_cells"]!=208482:
        raise AssertionError("Original 155 butterfly independent cohort drift")
    for name,orig in [("time_held",forecast_time),("spatial_held",forecast_space),
                     ("over_100km",forecast_distant)]:
        rows=orig["results"]
        recalculated=100*(rows["M3_both"]["species_equal_mse"]-
                          rows["M5_host_composition"]["species_equal_mse"])/rows["M3_both"]["species_equal_mse"]
        near(recalculated,f155["mse_reduction_percent"][name])
    f224=f4["panel_B_224_butterfly_source"]
    if f224["regions"]!=169 or f224["pairs"]!=14196:
        raise AssertionError("Original 224 GBIF observed beta source cohort changed")
    near(f224["original_H_after_O_MSE_improvement_pct"],gbif["provenance_assertions"]["original_H_after_O_percent"])
    for field,key in [("dataset_days_MSE_improvement_pct_median","primary_dataset_key_calendar_day"),
                      ("record_count_MSE_improvement_pct_median","secondary_indexed_candidate_record_count")]:
        near(f224[field],gbif[key]["H_after_O_MSE_improvement_median_percent"])
    fw=f4["panel_C_eight_unseen_worlds"]
    near(fw["original_H_after_O_world_equal_pct"],world["raw_original_world_equal_MSE_improvement_percent"])
    near(fw["dataset_day_median"],world["per_proxy"]["dataset_key_calendar_day"]["region_level1_equal_MSE_gain_percent_median"])
    near(fw["record_count_median"],world["per_proxy"]["indexed_candidate_record_count"]["region_level1_equal_MSE_gain_percent_median"])
    saout=f4["panel_D_south_america"]
    if saout["within_pairs"]!=190 or saout["cross_world_pairs"]!=2980:
        raise AssertionError("South American within and cross-world support shifted")
    day=sa["registered_H_after_O_world_transfer_under_original_thinning"]["datasetKey_calendar_day"]
    record=sa["registered_H_after_O_world_transfer_under_original_thinning"]["indexed_record_count"]
    for value_name,orig in [
        ("dataset_day_within_median",day["South_America_within_gain_median_pct"]),
        ("dataset_day_cross_median",day["South_America_to_other_gain_median_pct"]),
        ("record_count_within_median",record["South_America_within_gain_median_pct"]),
        ("record_count_cross_median",record["South_America_to_other_gain_median_pct"]),
    ]:
        near(saout[value_name],orig)

def average_ranks(values):
    order=sorted(range(len(values)),key=lambda i:(values[i],i))
    result=[0.0]*len(values)
    at=0
    while at<len(order):
        end=at+1
        while end<len(order) and values[order[end]]==values[order[at]]:
            end+=1
        rank=(at+1+end)/2
        for i in range(at,end):
            result[order[i]]=rank
        at=end
    return result

def pearson(a,b):
    ma=sum(a)/len(a)
    mb=sum(b)/len(b)
    num=sum((x-ma)*(y-mb) for x,y in zip(a,b))
    d1=sum((x-ma)**2 for x in a)
    d2=sum((y-mb)**2 for y in b)
    return num/math.sqrt(d1*d2)

def validate_enhanced_figures(receipt,source,repo_root,figurepaths):
    fig2=receipt["figures"]["fig2"]
    sea=fig2["marine_circuit"]
    source_strict=load_j(source/sea["source_report"])
    if sea["source_n_botanical_islands"]!=104 or sea["original_source_pair_count"]!=5355:
        raise AssertionError("Original marine cohort drift")
    orig=sea["source_original_one_island_heldout"]
    near(orig["introduced_MSE_gain_percent"],4.239)
    near(orig["native_MSE_gain_percent"],-.110)
    source_readme=(source/"results/WCVP126_circuit_strict_geographic_transfer_20261118/README.md").read_text(encoding="utf8")
    if "+4.239%" not in source_readme or "−0.110%" not in source_readme:
        raise AssertionError("Archived original marine model sign/magnitude not verified")
    strict=sea["strict_both_endpoints_unseen"]
    if strict["test_pairs"]!=1089 or strict["fixed_geography_folds"]!=5:
        raise AssertionError("Strict source geographic transfer support drift")
    for label,k in [("native","native"),("introduced","introduced")]:
        row=source_strict["uniform_circuit"][k]
        near(strict[label+"_MSE_gain_percent"],row["gain_pct_equal_fold"])
        for a,b in zip(strict[label+"_descriptive_interval"],row["descriptive_5fold_paired_bootstrap95"]):
            near(a,b)
        if strict[label+"_positive_folds"]!=row["circuit_folds_improved"]:
            raise AssertionError("Strict marine CV fold count drift")
    svg2=figurepaths[1].read_text(encoding="utf8")
    for token in ("+4.239%","-0.110%","+2.735%","-0.044%","Both islands unseen","5,355"):
        if token not in svg2:
            raise AssertionError("Missing original/strict positive or adverse sea plot: "+token)

    # Fig2D: source-backed newly recomputed original fixed239 isolation model.
    new239=load_j(repo_root/"results/isolation_port_resource_reanalysis_20261120/SOURCE239_PAPER_RECEIPT_V01.json")
    ip=fig2.get("isolation_panel") or {}
    if ip.get("source_receipt")!="results/isolation_port_resource_reanalysis_20261120/SOURCE239_PAPER_RECEIPT_V01.json" or new239["cohort"]["islands"]!=126:
        raise AssertionError("Fig2D current source/126 islands missing")
    if new239["cohort"]["native_island_coverage_units"]!=4829 or new239["cohort"]["contemporary_island_coverage_units"]!=10570:
        raise AssertionError("Fixed239 island resource denominators have drifted")
    for field,label in [("native_host_distance_standardized_beta","native_host239"),
                        ("contemporary_host_distance_standardized_beta","contemporary_host239")]:
        near(ip[field],new239["models"][label]["distance_beta"])
    near(ip["paired_distance_slope_shift_standardized_beta"],new239["models"]["paired_native_current_shift"]["distance_beta"])
    near(ip["paired_shift_nine_level1_t8_p"],new239["models"]["paired_native_current_shift"]["distance_nine_cluster_p"])
    near(ip["introduced_only_share_distance_logodds_beta"],new239["models"]["introduced_only_share_logodds"]["distance_beta"])
    near(ip["introduced_only_share_nine_level1_t8_p"],new239["models"]["introduced_only_share_logodds"]["distance_nine_cluster_p"])
    spatial=new239.get("spatial_sensitivity") or {}
    if spatial.get("source_workflow",{}).get("conclusion")!="success":
        raise AssertionError("Fig2D extra spatial sensitivity lacks successful source run")
    near(spatial["spatial_mem10"]["native_host"]["distance_beta"],-0.1330290218602384)
    near(spatial["spatial_mem10"]["introduced_only_share_logodds"]["distance_beta"],0.3948940113274428)
    if not (repo_root/"scripts/audit_island_port_spatial_239_20261120.py").is_file():
        raise AssertionError("Spatial source validator script missing from journal branch")
    for txt in ("Mainland isolation and candidate butterfly food","Original fixed 239-host panel",
                "55%","66%","72%","Nine-Level1 p = 0.014","p = 0.068","Both islands unseen"):
        if txt not in svg2: raise AssertionError("Fig2D source or uncertainty label lost: "+txt)

    fig3=receipt["figures"]["fig3"]
    f3=fig3["evidence_sources"]
    concentration=load_j(source/f3["plant_concentration"])
    p=concentration["plant_species"]
    plot=fig3["panelA_plants"]
    if (plot["total_introduced_host_contributors"]!=p["contributors"]
        or plot["contributors_for_half"]!=p["contributors_for_50pct"]
        or plot["total_introduced_host_contributors"]!=670
        or plot["contributors_for_half"]!=38):
        raise AssertionError("Plant source concentration audit drift")
    near(plot["concentration_Gini"],p["gini"])
    for k,v in p["top_k_share"].items():
        near(plot["top_k_cumulative_share"][k],v)
    if fig3["source_original_top10_port_plot_retained"] not in [
        "figures/integrated_island_port_butterfly_v03/Fig2_global_contributors_at_ports.svg"
    ] or not (repo_root/fig3["source_original_top10_port_plot_retained"]).is_file():
        raise AssertionError("Original supplemental top-ten port source graphic not retained")
    metrics_path=source/f3["all_butterfly_metrics"]
    region_path=source/f3["island_butterfly_gains"]
    # Both panels derive from original archived main, never from a rewritten species table.
    with metrics_path.open(encoding="utf8",newline="") as fh:
        metrics=list(csv.DictReader(fh))
    with region_path.open(encoding="utf8",newline="") as fh:
        regions=list(csv.DictReader(fh))
    ids={}
    for row in regions:
        if row["group"]=="island":
            if row["species"] in ids:
                raise AssertionError("Duplicate original island species")
            ids[row["species"]]=row
    if len(metrics)!=239 or len(ids)!=239:
        raise AssertionError("Original butterfly sample is not 239 species")
    points=[]
    for row in metrics:
        name=row["species"]
        if name not in ids:
            raise AssertionError("Original 239 butterfly species do not align across source CSVs")
        points.append({
            "species":name,
            "native":float(row["native_resource_units"]),
            "added":float(ids[name]["added"]),
            "family":float(row["host_family_count"]),
            "global_added":float(row["introduced_added_units"])
        })
    if sum(p["added"] for p in points)!=5741:
        raise AssertionError("Original island resource increments changed")
    q=fig3["panelB_butterflies"]
    rhos=pearson(average_ranks([p["native"] for p in points]),
                 average_ranks([p["added"] for p in points]))
    near(q["source_native_global_resources_vs_added_island_opportunities_Spearman_rho"],rhos,1e-10)
    ordered=sorted(points,key=lambda r:r["native"])
    top=ordered[-59:]
    bottom=ordered[:60]
    if (q["n_species"]!=239 or q["total_added_island_units"]!=5741
        or q["native_global_broadest_n"]!=59
        or q["native_global_narrowest_n"]!=60):
        raise AssertionError("Original resource quartile cohort labels drift")
    near(q["native_global_broadest_added"],sum(z["added"] for z in top))
    near(q["native_global_narrowest_added"],sum(z["added"] for z in bottom))
    near(q["native_global_broadest_share"],sum(z["added"] for z in top)/5741)
    near(q["native_global_narrowest_share"],sum(z["added"] for z in bottom)/5741)
    f=fig3["panelC_butterflies"]
    lone=[r for r in points if r["family"]==1]
    if (len(lone)!=f["one_documented_host_family_n"]
        or sum(r["global_added"]>0 for r in lone)!=f["one_documented_host_family_gaining_global_n"]
        or len(points)!=f["all_source_butterflies_n"]
        or sum(r["global_added"]>0 for r in points)!=f["all_source_butterflies_gaining_global_n"]
        or sum(r["added"]>0 for r in points)!=f["gaining_at_least_one_island_n"]):
        raise AssertionError("One-family consumer or overall butterfly gains drift")
    svg3=figurepaths[2].read_text(encoding="utf8")
    for token in ("670 accepted","38 supply half","Gini = 0.777","239 butterflies",
                  "Spearman","39.6%","11.2%","82/97","+46–47%"):
        if token not in svg3:
            raise AssertionError("New source numerical butterfly/plant graphic missing: "+token)
    if svg3.count("<circle ")<239:
        raise AssertionError("Expected actual 239 species plotted, not a summary-only cartoon")

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--repo-root",type=Path,default=Path("."))
    arg.add_argument("--source-root",type=Path,required=True)
    arg.add_argument("--frozen-v03",type=Path,required=True)
    arg.add_argument("--preview-dir",type=Path,default=None)
    args=arg.parse_args()
    root=args.repo_root.resolve()
    main=(root/MAIN).read_text(encoding="utf8")
    captions=(root/CAPTIONS).read_text(encoding="utf8")
    supp=(root/SUP_G).read_text(encoding="utf8")
    receipt=load_j(root/NUMBERS)
    old=args.frozen_v03.read_text(encoding="utf8")
    if main.count(ORIGINAL_PREFIX)!=1 or old.count(ORIGINAL_PREFIX)!=1:
        raise AssertionError("Original reference/supplement protected marker moved")
    # Full unedited historical bibliography, Supplement S.A-S.F and source details
    # remain at exactly their original bytes in a separate source archive.
    original_suffix=ORIGINAL_PREFIX+old.split(ORIGINAL_PREFIX,1)[1]
    archived=(root/LEGACY_ARCHIVE).read_text(encoding="utf8")
    if original_suffix!=archived:
        raise AssertionError("Unedited v0.3 full source archive is not byte-exact")
    original_references=original_suffix.split("## Supplementary results:",1)[0]
    current_references=main.split(ORIGINAL_PREFIX,1)[1].split("## Supplementary results:",1)[0]
    original_reference_lines=[x for x in original_references.splitlines()[1:] if x.startswith("- ")]
    current_reference_lines=[x for x in current_references.splitlines() if x.startswith("- ")]
    verified_literature_dois={
        "10.1111/geb.70292",  # Couto et al. 2026: global non-native Lepidoptera establishment
        "10.1111/geb.70195",  # Hatfield et al. 2026: actual species-level beta turnover
        "10.1073/pnas.0802290105",  # Sax & Gaines 2008
        "10.1016/s0169-5347(99)01679-1",  # McKinney & Lockwood 1999
        "10.1111/j.1466-822x.2006.00212.x",  # Hobbs et al. 2006
    }
    # MacArthur & Wilson book has no DOI; exact publication/author identity required.
    def source_reference_identity(textline):
        if ("National Geospatial-Intelligence Agency (2019)." in textline and
                "*World Port Index (Publication 150)*" in textline):
            return "wpi-2019-publication150"
        m=re.search(r"https://doi\.org/([^\s]+)",textline)
        return m.group(1).lower().rstrip(".,;") if m else textline.strip().lower()
    def unique_by_identity(lines):
        output={}
        for item in lines:
            key=source_reference_identity(item)
            if key not in output:
                output[key]=item
        return output
    original_by_id=unique_by_identity(original_reference_lines)
    current_by_id=unique_by_identity(current_reference_lines)
    # The archived pre-review bibliography must remain byte exact. This one
    # independently source-verified checklist correction is the ONLY allowed
    # in-place copyedit to an original reference; it is not source deletion.
    checklist_doi = "10.15468/ye7whj"
    permitted_checklist_source = "- INBO. *European butterfly checklist and distributions*, GBIF dataset. https://doi.org/10.15468/ye7whj. **Dataset-version citation to finalize before submission.**"
    permitted_checklist_current = "- Maes, D., Wiemers, M., Verovnik, R., Warren, M., Brosens, D. & Desmet, P. (2020). *National checklists and red lists for European butterflies* (Version 1.1) [checklist data set]. Research Institute for Nature and Forest (INBO). https://doi.org/10.15468/ye7whj"
    # DOIs and all title/journal/source text are protected. Only the lead-author-
    # preserving author lists of these 14 Crossref-verified abbreviations may
    # change, with the year and everything after it byte-identical. See the
    # source audit branch's 21-Nov 14/14 DOI/title/year/author metadata receipt.
    approved_author_list_dois = {
        "10.1111/ecog.06227": "Bach",
        "10.1038/s41467-021-27186-8": "Daru",
        "10.1093/aobpla/plz051": "Essl",
        "10.1038/s41597-021-00997-6": "Govaerts",
        "10.1111/een.12845": "Halsch",
        "10.1111/j.1466-822x.2006.00212.x": "Hobbs",
        "10.1038/s41559-023-02041-9": "Kawahara",
        "10.1111/jbi.14258": "Neu",
        "10.1111/ddi.70245": "Rashid",
        "10.1002/ajb2.16224": "Schmidt",
        "10.1038/s41597-022-01473-5": "Shirey",
        "10.1038/s41467-024-49079-2": "Suchan",
        "10.1007/s11252-017-0667-0": "Tordoni",
        "10.1038/nature14910": "van Kleunen",
    }
    for key, original_entry in original_by_id.items():
        if key not in current_by_id:
            raise AssertionError("Original unique reference lost: " + key)
        revised = current_by_id[key]
        if key == checklist_doi:
            if original_entry != permitted_checklist_source or revised != permitted_checklist_current:
                raise AssertionError("Official INBO checklist correction does not match the verified six-author version-1.1 citation")
        elif key == "wpi-2019-publication150":
            source_wpi = "- National Geospatial-Intelligence Agency (2019). *World Port Index (Publication 150)*, 27th edition. https://msi.nga.mil/Publications/WPI. **Note:** the exact archived source snapshot of the study’s port file must still be matched to the original version and date."
            amended_wpi = "- National Geospatial-Intelligence Agency (2019). *World Port Index (Publication 150)* [global port-location dataset]. https://msi.nga.mil/Publications/WPI"
            if original_entry != source_wpi or revised != amended_wpi:
                raise AssertionError("WPI source edition not correctly qualified against the 2019 public mirror")
        elif key in approved_author_list_dois:
            old_tail = re.search(r"\(\d{4}\)\. .*", original_entry)
            new_tail = re.search(r"\(\d{4}\)\. .*", revised)
            if not old_tail or not new_tail or old_tail.group() != new_tail.group():
                raise AssertionError("Title, year or publisher text changed in Crossref author-list-only edit: " + key)
            if ("et al." not in original_entry or "et al." in revised or
                not revised.startswith("- " + approved_author_list_dois[key] + ",")):
                raise AssertionError("Crossref-author-list correction/lead author not applied: " + key)
        elif revised != original_entry:
            raise AssertionError("Original unique reference altered without source-audited copyedit: " + key)
    added_references=[key for key in current_by_id if key not in original_by_id]
    if len(added_references)!=6 or not verified_literature_dois.issubset(set(added_references)):
        raise AssertionError("GEB exact published comparison and original theory references not all in sorted bibliography")
    if not any("MacArthur, R. H. & Wilson, E. O. (1967)" in current_by_id[k] for k in added_references):
        raise AssertionError("Original island biogeography book missing")
    if len(current_by_id)!=len(current_reference_lines):
        raise AssertionError("Duplicate citations remain in cleaned GEB bibliography")
    if len(current_reference_lines)!=39:
        raise AssertionError("Expected 39 distinct original+verified new GEB source references")
    if current_reference_lines != sorted(current_reference_lines,
        key=lambda row:re.sub(r"[^A-Za-z0-9]","",row[2:]).lower()):
        raise AssertionError("References are not in consistent first-author alphabetic order")
    old_sup="## Supplementary results:"+old.split("## Supplementary results:",1)[1]
    actual_sup="## Supplementary results:"+main.split("## Supplementary results:",1)[1]
    replacements=(
        ("### 3.2. Globally influential host plants in island port vegetation",
         "**Archived v0.3 result: globally influential host plants in island port vegetation.**"),
        ("### 3.3. Source-independent island-endemic larvae reveal local use but also incomplete coverage",
         "**Archived v0.3 result: source-independent island-endemic larvae and incomplete taxonomic coverage.**"),
        ("### 3.4. Independent butterfly-record prediction supports an information role for host geography",
         "**Archived v0.3 result: independent butterfly-record prediction and host-geography information.**"),
        ("### 3.6. Historical butterfly records constrain a port–host establishment interpretation",
         "**Archived v0.3 result: historical butterfly records constrain a port–host establishment interpretation.**"),
        ("#### port whole flora results","#### S.F.1. Port whole flora results"),
        ("#### endemic source selection","#### S.F.2. Endemic source selection"),
        ("#### predictive source results","#### S.F.3. Independent butterfly recording source results"),
        ("#### historical record safeguards","#### S.F.4. Historical record safeguards"),
        ("#### data availability concise","#### S.F.5. Historical data and code availability (source archive)"),
        ("## Data and code availability (draft; final rights and archive check required)",
         "**Archived v0.3 data-availability detail (supplementary source record).**"),
        ("### Source archive and reproducibility pointers (NOT independent extra studies)",
         "#### S.F.6. Source archive and reproducibility pointers (NOT independent extra studies)"),
        ("## v0.3 provenance and scientific claims","#### S.F.7. Original v0.3 provenance and scientific claims"),
        ("**conclusions.** ## 6. Conclusions","**Archived v0.3 conclusions.**")
    )
    expected_sup=old_sup
    for old_header,new_header in replacements:
        if expected_sup.count(old_header)!=1:
            raise AssertionError("Frozen old Supplement heading changed: "+old_header)
        expected_sup=expected_sup.replace(old_header,new_header)
    sf="### S.F. Full island–port–marine estimates, alternate comparisons and prior main-text detail"
    archive_note="**Provenance:** The entire earlier reference list and supplementary source, before editorial heading renumbering, is preserved verbatim in the [original unedited v0.3 References and Supplement S.A–S.F](SUPPLEMENTARY_LEGACY_REFERENCES_AND_APPENDICES_V03_UNEDITED.md). Only heading structure was regularized; source-positive and adverse numerical prose remains intact."
    if expected_sup.count(sf)!=1:
        raise AssertionError("Original S.F marker changed")
    expected_sup=expected_sup.replace(sf,sf+"\n\n"+archive_note)
    if expected_sup!=actual_sup:
        raise AssertionError("Historical Supplement content changed beyond approved headings and one provenance note")
    if main.count("## Data and code availability (draft; final rights and archive check required)")!=1:
        raise AssertionError("Main manuscript contains duplicate Data availability headings")
    if re.search(r"^### 3\.[2-6]\.",actual_sup,re.M):
        raise AssertionError("Historical Results section numbering still masquerades as a second main Results")
    new_suffix=main.split(ORIGINAL_PREFIX,1)[1]
    first=main.split(ORIGINAL_PREFIX,1)[0]
    abstract_text=first.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    abstract_plain=re.sub(r"[*#]+","",abstract_text)
    abstract_words=len(abstract_plain.strip().split())
    if abstract_words>300:
        raise AssertionError(f"GEB Research Article structured abstract exceeds 300 words: {abstract_words}")
    required_abstract_headers=("Aim","Location","Time period","Major taxa studied",
                               "Methods","Results","Main conclusions")
    for title in required_abstract_headers:
        if f"**{title}.**" not in abstract_text:
            raise AssertionError("GEB abstract missing official required section: "+title)
    for key in ("82/97","39.6%","11.2%","10.3-fold","4.24%",
                "−0.11%","+2.95","+54.9%","+118.9%","+40.9%"):
        if key not in abstract_text:
            raise AssertionError("Restored positive source numeric missing from abstract: "+key)
    body_text=first.split("## 1. Introduction",1)[1].split("## Data and code availability",1)[0]
    stripped_main=re.sub(r"!\[[^]]*\]\([^)]*\)","",body_text)
    stripped_main=re.sub(r"[#*\x60]+","",stripped_main)
    main_body_words=len(stripped_main.strip().split())
    if main_body_words>5000:
        raise AssertionError(f"GEB regular 5000-word research-body target exceeded: {main_body_words}")
    five_headings=(
        "### 4.1. Islands show larger relative resource gains",
        "### 4.2. Ports and the sea reveal two spatial layers",
        "### 4.3. Plants and butterflies receive unequal shares",
        "### 4.4. Independent butterfly records reveal a biological filter",
        "### 4.5. Conservation: manage introduced host plants with consumer interactions"
    )
    for title in five_headings:
        if first.count(title)!=1:
            raise AssertionError("Five ecologically positive Discussion sections missing or repeated: "+title)
    if len(re.findall(r"^### 4.[0-9]+.",first,re.M))!=5:
        raise AssertionError("Discussion does not contain the intended five implications")
    if not all(w in first for w in ("54.3%","38 (5.7%)","Gini coefficient was 0.777",
                                   "39.6%","11.2%","McKinney & Lockwood, 1999",
                                   "Hobbs et al., 2006","Sax & Gaines, 2008")):
        raise AssertionError("New ecological discussion and cited island/nonrandom resource literature incomplete")
    if "22 islands" not in first or "seawater-dispersal morphology" not in first:
        raise AssertionError("Original EuDiS 22-island positive trait estimate missing")
    if "54.3%" not in first:
        raise AssertionError("Source island contemporary resource share missing")
    if not main.startswith("# Ports, the sea and two layers of island floras reshape butterfly larval-host opportunities"):
        raise AssertionError("GEB island-port-sea title not installed")
    expected=[
        "### 3.1. Island resource gains were amplified",
        "### 3.2. Ports and introduced floras formed",
        "### 3.3. A small set of redistributed plant identities",
        "### 3.4. Host geography improves prediction of independently recorded butterflies",
    ]
    for line in expected:
        if first.count(line)!=1:
            raise AssertionError(f"Missing or duplicate main Results: {line}")
    if len(re.findall(r"^### 3\.[0-9]+\.",first,re.M))!=4:
        raise AssertionError("Main Results is not exactly the planned four sections")
    captions_h=re.findall(r"^## Figure ([1-4]) \|",captions,re.M)
    if captions_h!=["1","2","3","4"]:
        raise AssertionError(f"Figure captions order is not exactly 1–4: {captions_h}")
    embedded=re.findall(r"!\[Figure ([1-4])\.[^\]]*\]\(([^)]+)\)",first)
    if len(embedded)!=4 or [x[0] for x in embedded]!=["1","2","3","4"]:
        raise AssertionError(f"Embedded *main* figures should be 1–4 exactly, got {embedded}")
    figurepaths=[]
    for number,link in embedded:
        resolved=(root/REL/link).resolve()
        if not resolved.is_relative_to(root) or not resolved.exists() or resolved.suffix.lower()!=".svg":
            raise AssertionError(f"Missing/non-source SVG for main figure {number}: {resolved}")
        figurepaths.append(resolved)
        parsed=ET.parse(resolved)
        if not parsed.getroot().tag.endswith("svg"):
            raise AssertionError(f"Malformed SVG root: {resolved}")
        for child in parsed.iter():
            if child.tag.endswith("script") or child.tag.endswith("foreignObject"):
                raise AssertionError("Main manuscript SVG has unexpected active/foreign content")
            for key,value in child.attrib.items():
                if (key.endswith("href") or key.endswith("src")) and value.startswith(("http://","https://")):
                    raise AssertionError("Figure embeds external data instead of frozen original source")
    figfiles=receipt["figures"]
    if (figurepaths[0].name!=Path(figfiles["fig1"]["source_unchanged"]).name
        or figurepaths[1].name!=Path(figfiles["fig2"]["new_svg"]).name
        or figurepaths[2].name!=Path(figfiles["fig3"]["new_svg"]).name
        or figurepaths[3].name!=Path(figfiles["fig4"]["new_svg"]).name):
        raise AssertionError("Figure source receipt and manuscript are not the same four SVGs")
    # Fig. 3 is now rebuilt from original source 670 plants and 239 butterfly CSVs,
    # not copied from the original top-ten host species/port observational panel.
    # The old original plot remains present and linked in the supplementary record.
    validate_enhanced_figures(receipt,args.source_root.resolve(),root,figurepaths)
    original_full=root/receipt["figures"]["fig4"]["supplementary_original_full_190_2980_South_America_diagram"]
    if not original_full.is_file():
        raise AssertionError("Full South America four-subgroup source graphic lost")
    ET.parse(original_full)
    main_f4=figurepaths[3].read_text(encoding="utf8")
    for required in ("Cross-world pairs (n = 2,980)","-4.036%","-2.697%",
                     "0%","Within-SA: −0.496% / −0.126%"):
        if required not in main_f4:
            raise AssertionError("Compact South America figure lost its signed source boundary: "+required)
    tree=ET.parse(figurepaths[3])
    for elem in tree.iter():
        for key in ("x1","x2","cx"):
            if key in elem.attrib and not (0<=float(elem.attrib[key])<=1240):
                raise AssertionError("Figure 4 element outside source SVG canvas: "+str((elem.tag,key,elem.attrib[key])))
    if main_f4.count("Cross-world pairs (n = 2,980)")!=1:
        raise AssertionError("Source regional qualification unexpectedly duplicated")

    if "### S.G.6." not in supp or supp.count("## S.G.")!=0:
        # S.G.1–6 are separate headings within a standalone appendix;
        # check exact positive and boundary-result section markers instead.
        if not all(f"## S.G.{i}." in supp for i in range(1,7)):
            raise AssertionError("New independent recorded butterfly source supplement lacks expected sections")
    for literal in [
        "14,553", "+54.9%", "+118.9%", "+40.9%", "22.7%", "24.5%",
        "+0.2854", "−0.0187", "+0.0913%", "+0.215%",
        "+0.438%", "+0.695%", "+0.779%", "−4.036%", "−2.697%",
        "18 living endemic", "not a causal", "29", "source",
    ]:
        if literal not in first:
            raise AssertionError(f"Missing result or interpretation guardrail: {literal}")
    if not re.search(r"(?<!\d)\+?4\.12%(?!\d)",first):
        raise AssertionError("Original 155-butterfly host forecast +4.12% improvement is absent")
    for forbidden in ["../../results/pr179_source_evidence_20261119/README.md",
                      "../../results/gbif224_thinned_world_transfer_20261119/README_JA.md"]:
        if forbidden in first:
            raise AssertionError("Draft links to newer main-only evidence use unavailable relative paths")
    if "0/18" not in supp and "none of those 18" not in first:
        raise AssertionError("Independent endemic exclusion was lost")
    assert_source(receipt,args.source_root.resolve())
    if args.preview_dir:
        from cairosvg import svg2png
        args.preview_dir.mkdir(parents=True,exist_ok=True)
        for i,p in enumerate(figurepaths,1):
            svg2png(url=str(p),write_to=str(args.preview_dir/f"geb_v04_Fig{i}.png"),
                    output_width=1500,background_color="#ffffff")
    print("GEB_V04_SOURCE_AND_FOUR_FIGURES_PASS",json.dumps({
        "main_figures":[str(x.relative_to(root)) for x in figurepaths],
        "new_supplement":"Supplement S.G independent GBIF source",
        "archived_original_reference_and_supplement_chars":len(original_suffix),
        "approved_heading_cleanup_only":True,
        "new_verified_4_ecological_literature_references":4,
        "compact_Fig4D_with_frozen_full_supplement":True,
        "original_four_main_sections":len(expected),
        "main_chars":len(first),
        "GEB_structured_abstract_words":abstract_words,
        "GEB_main_body_words":main_body_words,
        "positive_island_port_sea_figures_verified":True,
        "Fig2C_exact_source_original_and_strict_marine_verified":True,
        "Fig3_670_plant_and_239_butterfly_points_verified":True,
        "SVG_visual_previews":bool(args.preview_dir),
        "main_source_commit":MAIN_SOURCE,
    },ensure_ascii=False))

if __name__=="__main__":
    main()
