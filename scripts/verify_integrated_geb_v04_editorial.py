#!/usr/bin/env python3
"""Fail-closed, source-backed audit of the four-figure GEB v0.4 editorial restructuring.

This verifies exactly four main figures, original v0.3 reference/supplement
survival, source values against immutable main evidence, SVG parseability,
and that supporting but non-independent positive results remain qualified.
Only renders PNG previews on request; does not fit biological models.
"""
from __future__ import annotations

import argparse
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
    # Exact original 18-Nov references, original supplementary positive and adverse data,
    # and v0.3 provenance are preserved in full, byte-for-byte after this boundary.
    original_suffix=old.split(ORIGINAL_PREFIX,1)[1]
    new_suffix=main.split(ORIGINAL_PREFIX,1)[1]
    if original_suffix!=new_suffix:
        raise AssertionError("Original v0.3 References / Supplement S.A–S.F source text was modified or lost")
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
    if "### 4.2. Ports mark where a shared introduced flora accumulates" not in first:
        raise AssertionError("Positive 4.2 port section heading not installed")
    if "### 4.3. Introduced host geography contributes to recorded butterflies" not in first:
        raise AssertionError("Supported positive 4.3 consumer heading not installed")
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
        "### 3.4. Independent butterfly observations contained",
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
        or figurepaths[2].name!=Path(figfiles["fig3"]["source_same_image_copied_with_new_sequential_filename"]).name
        or figurepaths[3].name!=Path(figfiles["fig4"]["new_svg"]).name):
        raise AssertionError("Figure source receipt and manuscript are not the same four SVGs")
    # New Figure 3 is the exactly same source plot, not an unsupported new numerical barplot.
    original_fig3=root/figfiles["fig3"]["original"]
    if figurepaths[2].read_bytes()!=original_fig3.read_bytes():
        raise AssertionError("Source-verifiable leading host plants were altered when renamed Figure 3")
    if "### S.G.6." not in supp or supp.count("## S.G.")!=0:
        # S.G.1–6 are separate headings within a standalone appendix;
        # check exact positive and boundary-result section markers instead.
        if not all(f"## S.G.{i}." in supp for i in range(1,7)):
            raise AssertionError("New independent recorded butterfly source supplement lacks expected sections")
    for literal in [
        "14,553", "+54.9%", "+118.9%", "+40.9%", "22.7%", "24.5%",
        "+0.2854", "−0.0187", "+4.12%", "+0.0913%", "+0.215%",
        "+0.438%", "+0.695%", "+0.779%", "−4.036%", "−2.697%",
        "18 living endemic", "not a causal", "29", "source",
    ]:
        if literal not in first:
            raise AssertionError(f"Missing result or interpretation guardrail: {literal}")
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
        "protected_original_reference_and_supplement_chars":len(new_suffix),
        "original_four_main_sections":len(expected),
        "main_chars":len(first),
        "GEB_structured_abstract_words":abstract_words,
        "GEB_main_body_words":main_body_words,
        "positive_island_port_sea_figures_verified":True,
        "SVG_visual_previews":bool(args.preview_dir),
        "main_source_commit":MAIN_SOURCE,
    },ensure_ascii=False))

if __name__=="__main__":
    main()
