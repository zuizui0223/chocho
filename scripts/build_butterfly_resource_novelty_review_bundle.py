#!/usr/bin/env python3
"""Make a deterministic source-audited ALTERNATIVE butterfly GEB review bundle.

This is not a publication deposit, not the official GEB submission bundle,
and does not modify PR #38 or any original manuscript and figures.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import zipfile
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parents[1]
TITLE = "When is an introduced host a new resource?"
EXPECTED_IMAGES = {
    "candidate_Fig1_accounting_and_structural_null": "figure_resource_novelty",
    "candidate_Fig2_species_region_source_turnover": "figure_added_resource_identity_turnover",
    "candidate_Fig3_diet_breadth_and_occurrence": "figure3_ecological_context",
    "supporting_FigS1_host_inventory_composition": "figure_bce_inventory_sensitivity",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(path:Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def review_conditions(text:str, receipt_dir:Path, sources:dict) -> dict:
    if not text.startswith("# "+TITLE):
        raise ValueError("not the approved alternative manuscript draft")
    if "PR #38" not in text or "not" not in text.lower():
        raise ValueError("original submission firewall absent")
    if any(text.count(h)!=1 for h in (
        "## Structured abstract","## 1. Introduction",
        "## 2. Methods","## 3. Results",
        "## 4. Discussion","## 5. Limitations",
        "## 6. Conclusions","## References",
        "## Data and Code Availability"
    )):
        raise ValueError("manuscript structure drift")
    abstract=text.split("## Structured abstract",1)[1].split("**Keywords:**",1)[0]
    main=text.split("## 1. Introduction",1)[1].split("## References",1)[0]
    refs=text.split("## References",1)[1].split("## Data and Code Availability",1)[0]
    wc=lambda x:len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",x))
    if wc(abstract)>300 or wc(main)>5000:
        raise ValueError("GEB Research Article length limit exceeded")
    if len([r for r in refs.splitlines() if r.startswith("- ")])>50:
        raise ValueError("too many references for GEB review package")

    inv=read_json(receipt_dir/"butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json")
    turnover=read_json(receipt_dir/"butterfly_bce_added_region_identity_turnover_v0.1.json")
    local=read_json(receipt_dir/"butterfly_spatial_null_longchain_qualified_v0.2.json")
    if local.get("status")!="QUALIFIED_EXPLORATORY_MCMC_SENSITIVITY":
        raise ValueError("unqualified source structural null")
    if any(not x["pass"] for x in local["qualification"]["checks"]):
        raise ValueError("numerical mixing qualification not met")
    if inv["original_HOSTS"]["candidate_incidences"]!=75346 or inv[
        "BCE_Clarke_partial_augmented"
    ]["candidate_incidences"]!=103349:
        raise ValueError("host incidence denominators changed")
    if turnover["affected_81"]["retained_in_both"]!=3225 or turnover[
        "affected_81"
    ]["union_of_either_source"]!=6445:
        raise ValueError("new-region source Jaccard changed")
    if [x["observed_absorption_minus_null"] for x in local["chains"]] != [271,276]:
        raise ValueError("spatial residual evidence changed")
    for token in ("75,346","103,349","3,225","6,445","0.500","0.36","not","post-hoc"):
        if token not in text:
            raise ValueError(f"key evidence or qualifier absent: {token}")
    def exists(key:str)->None:
        if not sources[key].is_dir():
            raise ValueError(f"source figure artifact missing: {key}")
    for key in ("combined","species_turnover","context"):
        exists(key)
    return {
        "abstract_words":wc(abstract),
        "main_text_words":wc(main),
        "references":len([r for r in refs.splitlines() if r.startswith("- ")]),
        "focal_butterflies":239,
        "bce_source_augmented_butterflies":81,
        "outside_native_occurrence_panel":32,
        "bce_source_agreement_jaccard_modified_species":3225/6445,
        "original_host_spatial_residual_incidences":[271,276],
        "claim_ceiling":"Post-hoc, mapped potential resources: no regional larval use or demographic response measured."
    }


def make_review_zip(
    manuscript:Path, receipt_dir:Path, sources:dict[str,Path],
    output_dir:Path, head_sha:str,
) -> dict:
    text=manuscript.read_text(encoding="utf-8")
    checks=review_conditions(text,receipt_dir,sources)
    stage=output_dir/"stage"
    stage.mkdir(parents=True,exist_ok=True)
    files:dict[str,bytes]={}
    files["article/ALTERNATIVE_RESEARCH_ARTICLE.md"]=manuscript.read_bytes()

    for alias, stem in EXPECTED_IMAGES.items():
        origin=sources[
            "context" if alias.startswith("candidate_Fig3") else
            "species_turnover" if alias.startswith("candidate_Fig2") else "combined"
        ]
        for ext in (".png",".svg"):
            p=origin/(stem+ext)
            if not p.is_file() or p.stat().st_size<=1500:
                raise ValueError(f"missing or empty source figure: {p}")
            files[f"figures/{alias}{ext}"]=p.read_bytes()

    numeric_names=(
        "butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json",
        "butterfly_bce_added_region_identity_turnover_v0.1.json",
        "butterfly_spatial_null_longchain_qualified_v0.2.json",
        "butterfly_plant_range_structural_null_two_chains_v0.2.json",
        "butterfly_climate_matched_absorption_v0.1.json",
        "butterfly_occurrence_overlap_null_v0.1.json",
        "butterfly_expansion_equivalence_v0.1.json",
    )
    for fn in numeric_names:
        p=receipt_dir/fn
        if not p.is_file():
            raise ValueError(f"source receipt missing: {fn}")
        files[f"evidence_receipts/{fn}"]=p.read_bytes()

    for source_dir, file, destination in (
        ("species_turnover","turnover_81_butterflies.csv","derived_data/turnover_81_butterflies.csv"),
        ("species_turnover","turnover_v0.1.json","derived_data/full_81_species_turnover.json"),
        ("context","figure3_context_source_values.json","derived_data/figure3_context_source_values.json"),
        ("combined","figure_resource_novelty_values.json","derived_data/figure1_null_source_values.json"),
        ("combined","figure_bce_inventory_sensitivity_values.json","derived_data/figure_S1_host_inventory_values.json"),
    ):
        p=sources[source_dir]/file
        if not p.is_file() or p.stat().st_size==0:
            raise ValueError(f"source artifact supporting data missing: {p}")
        files[destination]=p.read_bytes()

    readme = (
        "# Blinded alternative ecological review package, NOT the submitted GEB bundle\n\n"
        "This package presents a SINGLE alternative interpretation and full article draft "
        "of the already documented 239-butterfly host-resource analysis. "
        "It does NOT indicate a second study or editorial approval, and original "
        "submission PR #38 remains unchanged.\n\n"
        "The source-specific 81-butterfly Jaccard is a catalogue-comparison of "
        "butterfly×WGSRPD3 potential resource cells, NOT a time-series of butterfly "
        "colonisation, and is based on Europe-derived candidate host records. "
        "The original HOSTS-specific 1,500-km climate/geographic null residual "
        "is only 0.36 percentage points and has NOT been recalculated under BCE.\n\n"
        "The three main candidate figures are in figures/: candidate_Fig1, "
        "candidate_Fig2, candidate_Fig3. The first combined figure also carries "
        "the supporting structural-null panels and should be split/relaid out "
        "by the author before actual submission. Supporting FigS1 shows source "
        "budget compositions. Uncertainty and sample denominators are in the "
        "legends and evidence receipts. The manuscript is a Markdown author "
        "review document, not a journal-final typeset/line-numbered PDF.\n\n"
        "The source artifacts were verified from GitHub Actions: "
        "38989771672 (resource accounting), 38992305236 (turnover), "
        "38994098049 (ecological context). Each was generated using frozen source "
        "tables and scripts. This bundle does not create a stable archival DOI.\n"
    )
    files["README_AUTHOR_REVIEW.md"]=readme.encode("utf-8")
    if len(files)!=1+len(EXPECTED_IMAGES)*2+len(numeric_names)+5+1:
        raise ValueError("unexpected reviewer package file count")
    hashes={k:{"sha256":sha256(v),"size_bytes":len(v)} for k,v in sorted(files.items())}
    manifest={
        "schema":"butterfly_resource_novelty_alternative_author_review_bundle_v0.1",
        "status":"REPRODUCIBLE_AUTHOR_REVIEW_ONLY_NOT_SUBMITTED",
        "source_commit":head_sha,
        "checks":checks,
        "original_submission_PR":"https://github.com/zuizui0223/chocho/pull/38",
        "exploratory_PR":"https://github.com/zuizui0223/chocho/pull/43",
        "actions_source_runs":[38989771672,38992305236,38994098049],
        "contents":hashes,
    }
    files["manifest.json"]=(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n").encode()
    output_dir.mkdir(parents=True,exist_ok=True)
    output=output_dir/"butterfly_resource_novelty_alternative_review_v01.zip"
    with zipfile.ZipFile(output,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=8) as archive:
        for key,val in sorted(files.items()):
            info=zipfile.ZipInfo(key,date_time=(2026,10,17,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o644<<16
            archive.writestr(info,val,compress_type=zipfile.ZIP_DEFLATED,compresslevel=8)
    output_sha=sha256(output.read_bytes())
    (output_dir/"REVIEW_PACKAGE_MANIFEST.json").write_text(json.dumps({
        **manifest,
        "bundle_sha256":output_sha,
        "bundle_file":output.name,
    },indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return {"path":str(output),"bytes":output.stat().st_size,"sha256":output_sha,
            "file_count":len(files), "checks":checks}


def main()->None:
    p=argparse.ArgumentParser()
    for x in ("manuscript","receipt_dir","combined_figures_dir",
              "turnover_figures_dir","context_figures_dir","output_dir"):
        p.add_argument("--"+x.replace("_","-"),type=Path,required=True)
    p.add_argument("--head-sha",required=True)
    args=p.parse_args()
    res=make_review_zip(
        args.manuscript,args.receipt_dir,
        {"combined":args.combined_figures_dir,
         "species_turnover":args.turnover_figures_dir,
         "context":args.context_figures_dir},
        args.output_dir,args.head_sha
    )
    print(json.dumps({"status":"PASS_ALTERNATIVE_AUTHOR_REVIEW_PACKAGE",**res},indent=2))


if __name__=="__main__":
    main()
