#!/usr/bin/env python3
"""Build an auditable, anonymized optional GEB supplement for BCE/Clarke checks.

Original BCE pages are secondary published checklists, NOT independent larval
performance experiments. All source snapshots are explicit immutable outputs
of pinned GitHub Actions runs; the archive is not yet part of the paper PR.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE_FILES={
    "inputs/bce_original_pages_derived_Clarke_2024.json":"lexical/evidence_filtered_all92_exact_host_names_v01.json",
    "inputs/bce_accepted_species_host_edges.json":"wcvp/bce_host_accepted_taxon_gap_v01.json",
    "inputs/bce_wcvp_name_map.csv":"wcvp/wcvp_name_map.csv",
    "inputs/bce_candidate_names.csv":"wcvp/candidate_names.csv",
    "results/bce_geography_sensitivity.json":"geo/partial_europe_source_host_geography_v01.json",
    "results/bce_continent_355_region_null.json":"level1/bce_clarke_continent_fixed_margin_comparison_v01.json",
}
REPO_FILES=(
    "scripts/audit_bce_clarke_evidence_host_pilot.py",
    "scripts/audit_bce_clarke_full_exact_overlap.py",
    "scripts/audit_bce_clarke_taxonomic_crosswalk.py",
    "scripts/resolve_bce_clarke_plant_synonyms_wcvp.R",
    "scripts/build_wcvp_hosts_null_sidecars.R",
    "scripts/build_bce_clarke_additional_host_wcvp_ranges.R",
    "scripts/analyze_bce_clarke_added_hosts_geography.py",
    "scripts/analyze_butterfly_resource_homogenization.py",
    "scripts/analyze_butterfly_resource_homogenization_level1_null.py",
    "scripts/analyze_bce_clarke_homogenization_fixed_domain.py",
    "scripts/analyze_bce_clarke_homogenization_level1_null.py",
    "tests/test_bce_clarke_level1_null.py",
    "docs/exploratory/BCE_CLARKE_92_SPECIES_HOST_LINK_AUDIT_PROTOCOL_V01.json",
    "docs/exploratory/BCE_CLARKE_TAXONOMIC_HOST_CROSSWALK_PROTOCOL_V01.json",
    "docs/exploratory/BCE_CLARKE_ONE_DIRECTION_HOST_AUGMENTATION_PROTOCOL_V01.json",
    "docs/exploratory/BCE_CLARKE_LEVEL1_NULL_HOST_COMPLETENESS_PROTOCOL_V01.json",
    "docs/exploratory/BCE_CLARKE_HOMOGENIZATION_FIXED_DOMAIN_PROTOCOL_V01.json",
    "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "provenance/reviewer_defenses/results/butterfly_resource_homogenization_v0.1.json",
)
BANNED=("zuizui0223","zhang.ruiqi","ruiqi","zhang ruiqi")
STRIP_KEYS={
    "artifact_id","artifact_name","run_id","workflow_run_id","source_run_id",
    "source_workflow_run","crosswalk_workflow_run","prior_geography_run",
    "BCE_accepted_host_crosswalk_run","previous_geography",
}

README="""# Anonymized appendix: butterfly larval-host knowledge sensitivity

This optional supplementary review archive documents a POST-HOC source
robustness test on 239 butterfly species and a fixed 355-region analysis domain.
The extra European butterfly foodplant associations are derived from public
Butterfly Conservation Europe (BCE) pages which curate Clarke (2024),
DOI 10.1002/ece3.10834, evidence ranks 1–3. They are NOT independent original
field experiments, exhaustive global host-use, or adult-survival outcomes.
The original Clarke relational Dryad files were unavailable from automated
provider endpoints; the BCE-derived page names and their page SHA256s are frozen
in inputs/bce_original_pages_derived_Clarke_2024.json.

Source versions for reconstruction:
  globalbioticinteractions/HOSTS 808e0b869f9ec1adf8efff87cf6a395adda103e0
  matildabrown/rWCVPdata 65bed76bae9d644ccb6ad200c05f9f5071d89e05
  tdwg/wgsrpd 52da7828aba9d461dd133c27b3bd7a4407161f54

The package includes accepted-WCVP-ID butterfly--foodplant candidate pairs,
accepted-name mapping, 239-species frozen figure metrics, analysis protocols,
and source-specific geography and continent-constrained null result JSON.
It intentionally does not bundle huge provider HOSTS/WCVP/GIS source datasets,
but pins their exact Git commits above. The scripts build their sidecars from
those original provider checkouts. Re-executing source retrieval from the live
BCE website is NOT guaranteed to produce byte-identical data: use the frozen
source records here for scientific comparisons.

Verification:
  python -m unittest discover -s tests -p test_bce_clarke_level1_null.py -v
  python scripts/build_bce_clarke_anonymous_review_addendum.py --inputs-dir \
     <DIRECTORY_WITH_SOURCE_ARTIFACT_FILES> --output <REBUILT_ZIP>

Interpretation:
- Source-conditional +1,027 accepted-ID candidate butterfly-host pairs for 81
  European species are not global HOSTS false-negative probabilities.
- European host records extrapolated to all WCVP regions provide structural
  upper-bound geography, not proof of local larval performance.
- The effect is on potential resource opportunity and is NOT a new causal
  demonstration of homogenized butterfly communities or competition.
- Two conditional nulls, one with native constraints and a stricter one also
  preserving butterfly x WGSRPD Level1 additions, are supplementary analyses.
- p=0.002 is the resolution floor for 499 random draws with +1 correction.
- Raw host website snapshot HTML bytes are NOT present; individual page URL and
  SHA256 metadata accompany the extracted exact names.
"""

def sanitize(value):
    if isinstance(value,dict):
        return {k:sanitize(v) for k,v in value.items() if k not in STRIP_KEYS}
    if isinstance(value,list):
        return [sanitize(v) for v in value]
    if isinstance(value,str):
        return re.sub("zuizui0223","anonymous-user",value,flags=re.IGNORECASE)
    return value

def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def verify_source(source_dir):
    parsed={}
    for key,path in SOURCE_FILES.items():
        full=source_dir/path
        if not full.is_file():raise RuntimeError("Missing pinned source archive input "+path)
        if path.endswith(".json"):
            parsed[key]=json.loads(full.read_text(encoding="utf-8"))
    lexical=parsed["inputs/bce_original_pages_derived_Clarke_2024.json"]
    accepted=parsed["inputs/bce_accepted_species_host_edges.json"]
    geo=parsed["results/bce_geography_sensitivity.json"]
    lev=parsed["results/bce_continent_355_region_null.json"]
    if lexical["schema"]!="chocho_bce_clarke_all_exact_92_source_names_v0.1" or lexical["exact_bce_chocho_intersection"]!=92:
        raise RuntimeError("Unexpected 92-butterfly BCE lexical source")
    if accepted["schema"]!="chocho_bce_clarke92_wcvp_accepted_id_gap_v0.1" or accepted["verified_butterflies"]!=83:
        raise RuntimeError("Unexpected WCVP-reconciled host source")
    if accepted["pooled_BCE_only_taxon_resolutions"].get("ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY")!=1027:
        raise RuntimeError("Candidate host source count changed")
    if geo["schema"]!="chocho_bce_clarke_one_direction_european_augmentation_v0.1" or not geo["frozen_239_baseline_verified"]:
        raise RuntimeError("Unverified 239-butterfly geography sensitivity")
    g=geo["all_239_descriptive_partial_augmented"]
    if [g["baseline_native"],g["augmented_native"],g["baseline_introduced_added"],g["augmented_introduced_added"]] != [26530,29852,14553,14967]:
        raise RuntimeError("Geographic source totals drift")
    if lev["schema"]!="chocho_bce_clarke_level1_null_source_sensitivity_result_v0.1" or not lev["frozen_source_and_domain_verified"]:
        raise RuntimeError("Incomplete Level1 run artifact")
    if lev["shared_region_count"]!=355:
        raise RuntimeError("Changed spatial comparison domain")
    for k in ("original","BCE_augmented"):
        row=lev[k]["level1_fixed_margin_null"]
        if row["permutations"]!=499 or row["seed"]!=20261007:
            raise RuntimeError("Wrong Level1 null protocol for "+k)
        if not row["sampled_margins_verified"]:
            raise RuntimeError("Unverified row/region/continent margins "+k)
        if row["regional_observed_minus_null_median"]<=0 or row["regional_one_sided_p"]!=0.002:
            raise RuntimeError("Level1 original or augmented result does not match")
    return {"original_butterfly_species":239,"regional_comparison_units":355,
            "source_evidence_species":83,"candidate_accepted_host_edges":1027,
            "native_units":[g["baseline_native"],g["augmented_native"]],
            "introduced_added_units":[g["baseline_introduced_added"],g["augmented_introduced_added"]],
            "regional_excess":[lev[k]["level1_fixed_margin_null"]["regional_observed_minus_null_median"] for k in ("original","BCE_augmented")],
            "regional_null_p":[lev[k]["level1_fixed_margin_null"]["regional_one_sided_p"] for k in ("original","BCE_augmented")]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--inputs-dir",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    summary=verify_source(a.inputs_dir)
    with tempfile.TemporaryDirectory(prefix="bce-review-") as temp:
        stage=Path(temp)
        for dest,relative in SOURCE_FILES.items():
            source=a.inputs_dir/relative
            target=stage/dest
            target.parent.mkdir(parents=True,exist_ok=True)
            if source.suffix.lower()==".json":
                obj=json.loads(source.read_text(encoding="utf-8"))
                target.write_text(json.dumps(sanitize(obj),ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
            else:
                shutil.copyfile(source,target)
        for name in REPO_FILES:
            src=ROOT/name
            if not src.is_file():raise RuntimeError("Missing scientific method "+name)
            dst=stage/name
            dst.parent.mkdir(parents=True,exist_ok=True)
            if src.suffix==".json":
                o=json.loads(src.read_text(encoding="utf-8"))
                dst.write_text(json.dumps(sanitize(o),ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
            else:
                txt=src.read_text(encoding="utf-8")
                txt=re.sub("zuizui0223","anonymous-user",txt,flags=re.IGNORECASE)
                dst.write_text(txt,encoding="utf-8")
        (stage/"README.md").write_text(README,encoding="utf-8")
        (stage/"VALIDATED_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
        violations=[]
        for file in sorted(t for t in stage.rglob("*") if t.is_file()):
            text=file.read_text(encoding="utf-8",errors="replace").lower()
            for token in BANNED:
                if token in text:violations.append(str(file.relative_to(stage))+" contains "+token)
        if violations:raise RuntimeError("Double-anonymous identity scan failure: "+"; ".join(violations))
        files=sorted(q for q in stage.rglob("*") if q.is_file())
        entries=[f"{sha256(q)}  {q.relative_to(stage).as_posix()}" for q in files]
        (stage/"SHA256SUMS").write_text("\n".join(entries)+"\n",encoding="utf-8")
        a.output.parent.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(a.output,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
            for file in sorted(t for t in stage.rglob("*") if t.is_file()):
                info=zipfile.ZipInfo(file.relative_to(stage).as_posix(),date_time=(1980,1,1,0,0,0))
                info.compress_type=zipfile.ZIP_DEFLATED
                info.external_attr=0o100644<<16
                z.writestr(info,file.read_bytes())
    print(json.dumps({"status":"SUCCESS_ANON_SOURCE_PINNED_REVIEW_ADDENDUM",
                      "files":len(REPO_FILES)+len(SOURCE_FILES)+3,
                      "sha256":sha256(a.output),"validated_summary":summary},
                     indent=2,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
