#!/usr/bin/env python3
"""Assemble four verified, SOURCE-LIMITED GEB v0.5 replay components.

This builds an internal, screened and source-scoped **candidate**, not a legal
release, complete reviewer data/code access, anonymous deposit, or DOI.
Every component must have already passed independent extracted-archive replay.
No transformations silently upgrade prior stored-score replays to raw refits.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

COMPONENTS = {
    "fig1_fig3bc": (
        "GEB_V05_LIMITED_FIG1_FIG3BC_REPLAY_FIG3A_RECEIPT_ONLY_NOT_FULL_REVIEW.zip",
        "5f35415b81d26f87538381c53376f4ab8e4069a85069e7d510ddff70fdd9fcec", 7,
    ),
    "fig2abc_derived": (
        "GEB_V05_FIG2ABC_DERIVED_ONLY_NOT_FULL_REVIEW.zip",
        "0d34e188fb6ee2184a00f54853b61d49dfa085008679e6c537bae2e9386357d3", 11,
    ),
    "fig2d_fixed239": (
        "GEB_V05_SCREENED_LIMITED_FIXED239_NOT_FULL_REVIEW_ARCHIVE.zip",
        "2b4a598f9baf830f187baef194c3cbe2bd513e8a35052697bf5cad88cb186c43", 8,
    ),
    "fig4_scores_only": (
        "GEB_V05_FIG4_155_224_ORIGINAL_EVALUATION_SCORES_ONLY_NOT_MODEL_REPLAY.zip",
        "aa87dcf6b58f68d2533557e2940a1b294b372ddd3c096257dc6e0de26bf74d25", 5,
    ),
}
ROOT_RECEIPTS=Path("results/geb_v05_private_source_recovery_20261120")
SOURCE_RECEIPTS={
    "Fig3A_original_host_fractional_credit_rebuild.json":
        ROOT_RECEIPTS/"ORIGINAL_FIG3A_WCVP_HOSTS_CREDIT_REBUILT_SOURCE_QA_20261121.json",
    "GBIF224_original_raw_annual_rebuild.json":
        ROOT_RECEIPTS/"ORIGINAL_GBIF224_RAW_ANNUAL_REBUILD_20261121.json",
    "GBIF224_original_169_region_ecological_model_refit.json":
        ROOT_RECEIPTS/"GBIF224_ORIGINAL_FULL_MODEL_REFIT_20261121.json",
    "GBIF224_original_source_code_closure.json":
        ROOT_RECEIPTS/"ORIGINAL_GBIF224_SOURCE_CODE_CLOSURE_20261121.json",
    "GBIF155_original_source_and_missing_RF_panel_gate.json":
        ROOT_RECEIPTS/"ORIGINAL155_FORECAST_SOURCE_AND_PANEL_GATE_20261121.json",
}
BAD={
    "owner_account": re.compile(rb"zuizui0223",re.I),
    "corresponding_author_name": re.compile(rb"ZHANG[ _\.-]?(?:RUIQI|RACHEL)",re.I),
    "personal_email": re.compile(rb"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}",re.I),
    "private_home": re.compile(rb"/(?:home|Users)/[A-Z0-9_.-]+/",re.I),
}
# Names of original private inputs may be mentioned as documented absences,
# but their raw bytes, occurrence IDs, or private longitude/latitude are not
# permitted in the released technical candidate.
ROW_ID_HEAD=(b"gbifid",b"occurrenceid",b"datasetkey")
README="""# GEB v0.5 — internal four-component source replay candidate

**IMPORTANT: INTERNAL PREPARATION ONLY. NOT A SUBMISSION-READY REVIEW ARCHIVE.**
This source-candidate has no authorized anonymous reviewer URL, no actual
permanent archive DOI, and no completed source-by-source redistribution audit.
The archive intentionally contains no GBIF point occurrences, individual
records, named photographs, original private WPI GeoJSON, or private historical
butterfly–island–year source.

Goal: reproduce existing *source-qualified*, noncausal resource opportunity
numbers without inventing source data, historical shipping, or butterfly
colonization. The four components have different response definitions:

1. `fig1_fig3bc`: exact 239-species aggregate resource opportunity baseline
   for 126 insular and 237 mainland botanical regions, and source butterfly
   geography/rank/group statistics; Fig3A inside this component is only an
   old numerical source receipt, but the separate source receipt below records
   a later original-WCVP/HOSTS Fig3A independent re-computation.
2. `fig2abc_derived`: matched port and whole-flora HC3 regressions, and
   count-derived native/introduced Jaccard on 126 botanical islands. The
   separate 104-island marine analysis refits the one-island and strict
   both-islands-unseen electrical-resistance tests from derived inputs.
   This component includes an eight-appendix coverage register.
3. `fig2d_fixed239`: original fixed-239 potential-food and mainland-distance
   model in 126 botanical islands; MEM10 spatial adjustment for 123 relevant
   islands; exact-source coefficients on previously derived tables.
4. `fig4_scores_only`: fixed model-score summaries for the separate 155-
   and 224-butterfly response models using one documented GBIF acquisition. This component
   performs **NO raw GBIF refit**. A separate immutable provenance receipt
   records the successful original 224-ecology-model re-fit with 676 heldout
   model-fold checks and the successful 14,637,152 record original ZIP-to-annual
   source aggregation. That raw data and model prerequisites are not included
   here.

Original 155 recovered-input/model-refit verification is tracked separately.
The included 155 source receipt is a historical pre-recovery audit. The present
archive replays saved 155 evaluation scores and contains no recovered RF inputs.

**Scope exclusions:** the original point-record filtering, upstream plant/host
extraction, the original WPI conversion, full main and SI source replay, and
third-party redistribution clearance.
The original 225x51 and 1596-port GBIF facets remain attributable to saved
speciesKey query receipts, not a fabricated GBIF download DOI or verified
historical publisher list.

Underlying source citations (not synthetic DOIs):
- Natural History Museum HOSTS, doi:10.5519/havt50xw, CC0
- Kew World Checklist of Vascular Plants v13, doi:10.34885/0yex-xv26,
  CC BY 3.0, 21 May 2024
- GBIF original 224-butterfly download used for 155/224 separate models,
  doi:10.15468/dl.umqmsf; each of its 1,716 source datasets has its own
  acknowledgement/licence terms. This DOI does NOT apply to unrelated 225x51
  island facet queries or 1596 botanical-port facet queries.
- World Port Index Publication 150, 2019; the exact edition of the converted
  original study GeoJSON is not yet source independently verified.
- GSHHG coastline version 2.3.7, doi:10.5281/zenodo.7007502.
- TDWG WGSRPD3 pinned original polygons per source replay component.

To test the limited components after extraction with Python 3.12 and internet:
`bash RUN_LIMITED_COMPONENTS.sh`.
Each subcomponent keeps its own terms and caveats. This is an internal
reproducibility feasibility package, not a dryad deposit or double-anonymous
source archive until all licence/anonymity/access gates explicitly pass.
"""
RUN="""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
( cd components/fig1_fig3bc && bash REPRODUCE.sh )
( cd components/fig4_scores_only && bash REPRODUCE.sh )
( cd components/fig2d_fixed239 && bash REPRODUCE.sh )
( cd components/fig2abc_derived && python replay_fig2abc.py )
echo GEB_V05_LIMITED_FOUR_COMPONENT_REPLAY_PASS_NOT_FULL_PAPER
"""

def checksum(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()

def zipinfo(name:str)->ZipInfo:
    obj=ZipInfo(name,(2026,11,21,0,0,0))
    obj.external_attr=0o644<<16
    obj.compress_type=ZIP_DEFLATED
    return obj

def secure_name(name:str)->bool:
    p=Path(name)
    return bool(name and not p.is_absolute() and ".." not in p.parts and
                not name.startswith(("/",".\\")) and "\\" not in name)

def source_items(incoming:Path, audit:Path) -> tuple[dict[str,bytes],dict]:
    entries={"README_INTERNAL_NOT_FOR_JOURNAL_REVIEW.md":README.encode("utf-8"),
             "RUN_LIMITED_COMPONENTS.sh":RUN.encode("utf-8")}
    provenance={}
    for key, (archive_name, expected_sha, expected_count) in COMPONENTS.items():
        found=list(incoming.rglob(archive_name))
        if len(found)!=1:
            raise AssertionError(f"Expected exactly one {archive_name}; got {len(found)}")
        f=found[0]
        archive_sha = checksum(f.read_bytes())
        if archive_sha != expected_sha:
            raise AssertionError(f"Unexpected source component SHA256 for {archive_name}: {archive_sha}")
        with ZipFile(f) as z:
            if z.testzip() is not None:raise AssertionError("CRC verification failed")
            nested = [item for item in z.infolist() if not item.is_dir()]
            if len(nested) != expected_count:
                raise AssertionError("Unexpected scoped source size in " + archive_name)
            for item in nested:
                if not secure_name(item.filename):
                    raise AssertionError("Unsafe source name "+item.filename)
                name=f"components/{key}/{item.filename}"
                if name in entries:raise AssertionError("Duplicate source path "+name)
                content=z.read(item)
                if name.endswith(".csv"):
                    header=content.splitlines()[0].lower().replace(b"_",b"")
                    if any(v in header for v in ROW_ID_HEAD):
                        raise AssertionError("Original GBIF occurrence-level input included: "+name)
                entries[name]=content
        provenance[key]={"source_artifact_component_filename":archive_name,
                         "source_artifact_component_sha256":archive_sha,
                         "file_count":len(nested)}
    for key,relative in SOURCE_RECEIPTS.items():
        item=audit/relative
        if not item.is_file():raise AssertionError("Missing verified source receipt: "+str(relative))
        data=json.loads(item.read_text(encoding="utf-8"))
        if not isinstance(data,dict) or "status" not in data:
            raise AssertionError("Not a valid source QA receipt "+key)
        entries["source_qa_receipts/"+key]=item.read_bytes()
    return entries,provenance

def build(incoming:Path,audit:Path,out:Path)->dict:
    entries,provenance=source_items(incoming,audit)
    manifest={
        "schema":"GEB_v05_internal_4component_extracted_replay_candidate_v2",
        "status":"INTERNAL_ONLY_NOT_ANONYMITY_OR_LICENCE_APPROVED",
        "component_source":provenance,
        "source_qa_receipts":list(SOURCE_RECEIPTS),
        "known_main_figures_covered_by_limited_archives":["Fig1","Fig2A-C","Fig2D","Fig3B-C","Fig4_stored_scores"],
        "Fig3A_original_source_rebuilt_outside_composite":True,
        "original_GBIF224_full_model_rebuilt_outside_composite":True,
        "original_GBIF155_recovered_input_model_refit_verification":"tracked separately",
        "original_GBIF155_RF_training_panel_included":False,
        "full_SI_replayed":False,
        "third_party_original_data_included":False,
        "source_licensing_and_redistribution_clearance":False,
        "indirect_author_anonymity_certified":False,
        "reviewer_access_URL_exists":False,
        "permanent_archive_DOI_exists":False,
        "scientific_manuscript_source_modified":False,
        "limitations":"Scoped derived-input and saved-score replay. Full upstream butterfly/plant source processing, original WPI conversion, and full SI replay remain outside this archive."
    }
    entries["MANIFEST_INTERNAL_LIMITED_SCOPE.json"]=(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n").encode()
    for name,content in entries.items():
        if not secure_name(name):raise AssertionError("unsafe output path: "+name)
        for label,pattern in BAD.items():
            if pattern.search(name.encode()) or pattern.search(content):
                raise AssertionError(f"Direct-name/source scan flagged {name} : {label}")
    out.mkdir(parents=True,exist_ok=True)
    dest=out/"GEB_V05_INTERNAL_COMPOSITE_REPRO_NOT_FOR_JOURNAL_REVIEW.zip"
    with ZipFile(dest,"w") as z:
        for name,content in sorted(entries.items()):
            z.writestr(zipinfo(name),content)
    with ZipFile(dest) as z:
        if z.testzip() is not None or sorted(z.namelist())!=sorted(entries):
            raise AssertionError("Zipped archive content mismatch")
        for item in z.infolist():
            if item.date_time!=(2026,11,21,0,0,0):
                raise AssertionError("Non-deterministic ZIP date")
            if (item.external_attr>>16)&0o777!=0o644:
                raise AssertionError("Unexpected ZIP mode")
            content=z.read(item)
            if content != entries[item.filename]:
                raise AssertionError("Archived member differs from verified source: " + item.filename)
            for label,pattern in BAD.items():
                if pattern.search(item.filename.encode()) or pattern.search(content):
                    raise AssertionError(f"Post-archive source scan flagged {label}")
        z.extractall(out/"extracted_qa")
    report={
        "status":manifest["status"],
        "file_count":len(entries),
        "component_file_counts":{k:v["file_count"] for k,v in provenance.items()},
        "source_receipt_count":len(SOURCE_RECEIPTS),
        "zip_sha256":checksum(dest.read_bytes()),
        "zip_size":dest.stat().st_size,
        "no_source_identity_tokens_in_actual_ZIP":True,
        "real_Dryad_reviewer_deposit":False,
        "current_science_unchanged":True,
    }
    (out/"INTERNAL_COMPOSITE_BUILD_QA.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--incoming",type=Path,required=True)
    p.add_argument("--audit-root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    x=p.parse_args()
    print(json.dumps(build(x.incoming,x.audit_root,x.out),ensure_ascii=False))

if __name__=="__main__":
    main()
