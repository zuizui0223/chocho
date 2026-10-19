#!/usr/bin/env python3
"""Deterministic source-audited, potentially author-anonymous OFFLINE review ZIP.

All content is source-conditional candidate larval-resource evidence; the
archive is not a raw-provider redistribution permit or an official submission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

INPUTS=(
    "s1_resource_descriptors.csv",
    "anthropogenic_species_metrics.csv",
    "level1_lookup.json",
    "original/insect_host_accepted.csv",
    "original/native_extant_nondoubtful_wgsrpd3.csv",
    "original/extant_nondoubtful_wgsrpd3.csv",
    "bce/accepted_host_edges.json",
    "bce/additional_native_wgsrpd3.csv",
    "bce/additional_contemporary_wgsrpd3.csv",
    "unit_climate_crossfit.csv",
    "level3.geojson",
)
SCRIPTS=(
    "analyze_bce_resource_novelty_source_sensitivity.py",
    "summarize_bce_added_resource_identity_turnover.py",
    "analyze_butterfly_shared_host_structural_null.py",
    "analyze_butterfly_climate_matched_plant_range_null.py",
    "analyze_butterfly_spatially_local_plant_range_null.py",
    "audit_butterfly_spatial_null_longchain.py",
    "render_butterfly_resource_novelty_figure.py",
    "render_butterfly_resource_novelty_separated_figures.py",
    "render_bce_resource_novelty_sensitivity_figure.py",
    "render_bce_added_resource_identity_turnover_figure.py",
    "render_butterfly_resource_ecological_context_figure.py",
    "build_wcvp_hosts_null_sidecars.R",
    "build_bce_additional_wcvp_ranges.R",
)
EXPECTED_SHA={
    "inputs/s1_resource_descriptors.csv":
        "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc",
    "inputs/unit_climate_crossfit.csv":
        "0087117e0915a0847abaae8d11762b75f35b8eac94434d50913c2781bea727e2",
    "inputs/level3.geojson":
        "c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840",
    "inputs/bce/accepted_host_edges.json":
        "f46d18703fe5cc8d97699af6f084f652a704e46b0f9e8c68e32dd70decc32c53",
    "evidence/original_HOSTS_spatial_mcmc_draws.json":
        "22c60accba86b576aa6ddc9f46ca7dd4d2003e2e81f1ba333ef41da5dbabbb91",
}
SOURCE_METADATA_KEYS=(
    "source_workflow_run",
    "BCE_source_workflow_run",
    "BCE_sensitivity_workflow_run",
)
PRIVACY_RE={
    "source-owner identifier":r"(?i)zuizui0223",
    "submitter name in Roman script":r"(?i)\b(?:ZHANG[ ]+Ruiqi|Ruiqi[ ]+ZHANG)\b",
    "submitter name in Han characters":r"(?:張瑞琪|张瑞琪)",
    "email address":r"(?i)\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
    "account-specific local home":r"(?i)(?:/home/zuizui/|C:\\Users\\zuizui\\)",
    "private client tokens":r"(?i)file_0000000[0-9a-f]+",
}


def hash_bytes(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def scrub_script(name:str,data:bytes)->bytes:
    if not name.endswith(".py"):
        return data
    s=data.decode("utf-8")
    # Provenance-only source-run indices: removing these does not alter
    # species/host mapping, fitted margins, or any ecological statistic.
    for key in SOURCE_METADATA_KEYS:
        s=re.sub(r"(?m)^[ \t]*[\"']"+re.escape(key)+
                 r"[\"'][ \t]*:[ \t]*[0-9]+[ \t]*,[ \t]*\n","",s)
    return s.encode("utf-8")


def privacy_audit(members:dict[str,bytes])->dict:
    audited=0
    for path,data in members.items():
        if path.startswith("/") or ".." in Path(path).parts:
            raise ValueError("relative-path escape in source ZIP")
        for label,pattern in PRIVACY_RE.items():
            if re.search(pattern,path):
                raise ValueError("identifying ZIP member path ("+label+")")
        if not path.lower().endswith((".py",".r",".json",".csv",".geojson",".md")):
            continue
        if b"\x00" in data:
            raise ValueError("binary-like text input not audited: "+path)
        src=data.decode("utf-8")
        audited+=1
        for label,pattern in PRIVACY_RE.items():
            if re.search(pattern,src):
                raise ValueError("anonymity pattern "+label+" in member "+path)
        if path.startswith("code/") and path.endswith(".py"):
            for key in SOURCE_METADATA_KEYS:
                if re.search(r"(?m)^\s*[\"']"+key+r"[\"']\s*:",src):
                    raise ValueError("source-linked external workflow ID in "+path)
    return {
        "status":"PASS_STATIC_DIRECT_IDENTIFIER_SCAN",
        "text_members_examined":audited,
        "patterns_checked":sorted(PRIVACY_RE),
        "limitation":"A static text scan does not establish licence rights or guarantee complete human-identity removal.",
    }


def build(repo:Path,inputs:Path,qualified:Path,archive:Path)->dict:
    members={}
    for name in INPUTS:
        p=inputs/name
        if not p.is_file() or not p.stat().st_size:
            raise ValueError("missing source material: "+name)
        members["inputs/"+name]=p.read_bytes()
    for name in SCRIPTS:
        f=repo/"scripts"/name
        if not f.is_file():
            raise ValueError("missing source script: "+name)
        members["code/"+name]=scrub_script(name,f.read_bytes())
    members["RUN_ME.py"]=(repo/"scripts/anonymous_resource_novelty_offline_runner.py").read_bytes()
    members["README.md"]=(repo/"docs/exploratory/ANONYMOUS_RESOURCE_NOVELTY_ARCHIVE_README_20261017.md").read_bytes()
    members["evidence/original_HOSTS_spatial_mcmc_draws.json"]=qualified.read_bytes()
    for path,expected in EXPECTED_SHA.items():
        if hash_bytes(members[path]) != expected:
            raise ValueError("historical independent data hash drift: "+path)
    # Audit full archived MCMC draw arrays, not only a transcribed summary.
    original=json.loads(members["evidence/original_HOSTS_spatial_mcmc_draws.json"])
    if original.get("status")!="QUALIFIED_EXPLORATORY_MCMC_SENSITIVITY":
        raise ValueError("archived original HOSTS MCMC gate was not qualified")
    if len(original["null"]["chains"])!=2 or not all(
        x["pass"] for x in original["qualification"]["checks"]
    ):
        raise ValueError("missing qualified two-chain MCMC samples")
    if any(len(x["primary_null_draws"])!=99 for x in original["null"]["chains"]):
        raise ValueError("incomplete 99-draw original HOSTS null")

    privacy=privacy_audit(members)
    members["ANONYMITY_AUDIT.json"]=(
        json.dumps(privacy,indent=2,sort_keys=True)+"\n").encode("utf-8")
    checksums={key:hash_bytes(val) for key,val in sorted(members.items())}
    manifest={
        "schema":"butterfly_anonymous_offline_resource_package_v0.1",
        "status":"AUTHOR_REVIEW_ANONYMISED_SOURCE_CONDITIONAL_NOT_A_PUBLIC_DATA_DEPOSIT",
        "sha256":checksums,
        "inputs_count":len(INPUTS),
        "code_count":len(SCRIPTS)+1,
        "upstream_source_pins":{
            "HOSTS":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "WCVP":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
            "TDWG":"52da7828aba9d461dd133c27b3bd7a4407161f54",
        },
        "sources_conditioned_on_host_catalogue":True,
        "source_licence_approval_required_before_distribution":True,
        "disclaimer":"Local host use and larval recruitment are not inferred; external code/data licensing requires manual approval.",
    }
    members["SOURCE_MANIFEST.json"]=(
        json.dumps(manifest,indent=2,sort_keys=True)+"\n").encode()
    privacy_audit(members)
    archive.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=8) as z:
        for path,data in sorted(members.items()):
            info=zipfile.ZipInfo(path,date_time=(2026,10,17,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o644<<16
            z.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=8)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or set(z.namelist())!=set(members):
            raise ValueError("ZIP byte or source manifest corruption")
    return {
        "status":"PASS_ANONYMISED_SOURCE_ARCHIVE_ASSEMBLED_NOT_LICENSE_APPROVED",
        "name":archive.name,
        "bytes":archive.stat().st_size,
        "sha256":hash_bytes(archive.read_bytes()),
        "members":len(members),
        "anonymity_audit":privacy,
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--repo-root",type=Path,required=True)
    p.add_argument("--input-dir",type=Path,required=True)
    p.add_argument("--qualified-spatial-null-json",type=Path,required=True)
    p.add_argument("--output-zip",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(build(a.repo_root,a.input_dir,a.qualified_spatial_null_json,a.output_zip),indent=2))


if __name__=="__main__":
    main()
