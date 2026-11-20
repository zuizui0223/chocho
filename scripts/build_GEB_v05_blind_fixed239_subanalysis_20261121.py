#!/usr/bin/env python3
"""Build a limited author-identifier-screened reproduction package for GEB Fig2D.

This is NOT a complete anonymous journal review archive or source licence check.
Input data are existing derived aggregate region tables, NOT original WCVP,
HOSTS, WPI2019, GBIF occurrence level inputs or third-party media.
The generated package is only for controlled coauthor QA before rights approval.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

CODES = {
 "scripts/audit_island_port_isolation_239_20261120.py":"184c92d9b2c6fba7f15949bf1bb0945f8d249904",
 "scripts/audit_island_port_spatial_239_20261120.py":"fc00f7ae4955376a2fa3d24b5c0373a5621d79be",
}
TABLES = {
 "results/island-mainland-resource-q1/region_metrics.csv":(
   "data/region_metrics.csv","4aba89565354edc6c44d566523864cec66a30cbe"),
 "results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv":(
   "data/whole_flora_wpi_126_islands.csv","53b618fa51a1a418ccc6e4363f4c866c389c90cc"),
}
RECEIPT="results/isolation_port_resource_reanalysis_20261120/SOURCE239_PAPER_RECEIPT_V01.json"
URL="https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson"
SHA="c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840"
BLOCKED={
 "account_token":re.compile(rb"zuizui0223",re.I),
 "author_name":re.compile(rb"ZHANG[ _.-]?(?:RUI"+rb"QI|RAC"+rb"HEL)|\bru"+rb"iqi\b|\bre"+rb"iqi\b",re.I),
 "email":re.compile(rb"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}",re.I),
 "local_user_home":re.compile(rb"/(?:Users|home)/[A-Z0-9._-]+/",re.I),
 "raw_historical_private_data":re.compile(rb"global_original_species_island_year\.parquet|ports_2019\.geojson",re.I),
}
README="""# GEB botanical-island known-butterfly-host opportunity: limited reproduction

This **limited** author-identifier-screened source package can replay the
126 botanical-island / fixed-239 documented butterfly host-resource isolation
subanalysis (including 10 Moran eigenvectors). It cannot reproduce all paper
Figures 1–4, independently observed butterflies, original 1,596-port speciesKey
queries, true ship trajectories, individual larval feeding, or abundance.

The two already-derived aggregate data tables represent 126 botanical islands
and 369 WGSRPD3 regions. The fixed WCVP/HOSTS individual host-edge extraction
is not provided. No raw GBIF sighting IDs or private 2019 port source files
are included. This is **not a full original-source replay**.

CITE PRIMARY SOURCES:
- Natural History Museum HOSTS database, DOI 10.5519/havt50xw (larval hosts).
- Kew World Checklist of Vascular Plants v13, 21 May 2024,
  DOI 10.34885/0yex-xv26 (native and introduced plant distributions).
- NGA World Port Index 2019, Publication 150, 27th edition.
- TDWG World Geographical Scheme for Recording Plant Distributions (WGSRPD3),
  publicly retrieved from the source-pinned file in REPRODUCE.sh.

RERUN (Python 3.12 with internet): bash REPRODUCE.sh.
The runner fetches a public TDWG boundary file, verifies its SHA256, fits the
fixed-239 resource and spatial models and tests results against a numeric
source receipt with tolerance 1e-4.

SOURCE AND LEGAL LIMITS: These derived public-table contents require
source-specific licence/redistribution review and corresponding-author approval
before inclusion in an actual double-anonymous peer-review deposit. A direct
identity string scan cannot guarantee full indirect author anonymity. This ZIP
is not a citable stable repository; no DOI was minted. Original independent
155/224 butterfly records use a separate GBIF DOI 10.15468/dl.umqmsf, not the
unattributed 225x51 or 1596-port botanical facet sources.
"""
SHELL="""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python -m pip install --disable-pip-version-check \
  "numpy>=1.26,<3" "pandas>=2.2,<3" \
  "statsmodels>=0.14,<1" "scipy>=1.13,<2" \
  "shapely>=2,<3" "pyproj>=3.6,<4" "matplotlib>=3.9,<4"
mkdir -p external_geojson recomputed
curl --fail --silent --show-error --location --retry 4 \
  https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson \
  -o external_geojson/level3.geojson
python - <<'PY'
from pathlib import Path
import hashlib
p=Path("external_geojson/level3.geojson")
assert hashlib.sha256(p.read_bytes()).hexdigest()=="c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840"
print("TDWG_2026_SOURCE_SHA256_PASS")
PY
python scripts/audit_island_port_isolation_239_20261120.py \
  --port-flora data/whole_flora_wpi_126_islands.csv \
  --region-hosts data/region_metrics.csv \
  --tdwg-geojson external_geojson/level3.geojson \
  --out-dir recomputed
python scripts/audit_island_port_spatial_239_20261120.py \
  --island-csv recomputed/fixed239_source_islands_with_isolation.csv \
  --tdwg-geojson external_geojson/level3.geojson \
  --out-dir recomputed
python VERIFY_REPLAY.py
"""
VERIFY="""#!/usr/bin/env python3
import json
from pathlib import Path
E=json.loads(Path("EXPECTED_FIXED239_SUBMODEL.json").read_text(encoding="utf8"))
B=json.loads(Path("recomputed/SOURCE239_ISOLATION_PORT_REANALYSIS.json").read_text(encoding="utf8"))
S=json.loads(Path("recomputed/SOURCE239_MORAN_MEM10_SPATIAL_SENSITIVITY.json").read_text(encoding="utf8"))
def check(k,value):
    assert abs(float(value)-float(E[k]))<0.0001,(k,value,E[k])
M=B["models"]
check("port_intro",M["plant_introduced"]["ports"]["beta"])
check("native_distance",M["host_native_239"]["distance"]["beta"])
check("paired_distance",M["host_current_minus_native_z"]["distance"]["beta"])
check("intro_only_distance",M["introduced_only_share_logodds_239"]["distance"]["beta"])
for model,term,key in [
 ("host_native_239","distance","MEM10_native_distance"),
 ("host_current_minus_native_z","distance","MEM10_paired_distance"),
 ("introduced_only_share_logodds_239","distance","MEM10_intro_only_distance"),
 ("plant_introduced","ports","MEM10_port_intro")]:
    check(key,S["results"][model]["spatial_MEM10"][term]["beta"])
assert B["models"]["counts"]["islands"]==E["islands"]==126
assert S["active_resource_islands"]==E["active_islands"]==123
print("GEB_LIMITED_FIXED239_INDEPENDENT_REVIEW_SOURCE_REPLAY_PASS")
"""
def blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
def sha(data):
    return hashlib.sha256(data).hexdigest()
def pinned(path,required):
    data=path.read_bytes()
    if blob(data)!=required:
        raise AssertionError("Subanalysis original aggregate blob mismatch: "+str(path))
    return data
def info(name):
    x=ZipInfo(name,(2026,11,21,0,0,0))
    x.external_attr=0o644<<16
    x.compress_type=ZIP_DEFLATED
    return x
def build(root,out):
    out.mkdir(parents=True,exist_ok=True)
    original=json.loads((root/RECEIPT).read_text(encoding="utf8"))
    if original.get("schema")!="chocho_GEB_fixed239_isolation_port_paper_receipt_20261120_v01":
        raise AssertionError("Original numeric source receipt schema changed")
    m=original["models"];s=original["spatial_sensitivity"]["spatial_mem10"]
    exp={
      "schema":"limited_GEB_island239_host_port_replay_expected_v1",
      "islands":126,"active_islands":123,"butterflies":239,
      "port_intro":m["introduced_flora"]["ports_beta"],
      "native_distance":m["native_host239"]["distance_beta"],
      "paired_distance":m["paired_native_current_shift"]["distance_beta"],
      "intro_only_distance":m["introduced_only_share_logodds"]["distance_beta"],
      "MEM10_native_distance":s["native_host"]["distance_beta"],
      "MEM10_paired_distance":s["paired_current_minus_native"]["distance_beta"],
      "MEM10_intro_only_distance":s["introduced_only_share_logodds"]["distance_beta"],
      "MEM10_port_intro":s["introduced_whole_flora"]["port_beta"],
      "model_kind":"potential larval host geographic coverage not ecological occupancy",
    }
    data={
      "README.md":README.encode(),
      "REPRODUCE.sh":SHELL.encode(),
      "VERIFY_REPLAY.py":VERIFY.encode(),
      "EXPECTED_FIXED239_SUBMODEL.json":(json.dumps(exp,indent=2)+"\n").encode(),
    }
    for path,digest in CODES.items():data[path]=pinned(root/path,digest)
    for path,(rel,digest) in TABLES.items():data[rel]=pinned(root/path,digest)
    if len(data)!=8:
        raise AssertionError("Unexpected limited source set")
    for name,payload in data.items():
        if name.startswith(("/",".")) or ".." in Path(name).parts:
            raise AssertionError("Unsafe ZIP relative path")
        for key,rx in BLOCKED.items():
            if rx.search(name.encode()) or rx.search(payload):
                raise AssertionError("Author-identifier scan blocked source: "+name+" "+key)
        if name.endswith(".csv"):
            head=payload.splitlines()[0]
            if any(word in head for word in [b"gbifID",b"datasetKey",b"occurrenceID"]):
                raise AssertionError("Original occurrence record-level file mistakenly included")
    zipname="GEB_V05_SCREENED_LIMITED_FIXED239_NOT_FULL_REVIEW_ARCHIVE.zip"
    dest=out/zipname
    with ZipFile(dest,"w") as z:
        for name,payload in sorted(data.items()):
            z.writestr(info(name),payload)
    with ZipFile(dest) as z:
        if z.testzip() is not None or sorted(z.namelist())!=sorted(data):
            raise AssertionError("Source ZIP integrity failure")
        for item in z.infolist():
            if item.date_time!=(2026,11,21,0,0,0) or (item.external_attr>>16)&0o777!=0o644:
                raise AssertionError("ZIP metadata identity or archive mode drift")
            for key,rx in BLOCKED.items():
                if rx.search(item.filename.encode()) or rx.search(z.read(item)):
                    raise AssertionError("ZIP name/content after-extraction identity scan blocked "+key)
    receipt={
      "schema":"GEB_v05_limited_fixed239_anonymity_screen_QA_v1",
      "status":"ANONYMITY_SCREENED_LIMITED_SOURCE_SUBANALYSIS_NOT_COMPLETE_REVIEW_ARCHIVE",
      "source_aggregate_counts":{"botanical_islands":126,"regions":369,"candidate_butterflies":239},
      "included_source_files":sorted(data),
      "included_source_sha256":{name:sha(data[name]) for name in sorted(data)},
      "zip_sha256":sha(dest.read_bytes()),
      "source_identity_scan":"PASS_DIRECT_TEXT_AND_ZIP_HEADER_ONLY",
      "raw_source_WCVP_HOSTS_and_GBIF":"NOT_INCLUDED",
      "independent_butterfly_observation_cohorts":"NOT_REPLAYED",
      "third_party_derived_table_rights":"PENDING_CORRESPONDING_AUTHOR_APPROVAL",
      "complete_GEB_anonymous_review_repository":"NOT_CREATED",
      "DOI_minted":False,
      "public_geometry_url":URL,"public_geometry_sha256":SHA,
    }
    (out/"LIMITED_FIXED239_SCREEN_AUDIT.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GEB_LIMITED_FIXED239_EIGHT_SOURCE_FILES_ANON_IDENTIFIER_SCREEN_PASS",
          json.dumps({"zip_sha256":receipt["zip_sha256"],"complete_archive":False}))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--current-root",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    build(a.current_root,a.out_dir)
if __name__=="__main__":
    main()
