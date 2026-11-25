#!/usr/bin/env python3
"""Internal-only extension of the exact 223-file archive; not a reviewer release."""
from __future__ import annotations
import argparse, collections, csv, hashlib, io, json, shutil, zipfile
from pathlib import Path

HIER="77c76099bd492d6a9f5cf989196c84946e7d16aa"
OBS="2b73c8328235508f5ea3fa57b8c0874419cbd0b9"
PAPER="7e182df03dc3f031688159e4ae0ac2914acd0735"
EXPECTED={
  14992228143:"b419b7986865b8d3e82007f11bae6d8acb7c119d949b200f07d9cbee3b370d71",
  14997553306:"4e665dd376606782a84fbb5e2517cc32cd54abcedf0a0a44f3591d09b001ceb5",
  15001345576:"57507353d9bfb81689618b1876b41df29428b197c07444f0eef9f2fc55454682",
}
PREFIX="components/s36_s55_island_consumer/"
SOURCE_FILES={
 "hierarchy": {
 "scripts/audit_fixed239_host_taxonomic_resolution_20261125.py":"ecad5a8669b1f0a9136d35cbf8d52f46da1bc883",
 "scripts/audit_island_port_isolation_239_20261120.py":"184c92d9b2c6fba7f15949bf1bb0945f8d249904",
 "scripts/audit_island_port_spatial_239_20261120.py":"fc00f7ae4955376a2fa3d24b5c0373a5621d79be",
 ".github/workflows/audit-fixed239-host-taxonomic-resolution.yml":"4d1e74f805ec4ce12216f0c16e3c00779b7cf339",
 "results/three-gate-port-botanical-effort/botanical_effort_flora_island_inputs.csv":"23f5e5c7f7c875af41ccb5c6094cb4f48f1071fe",
 },
 "observation": {
 "scripts/audit_gbif224_island_observed_host_isolation_20261125.py":"b7f0242bae88fc01c1c8a8bba23b80d4ecc5d3d8",
 ".github/workflows/audit-gbif224-island-observed-host-isolation.yml":"a78e71452208acf9c5bb49bd2c6ae6ccc315cedb",
 }}
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def readzip(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        assert z.testzip() is None
        ns=[n for n in z.namelist() if not n.endswith("/")]
        assert len(ns)==len(set(ns))
        assert all(not Path(n).is_absolute() and ".." not in Path(n).parts and "\\" not in n for n in ns)
        return {n:z.read(n) for n in ns}
def find(files,name):
    ks=[k for k in files if Path(k).name==name]
    assert len(ks)==1,(name,ks)
    return files[ks[0]]
def table(b):return list(csv.DictReader(io.StringIO(b.decode("utf-8"))))
def js(b):return json.loads(b)
def main():
    p=argparse.ArgumentParser()
    for k in ("downloads","hierarchy","observation","out","extracted","replay"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    arc={}
    for aid,h in EXPECTED.items():
        b=(a.downloads/f"{aid}.zip").read_bytes()
        assert sha(b)==h,(aid,sha(b),h)
        arc[aid]=readzip(b)
    prior=arc[14992228143]
    inner=find(prior,"GEB_V05_223_INTERNAL_STAGING_NOT_FOR_RELEASE.zip")
    assert sha(inner)=="10561b8baf7d87e68df710a664c745e1c9d7a44047ae220ee8bbf3b6c953b84d"
    old=readzip(inner)
    rows=table(find(prior,"EXACT_223_FILE_SOURCE_DELIVERY_DECISIONS.csv"))
    reg=js(find(prior,"SOURCE_REGISTRY.json"))
    prior_summary=js(find(prior,"INTEGRATED_DELIVERY_SUMMARY.json"))
    assert len(old)==len(rows)==223 and len(reg)==26 and prior_summary["authorized_files"]==0
    assert set(old)=={r["file"] for r in rows}
    for r in rows:
        assert r["authorized_for_distribution"]=="False" and sha(old[r["file"]])==r["sha256"] and len(old[r["file"]])==int(r["bytes"])
    h=arc[14997553306]; o=arc[15001345576]
    assert {Path(k).name for k in h}=={"host_ids_outside_accepted_species_background.csv","island_taxonomic_resolution_counts.csv","source_gate.json","species_count_mismatch.csv","taxonomic_resolution_report.json"}
    assert {Path(k).name for k in o}=={"island_aggregates.csv","recorded_island_host_isolation.json"}
    gate=js(find(h,"source_gate.json")); hr=js(find(h,"taxonomic_resolution_report.json")); obs=js(find(o,"recorded_island_host_isolation.json"))
    assert gate["n_islands"]==126 and gate["n_butterflies"]==239 and gate["global_per_butterfly_source_count_gate"]==239
    assert gate["n_species_count_mismatches"]==0 and gate["native_species_total"]==[4829,4829] and gate["current_species_total"]==[10570,10570]
    assert hr["status"]=="EXACT_SPECIES_COUNT_GATE_PASSED_EXPLORATORY_TAXONOMIC_AUDIT"
    assert obs["status"]=="SOURCE_GATED_RECORDED_OPPORTUNITY_ASSOCIATIONS_NOT_OCCUPANCY"
    assert obs["original_species"]==224 and obs["frozen_resource_species"]==239
    assert len(table(find(h,"species_count_mismatch.csv")))==0
    assert len(table(find(h,"host_ids_outside_accepted_species_background.csv")))==37
    ht=table(find(h,"island_taxonomic_resolution_counts.csv"))
    ot=table(find(o,"island_aggregates.csv"))
    assert len(ht)==len(ot)==126 and {r["code"] for r in ht}=={r["code"] for r in ot}
    assert all(int(r["native_species_delta"])==int(r["current_species_delta"])==0 for r in ht)
    assert tuple(sum(int(r["all1_"+k]) for r in ot) for k in ("N","I","Z"))==(916,193,62)
    assert sum(int(r["all1_I"])>0 for r in ot)==69 and sum(int(r["all2_I"]) for r in ot)==143
    assert len(obs["monarch"]["recorded_introduced_only_regions"])==20 and len(obs["monarch"]["recurrently_recorded_introduced_only_regions"])==18
    new={}
    for name,data in h.items():new[PREFIX+"s36/"+Path(name).name]=data
    for name,data in o.items():new[PREFIX+"s55/"+Path(name).name]=data
    for label,files in SOURCE_FILES.items():
        root=getattr(a,label)
        for name,expected in files.items():
            b=(root/name).read_bytes()
            assert blob(b)==expected,(label,name,blob(b),expected)
            new[PREFIX+name]=b
    new[PREFIX+"replay.py"]=a.replay.read_bytes()
    new[PREFIX+"SOURCE_PINS.json"]=(json.dumps({
      "scope":"internal derived-panel statistical replay; original raw-source workflows separately retained",
      "hierarchy":{"head":HIER,"run":45560130341,"artifact":14997553306,"sha256":EXPECTED[14997553306]},
      "observation":{"head":OBS,"run":45567964403,"artifact":15001345576,"sha256":EXPECTED[15001345576]},
      "paper":{"head":PAPER,"export_run":45568373174,"export_artifact":15000204579},
      "tdwg_geometry":{"url":"https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson","sha256":gate["input_sha256"]["tdwg_geojson"]},
      "original_raw_inputs":"not redistributed; see original workflows and source checks",
      "source_file_git_blobs":SOURCE_FILES,
    },indent=2)+"\n").encode()
    new[PREFIX+"README.md"]=(
      "Internal-only S3.6/S5.5 source component. Original derived artifact members and source scripts preserved byte-for-byte. "
      "Run python replay.py --tdwg-geojson /path/to/pinned/level3.geojson from this component after installing numpy, pandas, scipy, statsmodels, shapely, pyproj, matplotlib. "
      "SOURCE_PINS.json gives the exact geometry URL and SHA256; no original raw annual occurrence, HOSTS or all-flora input files are redistributed. "
      "Replay covers S3.6 full and non-Fabaceae/Poaceae suites and range/effort diagnostics, and S5.5 all224 regional summaries and model variants; the without-monarch output remains source-verified archived evidence only. "
      "The original 223 archive members and their wrapper are unchanged; the old wrapper does not run this component. "
      "All new members are HOLD, including source code with internal provenance and the historical botanical-record effort table. "
      "Neither host taxon presence nor recorded butterfly occurrence identifies larval suitability, occupancy, establishment, or causation. "
      "No reviewer access, approval, or licence grant.\n").encode()
    assert not set(new)&set(old)
    allfiles=old|new
    for n,b in sorted(new.items()):
        # Deliberately conservative, not a licence determination. GBIF and HIST
        # cover the original butterfly archive and separate older plant indexing.
        sids={"AUTHOR","HOSTS","WCVP","TDWG","WPI","LEP"}
        if "/s55/" in n or "gbif224" in n or n.endswith(("replay.py","README.md","SOURCE_PINS.json")):sids.add("GBIF")
        if "/s36/" in n or "taxonomic_resolution" in n or "botanical_effort" in n or n.endswith(("replay.py","README.md","SOURCE_PINS.json")):sids.add("HIST")
        assert sids<=set(reg)
        role="code" if n.endswith((".py",".yml")) else "table" if n.endswith(".csv") else "metadata_or_summary"
        rows.append({"file":n,"sha256":sha(b),"bytes":len(b),"role":role,"source_ids":";".join(sorted(sids)),
          "mapping_basis":"conservative S3.6/S5.5 linked derived/source code dependency; inspect each source before release",
          "source_rights_states":";".join(sorted({reg[k]["source_state"] for k in sids})),
          "source_citations":" | ".join(reg[k]["citation_or_url"] for k in sorted(sids)),
          "source_terms":" | ".join(k+": "+reg[k]["source_terms"] for k in sorted(sids)),
          "delivery_decision":"HOLD_S36_S55_AUTHOR_SOURCE_AND_ANONYMITY_REVIEW",
          "proposed_route":"Internal only; original provenance retained; assess derivative rights and anonymized reviewer access separately",
          "author_approval_required":True,"authorized_for_distribution":False})
    total=len(allfiles)
    assert total==223+len(new)==len(rows)==len({r["file"] for r in rows})
    assert all(str(r["authorized_for_distribution"])=="False" for r in rows)
    a.out.mkdir(parents=True,exist_ok=False)
    dest=a.out/f"GEB_V05_{total}_INTERNAL_STAGING_NOT_FOR_RELEASE.zip"
    with zipfile.ZipFile(dest,"w") as z:
        for n,b in sorted(allfiles.items()):
            info=zipfile.ZipInfo(n,(2026,11,25,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,b)
    assert readzip(dest.read_bytes())==allfiles
    with zipfile.ZipFile(dest) as z:z.extractall(a.extracted)
    for r in rows:
        b=(a.extracted/r["file"]).read_bytes()
        assert sha(b)==r["sha256"] and len(b)==int(r["bytes"])
    with (a.out/f"EXACT_{total}_FILE_SOURCE_DELIVERY_DECISIONS.csv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for name,b in prior.items():
        if name not in ("EXACT_223_FILE_SOURCE_DELIVERY_DECISIONS.csv","GEB_V05_223_INTERNAL_STAGING_NOT_FOR_RELEASE.zip"):
            (a.out/("PREVIOUS_223_"+Path(name).name)).write_bytes(b)
    (a.out/"SOURCE_REGISTRY.json").write_bytes(find(prior,"SOURCE_REGISTRY.json"))
    summary={"status":"INTERNAL_INTEGRATION_BYTES_PASS_NOT_REVIEWER_RELEASE","files":total,
      "previous_members_byte_identical":223,"new_members":len(new),"source_families":len(reg),
      "file_roles":dict(collections.Counter(r["role"] for r in rows)),
      "delivery_decisions":dict(collections.Counter(r["delivery_decision"] for r in rows)),
      "authorized_files":0,"reviewer_access_created":False,"sha256":sha(dest.read_bytes()),"bytes":dest.stat().st_size,
      "paper_head":PAPER,"replay_pending":True}
    (a.out/"INTEGRATED_DELIVERY_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("INTEGRATED_S36_S55_BYTES_PASS",json.dumps(summary))
if __name__=="__main__":main()
