#!/usr/bin/env python3
"""Diagnose original HOSTS accepted butterfly edges against BSBI 2020 plant names.

This is source coverage/identification, NOT a butterfly host climate hindcast.
Original user-author HOSTS source stays on GitHub Actions transient filesystem.
Never export individual botanical hectad rows, UKBMS transect IDs/locations or
restricted user-author source edges.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,re,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from audit_bsbi_historical_hectad_source_20261029 import obtain
from audit_bsbi_nested_hectad_layout_20261029 import NESTED
from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK

AUTHOR_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
AUTHOR_INNER="dist/source/butterfly_resource_novelty_anonymous_OFFLINE_source_v01.zip"
AUTHOR_EDGES="inputs/original/insect_host_accepted.csv"
SCHEMA="chocho_BSBI_PA2020_original_HOSTS_27_species_name_matching_audit_20261029_v01"

def get_host_edges(path:Path):
    with zipfile.ZipFile(path) as outer:
        source=outer.read(AUTHOR_INNER)
    if hashlib.sha256(source).hexdigest()!=AUTHOR_SHA:
        raise ValueError("User-author original HOSTS spatial source SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(source)) as inner:
        if inner.testzip() is not None:raise ValueError("Corrupt original author HOSTS source")
        text=inner.read(AUTHOR_EDGES).decode("utf-8-sig")
    taxa=set(CROSSWALK)
    rows=csv.DictReader(io.StringIO(text))
    need={"insect_species","accepted_name","accepted_plant_name_id"}
    if not need<=set(rows.fieldnames or []):
        raise ValueError("Author HOSTS accepted plant/insect source schema drift")
    links=defaultdict(set)
    for row in rows:
        sp=row.get("insect_species","").strip()
        plant=row.get("accepted_name","").strip()
        if sp in taxa and plant:
            links[sp].add(plant)
    return links,{"original_SHA256":AUTHOR_SHA,
        "species_with_one_or_more_source_accepted_hosts":sum(bool(links.get(s)) for s in CROSSWALK),
        "source_edge_count":sum(len(x) for x in links.values())}

def canonical(n):
    """Explicit name-only exact matching; no unjustified botanical synonym inference."""
    return re.sub(r"\s+"," ",str(n or "").strip()).casefold()

def read_bsbi_names(source):
    with zipfile.ZipFile(io.BytesIO(source)) as z:
        name=z.read("PA2020_10by10kmSquareData/nameLookup.csv").decode("utf-8-sig")
        nested=z.read(NESTED)
    rows=list(csv.DictReader(io.StringIO(name)))
    if len(rows)!=3497:raise ValueError("Original 3497 BSBI name lookup changed")
    fields={"ddbidLookup","ddbidOrig","taxonName","canonical"}
    if not fields<=set(rows[0]):raise ValueError("BSBI original lookup/source fields drift")
    with zipfile.ZipFile(io.BytesIO(nested)) as z:
        files={Path(info.filename).stem for info in z.infolist() if info.filename.endswith(".csv")}
    if len(files)!=3497:raise ValueError("Source taxon to nested distribution file count mismatch")
    indices={}
    for row in rows:
        for field in ("taxonName","canonical"):
            nm=canonical(row.get(field,""))
            if nm:
                indices.setdefault(nm,set()).add(row["ddbidLookup"])
    if len(indices)<2800:raise ValueError("BSBI taxon source names cannot be resolved")
    # IDs are source-only table keys, not original geographical occurrences.
    pairs=[]
    for row in rows:
        key=row["ddbidLookup"].strip()
        orig=row["ddbidOrig"].strip()
        options=[key,orig,"2cd4p9h_"+key,"2cd4p9h_"+orig]
        exact=[x for x in options if x in files]
        pairs.append((row,exact))
    key_matches=sum(bool(q) for row,q in pairs)
    return indices,files,rows,{"original_taxon_lookup_rows":len(rows),
        "original_nested_taxon_csv_members":len(files),
        "filename_join_simple_key_matches":key_matches,
        "filename_join_key_examples":[{"type":"key","chars":len(str(rows[0]["ddbidLookup"]))},
                                      {"type":"original","chars":len(str(rows[0]["ddbidOrig"]))}],
        "name_lookup_taxa_names_unique":len(indices)}

def inspect(path:Path):
    links,source=get_host_edges(path)
    plantdata,url,errors=obtain()
    atlas,files,rows,atlas_qc=read_bsbi_names(plantdata)
    out={}
    for species in sorted(CROSSWALK):
        name_list=links.get(species,set())
        exact={nm for nm in name_list if canonical(nm) in atlas}
        out[species]={
            "original_author_HOSTS_resolved_named_host_taxa":len(name_list),
            "exact_matching_original_BSBI_taxon_names":len(exact),
            "unmatched_original_host_taxon_names":len(name_list)-len(exact),
            "possible_exact_plant_atlas_taxon_overlap_fraction":round(len(exact)/len(name_list),6) if name_list else None,
            "name_matching_is_not_accepted_local_host_use":True,
        }
    return {"schema":SCHEMA,"status":"SOURCE_VERIFIED_AUTHOR_HOSTS_TO_BSBI_TAXON_NAME_COVERAGE",
        "original_user_author_HOSTS":source,
        "original_BSBI":{**atlas_qc,"original_source_archive_sha256":hashlib.sha256(plantdata).hexdigest(),
                         "original_zenodo":"https://zenodo.org/records/10998003"},
        "species_coverage":out,
        "source_taxa_27":len(out),
        "source_accepted_host_plant_atlas_taxa_matching":sum(x["exact_matching_original_BSBI_taxon_names"] for x in out.values()),
        "no_original_AUTHOR_HOSTS_edge_or_BSBI_grid_rows_exported":True,
        "no_historical_site_host_geography_inferred":True,
        "known_source_gap_flag":"Original HOSTS Urtica dioica missing for Aglais urticae, Rumex hydrolapathum missing for Lycaena dispar; independent UKBMS sources document both.",
        "historical_cultivation_or_maintenance_measured":False}

def run(source,output):
    try:
        result=inspect(source)
    except (ValueError,OSError,zipfile.BadZipFile) as exc:
        result={"schema":SCHEMA,"status":"SOURCE_IDENTIFICATION_FAILED_CLOSED",
                "blocker":str(exc)[:320],"source_rows_exported":False}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"source_taxa_27":result.get("source_taxa_27"),
        "original_BSBI":result.get("original_BSBI"),
        "host_name_matches":result.get("source_accepted_host_plant_atlas_taxa_matching"),
        "blocker":result.get("blocker")},ensure_ascii=False))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--author-artifact-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    x=run(a.author_artifact_zip,a.output)
    if x["status"]!="SOURCE_VERIFIED_AUTHOR_HOSTS_TO_BSBI_TAXON_NAME_COVERAGE":
        raise SystemExit(2)
