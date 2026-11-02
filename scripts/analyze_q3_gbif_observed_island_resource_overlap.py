#!/usr/bin/env python3
"""Q3: independent GBIF occurrence evidence at exact added butterfly x WGSRPD3 island cells.

GBIF DOI download has a fixed 31-observed-species validation panel, NOT all 239
butterflies. Non-observation is not absence, uptake, establishment or breeding.
The WCVP contemporary layer is not a dated pre-observation prediction.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, re, zipfile
from collections import Counter, defaultdict
from pathlib import Path
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

KEY="0008693-260928105237408"
DOI="10.15468/dl.pp5nc9"
EXPECTED_RECORDS=53144
EXPECTED_API=53434
NAME_RE=re.compile(r"^\s*([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b")
# GBIF/Catalogue of Life: Historis acheronta (Fabricius, 1775) is a synonym
# of the accepted Coea acheronta. Source: https://www.gbif.org/taxon/PTQJ2
VERIFIED_SYNONYMS={"Coea acheronta":"Historis acheronta"}

def load_csv(path):
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h))

def dump_csv(path, rows):
    if not rows:raise RuntimeError("No rows to write")
    with path.open("w", newline="", encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def parse_archive(path):
    with zipfile.ZipFile(path) as zf:
        names=zf.namelist()
        files=[x for x in names if x.lower().endswith((".csv",".txt")) and
               not x.lower().endswith(("readme.txt","rights.txt"))]
        priorities=sorted(files,key=lambda x:(0 if "occurrence" in x.lower() else 1,-zf.getinfo(x).file_size))
        if not priorities:raise RuntimeError(f"No occurrence CSV/TXT inside GBIF DOI: {names[:30]}")
        for candidate in priorities:
            data=zf.read(candidate)
            head=data[:2048].decode("utf-8-sig",errors="replace").splitlines()[0]
            sep="\t" if "\t" in head else ","
            reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig")),delimiter=sep)
            fields=list(reader.fieldnames or [])
            canonical={x.strip().lower():x for x in fields}
            if all(x in canonical for x in ("decimallatitude","decimallongitude")):
                return list(reader),canonical,{"archive_file":candidate,"csv_delimiter":"tab" if sep=="\t" else "comma",
                                              "columns":fields,"zip_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
    raise RuntimeError("No GBIF occurrence file with decimalLatitude/decimalLongitude")

def get(row,keys,field):
    key=keys.get(field.lower())
    return "" if key is None else str(row.get(key) or "").strip()

def canon(s):
    match=NAME_RE.match(s)
    return None if not match else f"{match.group(1)} {match.group(2)}"

def geometry_lookup(path):
    data=json.loads(path.read_text(encoding="utf-8"))
    geoms=[];code=[]
    for feat in data["features"]:
        reg=str((feat.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        geom=feat.get("geometry")
        if geom and reg:
            poly=shape(geom)
            if not poly.is_empty:
                geoms.append(poly)
                code.append(reg)
    if len(set(code))!=369:raise RuntimeError(f"Wrong WGSRPD geometry codes: {len(set(code))}")
    return STRtree(geoms),geoms,code

def main():
    p=argparse.ArgumentParser()
    for name in ("gbif-download-zip","validation-species","q2-edges","q1-regions",
                 "q1-summary","wgsrpd3-geojson","out-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    old=json.loads(a.q1_summary.read_text(encoding="utf-8"))
    if old["panel_total"]["new_total"]!=14553 or old["panel_total"]["n_panel_mismatches"]!=0:
        raise RuntimeError("Q1 fixed opportunity reconstruction mismatch")
    panel=load_csv(a.validation_species)
    if len(panel)!=32 or len({r["species"] for r in panel})!=32:
        raise RuntimeError(f"Expected 32 fixed GEB validation candidates, found {len(panel)}")
    observed_species={r["species"] for r in panel if int(r["occurrence_records"])>0}
    if len(observed_species)!=31:raise RuntimeError("Expected 31 GBIF sampled species")
    all_species={r["species"] for r in panel}
    islands={r["area_code_l3"]:r for r in load_csv(a.q1_regions) if r["class"]=="island"}
    if len(islands)!=125:raise RuntimeError("Q1 strict island classification drift")
    added=defaultdict(set)
    for r in load_csv(a.q2_edges):
        added[r["species"]].add(r["island_code"])
    if sum(len(v) for v in added.values())!=5674:raise RuntimeError("Frozen 5674 island additions changed")
    records,columns,archive_info=parse_archive(a.gbif_download_zip)
    if len(records)!=EXPECTED_RECORDS:
        raise RuntimeError(f"GBIF DOI returned {len(records)} rows instead of archived {EXPECTED_RECORDS}; do not silently replace")
    geo,geoms,codes=geometry_lookup(a.wgsrpd3_geojson)
    outcome=defaultdict(lambda:{"n_records":0,"years":set(),"source_gbif_ids":set()})
    counts=Counter()
    missing_taxon=Counter()
    seen_keys=set()
    for rec in records:
        idv=get(rec,columns,"gbifID") or get(rec,columns,"gbif_key") or get(rec,columns,"key")
        if not idv:raise RuntimeError("GBIF DOI occurrence lacks unique ID")
        if idv in seen_keys:raise RuntimeError(f"Duplicate original occurrence GBIF id {idv}")
        seen_keys.add(idv)
        species=(canon(get(rec,columns,"species")) or
                 canon(get(rec,columns,"acceptedScientificName")) or
                 canon(get(rec,columns,"scientificName")))
        species=VERIFIED_SYNONYMS.get(species,species)
        if species not in all_species:
            counts["not_in_frozen_32_species"]+=1
            if species:missing_taxon[species]+=1
            continue
        if species not in observed_species:raise RuntimeError(f"GBIF zero-record species unexpectedly present: {species}")
        counts["matched_species_occurrences"]+=1
        lat,lon=get(rec,columns,"decimalLatitude"),get(rec,columns,"decimalLongitude")
        try:lat=float(lat);lon=float(lon)
        except (TypeError,ValueError):counts["missing_coords"]+=1;continue
        if not -90<=lat<=90 or not -180<=lon<=180:
            counts["invalid_coords"]+=1;continue
        pt=Point(lon,lat)
        matches={codes[int(i)] for i in geo.query(pt) if geoms[int(i)].covers(pt)}
        if not matches:
            counts["unmapped_geo"]+=1
            continue
        if len(matches)>1:
            counts["ambiguous_wgsrpd3"]+=1
            continue
        region=next(iter(matches))
        if region not in islands:
            counts["mainland_or_mixed_occurrences"]+=1
            continue
        try:year=int(float(get(rec,columns,"year")))
        except (ValueError,TypeError):year=None
        row=outcome[(species,region)]
        row["n_records"]+=1
        row["source_gbif_ids"].add(idv)
        if year and 2010<=year<=2026:row["years"].add(year)
        counts["island_occurrences"]+=1
    if len(seen_keys)!=EXPECTED_RECORDS:raise RuntimeError("GBIF DOI ID uniqueness mismatch")
    island_out=[]
    for region in sorted(islands):
        obs={sp:rec for (sp,code),rec in outcome.items() if code==region}
        new={sp for sp in observed_species if region in added[sp]}
        matched=set(obs)&new
        island_out.append({
          "area_code_l3":region,"region_name":islands[region]["region_name"],
          "validation_species_n":len(observed_species),
          "new_opportunity_cells_within_validation_species":len(new),
          "new_opportunity_cells_with_gbif_records":len(matched),
          "gbif_observed_species_cells_total":len(obs),
          "gbif_records_all_validation_species":sum(v["n_records"] for v in obs.values()),
          "gbif_records_at_added_opportunity":sum(obs[s]["n_records"] for s in matched),
          "years_with_any_observed_added_cell":
              len(set().union(*(obs[s]["years"] for s in matched))) if matched else 0
        })
    species_out=[]
    for r in panel:
        species=r["species"]
        observed_cells={region for sp,region in outcome if sp==species}
        added_cells=added[species]
        overlap=observed_cells&added_cells
        species_out.append({
         "species":species,"validation_sample_gbif_records":int(r["occurrence_records"]),
         "archival_gbif_records_islands":sum(v["n_records"] for (sp,code),v in outcome.items() if sp==species),
         "added_island_opportunity_cells":len(added_cells),
         "observed_island_units":len(observed_cells),
         "observed_at_new_opportunity_island_units":len(overlap),
         "added_opportunity_observed_fraction_not_realization_rate":
             len(overlap)/len(added_cells) if added_cells else "",
         "other_observed_island_units":len(observed_cells-added_cells),
        })
    edges=[]
    for sp,region in sorted(outcome):
        v=outcome[(sp,region)]
        edges.append({
          "species":sp,"area_code_l3":region,
          "observed_butterfly":1,
          "introduced_host_resource_added":int(region in added[sp]),
          "observed_records":v["n_records"],
          "first_indexed_year_in_DOI_subset":min(v["years"]) if v["years"] else "",
          "last_indexed_year_in_DOI_subset":max(v["years"]) if v["years"] else "",
          "indexed_2010_2017":int(any(2010<=y<=2017 for y in v["years"])),
          "indexed_2018_2026":int(any(2018<=y<=2026 for y in v["years"])),
        })
    dump_csv(a.out_dir/"island_gbif_observation_summary.csv",island_out)
    dump_csv(a.out_dir/"species_observation_summary.csv",species_out)
    if edges:dump_csv(a.out_dir/"observed_species_island_cells.csv",edges)
    N=len(observed_species)
    added_subset=sum(len(added[x]) for x in observed_species)
    matched=sum(r["observed_at_new_opportunity_island_units"] for r in species_out)
    report={
     "schema":"chocho_q3_gbif_doi_31species_island_opportunity_observation_v0.1",
     "status":"INDEPENDENT_GBIF_DOI_RECORD_OVERLAP_DESCRIPTIVE_ONLY",
     "gbif":{"doi":DOI,"download_key":KEY,"archive_total_records_expected":EXPECTED_RECORDS,
             "verified_taxonomic_synonyms":VERIFIED_SYNONYMS,
             "original_api_records":EXPECTED_API,"not_returned_at_doi_rebuild":EXPECTED_API-EXPECTED_RECORDS,
             **archive_info},
     "validation_panel":{
         "original_GEB_species":239,"original_added_island_resource_cells":5674,
         "gbif_validation_species":len(panel),"gbif_species_with_original_records":N,
         "sampled_species_added_island_opportunity_cells":added_subset,
         "recorded_butterfly_in_new_opportunity_cells":matched,
         "fraction_sampled_added_opportunity_cells_with_any_gbif_record":
             matched/added_subset if added_subset else None,
         "sampled_islands_with_any_butterfly_record":sum(x["gbif_observed_species_cells_total"]>0 for x in island_out),
         "sampled_islands_with_record_at_added_opportunity":sum(x["new_opportunity_cells_with_gbif_records"]>0 for x in island_out),
         "island_butterfly_species_cells_with_any_gbif_record":len(edges),
         "total_island_gbif_records":counts["island_occurrences"],
     },
     "audit_counts":dict(counts),
     "unmatched_species_top10":missing_taxon.most_common(10),
     "inference_limits":[
        "Archive is an independent historic GBIF *sample* for 31 observed butterfly species only, NOT the complete 239 species or full GBIF inventory.",
        "A species x island cell lacking GBIF records is unobserved in this downloaded subset and cannot be treated as butterfly absence or unrealized opportunity.",
        "Every observed overlap is an index record, NOT evidence of breeding, larval feeding, establishment, or colonial history.",
        "WCVP contemporary introduced-host ranges are not time-indexed; no prospective prediction or causal arrival-before-host ordering is tested.",
        "Island classification is a coarse WGSRPD3 region, which may comprise archipelagos, not individual local ports or island populations.",
        "The GBIF DOI reconstructed fewer records than original API; not all 53,434 original IDs were returned.",
        "The validation species were selected for a prior GEB ecological panel and do not represent the unobserved 208 GEB species."
     ]
    }
    (a.out_dir/"q3_archival_gbif_overlap_audit.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("audit_counts","unmatched_species_top10","inference_limits")},indent=2,ensure_ascii=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
