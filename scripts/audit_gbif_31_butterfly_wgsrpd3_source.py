#!/usr/bin/env python3
"""Independent GBIF observation-source feasibility at the SAME TDWG regions.

Uses a previously minted, citable GBIF download DOI for 31 butterflies:
https://doi.org/10.15468/dl.pp5nc9, key 0008693-260928105237408.
Archival GBIF download returned 53,144 of originally requested 53,434 IDs.
This is *not* the 239-butterfly global potential-resource cohort.

The purpose is to audit factual, independently observed butterfly assemblage
coverage BEFORE testing any host-predictor association. No GBIF non-record
means true absence. Raw GPS points, GBIF IDs, and observer/provider identity
are NOT committed back to GitHub.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

from shapely.geometry import Point,shape
from shapely.strtree import STRtree

SOURCE_KEY="0008693-260928105237408"
SOURCE_DOI="10.15468/dl.pp5nc9"
EXPECTED_RETURNED=53144
SOURCE=Path("_private_observations/gbif_31_2026.zip")
GEOM=Path("_external/wgsrpd/geojson/level3.geojson")
REGIONS=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/lenzner_wgsrpd3_colonial_area_crosswalk.csv")
PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/gbif_31_observed_region_source_gate.json")
OUT_REG=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/gbif_31_observed_region_counts.csv")


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def find_occurrence_member(z):
    found=[]
    diagnostics=[]
    for x in z.infolist():
        if x.is_dir() or x.file_size<1000 or not x.filename.lower().endswith((".txt",".tsv",".csv")):
            continue
        with z.open(x.filename) as inp:
            head=inp.readline().decode("utf-8-sig",errors="replace").rstrip("\r\n")
        delimiter="\t" if head.count("\t")>=head.count(",") else ","
        columns=[v.strip().strip('"') for v in head.split(delimiter)]
        diagnostics.append({"filename":x.filename,"bytes":x.file_size,
                            "fields":len(columns),"delimiter":"tab" if delimiter=="\t" else "comma"})
        if ({"decimallatitude","decimallongitude"} <= {s.casefold() for s in columns}
                and ({"gbifid","species","scientificname"} & {s.casefold() for s in columns})):
            found.append((x.filename,delimiter,columns))
    if len(found)!=1:
        raise RuntimeError(f"Expected one occurrence table; found {len(found)} candidates: {diagnostics[:15]}")
    return found[0], diagnostics


def field(row, names):
    for name in names:
        value=row.get(name)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def normalize_species(row):
    # GBIF 'species' is normally a canonical binomial. Rejected ambiguous
    # scientificName authors/infra-specific strings are not silently truncated.
    name=field(row,("species","acceptedscientificname","scientificname"))
    if not name:
        genus=field(row,("genus",))
        epithet=field(row,("specificepithet",))
        if genus and epithet:
            name=f"{genus} {epithet}"
    values=name.split()
    return (" ".join(values[:2]).casefold() if len(values)>=2 else None)


def geometry():
    o=json.loads(GEOM.read_text(encoding="utf8"))
    shapes=[]
    codes=[]
    for f in o["features"]:
        name=str(f["properties"].get("LEVEL3_COD","")).strip()
        if not name:
            continue
        shapes.append(shape(f["geometry"]))
        codes.append(name)
    if len(shapes)!=369 or len(set(codes))!=369:
        raise RuntimeError("Expected exact 369 official Level3 polygons")
    return STRtree(shapes),shapes,codes


def audit(source=SOURCE):
    if not source.exists():
        raise RuntimeError("GBIF downloaded original ZIP is missing")
    recordset=read_csv(PANEL)
    panel={r["species"].strip().casefold() for r in recordset}
    if len(recordset)!=239 or len(panel)!=239:
        raise RuntimeError("Fixed HOSTS butterfly panel has drifted")
    historical=read_csv(REGIONS)
    eligible={r["wgsrpd3_code"] for r in historical
              if float(r["glonaf_covered_area_fraction_unclipped"])>=.8
              and float(r["glonaf_covered_area_fraction_unclipped"])<=1.05}
    if len(eligible)!=284:
        raise RuntimeError(f"Original colonial matched region set changed: {len(eligible)}")
    tree,shapes,codes=geometry()
    totals=Counter()
    per_region=defaultdict(lambda:{"records":0,"species":set(),"dates":set(),"providers":set(),
                                   "provider_days":set(),"records_uncertainty_<=1000m":0})
    observed_species=set()
    seen=set()
    with zipfile.ZipFile(source) as z:
        if z.testzip() is not None:
            raise RuntimeError("GBIF archive contains corrupt ZIP members")
        (member,delimiter,columns),members=find_occurrence_member(z)
        names={n.casefold():n for n in columns}
        required={"decimallatitude","decimallongitude"}
        if not required <= set(names):
            raise RuntimeError("Unrecognized decimal coordinate column schema")
        with z.open(member) as src:
            reader=csv.DictReader(io.TextIOWrapper(src,encoding="utf-8-sig",newline=""),delimiter=delimiter)
            for r in reader:
                totals["source_occurrence_rows"]+=1
                row={k.casefold():v for k,v in r.items() if k is not None}
                gbifkey=field(row,("gbifid","gbifkey","key"))
                if gbifkey:
                    if gbifkey in seen:
                        totals["duplicate_gbif_key"]+=1
                    seen.add(gbifkey)
                name=normalize_species(row)
                if not name:
                    totals["no_taxonomically_resolvable_species"]+=1
                    continue
                if name not in panel:
                    totals["not_in_frozen_239_panel"]+=1
                    continue
                observed_species.add(name)
                totals["focal_239_occurrence_rows"]+=1
                try:
                    lat=float(field(row,("decimallatitude",)))
                    lon=float(field(row,("decimallongitude",)))
                except ValueError:
                    totals["missing_coordinates"]+=1
                    continue
                if not (math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180):
                    totals["bad_coordinates"]+=1
                    continue
                point=Point(lon,lat)
                covered=[codes[int(k)] for k in tree.query(point)
                         if shapes[int(k)].covers(point)]
                if not covered:
                    totals["not_in_official_wgsrpd3_polygon"]+=1
                    continue
                if len(covered)>1:
                    totals["multiple_polygons_cover_point"]+=1
                code=sorted(covered)[0]
                totals["matched_official_region_occurrences"]+=1
                if code not in eligible:
                    totals["outside_matched_284_regions"]+=1
                    continue
                year=field(row,("year",))
                date=field(row,("eventdate",))
                dataset=field(row,("datasetkey",))
                region=per_region[code]
                region["records"]+=1
                region["species"].add(name)
                if date and len(date)>=10 and date[4]=="-":
                    region["dates"].add(date[:10])
                    if dataset:
                        region["provider_days"].add((dataset,date[:10]))
                if dataset:
                    region["providers"].add(dataset)
                uncertainty=field(row,("coordinateuncertaintyinmeters",))
                try:
                    if uncertainty and 0<=float(uncertainty)<=1000:
                        region["records_uncertainty_<=1000m"]+=1
                except ValueError:
                    pass
    if totals["source_occurrence_rows"] != EXPECTED_RETURNED:
        raise RuntimeError("DOI archival GBIF returned count changed: "+
                           repr(dict(totals)))
    if totals["duplicate_gbif_key"]>0:
        raise RuntimeError("Duplicate GBIF keys in DOI archive")
    region_counts=[]
    for code in sorted(eligible):
        r=per_region.get(code)
        region_counts.append({
            "wgsrpd3_code":code,
            "indexed_butterfly_species":len(r["species"]) if r else 0,
            "gbif_indexed_records":r["records"] if r else 0,
            "gbif_observed_days":len(r["dates"]) if r else 0,
            "gbif_dataset_keys":len(r["providers"]) if r else 0,
            "gbif_provider_days":len(r["provider_days"]) if r else 0,
            "gbif_records_with_uncertainty_le1000m":r["records_uncertainty_<=1000m"] if r else 0,
        })
    with_data=[r for r in region_counts if r["gbif_indexed_records"]>0]
    has_5=[r for r in region_counts if r["indexed_butterfly_species"]>=5
           and r["gbif_provider_days"]>=5]
    count_species=Counter()
    for r in per_region.values():
        count_species.update(r["species"])
    result={
        "status":"INDEPENDENT_GBIF_BUTTERFLY_SOURCE_COVERAGE_NOT_HOST_REALIZATION_EFFECT",
        "source":{
            "doi":SOURCE_DOI,"gbif_download_key":SOURCE_KEY,
            "archive_sha256":sha256(source),
            "public_download_url":"https://api.gbif.org/v1/occurrence/download/request/"+SOURCE_KEY+".zip",
            "requested_gbif_ids":53434,"returned_archival_download_rows":EXPECTED_RETURNED,
            "missing_from_exact_ID_archival_export":290,
            "frozen_butterfly_panel_count":239,
            "source_original_selected_species":31,
            "geography":"official TDWG 369 level3, GloNAF 284-region matched subset",
        },
        "archive_member_diagnostics":members,
        "column_names":columns,
        "source_row_quality":dict(totals),
        "coverage":{
            "source_taxa_matching_frozen_239":len(observed_species),
            "matched_regions_with_any_observed_butterfly":len(with_data),
            "matched_regions_with_5plus_species_and_5plus_provider_days":len(has_5),
            "all_matched_regions":len(eligible),
            "species_recorded_in_5plus_matched_regions":sum(x>=5 for x in count_species.values()),
            "species_recorded_in_10plus_matched_regions":sum(x>=10 for x in count_species.values()),
            "median_observed_focal_species_per_nonempty_region":(
                sorted(r["indexed_butterfly_species"] for r in with_data)[len(with_data)//2]
                if with_data else None),
        },
        "pilot_gate":"REGION_BETA_FEASIBLE_FOR_LIMITED_31_SPECIES" if
          (len(has_5)>=30 and len(observed_species)>=20) else
          "INSUFFICIENT_REGIONAL_SUPPORT_FOR_INDEPENDENT_31_SPECIES_BETA",
        "limits":[
            "Only 31 selected butterflies in the frozen GBIF DOI, not the 239-focal or 225-species global island forecast.",
            "GBIF first-index years and non-detections are not true arrival years or species absences.",
            "Same-region presence reflects heterogeneous observer and provider-day reporting; no absence correction yet.",
            "GloNAF and WCVP host source are current retrospective distributions, not true historical introductions.",
            "A good pilot feasibility result would authorize a separate observational model, not prove plant-mediated butterfly colonization.",
        ],
    }
    return result,region_counts


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,default=SOURCE)
    args=p.parse_args()
    output,regional=audit(args.zip)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with OUT_REG.open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(regional[0]))
        writer.writeheader()
        writer.writerows(regional)
    print(json.dumps({k:v for k,v in output.items()
                      if k not in ("archive_member_diagnostics","column_names")},
                     indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
