#!/usr/bin/env python3
"""Stream global GBIF butterfly occurrences into species × WGSRPD3 × year.

Source: GBIF DOI 10.15468/dl.umqmsf, archive key
0013264-261023074918325; 14,637,152 records, ~1.42GB ZIP.
Original candidate set: frozen 225 butterflies in chocho's historical
forecast_all_225_disposition.csv. The 239 potential-resource study panel is
a DIFFERENT panel: intersections are audited, never assumed to be 224/239.

This program makes source-derived, non-coordinate aggregates ONLY:
 - GBIF-indexed positive species/region/year records and dataset × date counts
 - regional recording opportunity by year, among the same 225 candidate species
 - species/region positive multi-year summaries
 - identity/geometry/count/quality source receipt.
A missing positive row is NOT biological absence. No model fitting or inferred
colonial or shipping effects. The archived ZIP and individual source GBIF
IDs, provider IDs, dates, and coordinates must stay outside published GitHub.

Requires Python >=3.11, pandas >=2, numpy >=1.26, shapely >=2.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import date
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from shapely.geometry import shape
from shapely import points as shapely_points
from shapely.strtree import STRtree

DOI="10.15468/dl.umqmsf"
GBIF_KEY="0013264-261023074918325"
EXPECTED_RECORDS=14637152
EXPECTED_ZIP_BYTES=1415705617
START_YEAR=2001
END_YEAR=2024
SOURCE_ZIP=Path("_private_observations/gbif_global_224_2026.zip")
TDWG=Path("_external/wgsrpd/geojson/level3.geojson")
LOOKUP=Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json")
CANDIDATES=Path("provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_all_225_disposition.csv")
RESOURCE_PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
OUTPUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/gbif_global224_annual_v1")
LOCK=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/gbif_224_annual_wgsrpd3_source_lock_v1.json")

REQUIRED={"decimallatitude","decimallongitude","species","year"}
OPTIONAL={"specieskey","taxonkey","eventdate","datasetkey","basisofrecord",
          "coordinateuncertaintyinmeters","gbifid","month","day",
          "occurrencestatus","issue"}


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<21),b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def candidate_panel(path=CANDIDATES,panel=RESOURCE_PANEL):
    rows=read_csv(path)
    if len(rows)!=225:
        raise RuntimeError(f"Frozen candidate list size changed: {len(rows)}")
    by_key={}
    by_name={}
    original_disposition=Counter()
    for r in rows:
        name=r["species"].strip()
        key=str(r["taxon_key"]).strip()
        if not key or not name or key in by_key or name.casefold() in by_name:
            raise RuntimeError("Missing or duplicate GBIF candidate species identity")
        by_key[key]=name
        by_name[name.casefold()]=name
        original_disposition[r["stage_disposition"]]+=1
    resource={r["species"].strip().casefold() for r in read_csv(panel)}
    if len(resource)!=239:
        raise RuntimeError("GEB resource 239 species panel changed")
    return by_key,by_name,resource,dict(original_disposition)


def official_geography(path=TDWG,lookup=LOOKUP):
    frozen=json.loads(Path(lookup).read_text(encoding="utf8"))
    if frozen.get("source_commit")!="52da7828aba9d461dd133c27b3bd7a4407161f54":
        raise RuntimeError("Official TDWG source SHA changed")
    expected=set(frozen["units"])
    if len(expected)!=369:
        raise RuntimeError("Frozen WGSRPD3 geography is not 369 units")
    raw=json.loads(Path(path).read_text(encoding="utf8"))
    items=[]
    for f in raw["features"]:
        c=str(f["properties"].get("LEVEL3_COD") or "").strip()
        if c:
            items.append((c,shape(f["geometry"])))
    items.sort(key=lambda x:x[0])
    codes=[x[0] for x in items]
    polygons=[x[1] for x in items]
    if len(codes)!=369 or set(codes)!=expected or len(set(codes))!=369:
        raise RuntimeError("Official TDWG polygon codes differ from source lock")
    if not all(not g.is_empty and g.is_valid for g in polygons):
        raise RuntimeError("Official polygon is empty/invalid")
    return STRtree(polygons),codes


def assign_xy_to_369(lons,lats,tree):
    """Vectorized TDWG lookup with alphabetically deterministic overlap ties.

    Shapely 2 STRtree.query returns [point_index, polygon_index] pairs.
    Official polygons are stored sorted by region code. For overlapping
    polygons or shared borders, assign lowest sorted code, count ambiguity.
    A missing polygon hit is -1, NEVER a biological absence.
    """
    lons=np.asarray(lons,float)
    lats=np.asarray(lats,float)
    if lons.shape!=lats.shape:
        raise ValueError("Coordinate dimensions differ")
    out=np.full(len(lons),-1,np.int32)
    valid=(np.isfinite(lons)&np.isfinite(lats)&
           (lons>=-180)&(lons<=180)&(lats>=-90)&(lats<=90))
    if not valid.any():
        return out,0
    original=np.flatnonzero(valid)
    pts=shapely_points(lons[valid],lats[valid])
    matched=tree.query(pts,predicate="intersects")
    if matched.shape[1]==0:
        return out,0
    # Second lexicographic key is polygon index, already alphabetical code.
    order=np.lexsort((matched[1],matched[0]))
    pairs=matched[:,order]
    unique,first,cnt=np.unique(pairs[0],return_index=True,return_counts=True)
    out[original[unique]]=pairs[1,first].astype(np.int32)
    ambiguous=int(np.sum(cnt>1))
    return out,ambiguous


def archive_table(z):
    found=[]
    found_info=[]
    for f in z.infolist():
        if f.is_dir() or not f.filename.lower().endswith((".csv",".tsv",".txt")):
            continue
        if f.file_size<1000:
            continue
        with z.open(f.filename) as src:
            header=src.readline().decode("utf-8-sig",errors="replace").rstrip("\r\n")
        sep="\t" if header.count("\t")>header.count(",") else ","
        columns=next(csv.reader([header],delimiter=sep))
        norm={name.strip().casefold():name for name in columns}
        found_info.append({"file":f.filename,"bytes_uncompressed":f.file_size,
                           "fields":len(columns),"delimiter":"tab" if sep=="\t" else "comma"})
        if REQUIRED <= set(norm) and ("specieskey" in norm or "taxonkey" in norm):
            found.append((f.filename,sep,norm))
    if len(found)!=1:
        raise RuntimeError("Expected one uniquely identified GBIF coordinate occurrence table: "+
                           repr(found_info[:10]))
    return found[0],found_info


def accepted_iso_date(value,source_year,cache):
    """Accept only explicit ISO calendar day matched to supplied year."""
    x=str(value).strip()[:10]
    if x in cache:
        yy=cache[x]
    else:
        yy=None
        if len(x)==10 and x[4]=="-" and x[7]=="-":
            try:
                yy=date.fromisoformat(x).year
            except ValueError:
                pass
        if len(cache)<50000:
            cache[x]=yy
    return x if yy==source_year else None


def process_chunks(file_obj,separator,source_columns,tree,codes,
                   key_to_name,name_to_name,chunk_size,
                   max_source_rows=None):
    original=[v for k,v in source_columns.items()
              if k in REQUIRED|OPTIONAL]
    if len(set(original))!=len(original):
        raise RuntimeError("Duplicate GBIF source column mapping")
    reader=pd.read_csv(file_obj,sep=separator,usecols=original,dtype=str,
                       chunksize=chunk_size,keep_default_na=False,
                       encoding="utf-8-sig",on_bad_lines="error",
                       low_memory=False)
    n_source=Counter()
    record_cells=Counter()
    precise_cells=Counter()
    annual_species_region=defaultdict(set)
    years_by_s_r=defaultdict(set)
    provider_days=defaultdict(set)
    region_days=defaultdict(set)
    providers=defaultdict(set)
    region_providers=defaultdict(set)
    # DatasetKey/date pair is only an ephemeral integer mapping. Neither ID
    # nor event date leaves this process; aggregated cardinalities do.
    provider_date_index={}
    date_cache={}
    original_taxa_seen=set()
    assigned_species=set()
    orig_years=Counter()
    basis_counts=Counter()
    chunks=0
    for df in reader:
        chunks+=1
        df.columns=[str(c).strip().casefold() for c in df.columns]
        length=len(df)
        if max_source_rows is not None and n_source["rows"]+length>max_source_rows:
            df=df.iloc[:max_source_rows-n_source["rows"]].copy()
            length=len(df)
        if length<=0:
            break
        n_source["rows"]+=length

        keycol="specieskey" if "specieskey" in df else "taxonkey"
        if "specieskey" not in df and "taxonkey" not in df:
            raise RuntimeError("GBIF speciesKey/taxonKey missing")
        bykey=df[keycol].map(key_to_name).fillna("")
        byname=df["species"].str.strip().str.casefold().map(name_to_name).fillna("")
        disagreement=(bykey.ne("")&byname.ne("")&bykey.ne(byname))
        taxon=bykey.where(bykey.ne(""),byname)
        n_source["candidate_key_match"]+=int(bykey.ne("").sum())
        n_source["candidate_name_only_match"]+=int((bykey.eq("")&byname.ne("")).sum())
        n_source["candidate_key_and_name_disagree"]+=int(disagreement.sum())
        n_source["not_a_frozen_candidate"]+=int(taxon.eq("").sum())
        taxon=taxon.where(~disagreement,"")
        years=pd.to_numeric(df["year"],errors="coerce").to_numpy(float)
        lat=pd.to_numeric(df["decimallatitude"],errors="coerce").to_numpy(float)
        lon=pd.to_numeric(df["decimallongitude"],errors="coerce").to_numpy(float)
        finite_coords=(np.isfinite(lat)&np.isfinite(lon)&
                       (lat>=-90)&(lat<=90)&(lon>=-180)&(lon<=180))
        good_year=(np.isfinite(years)&(years>=START_YEAR)&(years<=END_YEAR)&
                   (years==np.floor(years)))
        n_source["no_valid_wgs84_coordinate"]+=int((~finite_coords).sum())
        n_source["out_of_2001_2024_year_scope"]+=int((~good_year).sum())
        valid_taxa=taxon.ne("").to_numpy(bool)
        n_source["candidate_indexed_with_valid_year_and_point"]+=int(
            (valid_taxa&finite_coords&good_year).sum())
        if "basisofrecord" in df:
            basis_counts.update(df["basisofrecord"].value_counts().to_dict())
        if "occurrencestatus" in df:
            nonpresent=(df["occurrencestatus"].str.upper().ne("PRESENT") &
                        df["occurrencestatus"].ne(""))
            n_source["source_occurrence_status_nonpresent"]+=int(nonpresent.sum())
        else:
            nonpresent=pd.Series(np.zeros(length,bool))
        allowed=(valid_taxa&finite_coords&good_year&~nonpresent.to_numpy(bool))
        if not allowed.any():
            if max_source_rows is not None and n_source["rows"]>=max_source_rows:break
            continue
        local=np.flatnonzero(allowed)
        local_codes,overlap=assign_xy_to_369(lon[local],lat[local],tree)
        n_source["ambiguous_multiple_level3_polygon_hits"]+=overlap
        assigned=(local_codes>=0)
        n_source["candidate_valid_rows_outside_official_369_regions"]+=int((~assigned).sum())
        selected=local[assigned]
        regions=local_codes[assigned]
        n_source["candidate_records_mapped_to_369"]+=len(selected)
        if not len(selected):
            if max_source_rows is not None and n_source["rows"]>=max_source_rows:break
            continue
        names=taxon.to_numpy(object)[selected]
        yrs=years[selected].astype(np.int16)
        dset=(df["datasetkey"].to_numpy(object)[selected]
              if "datasetkey" in df else np.array([""]*len(selected),object))
        dates=(df["eventdate"].to_numpy(object)[selected]
               if "eventdate" in df else np.array([""]*len(selected),object))
        unc=(pd.to_numeric(df["coordinateuncertaintyinmeters"],errors="coerce").to_numpy(float)[selected]
             if "coordinateuncertaintyinmeters" in df
             else np.full(len(selected),np.nan))
        highprecision=(np.isfinite(unc)&(unc>=0)&(unc<=1000))
        n_source["candidate_mapped_records_reported_uncertainty_le1000m"]+=int(highprecision.sum())
        n_source["candidate_mapped_records_without_uncertainty"]+=int((~np.isfinite(unc)).sum())
        for k,(species,regionidx,year) in enumerate(zip(names,regions,yrs)):
            species=str(species)
            yy=int(year)
            area=codes[int(regionidx)]
            original_taxa_seen.add(species)
            assigned_species.add(species)
            orig_years[yy]+=1
            cell=(species,area,yy)
            record_cells[cell]+=1
            if bool(highprecision[k]):
                precise_cells[cell]+=1
            years_by_s_r[(species,area)].add(yy)
            annual_species_region[(area,yy)].add(species)
            dataset=str(dset[k]).strip()
            if dataset:
                providers[cell].add(dataset)
                region_providers[(area,yy)].add(dataset)
            dt=accepted_iso_date(dates[k],yy,date_cache)
            if dataset and dt:
                # A provider/date fingerprint is an ephemeral exact integer,
                # not a probabilistic hash with collisions.
                dk=(dataset,dt)
                pkey=provider_date_index.setdefault(dk,len(provider_date_index))
                provider_days[cell].add(pkey)
                region_days[(area,yy)].add(pkey)
            else:
                n_source["mapped_records_without_valid_dataset_date"]+=1
        if max_source_rows is not None and n_source["rows"]>=max_source_rows:
            break
        if chunks%10==0:
            print(json.dumps({"chunks":chunks,"source_rows":n_source["rows"],
                              "mapped":n_source["candidate_records_mapped_to_369"],
                              "source_positive_species":len(assigned_species),
                              "cells":len(record_cells),
                              "provider_date_groups":len(provider_date_index)}),flush=True)
    return {
        "source_counts":dict(n_source),
        "original_taxa_seen":sorted(original_taxa_seen),
        "years":dict(sorted(orig_years.items())),
        "basis_of_record":dict(basis_counts),
        "record_cells":record_cells,
        "precision_cells":precise_cells,
        "year_species":annual_species_region,
        "species_region_years":years_by_s_r,
        "provider_days":provider_days,
        "region_days":region_days,
        "providers":providers,
        "region_providers":region_providers,
        "total_unique_provider_date_keys":len(provider_date_index),
        "chunks":chunks,
    }


def save_outputs(state,bykey,resource,src_sha,meta,output=OUTPUT,
                 full_run=True):
    output.parent.mkdir(parents=True,exist_ok=True)
    species_to_key={v:k for k,v in bykey.items()}
    cells=state["record_cells"]
    head=["butterfly_species","frozen_gbif_taxon_key","wgsrpd3_code","year",
          "indexed_gbif_records","gbif_datasetkey_calendar_days",
          "distinct_gbif_dataset_keys","reported_uncertainty_le1000m_records"]
    file1=Path(str(output)+".species_region_year.csv.gz")
    with gzip.open(file1,"wt",encoding="utf8",newline="") as f:
        w=csv.writer(f);w.writerow(head)
        for (species,region,year),n in sorted(cells.items()):
            k=(species,region,year)
            w.writerow([species,species_to_key[species],region,year,n,
                        len(state["provider_days"].get(k,set())),
                        len(state["providers"].get(k,set())),
                        state["precision_cells"].get(k,0)])
    regionyear=Counter()
    for (_,region,year),n in cells.items():
        regionyear[(region,year)]+=n
    file2=Path(str(output)+".region_year_effort.csv")
    with file2.open("w",encoding="utf8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["wgsrpd3_code","year","indexed_candidate_butterfly_records",
                    "indexed_candidate_butterfly_species",
                    "distinct_gbif_datasetkey_calendar_days",
                    "distinct_gbif_dataset_keys"])
        for (region,year),n in sorted(regionyear.items()):
            w.writerow([region,year,n,
                        len(state["year_species"].get((region,year),set())),
                        len(state["region_days"].get((region,year),set())),
                        len(state["region_providers"].get((region,year),set()))])
    speciesregions=defaultdict(list)
    for (species,region,year),n in cells.items():
        speciesregions[(species,region)].append((year,n))
    file3=Path(str(output)+".species_region.csv.gz")
    with gzip.open(file3,"wt",encoding="utf8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["butterfly_species","frozen_gbif_taxon_key","wgsrpd3_code",
                    "gbif_indexed_records","gbif_provider_days_across_years",
                    "distinct_recorded_years","earliest_indexed_record_year",
                    "latest_indexed_record_year",
                    "records_2001_2010","records_2011_2017","records_2018_2024",
                    "years_2001_2010","years_2011_2017","years_2018_2024"])
        for (species,region),v in sorted(speciesregions.items()):
            byyr=dict(v)
            grouped=[[(year,n) for year,n in v if a<=year<=b]
                     for a,b in ((2001,2010),(2011,2017),(2018,2024))]
            w.writerow([species,species_to_key[species],region,sum(byyr.values()),
                        sum(len(state["provider_days"].get((species,region,y),set()))
                            for y in byyr),
                        len(byyr),min(byyr),max(byyr),
                        *[sum(n for year,n in group) for group in grouped],
                        *[len(group) for group in grouped]])
    source=state["source_counts"]
    positive_species=set(s for s,_,_ in cells)
    all_candidates=set(bykey.values())
    intersect={s for s in positive_species if s.casefold() in resource}
    sums=Counter(n for n in [len(s) for s in state["year_species"].values()])
    receipt={
        "schema":"chocho_global_GBIF_224_butterfly_369_region_year_effort_v1",
        "status":"SOURCE_RECONSTRUCTED_INDEPENDENT_INDEXED_OCCURRENCE_GEOGRAPHY_NOT_OCCUPANCY",
        "source":{
            "gbif_doi":DOI,"gbif_download_key":GBIF_KEY,
            "GBIF_archive_sha256":src_sha,
            "GBIF_zip_size_bytes":meta["zip_size_bytes"],
            "download_expected_rows":EXPECTED_RECORDS,
            "frozen_225_candidate_csv_git_blob_sha":"d4abc2edd762ef1003fff7c84c096c38883f0c20",
            "official_TDWG_git_commit":"52da7828aba9d461dd133c27b3bd7a4407161f54",
            "method_lock_sha256":sha256(LOCK),
            "ZIP_occurrence_member":meta["occurrence_member"],
            "ZIP_occurrence_member_uncompressed_bytes":meta["uncompressed_bytes"],
        },
        "records":source,
        "year_counts_mapped":state["years"],
        "source_basis_of_record":state["basis_of_record"],
        "coverage":{
            "original_candidate_species":len(all_candidates),
            "candidate_species_matched_to_a_369_region":len(positive_species),
            "candidate_species_with_no_369_region_record":sorted(all_candidates-positive_species),
            "number_of_frozen239_resource_species_among_224_indexed":len(intersect),
            "indexed_but_not_in_frozen239_panel_species":sorted(positive_species-intersect),
            "positive_species_region_year_cells":len(cells),
            "positive_species_region_cells":len(speciesregions),
            "region_year_cells_with_any_candidate_record":len(regionyear),
            "official_level3_region_codes":369,
            "regions_with_any_record":len({r for r,_ in regionyear}),
            "years_with_any_record":len(state["years"]),
            "total_distinct_source_dataset_calendar_day_ids_included":state["total_unique_provider_date_keys"],
            "max_candidate_species_indexed_in_single_region_year":
                max((len(v) for v in state["year_species"].values()),default=0),
            "max_records_in_single_species_region_year_cell":max(cells.values(),default=0)
        },
        "files":{
            file1.name:{"sha256":sha256(file1),"rows":len(cells)},
            file2.name:{"sha256":sha256(file2),"rows":len(regionyear)},
            file3.name:{"sha256":sha256(file3),"rows":len(speciesregions)}
        },
        "sampled_mode":not full_run,
        "limitations":[
            "Record-based GBIF catalogued detections are not occupancy, historical arrival, local host use or naturalization.",
            "Candidate source is 225 targeted butterflies, not general insect observation effort or all Lepidoptera.",
            "Provider-day is unique datasetKey×ISO day, not a verified biological survey visit.",
            "No species-absence cells were generated. Sources without GBIF records may have poor recording coverage.",
            "Date fields inconsistent with GBIF year are not counted as provider-days.",
            "No full distinct-GBIF-id audit; archive row count and source provider identity are reported.",
            "Original GBIF source and any individual latitude/longitude, datasetKey, GBIF ID, day never committed.",
            "Yearly GBIF observations are not historical host-plant occurrence histories.",
            "Subsequent host-specificity requires matched source, family/geographic negative controls and region-held-out validation.",
        ]
    }
    receipt_path=Path(str(output)+".audit.json")
    receipt_path.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",
                            encoding="utf8")
    return receipt


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-zip",type=Path,default=SOURCE_ZIP)
    p.add_argument("--geojson",type=Path,default=TDWG)
    p.add_argument("--candidate-csv",type=Path,default=CANDIDATES)
    p.add_argument("--resource-panel",type=Path,default=RESOURCE_PANEL)
    p.add_argument("--output-prefix",type=Path,default=OUTPUT)
    p.add_argument("--chunk-size",type=int,default=60000)
    p.add_argument("--max-source-rows",type=int)
    args=p.parse_args()
    if args.chunk_size<1000 or args.chunk_size>250000:
        p.error("chunk-size 1000–250000")
    src=args.source_zip
    size=src.stat().st_size
    if args.max_source_rows is None and size!=EXPECTED_ZIP_BYTES:
        raise RuntimeError("Full GBIF original ZIP size changed; refuse different source")
    bykey,byname,resource,_source_disposition=candidate_panel(
        args.candidate_csv,args.resource_panel)
    tree,codes=official_geography(args.geojson)
    checksum=sha256(src)
    with zipfile.ZipFile(src) as z:
        (member,sep,cols),audit=archive_table(z)
        zipinfo=z.getinfo(member)
        print(json.dumps({
            "source_sha256":checksum,"archive_bytes":size,
            "source_member":member,"uncompressed_bytes":zipinfo.file_size,
            "source_columns":list(cols),"diagnostics":audit
        }),flush=True)
        with z.open(member) as source:
            state=process_chunks(source,sep,cols,tree,codes,bykey,byname,
                                 args.chunk_size,args.max_source_rows)
    if args.max_source_rows is None and state["source_counts"]["rows"]!=EXPECTED_RECORDS:
        raise RuntimeError("Full GBIF source row-count integrity failed: "+
                           repr(state["source_counts"]["rows"]))
    receipt=save_outputs(state,bykey,resource,checksum,
        {"zip_size_bytes":size,"occurrence_member":member,
         "uncompressed_bytes":zipinfo.file_size},args.output_prefix,
        full_run=args.max_source_rows is None)
    print(json.dumps({
        "status":receipt["status"],
        "records":receipt["records"],
        "coverage":receipt["coverage"],
        "outputs":receipt["files"],
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
