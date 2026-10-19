#!/usr/bin/env python3
"""Audit modern iNaturalist host observations against source-pinned WCVP native TDWG L3.

This does NOT infer historical plant occupancy, dispersal corridors, population
genetics or that a record outside WCVP's native map is necessarily introduced.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from shapely.geometry import Point,shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree
from pyproj import Geod

RDATA_SOURCE="65bed76bae9d644ccb6ad200c05f9f5071d89e05"
TTF_SOURCE="187e15c6a7b777c119b97b8a43a0d92783ef7c5c"
SOURCE_SUMMARY_SHA="cffb1144a61151c9da78737b1c7a84034398fb50d03ea5ecb46642b317d93634"
MANIFEST_SHA="a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
GEOD=Geod(ellps="WGS84")
NATIVE="WITHIN_DOCUMENTED_NATIVE_L3"
OUTSIDE="OUTSIDE_DOCUMENTED_NATIVE_L3_NOT_PROVEN_INTRODUCED"
BOUNDARY="REGIONAL_BOUNDARY_UNCERTAIN"
UNMAPPED="UNMAPPED_TO_TDWG_L3"
NO_NATIVE="NO_DOCUMENTED_NATIVE_L3"

def sha256_path(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def canonical(data):
    return json.dumps(data,ensure_ascii=False,sort_keys=True,separators=(",",":"),allow_nan=False).encode()

def native_table(path):
    by_id=defaultdict(set)
    with Path(path).open(encoding="utf8",newline="") as f:
        for r in csv.DictReader(f):
            k=str(int(r["accepted_plant_name_id"]));code=r["area_code_l3"].strip()
            if not code:raise ValueError("empty WCVP region code")
            by_id[k].add(code)
    if len(by_id)!=517:raise ValueError("expected exactly 517 documented native WCVP host taxa")
    return by_id

def load_regions(path):
    payload=json.loads(Path(path).read_text(encoding="utf8"))
    features=payload["features"]
    geoms=[];codes=[]
    for f in features:
        code=str(f["properties"].get("area_code_l3","")).strip()
        geom=shape(f["geometry"])
        if not code or geom.is_empty or not geom.is_valid:
            raise ValueError("empty or invalid frozen WCVP TDWG polygon")
        codes.append(code);geoms.append(geom)
    if len(codes)!=369 or len(set(codes))!=len(codes):
        raise ValueError("original 369 TDWG region identities changed")
    return geoms,codes

def make_matcher(geoms,codes):
    tree=STRtree(geoms)
    old_shapely_lookup={id(g):i for i,g in enumerate(geoms)}
    def matching_regions(lon,lat):
        point=Point(lon,lat)
        matched=[]
        for candidate in tree.query(point):
            i=int(candidate) if isinstance(candidate,(int,np.integer)) else old_shapely_lookup[id(candidate)]
            if geoms[i].covers(point):matched.append(i)
        return point,sorted(set(matched))
    return matching_regions

def boundary_distance_m(point,geom):
    border=geom.boundary
    near=nearest_points(point,border)[1]
    _,_,meters=GEOD.inv(point.x,point.y,near.x,near.y)
    return abs(float(meters))

def classify_observation(*,lat,lon,uncertainty_m,known_native_codes,matcher,geoms,codes):
    if not np.isfinite([lat,lon,uncertainty_m]).all() or not -90<=lat<=90 or not -180<=lon<=180 or not 0<=uncertainty_m<=1000:
        raise ValueError("invalid quality-accepted coordinates/uncertainty")
    if not known_native_codes:return NO_NATIVE, ""
    point,ids=matcher(lon,lat)
    if not ids:return UNMAPPED,""
    matches=[codes[i] for i in ids]
    near=any(boundary_distance_m(point,geoms[i])<=uncertainty_m+1e-6 for i in ids)
    if len(ids)>1 or near:return BOUNDARY,"|".join(sorted(matches))
    code=matches[0]
    return (NATIVE if code in known_native_codes else OUTSIDE),code

def iter_host_checkpoints(directory):
    root=Path(directory)
    if sha256_path(root/"summary.json")!=SOURCE_SUMMARY_SHA:
        raise ValueError("v0.2 source summary bytes changed")
    summary=json.loads((root/"summary.json").read_text())
    binding=summary["binding"]
    if (summary["attempted_taxa"],summary["pending_taxa"],summary["quality_records"])!=(518,0,296697):
        raise ValueError("original host source completeness mismatch")
    if binding["manifest_sha256"]!=MANIFEST_SHA or not binding["not_equivalent_to_v01"]:
        raise ValueError("host source binding mismatch")
    checkpoints=sorted((root/"checkpoints").glob("*.json"))
    if len(checkpoints)!=518:raise ValueError("518 immutable host checkpoints not available")
    for path in checkpoints:
        obj=json.loads(path.read_text())
        result=obj["result"]
        key=str(int(result["wcvp_id"]))
        if obj["binding"]!=binding or path.stem!=key:
            raise ValueError("source WCVP ID/binding drift")
        if hashlib.sha256(canonical(result)).hexdigest()!=obj["result_sha256"]:
            raise ValueError("source checkpoint digest mismatch")
        if result["status"]=="OBSERVATIONS_ACQUIRED":
            if result["quality_records"]!=len(result["points"]):
                raise ValueError("source count mismatch")
        elif result["points"] or result["quality_records"] is not None:
            raise ValueError("unmatched taxonomy cannot be treated as 0 observations")
        yield result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native",required=True,type=Path)
    ap.add_argument("--geojson",required=True,type=Path)
    ap.add_argument("--host-snapshot",required=True,type=Path)
    ap.add_argument("--host-manifest",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    if sha256_path(args.host_manifest)!=MANIFEST_SHA:
        raise ValueError("TTF frozen 1037 host associations changed")
    native=native_table(args.native)
    geoms,codes=load_regions(args.geojson)
    if set().union(*native.values())-set(codes):
        raise ValueError("WCVP native code missing matching region polygon")
    matcher=make_matcher(geoms,codes)
    args.output.mkdir(parents=True,exist_ok=True)
    species=[];global_counts=Counter();source_n=0
    gzpath=args.output/"native_supported_observations.csv.gz"
    fields=["wcvp_id","host_name","gbif_id","latitude","longitude",
            "uncertainty_m","event_date","native_wgsrpd3_code","geography_status"]
    with gzip.open(gzpath,"wt",newline="",encoding="utf8",compresslevel=6) as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
        for host in iter_host_checkpoints(args.host_snapshot):
            key=str(int(host["wcvp_id"]));counts=Counter()
            expected_native=native.get(key,set())
            for rec in host["points"]:
                state,code=classify_observation(
                    lat=float(rec["latitude"]),lon=float(rec["longitude"]),
                    uncertainty_m=float(rec["uncertainty_m"]),
                    known_native_codes=expected_native,
                    matcher=matcher,geoms=geoms,codes=codes)
                counts[state]+=1;global_counts[state]+=1;source_n+=1
                if state==NATIVE:
                    writer.writerow({"wcvp_id":key,"host_name":host["name"],
                        "gbif_id":rec["gbif_id"],
                        "latitude":rec["latitude"],"longitude":rec["longitude"],
                        "uncertainty_m":rec["uncertainty_m"],
                        "event_date":rec.get("event_date") or "",
                        "native_wgsrpd3_code":code,"geography_status":state})
            species.append({"wcvp_id":key,"host_name":host["name"],
                "taxon_status":host["status"],
                "native_L3_count":len(expected_native),
                "quality_observations":len(host["points"]),
                "within_documented_native_L3":counts[NATIVE],
                "outside_documented_native_L3":counts[OUTSIDE],
                "boundary_uncertain":counts[BOUNDARY],
                "unmapped_L3":counts[UNMAPPED],
                "no_native_region_source":counts[NO_NATIVE]})
    if source_n!=296697 or len(species)!=518:raise ValueError("complete source records not reconciled")
    assert sum(global_counts.values())==source_n
    spath=args.output/"host_native_observation_audit.csv"
    with spath.open("w",newline="",encoding="utf8") as stream:
        w=csv.DictWriter(stream,fieldnames=list(species[0]));w.writeheader()
        w.writerows(sorted(species,key=lambda x:int(x["wcvp_id"])))
    summary={"schema":"chocho_TTF_WCVP_native_occurrence_audit_v0.1",
        "status":"DESCRIPTIVE_MODERN_NATIVE_L3_CONSISTENCY_NOT_LGM_OCCUPANCY",
        "all_host_taxa":518,"with_documented_native_L3":517,
        "source_quality_observations":source_n,
        "source_status_counts":dict(global_counts),
        "plants_with_native_observations":sum(x["within_documented_native_L3"]>0 for x in species),
        "plants_with_at_least_30_native_observations":sum(x["within_documented_native_L3"]>=30 for x in species),
        "original_WCVP_2024_snapshot":RDATA_SOURCE,
        "TTF_source_commit":TTF_SOURCE,
        "source_sha256":{p.name:sha256_path(p) for p in (args.native,args.geojson,args.host_manifest,
                                                  args.host_snapshot/"summary.json")},
        "inferential_limits":[
            "WDWSRPD L3 polygons are present-day broad native regional attribution, not individual-site evidence or LGM occupancy",
            "Outside a documented native L3 region is NOT proof of introduced status, just undocumented native at coarse resolution",
            "Input GBIF iNaturalist research-grade observations are nonrandom, and 2010-2026 year-bin retrieval is limited to 300 per bin",
            "No genetic response or inferred LGM migration path used"]}
    (args.output/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()
