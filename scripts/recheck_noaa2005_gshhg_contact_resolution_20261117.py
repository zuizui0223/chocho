#!/usr/bin/env python3
"""Independent fixed-source audit of the 19 actual and 10 reversed GSHHG coast hits.

Only reuses SHA-pinned exact-vector matching outcome from original Actions run.
No newly selected islands, imagined long-distance landfalls or biological rates.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED={
 "OBSERVED_NOAA_FLOW":{"SAME_ORIGIN_GSHHG_POLYGON":13,
     "OTHER_GSHHG_POLYGON_IN_ORIGIN_TDWG3_REGION":3,
     "OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION":3},
 "MATHEMATICAL_REVERSED_FLOW":{"SAME_ORIGIN_GSHHG_POLYGON":3,
     "OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION":7}
}
THRESHOLDS=(0.001,0.5,1.,2.,4.,8.)

def find(root,name):
    z=list(root.rglob(name))
    if len(z)!=1:raise RuntimeError(f"Expected exact one source artifact {name}, got {len(z)}")
    return z[0]

def read(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))

def haversine_km(la1,lo1,la2,lo2):
    k=math.pi/180
    a=math.sin((la2-la1)*k/2)**2+math.cos(la1*k)*math.cos(la2*k)*math.sin((lo2-lo1)*k/2)**2
    return 12742.0176*math.atan2(math.sqrt(max(a,0)),math.sqrt(max(1-a,0)))

def median(v):
    return round(float(statistics.median(v)),4) if v else None

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    rows=read(find(a.original_root,"original_noaa2005_gshhg_land_contact_identity.csv"))
    qa=json.loads(find(a.original_root,"original_noaa2005_gshhg_land_contact_identity_audit.json").read_text())
    if qa["source_particle_CI"]!=43902616033 or qa["source_recontact_CI"]!=43903614189 or len(rows)!=29:
        raise RuntimeError("Original landmask-contact source population differs")
    original=qa["original_WPI_coastal_anchor_polygon_gates"]
    actual=defaultdict(list)
    for r in rows:
        actual[r["scenario"]].append(r)
        if r["is_external_destination_verified"]!="False":
            raise RuntimeError("An independently verified biological arrival was invented")
    if {k:len(v) for k,v in actual.items()} != {"OBSERVED_NOAA_FLOW":19,"MATHEMATICAL_REVERSED_FLOW":10}:
        raise RuntimeError("Source scenario counts changed")
    summary={}
    details=[]
    for scenario,rs in actual.items():
        types=Counter(z["coast_identity_class"] for z in rs)
        if types!=EXPECTED[scenario]:
            raise RuntimeError(f"Source classification drift {scenario}: {types}")
        geo=Counter(z["nearest_coast_method"] for z in rs)
        source=Counter(z["source"] for z in rs)
        offsets=[]
        for z in rs:
            orig=original[z["source"]]["port_point"]
            dist=haversine_km(float(z["contact_lat"]),float(z["contact_lon"]),
                              float(orig[1]),float(orig[0]))
            vdist=float(z["nearest_coast_distance_km"])
            offsets.append(vdist)
            details.append({
              "scenario":scenario,"source":z["source"],
              "failure_day":int(z["failure_day"]),
              "original_port_distance_km":round(dist,3),
              "GSHHG_vector_coast_distance_km":vdist,
              "vector_coast_point_inside_polygon":z["nearest_coast_method"]=="VECTOR_POLYGON_CONTAINS_POINT",
              "exactly_same_GSHHG_source_polygon":
                    z["coast_identity_class"]=="SAME_ORIGIN_GSHHG_POLYGON",
              "origin_TDWG3_exact_region":z["coast_identity_class"]=="OTHER_GSHHG_POLYGON_IN_ORIGIN_TDWG3_REGION",
              "region_unresolved":z["coast_identity_class"]=="OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION",
              "native_gshhg_polygon_id":z["contact_GSHHG_native_id"],
            })
        mask_inside=sum(z["nearest_coast_method"]=="VECTOR_POLYGON_CONTAINS_POINT" for z in rs)
        summary[scenario]={
         "original_contacts":len(rs),
         "class_counts":dict(types),
         "by_source":dict(source),
         "point_inside_original_GSHHG_vector_polygon":mask_inside,
         "point_nearest_GSHHG_vector_but_outside":len(rs)-mask_inside,
         "contact_point_vector_distance_median_km":median(offsets),
         "contact_point_vector_within_km":{
              str(k):sum(z<=k for z in offsets) for k in THRESHOLDS},
         "source_WPI_port_distance_median_km":median(
              [z["original_port_distance_km"] for z in details if z["scenario"]==scenario]),
         "class_by_source":{
            src:dict(Counter(z["coast_identity_class"] for z in rs if z["source"]==src))
            for src in ("AZO","CNY","CVI")
         },
         "original_ambiguous_three_regions":[{
           "source":z["source"],"day":int(z["failure_day"]),
           "coast_km":float(z["nearest_coast_distance_km"]),
           "gshhg_poly_index":int(z["contact_GSHHG_polygon_index"]),
           "TDWG3_region":z["contact_TDWG3_regions"],
           "port_name":z["source_WPI_name"]
         } for z in rs if z["coast_identity_class"]=="OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION"]
        }
    audit={
      "schema":"chocho_original_GSHHG_vector_contact_29_identifiability_review_v1",
      "status":"VERIFIED_SOURCE_LAND_CONTACT_IDENTIFIABILITY_NOT_LARVAL_OR_SEED_LANDFALL",
      "parent_geometric_CI":43906815380,
      "parent_trajectory_CI":43902616033,
      "source_original_number_tracks":81,
      "original_GSHHG_archive_SHA256":qa["GSHHG_original_archive_sha256"],
      "vector_distance_thresholds_km":THRESHOLDS,
      "source_groups":summary,
      "stratified_contact_rows":details,
      "inference_guards":[
       "Exact GSHHG polygon index match supports physical recontact of the same mapped land polygon, not a biologically viable seed return.",
       "A raster collision close to but outside vector polygon is a mask-resolution mismatch; classification is geographically uncertain.",
       "TDWG3 source region is an archipelago, so other vector polygon within source region is not necessarily transport to a single geographically distinct island.",
       "A 50km port-distance negative proxy is NOT proof of absent remote island coastline encounters.",
       "Do not use mathematical reversed current as independent placebo or conclude arrival probability from 27 simulated source releases.",
       "Original EuDiS/WCVP plant trait/host geography, actual historical ship OD, seed viability and butterflies are not in this stage."
      ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("GSHHG_CONTACT_IDENTIFIABILITY_DIAGNOSTIC",json.dumps({
        "source":43906815380,
        "actual":{k:v for k,v in summary["OBSERVED_NOAA_FLOW"].items() if k!="original_ambiguous_three_regions"},
        "original_ambiguous_three_regions":summary["OBSERVED_NOAA_FLOW"]["original_ambiguous_three_regions"],
        "reversed_counts":summary["MATHEMATICAL_REVERSED_FLOW"]["class_counts"],
        "status":audit["status"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
