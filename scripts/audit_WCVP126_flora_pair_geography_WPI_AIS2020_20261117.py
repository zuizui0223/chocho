#!/usr/bin/env python3
"""Strict 126-island WCVP pair geography and external 2020 vessel-source crosswalk.

This is only a source and identity audit, not a 2005 ocean-current model.
The 2020 name+position-verified 194-WPI marine subset is not a sample of all
worldwide shipping routes, and an unobserved ship link is UNKNOWN COVERAGE.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
from pyproj import Geod
from shapely.geometry import Point,shape
from shapely.strtree import STRtree

BLOBS={
 "wpi2019":"0c8b5de466c68bbd077dabce1f5511c2e84d42b8",
 "ais_direct":"b165c00d70c8ca3c1897745513c850247e44bec3",
 "ais_verified":"82fbbdeda2e5d350654802ead129c2743082ea27",
 "ais_receipt":"ed0b74df50ed6fe3ec02735fd43508e625691c83",
 "original_WCVP_port":"53b618fa51a1a418ccc6e4363f4c866c389c90cc"
}
TDWG_SHA="52da7828aba9d461dd133c27b3bd7a4407161f54"
G=Geod(ellps="WGS84")

def sha(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def read(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def numeric(s):
    return float(s) if s not in (None,"") else math.nan

def main(a):
    for k,p in (("wpi2019",a.wpi2019),("ais_direct",a.ais_direct),
      ("ais_verified",a.ais_verified),("ais_receipt",a.ais_receipt),
      ("original_WCVP_port",a.original_flora_port)):
        if sha(p)!=BLOBS[k]:raise RuntimeError(f"Source Git-blob mismatch {k}")
    isles=read(a.original_island_counts)
    pairs=read(a.pair_wcvp)
    if len(isles)!=126 or len(pairs)!=7875 or len({r["code"] for r in isles})!=126:
        raise RuntimeError("Original 126 WCVP geographic region pair cohort mismatch")
    codes=sorted(r["code"] for r in isles)
    orig={r["code"]:r for r in isles}
    old={r["code"]:r for r in read(a.original_flora_port)}
    if set(codes)!=set(old):raise RuntimeError("Frozen WCVP and island ports 126 codes misaligned")
    for code,r in orig.items():
        if (int(r["native_WCVP_accepted"])!=int(old[code]["native"]) or
            int(r["alien_WCVP_accepted"])!=int(old[code]["introduced"]) or
            int(r["WPI_2023_ports"])!=int(old[code]["wpi_port_count"])):
            raise RuntimeError("Accepted WCVP count or WPI ports drift "+code)
    d=json.loads(a.tdwg3.read_text())
    region_poly={}
    for f in d["features"]:
        code=str(f.get("properties",{}).get("LEVEL3_COD","")).strip()
        if code not in orig:continue
        if code in region_poly:raise RuntimeError("Duplicate original TDWG3")
        region_poly[code]=shape(f["geometry"])
    if set(region_poly)!=set(codes):raise RuntimeError("Pinned TDWG3 does not include original 126 islands")
    geoms=[region_poly[c] for c in codes];spatial=STRtree(geoms)
    reps={c:region_poly[c].representative_point() for c in codes}
    if any(not region_poly[c].covers(reps[c]) for c in codes):
        raise RuntimeError("TDWG island geometrical representative not inside region")
    wpi=json.loads(a.wpi2019.read_text())
    feat=wpi["features"]
    if len(feat)!=3630:raise RuntimeError("Original 2019 WPI source missing its full 3630 ports")
    port_pt={}
    matched=defaultdict(list)
    ambiguous=[]
    without=[]
    for f in feat:
        pid=int(f["properties"]["OBJECTID"])
        if pid in port_pt:raise RuntimeError("Duplicated original WPI OBJECTID")
        pt=shape(f["geometry"]);port_pt[pid]=pt
        index=spatial.query(pt)
        lands=[codes[int(i)] for i in index if geoms[int(i)].covers(pt)]
        if len(lands)==1:matched[lands[0]].append(pid)
        elif len(lands)>1:ambiguous.append(pid)
        else:without.append(pid)
    if len(port_pt)!=3630:raise RuntimeError("2019 WPI port ID uniqueness failed")
    strict_port={pid:code for code,ids in matched.items() for pid in ids}
    cover=read(a.ais_verified); links=read(a.ais_direct)
    receipt=json.loads(a.ais_receipt.read_text())
    if len(cover)!=194 or len(links)!=94 or receipt["distinct_WPI_origin_destination_port_pairs"]!=94:
        raise RuntimeError("Original full AIS2020 matched link cohort changed")
    valid={int(x["original_WPI_port_id"]) for x in cover
           if x["verified_external_2020_endpoint"]=="True"}
    if len(valid)!=46:raise RuntimeError("Original 46 confirmed 2020 WPI original port identities changed")
    if not valid.issubset(port_pt):raise RuntimeError("AIS2020 verified source not inside 2019 WPI3630")
    valid_region=defaultdict(set)
    for pid in valid:
        if pid in strict_port:valid_region[strict_port[pid]].add(pid)
    idx=defaultdict(int)
    qualified_source=0
    original_region_different=0
    remap_same=0
    unmatched_direct=0
    for r in links:
        i=int(r["origin_WPI_port_id"]);j=int(r["dest_WPI_port_id"])
        if i not in valid or j not in valid:raise RuntimeError("Unmatched port silently included in verified direct AIS source")
        n=int(r["source_observed_2020_triplet_count"])
        if n<=0:raise RuntimeError("Original 2020 OD sample has nonpositive record count")
        if i in strict_port and j in strict_port:
            qualified_source+=n
            a_code=strict_port[i];b_code=strict_port[j]
            if a_code!=b_code:
                idx[(a_code,b_code)]+=n
                original_region_different+=n
            else:
                remap_same+=n
        else:unmatched_direct+=n
    if sum(int(r["source_observed_2020_triplet_count"]) for r in links)!=442:
        raise RuntimeError("Original 442 source-verified WPI2020 matched directed links changed")
    zonepair={}
    for r in pairs:
        ca,cb=r["source_region"],r["destination_region"]
        key=tuple(sorted([ca,cb]))
        if key in zonepair:raise RuntimeError("Repeated 126 island unordered pair")
        if ca==cb or ca not in orig or cb not in orig:raise RuntimeError("Original region key mismatch")
        one,two=reps[ca],reps[cb]
        _,_,distance=G.inv(one.x,one.y,two.x,two.y)
        direct_ab=idx.get((ca,cb),0);direct_ba=idx.get((cb,ca),0)
        source_ready=len(valid_region[ca])>0 and len(valid_region[cb])>0
        if (direct_ab or direct_ba) and not source_ready:
            raise RuntimeError("Verified ship OD without confirmed WPI endpoint coverage on both source regions")
        zonepair[key]=dict(r,
            **{
            "centroid_geodesic_km":round(distance/1000,5),
            "geographic_distance_reference":"one shape.representative_point per official TDWG3 region",
            "same_plant_level1":int(orig[ca]["level1"]==orig[cb]["level1"]),
            "source_representative_lat":round(one.y,5),
            "source_representative_lon":round(one.x,5),
            "dest_representative_lat":round(two.y,5),
            "dest_representative_lon":round(two.x,5),
            "WPI_2019_strict_inside_TDWG_source_port_count":len(matched[ca]),
            "WPI_2019_strict_inside_TDWG_dest_port_count":len(matched[cb]),
            "AIS2020_verified_source_ports_in_region_A":len(valid_region[ca]),
            "AIS2020_verified_source_ports_in_region_B":len(valid_region[cb]),
            "AIS2020_direct_WPI_source_trips_A_to_B":direct_ab,
            "AIS2020_direct_WPI_source_trips_B_to_A":direct_ba,
            "AIS2020_symmetric_direct_source_rows":direct_ab+direct_ba,
            "AIS2020_both_regions_have_verified_port_endpoint":source_ready,
            "AIS2020_not_observed_is_not_actual_absence":True,
            "AIS2020_year_after_WCVP_unknown_introduction_year":True,
            "modern_port_pair_proxy_log1p_product":
                round(math.log1p(int(old[ca]["wpi_port_count"]))*
                      math.log1p(int(old[cb]["wpi_port_count"])),9),
            })
    if len(zonepair)!=7875:raise RuntimeError("Original all island pairs not exactly preserved")
    out=list(zonepair.values())
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"WCVP126_flora_similarity_strict_WPI2020_geographic_pairs.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),extrasaction="raise");w.writeheader();w.writerows(out)
    botht=sum(r["AIS2020_both_regions_have_verified_port_endpoint"] for r in out)
    vship=[r for r in out if r["AIS2020_symmetric_direct_source_rows"]>0]
    per_site=[{
        "region":c,"TDWG_Level1":orig[c]["level1"],
        "representative_lat":reps[c].y,"representative_lon":reps[c].x,
        "original_2023_WPI_indexed_port_count":int(old[c]["wpi_port_count"]),
        "2019_WPI_strict_matched_port_count":len(matched[c]),
        "verified_2020_AIS_port_IDs_strict_126":sorted(valid_region[c])
    } for c in codes]
    (a.out_dir/"strict126_TDWG_region_geography_WPI_AIS_coverage.json").write_text(json.dumps(per_site,indent=2)+"\n")
    audit={
        "schema":"chocho_exact_WCVP126_pair_source_geography_WPI_ports_AIS2020_gap_v1",
        "status":"SOURCE_7875_WCVP_PAIR_GEOGRAPHY_AND_VERIFIED_2020_AIS_COVERAGE_NOT_CAUSAL_TRANSPORT",
        "source_Git_blobs":BLOBS,
        "source_TDWG_WGSRPD3_commit":TDWG_SHA,
        "N_WCVP_native_alien_pairs":7875,
        "original_WPI2019_port_records":len(feat),
        "strict_WPI_inside_126_regions":len(strict_port),
        "WPI_outside_all_126_regions":len(without),
        "WPI_in_ambiguous_TDWG_polygons":len(ambiguous),
        "126_regions_with_one_or_more_strict_WPI2019_ports":sum(bool(matched[c]) for c in codes),
        "126_regions_with_one_or_more_source_verified_AIS2020_ports":sum(bool(valid_region[c]) for c in codes),
        "source_verified_2020_WPI194_port_IDs":len(valid),
        "source_verified_2020_WPI_link_rows":len(links),
        "source_verified_2020_port_pair_triplet_rows":442,
        "mapped_2020_WPI_port_pair_triplet_rows_to_126_TDWG":qualified_source,
        "unmapped_original_WPI_triplets":unmatched_direct,
        "mapped_2020_triplet_rows_to_same_TDWG3_region":remap_same,
        "mapped_2020_triplet_rows_to_different_TDWG3_regions":original_region_different,
        "strict126_unordered_plant_pairs_with_both_verified_2020_AIS_port_endpoints":botht,
        "strict126_unordered_plant_pairs_with_one_or_more_actual_source_2020_port_OD":len(vship),
        "verified_port_OD_unordered_pair_details":[
           {"regions":[r["source_region"],r["destination_region"]],
            "AIS2020_OBSERVED_2020_A_to_B":r["AIS2020_direct_WPI_source_trips_A_to_B"],
            "AIS2020_OBSERVED_2020_B_to_A":r["AIS2020_direct_WPI_source_trips_B_to_A"],
            "both_2020_port_verified":r["AIS2020_both_regions_have_verified_port_endpoint"]}
           for r in vship
        ],
        "not_2005_OSCAR_particle_access_matrix":True,
        "no_shipping_minus_current_causal_comparison_fitted":True,
        "limits":[
          "Region-to-region geodesic point distance is not marine shortest-water path or ocean current directed accessibility.",
          "TDWG3 code may include archipelagos; a single representative point is a geographic approximation, not one island centroid.",
          "Strict original WPI2019 coastal points need to fall INSIDE official TDWG polygons; ports just outside may be real but not forced.",
          "No strict matching of 2019 WPI ports to 2023 index is assumed because provider snapshots have different dates.",
          "Original AIS2020 WPI ports were matched by external exact name and <=25km against original WPI; most original port identities and shipping exposures unknown.",
          "A 2020 AIS observed port OD episode is not a seed/cargo voyage, and 2020 cannot predict pre-2020 plant introduction.",
          "No observed OD for a region pair does not mean absent ships, especially when most ports are unmatched.",
          "Native/alien WCVP similarity is co-reported contemporary accepted taxon distributions, not demonstrated dispersal event."
        ]
    }
    (a.out_dir/"WCVP126_geodesic_port_shipment_source_coverage_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n")
    print("WCVP126_TDWG_WPI2020_PAIR_COVERAGE",json.dumps({
        "status":audit["status"],"n_pairs":len(out),
        "WPI2019_strict_ports":len(strict_port),
        "AIS2020_verified_regions":audit["126_regions_with_one_or_more_source_verified_AIS2020_ports"],
        "AIS2020_source_unordered_OD_pairs":len(vship),
        "AIS2020_both_endpoints_covered_pairs":botht,
        "AIS2020_triplets_mapped":qualified_source,
        "original_AIS2020_triplets":442
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    for k in ("original-island-counts","pair-wcvp","original-flora-port","wpi2019",
              "ais-direct","ais-verified","ais-receipt","tdwg3","out-dir"):
        ap.add_argument("--"+k,required=True,type=Path)
    main(ap.parse_args())
