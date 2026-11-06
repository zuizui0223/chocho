#!/usr/bin/env python3
"""Recover the last two original 51-island eligible port groups *exactly*
from source-verified contemporary WPI OBJECTIDs and the original fixed
1,275 geodesic-distance matrix, with 49 other island locations independently
recovered to machine precision.

This is archival SOURCE RECONSTRUCTION from pre-existing frozen geographic
distances (not selection on plant or butterfly response). Do not accept merely
the best-fitting/closest subset: require unique extra WPI port IDs yielding
all original distances within a specified strict numeric tolerance.

Inputs: successful source-pinned workflow 42173280157 artifact, official NOAA
WPI 3630 source workflow 42172382077 artifact, original JBI 51-pair CSV.
"""
from __future__ import annotations
import argparse,itertools,hashlib,json,math,statistics
from pathlib import Path
import numpy as np
import pandas as pd

R=6371.0088
TOL_KM=1e-6  # Far below 0.01 km metadata acceptance, near floating-point noise.
KNOWN_MIN=49
N=51

def mean_xyz(port_ids,coords):
    ll=np.radians(np.array([coords[k] for k in port_ids],dtype=float))
    x,y=ll[:,0],ll[:,1]
    xyz=np.array([np.cos(y)*np.cos(x),np.cos(y)*np.sin(x),np.sin(y)]).mean(axis=1)
    norm=np.linalg.norm(xyz)
    if norm<.01:raise ValueError("Source port vectors have no valid direction")
    return xyz/norm

def angular_distance(a,b):
    return float(R*np.arccos(np.clip(a@b,-1,1)))

def invert_xyz(centroid):
    x,y,z=map(float,centroid)
    return {"longitude":float(np.degrees(np.arctan2(y,x))),
            "latitude":float(np.degrees(np.arctan2(z,np.hypot(x,y))))}

def candidates_one_missing(base,raw,n_expected,coords,source_distance,anchors,tolerance=TOL_KM):
    """Enumerate all additions (usually exactly one); don't select approximate."""
    b=sorted(set(base),key=int)
    pool=sorted(set(raw)-set(b),key=int)
    missing=n_expected-len(b)
    if missing<=0:raise ValueError("Source candidate island is not missing a port")
    if len(pool)<missing:raise ValueError("No candidate port ID available to recover")
    details=[]
    for added in itertools.combinations(pool,missing):
        ids=b+list(added)
        centroid=mean_xyz(ids,coords)
        errors=[abs(angular_distance(centroid,pos)-source_distance[k])
                for k,pos in anchors.items()]
        row={"extra_port_ids":list(added),"full_port_ids":sorted(ids,key=int),
             "max_source_distance_error_km":max(errors),
             "median_source_distance_error_km":float(np.median(errors)),
             "mean_source_distance_error_km":float(np.mean(errors)),
             "n_exact_anchors":sum(e<=tolerance for e in errors)}
        details.append((row,centroid))
    details.sort(key=lambda x:x[0]["max_source_distance_error_km"])
    successful=[x for x in details if x[0]["max_source_distance_error_km"]<=tolerance]
    return {
      "candidate_count":len(details),
      "candidate_ids_not_passing_single_landmass":pool,
      "n_candidates_uniquely_satisfying_all_known_distances":len(successful),
      "best_6": [x[0] for x in details[:6]],
      "winner":successful[0][0] if len(successful)==1 else None,
      "center":successful[0][1] if len(successful)==1 else None
    }


def run(probe_path,port_path,source_pair_path,source_island_path):
    probe=json.loads(probe_path.read_text(encoding="utf-8"))
    if probe["schema"]!="chocho_WPI3630_GSHHG237_51island_provenance_reconstruction_v1":
        raise ValueError("Wrong original GSHHG/WPI receipt source")
    g=probe["source_mapping_methods"]["geodesic_5km_intersects_exactly_one_global_landmass"]
    raw=probe["source_mapping_methods"]["geodesic_5km_buffer_intersects_target_landmass"]
    if g["count_matched_islands"]!=49 or g["n_port_count_mismatch"]!=2:
        raise ValueError("Exactly 49 independently source-recovered port groups required")
    failed=set(g["port_count_mismatch"])
    if failed!={"12","26"}:
        raise ValueError("Unexpected historical 2 remaining mismatch ids "+str(failed))
    if g["distances_both_port_counts_matched"]["n_pairs"]!=1176 or \
       g["distances_both_port_counts_matched"]["max_error_km"]>TOL_KM:
        raise ValueError("49 native source reconstructed island pair distances are not exact")

    obj=json.loads(port_path.read_text(encoding="utf-8"))
    rows=obj["features"]
    if len(rows)!=3630:raise ValueError("NOAA source must retain all WPI 2019 ports")
    coords={str(int(v["properties"]["OBJECTID"])):v["geometry"]["coordinates"] for v in rows}
    if len(coords)!=3630:raise ValueError("Repeated or incomplete original WPI port IDs")
    orig=pd.read_csv(source_pair_path,dtype={"island_a":str,"island_b":str})
    islands=pd.read_csv(source_island_path,dtype={"landmass_id":str})
    valid=islands[islands.eligible.astype(str).str.casefold().eq("true")]
    counts=dict(zip(valid.landmass_id,valid.ports.astype(int)))
    if len(counts)!=N or len(orig)!=1275:raise ValueError("Archived 51 island input size changed")
    ids=sorted(counts,key=int)
    lookup={}
    for r in orig.itertuples(index=False):
        lookup[(r.island_a,r.island_b)]=float(r.port_distance_km)
        lookup[(r.island_b,r.island_a)]=float(r.port_distance_km)

    existing={id:g["matched_port_id_lists"][id] for id in ids}
    known={id:mean_xyz(existing[id],coords) for id in ids if id not in failed}
    if len(known)!=49:raise ValueError("Missing one of 49 externally source-recovered coordinates")
    strict_errors=[
       abs(angular_distance(known[a],known[b])-lookup[(a,b)])
       for a,b in itertools.combinations(known,2)]
    if len(strict_errors)!=1176 or max(strict_errors)>TOL_KM:
        raise ValueError("49 known recovered WPI centroids fail independent pair distances")

    recovery={}
    for id in sorted(failed,key=int):
        anchors={k:xyz for k,xyz in known.items()}
        archived={k:lookup[(id,k)] for k in known}
        results=candidates_one_missing(existing[id],raw["matched_port_id_lists"][id],
                                       counts[id],coords,archived,anchors)
        recovery[id]={k:v for k,v in results.items() if k!="center"}
        if results["winner"] is not None:
            existing[id]=results["winner"]["full_port_ids"]
            known[id]=results["center"]

    source_status="NO_UNIQUE_SOURCE_RECOVERY"
    max_error=None
    if len(known)==N:
        all_err=[abs(angular_distance(known[a],known[b])-lookup[(a,b)])
                 for a,b in itertools.combinations(ids,2)]
        max_error=max(all_err)
        if len(all_err)!=1275:raise ValueError("Incomplete full archival check")
        if max_error<=TOL_KM:
            source_status="ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES"
        else:source_status="SOURCE_PORT_CANDIDATE_MATCH_INCONSISTENT_ON_1275_PAIRS"
    result={
      "schema":"chocho_original_WPI2019_GSHHG237_island51_unique_port_reconstruction_v1",
      "status":source_status,
      "archive_original_49_source_recovered_distance_max_km":max(strict_errors),
      "archive_original_1275_recovered_distance_max_km":max_error,
      "source_original_islands":N,"source_original_pairs":len(orig),
      "source_missing_ids":sorted(failed,key=int),
      "source_missing_port_candidate_diagnostics":recovery,
      "previously_known_source_islands":49,
      "source_WPI_port_checksum":hashlib.sha256(port_path.read_bytes()).hexdigest(),
      "source_GSHHG_portmap_receipt_checksum":hashlib.sha256(probe_path.read_bytes()).hexdigest(),
      "source_pair_sha256":hashlib.sha256(source_pair_path.read_bytes()).hexdigest(),
      "source_semantics":"NOAA NGA2019 ports + GSHHG237 sibling IDs; original complete-port status and archived pairwise geodesics, not inferred from biological response",
      "data_availability_constraint":"Source port subset recovered by exact geographical inverse-matching; raw original mapping table remains unarchived.",
      "limitation":["Port subset reconciliation is an archival inverse problem using the old distance matrix. It is valid only if unique under exact verification.",
                    "This proves source centroid consistency, not ship connectivity, cargo movement or port-mediated plant transport.",
                    "WPI2019 coordinates, GSHHG237 shapes and 2001–2024 plant record completeness are not historical time-varying navigation data."]
    }
    if source_status.startswith("ALL_51_"):
        result["source_51_port_groups"]={id:existing[id] for id in ids}
        result["verified_source_51_mean_port_points"]={
            id:{**invert_xyz(known[id]),"n_complete_single_landmass_ports":len(existing[id])}
            for id in ids}
        result["total_51_original_port_assignments"]=sum(len(x) for x in existing.values())
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("probe","ports","pairs","islands","out"):p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    r=run(a.probe,a.ports,a.pairs,a.islands)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ISLAND51_EXACT_NOAA_SOURCE_PORT_INVERSE_RECONCILIATION",json.dumps({
        "status":r["status"],"missing_port_ids":r["source_missing_port_candidate_diagnostics"],
        "all_1275_source_distance_error_km":r["archive_original_1275_recovered_distance_max_km"],
        "total_unique_source_assigned_ports":r.get("total_51_original_port_assignments")
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
