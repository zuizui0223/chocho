#!/usr/bin/env python3
"""Post-result audit: actual path lengths and fixed recipient wet anchors.

This diagnostic is necessary-reach geometry, not observed landing, not a
complete coastline test, and not a new confirmatory ecological regression.
"""
import argparse, collections, hashlib, json, math
from pathlib import Path
import numpy as np
from pyproj import Geod
G=Geod(ellps="WGS84")
THRESH=50.
def dist(y,x,ys,xs):
    _,_,m=G.inv(np.full(len(ys),x,dtype=float),np.full(len(ys),y,dtype=float),xs,ys)
    return np.abs(m)/1000
def q(v):
    return np.quantile(v,[0,.25,.5,.75,1]).tolist() if len(v) else []
def main():
    p=argparse.ArgumentParser()
    for name in ("physical","tracks","anchors","preflight","robustness","out"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    phy=json.loads(a.physical.read_text());tracks=json.loads(a.tracks.read_text());pre=json.loads(a.preflight.read_text())
    rob=json.loads(a.robustness.read_text())
    assert phy["launches"]==len(tracks)==306 and phy["day30_valid"]==52 and phy["source_regions_valid30"]==12
    assert phy["technical_gate_pass"] and pre["all_306_anchor_date_queries_verified"]
    assert rob["original_model_coefficients_reproduced"] is True
    expected={(s["source"],d["date"],i+1):v for s in pre["region_launches"] for d in s["dates"] for i,v in enumerate(d["starts"])}
    actual={(r["source"],r["date"],r["start_index"]):r["start"] for r in tracks}
    assert len(expected)==len(actual)==306 and actual==expected
    ac=json.loads(a.anchors.read_text())["coastal_anchors_by_TDWG_region"]
    assert len(ac)==126 and sum(bool(v) for v in ac.values())==104
    dest={s:np.array([(float(x["snapped_wet_latitude"]),float(x["snapped_wet_longitude"])) for x in v]) for s,v in ac.items() if v}
    assert all(np.isfinite(v).all() for v in dest.values())
    allrows=[];per=[]
    for r in tracks:
        s=r["source"];start=r["start"];L=float(r["path_length_km"])
        assert L>=0 and r["net_displacement_km"]<=L+1e-6
        assert s in dest and r["elapsed_hours"]<=720
        rows=[]
        for d,coords in dest.items():
            if d==s:continue
            d0=float(dist(start["lat"],start["lon"],coords[:,0],coords[:,1]).min())
            d1=float(dist(r["end_lat"],r["end_lon"],coords[:,0],coords[:,1]).min())
            # Necessary distance bound for approaching any of these specific
            # wet anchor points within 50km over the *computed* path. Does not
            # describe the entire coast or uncomputed future path.
            eligible=d0<=L+THRESH+1e-6
            endpoint_approach=d0>THRESH and d1<=THRESH
            assert not endpoint_approach or eligible
            rows.append({"source":s,"destination":d,"date":r["date"],"start_index":r["start_index"],"day30_valid":r["day30_valid"],"status":r["status"],"elapsed_hours":r["elapsed_hours"],"path_length_km":L,"initial_anchor_distance_km":d0,"final_anchor_distance_km":d1,"initial_within50":d0<=THRESH,"path_length_compatible_with_fixed_anchor50":eligible,"new_endpoint_within50":endpoint_approach,"endpoint_distance_reduction_km":d0-d1})
        assert len(rows)==103
        allrows.extend(rows)
        per.append({"source":s,"date":r["date"],"start_index":r["start_index"],"valid30":r["day30_valid"],"path_length_km":L,"net_displacement_km":r["net_displacement_km"],"minimum_other_fixed_anchor_km_start":min(x["initial_anchor_distance_km"] for x in rows),"minimum_other_fixed_anchor_km_end":min(x["final_anchor_distance_km"] for x in rows),"candidate_other_recipients":sum(x["path_length_compatible_with_fixed_anchor50"] and not x["initial_within50"] for x in rows)})
    assert len(allrows)==306*103
    def describe(rows):
        source_set={r["source"] for r in rows}
        reach=[r for r in rows if r["path_length_compatible_with_fixed_anchor50"] and not r["initial_within50"]]
        initial=[r for r in rows if r["initial_within50"]]
        contact=[r for r in rows if r["new_endpoint_within50"]]
        dyads={(r["source"],r["destination"]) for r in reach}
        rs={r["source"] for r in reach};rd={r["destination"] for r in reach}
        return {"tracks":len({(r["source"],r["date"],r["start_index"]) for r in rows}),"source_regions":len(source_set),"track_recipient_combinations":len(rows),"initially_within50_combinations":len(initial),"length_compatible_noninitial_combinations":len(reach),"length_compatible_directed_dyads":len(dyads),"length_compatible_source_regions":len(rs),"length_compatible_recipient_regions":len(rd),"length_compatible_source_codes":sorted(rs),"new_endpoint_within50_combinations":len(contact),"new_endpoint_within50_dyads":len({(r["source"],r["destination"]) for r in contact}),"length_compatible_fraction_among_noninitial":len(reach)/max(1,len(rows)-len(initial)),"source_gate_for_ecological_fit_pass":len(rs)>=10 and len(rd)>=20 and len(dyads)>=100,"contact_examples":[{k:r[k] for k in ("source","destination","date","start_index","day30_valid","initial_anchor_distance_km","final_anchor_distance_km","path_length_km")} for r in contact],"dyad_list":[list(x) for x in sorted(dyads)]}
    valid=[r for r in allrows if r["day30_valid"]]
    truncated=[r for r in allrows if not r["day30_valid"]]
    complete=describe(valid)
    assert complete["tracks"]==52 and complete["source_regions"]==12
    per_source={}
    for s in sorted({r["source"] for r in allrows}):
        pp=[r for r in per if r["source"]==s and r["valid30"]]
        d=describe([r for r in valid if r["source"]==s])
        d.update({"day30_path_length_km_quantiles":q([r["path_length_km"] for r in pp]),"day30_net_displacement_km_quantiles":q([r["net_displacement_km"] for r in pp])})
        per_source[s]=d
    # Preserve the original global positive slope and identify its source
    # leverage, without interpreting distant endpoint projections as routes.
    b=rob["diagnostics"]["primary_fixed22_50km_hycom"]["all_pairs"]["results"]
    japan={k:{"full_beta":v["beta"],"without_JAP_beta":v["leave_one_source"]["JAP"],"JAP_leverage":v.get("source_leverage",{}).get("JAP")} for k,v in b.items()}
    result={"status":"FIXED_ANCHOR_PATH_LENGTH_COMPATIBILITY_AUDIT_COMPLETE","complete_day30":complete,"truncated_computed_portion_only":describe(truncated),"all_computed_paths":describe(allrows),"complete_path_length_km_quantiles":q([r["path_length_km"] for r in per if r["valid30"]]),"complete_net_displacement_km_quantiles":q([r["net_displacement_km"] for r in per if r["valid30"]]),"per_source_complete":per_source,"original_JAP_leverage":japan,"regression_fitted":False,"limits":["Distance to fixed original wet anchor points, not all coastline pixels or true shore landing.","Path length is a necessary geometric bound for reaching those points within 50km, not a sufficient water route or observed arrival.","End position only; intermediate near passes cannot be recovered from endpoint-only archives.","For truncated trajectories the bound applies only to the computed prefix, not the missing remaining days.","All outcomes post-result exploratory, original positive global association and adverse sensitivities retained."]}
    (a.out/"FIXED_ANCHOR_REACH_RESULT.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FIXED_ANCHOR_REACH_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
