#!/usr/bin/env python3
"""Native-resolution OSCAR source-coverage and fixed-30d particle diagnostic.

New physical design after the coarse-grid no-fit outcome. Frozen plant support,
same 50-km anchor and 30-day thresholds, six prespecified 2005 release dates,
up to three spatially separated source starts per date. No observed arrivals.
"""
import argparse, hashlib, io, json, math, time, zipfile
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request,urlopen
import sys
import netCDF4
import numpy as np
import pandas as pd
from pyproj import Geod

G=Geod(ellps="WGS84")
BASE="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.nc?"
RELEASES=[907,913,919,925,931,937]
FIRST=907;LAST=945
HORIZON=30
MAX_OFFSET=50.
MIN_SEP=15.
NSTART=3
GRAPH_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
STATUS="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE"
def sha(b):return hashlib.sha256(b).hexdigest()
def get(start,stop,step=1):
    # <=3 full native grid time slices per request, exact original 15-m field.
    expr=",".join(f"{v}[{start}:{step}:{stop}][0][0:1:480][0:1:1200]" for v in ("u","v"))
    url=BASE+quote(expr,safe="[],:")
    last=None
    for k in range(3):
        try:
            with urlopen(Request(url,headers={"User-Agent":"chocho-native-OSCAR-coast-source-audit/1.0"}),timeout=150) as r:
                b=r.read(32_000_001)
                assert r.status==200 and len(b)<=32_000_000
            assert b.startswith(b"CDF") or b.startswith(b"\x89HDF\r\n\x1a\n")
            with netCDF4.Dataset("source.nc",memory=b) as d:
                def arr(v):
                    a=d[v][:];return np.asarray(a.filled(np.nan) if np.ma.isMaskedArray(a) else a)
                lat=arr("latitude").astype(float);lon=arr("longitude").astype(float)
                u=arr("u")[:,0].astype(np.float32);v=arr("v")[:,0].astype(np.float32)
                t=netCDF4.num2date(d["time"][:],d["time"].units,calendar=getattr(d["time"],"calendar","standard"))
                T=np.array([datetime.fromisoformat(str(x).replace(" ","T")[:19]).replace(tzinfo=timezone.utc).timestamp() for x in t])
            assert u.shape==v.shape==(len(T),481,1201)
            assert len(T)==len(range(start,stop+1,step))
            assert np.allclose(lat[[0,-1]],[80,-80]) and np.allclose(lon[[0,-1]],[20,420])
            assert np.allclose(np.diff(lat),-1/3,atol=2e-5) and np.allclose(np.diff(lon),1/3,atol=4e-5)
            receipt={"indices":list(range(start,stop+1,step)),"sha256":sha(b),"bytes":len(b),"valid_uv":int((np.isfinite(u)&np.isfinite(v)).sum()),"timestamps":[datetime.fromtimestamp(x,timezone.utc).isoformat() for x in T]}
            print("NATIVE_SOURCE_CHUNK",json.dumps(receipt),flush=True)
            return lat,lon,T,u,v,receipt
        except Exception as e:
            last=e
            if k<2:time.sleep(2*(k+1))
    raise RuntimeError(f"Native OSCAR source unavailable for {start}:{step}:{stop}: {last}")
def combine(chunks):
    lat,lon=chunks[0][:2]
    for c in chunks:assert np.array_equal(lat,c[0]) and np.array_equal(lon,c[1])
    T=np.concatenate([c[2] for c in chunks]);u=np.concatenate([c[3] for c in chunks]);v=np.concatenate([c[4] for c in chunks])
    idx=[i for c in chunks for i in c[5]["indices"]]
    assert len(set(idx))==len(idx) and np.all(np.diff(T)>0)
    return lat,lon,T,u,v,idx
def wet(mask,lat,lon):
    r=int(round((85-lat)/.0625-.5));c=int(round((((lon+180)%360)/.0625)-.5))%mask.shape[1]
    # longitude axis of the original fine mask starts at -180, not 0.
    return 0<=r<mask.shape[0] and bool(mask[r,c]>=0)
def gc(x,y):
    return abs(G.inv(x[1],x[0],y[1],y[0])[2])/1000
def launches(lat,lon,u,v,anchor_rows,mask):
    # Preserve all original anchors and the original finite 2x2 interpolation
    # gate; choose starts only by source geometry and current-data availability.
    out={}
    for anchor in anchor_rows:
        a=float(anchor["snapped_wet_latitude"]);b=float(anchor["snapped_wet_longitude"])
        if abs(a)>80:continue
        L=(b-20)%360+20
        rr=np.where(abs(lat-a)<=.6)[0]
        delta=.6/max(.1,math.cos(math.radians(a)))+.35
        cc=np.where(abs(lon-L)<=delta)[0]
        for r in rr:
            for c in cc:
                r=int(r);c=int(c)
                if r>=len(lat)-1 or c>=len(lon)-1:continue
                x=float(lat[r]);y=float((lon[c]+180)%360-180)
                dist=gc((a,b),(x,y))
                if dist>MAX_OFFSET or not wet(mask,x,y):continue
                if not (np.isfinite(u[r:r+2,c:c+2]).all() and np.isfinite(v[r:r+2,c:c+2]).all()):continue
                key=(round(x,6),round(y,6))
                z={"row":r,"col":c,"lat":x,"lon":y,"offset_km":dist}
                if key not in out or dist<out[key]["offset_km"]:out[key]=z
    result=[]
    for z in sorted(out.values(),key=lambda z:(z["offset_km"],z["row"],z["col"])):
        if all(gc((z["lat"],z["lon"]),(x["lat"],x["lon"]))>=MIN_SEP for x in result):
            result.append(z)
            if len(result)==NSTART:break
    return result
def prior_losses(tracks,prep,anchors,sources):
    per=[]
    for s in sorted(sources):
        ts=[q for q in tracks if q["source"]==s]
        p=prep[s];d=p["nearest_valid_offshore_source_km_per_season"]
        z=[]
        for q in ts:
            off=d[[907,919,931].index(q["original_NOAA_source_time_index"])]
            z.append((q,off))
        per.append({"source":s,"original_anchors":len(anchors[s]),"coarse_valid_starts":len(ts),"coarse_starts_le50":sum(v is not None and v<=50 for v in d),"coarse_complete30_le50":sum(q["complete_elapsed_days"]>=30 and off<=50 for q,off in z),"coarse_complete30_le200":sum(q["complete_elapsed_days"]>=30 and off<=200 for q,off in z),"coarse_stop_reasons":dict(Counter(q["status"] for q in ts))})
    return per
def integrate_one(T,lat,lon,U,V,mask,q,idx,i0,helpers):
    current=(q["lat"],q["lon"]);start=current;sec=float(T[i0]);status="COMPLETED_30_DAYS";done=0
    for day in range(1,HORIZON+1):
        now=sec+(day-1)*86400.
        u=helpers.uv_interpolate(T,lat,lon,U,V,now,*current)
        if u is None:status="MISSING_INITIAL_INTERPOLATION";break
        half=helpers.integrate_vector(*current,*u,43200.)
        v=helpers.uv_interpolate(T,lat,lon,U,V,now+43200.,*half)
        if v is None:status="MISSING_MIDPOINT_INTERPOLATION";break
        dest=helpers.integrate_vector(*current,*v,86400.)
        # Same 2km segment resolution as original, periodic longitude indexing.
        hit=False
        n=max(1,int(math.ceil(gc(current,dest)/2.)))
        for j in range(1,n+1):
            w=j/n;la=current[0]+w*(dest[0]-current[0]);lo=(current[1]+w*((dest[1]-current[1]+180)%360-180)+180)%360-180
            if not wet(mask,la,lo):hit=True;break
        if hit:status="GSHHG_LAND_CONTACT_BEFORE_DAY30";break
        current=dest;done=day
    return {"release_index":idx,"start":q,"day30_valid":done==30,"elapsed_days":done,"status":status,"end_lat":current[0],"end_lon":current[1]}
def main():
    p=argparse.ArgumentParser()
    for name in ("traits","paired","tracks","preflight","anchors","graph","out"):p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--original-helper-dir",type=Path)
    p.add_argument("--stage",choices=["preflight","integrate"],required=True)
    p.add_argument("--native-preflight",type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    assert sha(a.graph.read_bytes())==GRAPH_SHA
    with np.load(a.graph) as g:mask=g["fine_cell_to_node"]
    assert mask.shape==(2720,5760)
    ac=json.loads(a.anchors.read_text());anchors=ac["coastal_anchors_by_TDWG_region"]
    assert len(anchors)==126 and sum(bool(x) for x in anchors.values())==104
    pre=json.loads(a.preflight.read_text());prep={q["region"]:q for q in pre["126_source_region_release_eligibility"]}
    tracks=json.loads(a.tracks.read_text());assert len(tracks)==280
    t=pd.read_csv(a.traits);fixed=pd.read_csv(a.paired)
    fixed22=set(fixed.loc[fixed["class"]=="island","code"]);assert len(fixed22)==22
    supported=set(t.loc[t.eligible.astype(str).str.lower().isin(["true","1"]),"source"])
    assert len(supported)==17
    sources=supported
    primary_sources=supported&fixed22
    # Keep the original frozen22 primary distinct from the wider 126 trait-supported sensitivity.
    old=prior_losses(tracks,prep,anchors,sources)
    assert sum(r["coarse_complete30_le50"]>0 for r in old)==4
    assert sum(r["coarse_complete30_le200"]>0 for r in old)==7
    if a.stage=="preflight":
        # Same six dates later used for integration; fetch three at once.
        chunks=[get(907,919,6),get(925,937,6)]
        lat,lon,T,U,V,idx=combine(chunks)
        assert idx==RELEASES
        for i,original_idx in enumerate([907,919,931]):
            assert datetime.fromtimestamp(T[idx.index(original_idx)],timezone.utc).strftime("%Y-%m-%d")==pre["noaa_actual_source_release_dates"][i][:10]
        rows=[]
        for r in old:
            src=r["source"];seasons=[]
            for k,ii in enumerate(idx):
                z=launches(lat,lon,U[k],V[k],anchors[src],mask)
                seasons.append({"release_index":ii,"n_starts":len(z),"starts":z})
            rows.append({**r,"in_fixed22":src in primary_sources,"native_release_dates_with_starts":sum(q["n_starts"]>0 for q in seasons),"native_total_starts":sum(q["n_starts"] for q in seasons),"native_seasons":seasons})
        n=sum(x["native_total_starts"]>0 for x in rows)
        result={"status":"NATIVE_GRID_LAUNCH_COVERAGE_GATE_ONLY","sources":len(sources),"primary_fixed22_trait_supported_sources":len(primary_sources),"coarse_le50_day30_sources":4,"coarse_le200_day30_sources":7,"native_sources_with_launch_le50":n,"native_fixed22_sources_with_launch_le50":sum(x["in_fixed22"] and x["native_total_starts"]>0 for x in rows),"native_sources_with_ge2_dates":sum(x["native_release_dates_with_starts"]>=2 for x in rows),"native_total_starts":sum(x["native_total_starts"] for x in rows),"minimum_sources_for_followup":10,"launch_gate_pass":n>=10,"source_rows":rows,"source_chunks":[q[5] for q in chunks],"no_new_regression":True,"notes":["Native 1/3-degree actual 2005 OSCAR data, same fixed coast anchors, <=50km; six dates, <=3 starts separated >=15km.","Preflight only requires finite launch grid; successful 30day tracking is not yet established.","No new plant records or source data distributed."]}
        (a.out/"NATIVE_PREFLIGHT.json").write_text(json.dumps(result,indent=2)+"\n")
        print("NATIVE_PREFLIGHT_SUMMARY",json.dumps({k:v for k,v in result.items() if k not in ("source_rows","source_chunks")}),flush=True)
        print("NATIVE_SOURCE_COVERAGE",json.dumps([{k:v for k,v in r.items() if k!="native_seasons"} for r in rows]),flush=True)
        return
    initial=json.loads(a.native_preflight.read_text());assert initial["launch_gate_pass"] and initial["native_sources_with_launch_le50"]>=10
    sys.path.insert(0,str(a.original_helper_dir))
    import run_atlantic_noaa_2005_midpoint_particles_20261117 as helpers
    ranges=[(i,min(i+2,LAST)) for i in range(FIRST,LAST+1,3)]
    chunks=[get(lo,hi) for lo,hi in ranges] # netCDF4/HDF access is serialized
    lat,lon,T,U,V,idx=combine(chunks)
    assert idx==list(range(FIRST,LAST+1)) and T[-1]>=T[idx.index(RELEASES[-1])]+HORIZON*86400.
    # Verify every prespecified launch against the separately retrieved preflight.
    rows=[]
    for r in initial["source_rows"]:
        s=r["source"]
        for season in r["native_seasons"]:
            ii=season["release_index"];k=idx.index(ii)
            z=launches(lat,lon,U[k],V[k],anchors[s],mask)
            assert z==season["starts"],(s,ii,"launch selection drift")
            for q in z:
                v=integrate_one(T,lat,lon,U,V,mask,q,ii,k,helpers)
                v["source"]=s;rows.append(v)
    assert len(rows)==initial["native_total_starts"]
    good=[r for r in rows if r["day30_valid"]]
    counts=Counter(q["source"] for q in good)
    report={"status":"NATIVE_GRID_30DAY_PHYSICAL_COVERAGE_DIAGNOSTIC","source_release_dates":RELEASES,"source_chunks":[q[5] for q in chunks],"launches":len(rows),"day30_valid_tracks":len(good),"source_regions_valid30":len(counts),"fixed22_source_regions_valid30":sum(s in primary_sources for s in counts),"source_regions_valid30_ge2_dates":sum(len({q["release_index"] for q in good if q["source"]==s})>=2 for s in counts),"source_valid_track_counts":dict(counts),"termination_states":dict(Counter(q["status"] for q in rows)),"minimum_sources":10,"physical_gate_pass":len(counts)>=10,"no_biological_transport_observed":True}
    (a.out/"NATIVE_TRACKS_INTERNAL.json").write_text(json.dumps(rows,indent=2)+"\n")
    (a.out/"NATIVE_PHYSICAL_RESULT.json").write_text(json.dumps(report,indent=2)+"\n")
    print("NATIVE_PHYSICAL_SUMMARY",json.dumps({k:v for k,v in report.items() if k!="source_chunks"}),flush=True)
if __name__=="__main__":main()
