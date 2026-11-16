#!/usr/bin/env python3
"""Date-locked 2001–2010 seasonal historical OSCAR flow alignment.

Uses only real legacy NOAA offshore u/v samples along three original
geodesic WPI/TDWG3-verified routes. February and August 15 (yearly) are
fixed BEFORE looking at u/v; represents two seasonal snapshots, NOT yearly
current integral, seeds, ships, or biotic dispersal.
"""
from __future__ import annotations
import argparse,csv,io,json,math,random,time
import urllib.parse,urllib.request,urllib.error
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from statistics import mean,median
from analyze_atlantic_noaa_legacy_2005_current_alignment_20261117 import (
    ERDDAP, LEG_PAIRS, FRACTIONS, choose_points, arc_and_tangent, noaa_indices,
    fetch_sample
)

YEARS=tuple(range(2001,2011))
SEASONS=("02-15","08-15")
AXIS_START,AXIS_END=0,1878

def read_url(url,maxbytes=500000):
    last=None
    for trial in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"chocho-noaa-oscar-date-lock/2026.11",
                   "Accept":"text/csv"})
            with urllib.request.urlopen(req,timeout=55) as rsp:
                if rsp.status!=200:raise RuntimeError(f"Source response {rsp.status}")
                d=rsp.read(maxbytes+1)
                if len(d)>maxbytes:raise RuntimeError("Unexpected response larger than source bound")
                return d.decode("utf-8-sig")
        except (urllib.error.URLError,RuntimeError,TimeoutError) as e:
            last=type(e).__name__
            if trial<3:time.sleep(2*(trial+1))
    raise RuntimeError("NOAA source query exhausted retries "+str(last))

def download_source_time_axis():
    # ERDDAP accepts its physical grid dimension as a single selected variable.
    attempts=("time[0:1:1878]","u[0:1:1878][0][150][390]")
    for expr in attempts:
        try:
            txt=read_url(ERDDAP+urllib.parse.quote(expr,safe="[],:"))
            rr=list(csv.reader(io.StringIO(txt)))
            if len(rr)<1850 or "time" not in rr[0]:continue
            idx=rr[0].index("time")
            start=2 if len(rr)>1 and "UTC" in str(rr[1]) else 1
            times=[]
            for rec in rr[start:]:
                if len(rec)>idx and rec[idx][:4].isdigit():
                    times.append(datetime.fromisoformat(rec[idx].replace("Z","+00:00")))
            if len(times)!=AXIS_END+1:continue
            if times[0].year!=1992 or times[-1].year!=2018:continue
            if len(set(times))!=len(times):raise RuntimeError("Nonunique NOAA source time axis")
            return times,expr
        except (urllib.error.URLError,ValueError,RuntimeError,TimeoutError):
            continue
    raise RuntimeError("Original NOAA time axis not retrievable; cannot assign seasonal query index")

def choose_historic_season_time_indices(axis):
    results=[]
    for year in YEARS:
        for m in SEASONS:
            target=datetime.fromisoformat(f"{year}-{m}T00:00:00+00:00")
            k=min(range(len(axis)),key=lambda i:abs((axis[i]-target).total_seconds()))
            dev=abs((axis[k]-target).total_seconds())/86400
            if dev>4:raise RuntimeError(f"Nearest OSCAR date too distant {target},{axis[k]},{dev}")
            if axis[k].year!=year:raise RuntimeError("Near source index crossed year boundary")
            results.append({"target":target.isoformat(),"year":year,"season":m,
                            "index":k,"source_time":axis[k].isoformat(),
                            "date_offset_days":dev})
    if len({z["index"] for z in results})!=20:raise RuntimeError("Duplicate source time index")
    return results

def read_source_pair(t,lat,lon,expected):
    # Same fixed finite-source schema as 2005 pilot, but permits all source years.
    dims=f"[{t}][0][{lat}][{lon}]"
    expr="u"+dims+",v"+dims
    txt=read_url(ERDDAP+urllib.parse.quote(expr,safe="[],:"),maxbytes=8192)
    rr=list(csv.reader(io.StringIO(txt)))
    if len(rr)!=3 or not all(x in rr[0] for x in ("time","depth","latitude","longitude","u","v")):
        raise RuntimeError(f"Unexpected NOAA CSV source schema {rr[:2]}")
    row=dict(zip(rr[0],rr[2]))
    if row["time"]!=expected.replace("+00:00","Z"):
        raise RuntimeError(f"Time axis mismatch: {row['time']} vs {expected}")
    try:u,v=float(row["u"]),float(row["v"])
    except (ValueError,TypeError):u,v=float("nan"),float("nan")
    if not all(math.isfinite(z) for z in (u,v)) or abs(u)>5 or abs(v)>5:
        return {"status":"MISSING_UV","date":row["time"]}
    return {"status":"FINITE_SOURCE_UV","date":row["time"],
            "u":u,"v":v,"lat":float(row["latitude"]),"lon":float(row["longitude"])}

def one(q):
    p,t=q
    try:
        x=read_source_pair(t["index"],p["ilat"],p["ilon"],t["source_time"])
        out={"status":x["status"],"date":x["date"]}
        if x["status"]=="FINITE_SOURCE_UV":
            a=x["u"]*p["east_component"]+x["v"]*p["north_component"]
            b=-x["u"]*p["north_component"]+x["v"]*p["east_component"]
            if abs(x["lat"]-p["latitude"])>0.18 or abs(((x["lon"]-p["longitude"]+180)%360)-180)>0.18:
                raise RuntimeError("NOAA result does not match predeclared real coordinate")
            out.update({"along_m_s":a,"cross_m_s":b,"u_m_s":x["u"],"v_m_s":x["v"],
                        "source_lat":x["lat"],"source_lon":x["lon"]})
        return out
    except (urllib.error.URLError,RuntimeError,TimeoutError) as err:
        return {"status":"SOURCE_QUERY_FAILED","error_type":type(err).__name__}

def median_bootstrap_year_ci(annual):
    if len(annual)<5:return None
    r=random.Random(20261117)
    vals=[]
    for _ in range(4999):
        ys=r.choices(annual,k=len(annual))
        vals.append(mean(ys))
    vals.sort()
    return [vals[124],vals[4873]]

def main(a):
    source_times,axis_method=download_source_time_axis()
    samples=choose_historic_season_time_indices(source_times)
    anchors=choose_points(a.strict_source_ports)
    # 2005 source byte/units replay stays pinned to older accepted implementation.
    original=fetch_sample(930,150,390)
    if original.get("status")!="SOURCE_UV_FINITE":raise RuntimeError("Old 2005 source replay failed")
    pts=[]
    for start,end in LEG_PAIRS:
        one=anchors[start];two=anchors[end]
        for f in FRACTIONS:
            p=arc_and_tangent((one["longitude"],one["latitude"]),
                  (two["longitude"],two["latitude"]),f)
            ilat,ilon=noaa_indices(p)
            pts.append({"source":start,"destination":end,"fraction":f,**p,"ilat":ilat,"ilon":ilon})
    jobs=[(p,t) for p in pts for t in samples]
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(one,q):q for q in jobs}
        for future in as_completed(futures):
            p,t=futures[future]
            obs=future.result()
            results.append({"source":p["source"],"destination":p["destination"],
              "fraction":p["fraction"],"transect_lat":p["latitude"],
              "transect_lon":p["longitude"],"transit_direction_east":p["east_component"],
              "transit_direction_north":p["north_component"],"distance_km":p["ocean_distance_km"],
              "year":t["year"],"season":t["season"],"time_index":t["index"],
              "source_time_expected":t["source_time"],"lat_index":p["ilat"],
              "lon_index":p["ilon"],**obs})
    results.sort(key=lambda r:(r["source"],r["destination"],r["year"],r["season"],r["fraction"]))
    by_route=defaultdict(list)
    for r in results:by_route[(r["source"],r["destination"])].append(r)
    routes=[]
    for src,dst in LEG_PAIRS:
        group=by_route[(src,dst)]
        slices=[]
        for year in YEARS:
            for season in SEASONS:
                sub=[x for x in group if x["year"]==year and x["season"]==season]
                finite=[x for x in sub if x["status"]=="FINITE_SOURCE_UV"]
                slices.append({"year":year,"season":season,"finite_points":len(finite),
                    "mean_along_m_s":mean(x["along_m_s"] for x in finite) if len(finite)>=2 else None})
        usable=[x["mean_along_m_s"] for x in slices if x["mean_along_m_s"] is not None]
        annual=[]
        for year in YEARS:
            vals=[x["mean_along_m_s"] for x in slices if x["year"]==year and
                  x["mean_along_m_s"] is not None]
            if len(vals)==2:annual.append(mean(vals))
        routes.append({"source":src,"destination":dst,
            "source_WPI_port":anchors[src]["name"],"destination_WPI_port":anchors[dst]["name"],
            "geodesic_distance_km":round(group[0]["distance_km"],2),
            "finite_source_cells":sum(x["status"]=="FINITE_SOURCE_UV" for x in group),
            "total_source_cells":len(group),"usable_time_slices":len(usable),
            "nominal_time_slices":20,"complete_years":len(annual),
            "2001_2010_two_seasons_mean_along_m_s":mean(usable) if usable else None,
            "median_along_season_m_s":median(usable) if usable else None,
            "positive_slices":sum(x>0 for x in usable),"negative_slices":sum(x<0 for x in usable),
            "year_block_bootstrap95_m_s":median_bootstrap_year_ci(annual),
            "sampled_seasons":slices,
            "estimand":"20_FIXED_FEB_AUG_SOURCE_CURRENT_SNAPSHOTS_NOT_ANNUAL_ADVECTION_OR_SEED_TRANSPORT"})
    output={
      "status":"NOAA_LEGACY_2001_2010_FIXED_SEASON_PHYSICAL_CORRIDOR_ONLY",
      "noaa_dataset":"yearly_336c_0b32_9cd3",
      "source_geometry_run":43895414730,
      "2005_original_source_run":43896071940,
      "source_axis_method":axis_method,
      "time_axis_length":len(source_times),"study_time_indices":samples,
      "strict_TDWG3_WPI_anchors":anchors,"nominal_queries":len(jobs),
      "finite_queries":sum(x["status"]=="FINITE_SOURCE_UV" for x in results),
      "routes":routes,"important_limits":[
         "Physical data are old OSCAR at 15 m on a 1/3-degree grid, 5-day, not NASA OSCAR v2 daily quarter-degree.",
         "Sampling mid-February and mid-August in 2001–2010 is not time-resolved particle tracking.",
         "Potential alignment along a great circle is not a modeled propagule arrival probability or coastline flux.",
         "The 3-point coastal port anchor is geographic only, not observed historical shipping OD.",
         "Year bootstrap is descriptive, not an independent current vs seed occurrence statistical test.",
         "No plant source-population history, propagule buoyancy/duration, observed seed delivery or butterfly outcome entered this model."
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"noaa_2001_2010_two_season_alignment.json").write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    with (a.out_dir/"noaa_2001_2010_source_uv_cells.csv").open("w",newline="",encoding="utf8") as f:
        keys=list(results[0].keys())
        for z in results:
            for key in z:
                if key not in keys:keys.append(key)
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(results)
    print("OSCAR_2001_2010_DECADAL_PILOT",json.dumps({"status":output["status"],
         "source_axis_length":len(source_times),"date_method":axis_method,
         "queries":len(jobs),"finite":output["finite_queries"],
         "routes":[{k:v for k,v in r.items() if k!="sampled_seasons"} for r in routes]
         },ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--strict-source-ports",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    main(p.parse_args())
