#!/usr/bin/env python3
"""Physical source-only pilot: historic OSCAR signed current alignment among three
predeclared Atlantic island regions, with exact WPI-in-TDWG3 anchors.

NOT a seed trajectory, propagule path, butterfly dispersal, shipping route,
invasion origin, or establishment model. Reverse signs are mathematical.
"""
from __future__ import annotations
import argparse, csv, io, json, math, random, time
import urllib.parse, urllib.request, urllib.error
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from statistics import mean,median

ERDDAP="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.csv?"
DATASET="yearly_336c_0b32_9cd3"
WPI_GEOJSON_SHA256="62159ee5ce4a31aff90ee46efdb3cdf11a9b368472a147392eb742148f340d83"
# All fixed BEFORE looking at currents or butterflies.
REGIONS=("AZO","CNY","CVI")
LEG_PAIRS=(("AZO","CNY"),("CNY","CVI"),("AZO","CVI"))
PREFERRED_IDS={"AZO":("2395","2370","3584","3588"),
               "CNY":("428","429","427","431"),
               "CVI":("396",)}
TIME_INDICES=(888,900,912,924,936,948)
FRACTIONS=(0.25,0.5,0.75)
R_KM=6371.0088

def rows(path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def vec(lon,lat):
    lo=math.radians(lon);la=math.radians(lat)
    return (math.cos(la)*math.cos(lo),math.cos(la)*math.sin(lo),math.sin(la))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def norm(a):return math.sqrt(dot(a,a))
def normalized(a):
    n=norm(a)
    if n<1e-12:raise ValueError("zero direction vector")
    return tuple(x/n for x in a)
def arc_and_tangent(a,b,f):
    x=vec(*a);y=vec(*b)
    omega=math.acos(max(-1,min(1,dot(x,y))))
    if not 1e-6<omega<math.pi-1e-6:raise ValueError("Invalid geodesic distance")
    # Spherical linear interpolation of great-circle coordinates.
    p=tuple((math.sin((1-f)*omega)*x[i]+math.sin(f*omega)*y[i])/math.sin(omega)
            for i in range(3))
    p=normalized(p)
    lon=math.degrees(math.atan2(p[1],p[0]))
    lat=math.degrees(math.atan2(p[2],math.hypot(p[0],p[1])))
    # Tangent pointing from current point toward island B, not a fixed-angle bearing.
    t=normalized(tuple(y[i]-dot(y,p)*p[i] for i in range(3)))
    lo=math.radians(lon);la=math.radians(lat)
    east=(-math.sin(lo),math.cos(lo),0.0)
    north=(-math.sin(la)*math.cos(lo),-math.sin(la)*math.sin(lo),math.cos(la))
    te,tn=dot(t,east),dot(t,north)
    if abs(math.hypot(te,tn)-1)>1e-7:raise ValueError("Local directional basis invalid")
    return {"longitude":lon,"latitude":lat,
      "east_component":te,"north_component":tn,"ocean_distance_km":omega*R_KM}

def noaa_indices(point):
    lat=point["latitude"];lon=point["longitude"]
    if not -80<=lat<=80:raise ValueError("Out of legacy OSCAR latitude range")
    lati=round((80-lat)*3)
    lon360=lon if lon>=20 else lon+360
    loni=round((lon360-20)*3)
    if not 0<=lati<=480 or not 0<=loni<=1200:raise ValueError("Out of legacy OSCAR grid")
    return lati,loni

def fetch_sample(t,lat,lon):
    dims=f"[{t}][0][{lat}][{lon}]"
    expr="u"+dims+",v"+dims
    url=ERDDAP+urllib.parse.quote(expr,safe="[],")
    last_error=None
    for trial in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"chocho-historic-oscar-source-only/2026.11"})
            with urllib.request.urlopen(req,timeout=50) as resp:
                if resp.status!=200:raise RuntimeError(f"Unexpected status {resp.status}")
                text=resp.read(8192).decode("utf-8-sig")
            rr=list(csv.reader(io.StringIO(text)))
            if len(rr)!=3:raise RuntimeError(f"Unexpected ERDDAP CSV rows {len(rr)}")
            fields=rr[0]
            if not all(f in fields for f in ("time","depth","latitude","longitude","u","v")):
                raise RuntimeError(f"Missing source velocity columns: {fields}")
            row=dict(zip(fields,rr[2]))
            year=int(row["time"][:4])
            if year!=2005:raise RuntimeError("Fixed 2005 time index did not return 2005")
            def number(s):
                try:return float(s)
                except:return float("nan")
            u,v=number(row["u"]),number(row["v"])
            out={"timestamp":row["time"],"longitude":float(row["longitude"]),
                 "latitude":float(row["latitude"]),"u_m_s":u,"v_m_s":v,
                 "source_time_index":t,"source_lat_index":lat,"source_lon_index":lon}
            if not all(math.isfinite(z) for z in [u,v]) or abs(u)>5 or abs(v)>5:
                out["status"]="MISSING_OR_NONFINITE_VELOCITY"
            else:out["status"]="SOURCE_UV_FINITE"
            return out
        except (urllib.error.URLError,urllib.error.HTTPError,RuntimeError,TimeoutError) as exc:
            last_error=type(exc).__name__
            if trial<3: time.sleep((trial+1)*2)
    return {"status":"NOAA_QUERY_FAILED","error_type":last_error,
            "source_time_index":t,"source_lat_index":lat,"source_lon_index":lon}

def choose_points(source_csv):
    by_code=defaultdict(list)
    for r in rows(source_csv):
        if r["WGSRPD3"] in REGIONS and r["original_WPI_2019_coordinate_inside_exact_TDWG3_polygon"]=="1":
            by_code[r["WGSRPD3"]].append(r)
    chosen={}
    for code in REGIONS:
        choices=by_code.get(code,[])
        if not choices:raise RuntimeError(f"No strictly polygon-verified WPI port for {code}")
        ids={r["WPI_source_OBJECTID"]:r for r in choices}
        prioritized=next((pid for pid in PREFERRED_IDS[code] if pid in ids),None)
        if prioritized:
            r=ids[prioritized];rule="PREDECLARED_HARBOR_PRIORITY"
        else:
            r=sorted(choices,key=lambda q:int(q["WPI_source_OBJECTID"]))[0]
            rule="LOWEST_WPI_ID_FALLBACK_NOT_BIOLOGICAL_SELECTION"
        chosen[code]={
            "WGSRPD3":code,"WPI_id":r["WPI_source_OBJECTID"],
            "name":r["port_name"],"country":r["country"],
            "longitude":float(r["longitude"]),"latitude":float(r["latitude"]),
            "anchor_selection_rule":rule,"verified_in_polygons":True,
            "other_eligible_ports":len(choices)-1}
    return chosen

def run(a):
    anchors=choose_points(a.strict_source_ports)
    # NOAA 2005 2-variable payload must reproduce a separately archived source readout.
    replay=fetch_sample(930,150,390)
    if replay.get("status")!="SOURCE_UV_FINITE":
        raise RuntimeError(f"NOAA two-variable payload source control failed: {replay}")
    if abs(replay["u_m_s"]-.11842477321624756)>1e-6 or abs(replay["v_m_s"]+.058310359716415405)>1e-6:
        raise RuntimeError(f"NOAA source replay does not match prior independent u/v: {replay}")

    grid=[]
    for source,dest in LEG_PAIRS:
        one=anchors[source];two=anchors[dest]
        for fraction in FRACTIONS:
            p=arc_and_tangent((one["longitude"],one["latitude"]),
                (two["longitude"],two["latitude"]),fraction)
            lati,loni=noaa_indices(p)
            grid.append({"source":source,"destination":dest,"fraction":fraction,
                  **p,"source_lat_index":lati,"source_lon_index":loni})
    if len(grid)!=9:raise ValueError("Not the original 3*3 transect design")
    jobs=[(p,t) for p in grid for t in TIME_INDICES]
    samples=[]
    with ThreadPoolExecutor(max_workers=3) as executor:
        fut={executor.submit(fetch_sample,t,p["source_lat_index"],
                              p["source_lon_index"]):(p,t) for p,t in jobs}
        for f in as_completed(fut):
            p,t=fut[f]
            o=f.result()
            point={**p,**o,"study_time_index":t}
            if o["status"]=="SOURCE_UV_FINITE":
                # The same local point and source grid must be close enough.
                if abs(o["latitude"]-p["latitude"])>0.18 or abs(((o["longitude"]-p["longitude"]+180)%360)-180)>0.18:
                    raise RuntimeError("NOAA index does not return intended real lat/lon coordinate")
                along=o["u_m_s"]*p["east_component"]+o["v_m_s"]*p["north_component"]
                cross=-o["u_m_s"]*p["north_component"]+o["v_m_s"]*p["east_component"]
                point.update({"aligned_m_s":along,"cross_m_s":cross,
                    "opposite_direction_aligned_m_s":-along})
            samples.append(point)
    samples.sort(key=lambda p:(p["source"],p["destination"],p["study_time_index"],p["fraction"]))

    groups=[]
    for src,dest in LEG_PAIRS:
        cells=[r for r in samples if r["source"]==src and r["destination"]==dest]
        good=[r for r in cells if r["status"]=="SOURCE_UV_FINITE"]
        times=sorted({r["timestamp"] for r in good})
        by_time=defaultdict(list)
        for r in good:by_time[r["timestamp"]].append(r["aligned_m_s"])
        complete_times=sorted(t for t,z in by_time.items() if len(z)>=2)
        time_means=[mean(by_time[t]) for t in complete_times]
        groups.append({
           "source":src,"destination":dest,
           "anchor_port_a":anchors[src]["name"],"anchor_port_b":anchors[dest]["name"],
           "distance_km":round(cells[0]["ocean_distance_km"],2),
           "source_cells":len(cells),"finite_cells":len(good),
           "available_months":sorted({t[:7] for t in times}),
           "complete_time_slices":len(complete_times),
           "mean_along_current_m_s":mean(time_means) if time_means else None,
           "median_along_current_m_s":median(time_means) if time_means else None,
           "positive_time_slices":sum(x>0 for x in time_means),
           "negative_time_slices":sum(x<0 for x in time_means),
           "within_season_directional_reversal":min(time_means,default=0)<0<max(time_means,default=0),
           "time_slice_values":[{"time":t,"n_points":len(by_time[t]),
                   "mean_along_m_s":mean(by_time[t])} for t in complete_times],
           "inference":"ALONG_CORRIDOR_OCEAN_CURRENT_ALIGNMENT_NOT_SEED_TRANSFER_OR_BUTTERFLY_ARRIVAL",
        })
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"source_2005_signed_uv_corridor_cells.csv").open("w",newline="",encoding="utf8") as f:
        colnames=list(samples[0])
        w=csv.DictWriter(f,fieldnames=colnames,extrasaction="ignore")
        w.writeheader();w.writerows(samples)
    out={
        "status":"NOAA_OSCAR_LEGACY_DIRECTIONAL_FLOW_ALIGNMENT_ONLY",
        "source_dataset":DATASET,
        "source_resolution":"legacy 5-day 1/3-degree 15-m, NOT NASA OSCAR v2.0",
        "source_year":2005,"time_indices":list(TIME_INDICES),
        "fractions":list(FRACTIONS),"strict_source_geometry_run":43895414730,
        "original_eudis_source_run":43892007427,
        "wpi_source_geojson_sha256":WPI_GEOJSON_SHA256,
        "anchor_selection_before_seeing_current":True,
        "anchors":anchors,"independent_source_sample_replay":replay,
        "routes":groups,"total_source_queries":len(jobs),
        "total_source_finite":sum(r["status"]=="SOURCE_UV_FINITE" for r in samples),
        "key_limits":[
          "WPI physical point is only a coastal geographic anchor, not a ship route.",
          "No actual current trajectory, drift duration, seed survival or water crossing simulated.",
          "Signed source->destination versus reverse is deterministically negative, not an independent null test.",
          "2005 limited time-slices do not establish long-term prevailing currents.",
          "No plant species-specific source availability or biological outcome is tested.",
          "Original WCVP native/alien and EuDiS morphological trait distributions remain separate, biased observational material.",
        ]
    }
    (a.out_dir/"corridor_current_alignment_2005.json").write_text(
      json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print("NOAA_DIRECTIONAL_CORRIDOR_ALIGN_SOURCE",json.dumps({
       "status":out["status"],"source_points":len(grid),"total_queries":len(jobs),
       "total_finite":out["total_source_finite"],
       "anchors":anchors,"routes":groups},ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--strict-source-ports",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    run(ap.parse_args())
