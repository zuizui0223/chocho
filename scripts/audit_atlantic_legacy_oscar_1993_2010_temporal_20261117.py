#!/usr/bin/env python3
"""Post-pilot temporal check: 1993–2010 legacy OSCAR u/v on the same 9 frozen
Atlantic transect grid cells. Read actual time-series source values, not
Circuitscape electrical currents. This is a physical alignment diagnostic,
NOT transport trajectories, vessel traffic, or species colonization.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
import random
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from analyze_atlantic_noaa_legacy_2005_current_alignment_20261117 import (
    REGIONS, LEG_PAIRS, FRACTIONS, R_KM, arc_and_tangent, choose_points, noaa_indices,
)

ERDDAP = "https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.csv?"
TIME_INDEX_START = 40
TIME_INDEX_END = 1450
YEARS = tuple(range(1993,2011))
MIN_YEARS = 15
MIN_3PT_SLICES_PER_YEAR = 25


def fetch_timeseries(lat_index: int, lon_index: int) -> list[dict]:
    dims = f"[{TIME_INDEX_START}:1:{TIME_INDEX_END}][0][{lat_index}][{lon_index}]"
    expr = "u" + dims + ",v" + dims
    url = ERDDAP + urllib.parse.quote(expr, safe="[],:")
    err = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent":"chocho-historical-OSCAR-directional-source-audit/2026.11",
                "Accept":"text/csv"
            })
            with urllib.request.urlopen(req, timeout=110) as rsp:
                raw = rsp.read(900_001)
                if rsp.status != 200 or len(raw) > 900_000:
                    raise RuntimeError("NOAA returned an unexpected response size/status")
            lines = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
            if len(lines) < 1000:
                raise RuntimeError(f"Too few source time rows: {len(lines)}")
            required={"time","depth","latitude","longitude","u","v"}
            if not required.issubset(lines[0]):
                raise RuntimeError(f"Missing NOAA fields: {set(lines[0])}")
            results=[]
            last_time=""
            for r in lines:
                ts=r["time"]
                if not ts.endswith("Z") or len(ts)<10:
                    continue # skip units and non-datetime CSV rows
                if ts<=last_time:
                    raise RuntimeError("Unsorted/duplicate timestamps in NOAA response")
                last_time=ts
                year=int(ts[:4])
                if not YEARS[0]<=year<=YEARS[-1]:
                    continue
                la=float(r["latitude"]);lo=float(r["longitude"])
                if abs(la-(80-lat_index/3))>0.20 or abs((lo-(20+lon_index/3)+180)%360-180)>0.20:
                    raise RuntimeError("NOAA returned unexpected grid coordinates")
                try:u=float(r["u"]);v=float(r["v"])
                except (ValueError,TypeError):u=float("nan");v=float("nan")
                if math.isfinite(u) and math.isfinite(v) and abs(u)<=5 and abs(v)<=5:
                    results.append({"time":ts,"year":year,"u":u,"v":v,"valid":True})
                else:
                    results.append({"time":ts,"year":year,"valid":False})
            if len(results)<1000:
                raise RuntimeError(f"Insufficient historical 1993–2010 source samples {len(results)}")
            return results
        except (urllib.error.URLError,TimeoutError,OSError,ValueError,RuntimeError) as exc:
            err=type(exc).__name__ + ": " + str(exc)[:180]
            if attempt<3:time.sleep(2*(attempt+1))
    raise RuntimeError(f"NOAA source unavailable after bounded retries: {err}")


def quantile(xs,q):
    values=sorted(xs)
    if not values:return None
    pos=(len(values)-1)*q
    a=int(math.floor(pos));b=int(math.ceil(pos))
    return values[a]+(values[b]-values[a])*(pos-a)


def yearly_bootstrap(values, iterations=3999):
    rng=random.Random(20261117)
    if len(values)<MIN_YEARS:return None
    draws=[]
    for _ in range(iterations):
        s=[values[rng.randrange(len(values))] for _ in values]
        draws.append(statistics.mean(s))
    return [quantile(draws,.025),quantile(draws,.975)]


def run(a):
    chosen=choose_points(a.strict_source_ports)
    positions=[]
    for src,dst in LEG_PAIRS:
        p=chosen[src];q=chosen[dst]
        for fraction in FRACTIONS:
            xy=arc_and_tangent((p["longitude"],p["latitude"]),
                              (q["longitude"],q["latitude"]),fraction)
            lai,loi=noaa_indices(xy)
            positions.append({"source":src,"destination":dst,"fraction":fraction,
                              **xy,"lat_index":lai,"lon_index":loi})
    if len(positions)!=9:raise RuntimeError("Original corridor geometry changed")
    # Nine independent grid points may resolve to duplicates; reuse exact source grid.
    queries=sorted({(p["lat_index"],p["lon_index"]) for p in positions})
    cache={}
    for lat,lon in queries:
        cache[(lat,lon)]=fetch_timeseries(lat,lon)
    cells=[]
    for p in positions:
        r=cache[(p["lat_index"],p["lon_index"])]
        for x in r:
            row={"route":p["source"]+"->"+p["destination"],
                 "source":p["source"],"destination":p["destination"],
                 "fraction":p["fraction"],"year":x["year"],"timestamp":x["time"],
                 "lat_index":p["lat_index"],"lon_index":p["lon_index"],
                 "latitude":p["latitude"],"longitude":p["longitude"],
                 "status":"NOAA_SOURCE_VALID" if x["valid"] else "NOAA_SOURCE_MISSING",
                 "u_m_s":x.get("u"),"v_m_s":x.get("v"),
                 "aligned_m_s":(x["u"]*p["east_component"]+
                                x["v"]*p["north_component"]) if x["valid"] else None,
                 "cross_m_s":(-x["u"]*p["north_component"]+
                               x["v"]*p["east_component"]) if x["valid"] else None}
            cells.append(row)

    annual=[]
    routes=[]
    for src,dst in LEG_PAIRS:
        candidates=[x for x in cells if x["source"]==src and x["destination"]==dst]
        groups=defaultdict(list)
        for x in candidates:
            if x["aligned_m_s"] is not None:
                groups[x["timestamp"]].append(x)
        by_year=defaultdict(list)
        for ts,g in sorted(groups.items()):
            if len(g)!=3:continue
            by_year[int(ts[:4])].append(statistics.mean(x["aligned_m_s"] for x in g))
        for yr in YEARS:
            x=by_year.get(yr,[])
            annual.append({"source":src,"destination":dst,"year":yr,
                           "valid_3pt_samples":len(x),
                           "mean_aligned_m_s":statistics.mean(x) if x else "",
                           "positive_sample_fraction":sum(v>0 for v in x)/len(x) if x else ""})
        valid_years={yr:vs for yr,vs in by_year.items() if len(vs)>=MIN_3PT_SLICES_PER_YEAR}
        if len(valid_years)<MIN_YEARS:
            raise RuntimeError(f"{src}->{dst} has insufficient distinct years {len(valid_years)}")
        year_means=[statistics.mean(vs) for yr,vs in sorted(valid_years.items())]
        all_slices=[v for vs in valid_years.values() for v in vs]
        matching=[x for x in annual if x["source"]==src and x["destination"]==dst and
                   x["valid_3pt_samples"]>=MIN_3PT_SLICES_PER_YEAR]
        route={
           "source":src,"destination":dst,
           "source_WPI_port_id":chosen[src]["WPI_id"],
           "destination_WPI_port_id":chosen[dst]["WPI_id"],
           "distance_km":round(next(p["ocean_distance_km"] for p in positions if p["source"]==src and p["destination"]==dst),2),
           "eligible_year_count":len(valid_years),"complete_3point_time_slices":len(all_slices),
           "mean_of_annual_signed_m_s":statistics.mean(year_means),
           "annual_mean_median_m_s":statistics.median(year_means),
           "annual_year_block_bootstrap_descriptive95_m_s":yearly_bootstrap(year_means),
           "annual_means_positive_years":sum(x>0 for x in year_means),
           "annual_means_negative_years":sum(x<0 for x in year_means),
           "source_fraction_positive_time_slices":sum(x>0 for x in all_slices)/len(all_slices),
           "mean_within_year_fraction_positive":statistics.mean(float(x["positive_sample_fraction"]) for x in matching),
           "interpretation":"SIGNED_OCEAN_CURRENT_ALIGNMENT_NOT_PARTICLE_ARRIVAL"
        }
        routes.append(route)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    for p,records in (
        (a.out_dir/"legacy_1993_2010_source_timeseries_cells.csv",cells),
        (a.out_dir/"legacy_1993_2010_annual_corridor_alignment.csv",annual),
    ):
        with p.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=list(records[0]))
            writer.writeheader();writer.writerows(records)
    audit={
       "status":"HISTORICAL_SIGNED_CURRENT_ALIGNMENT_NOT_SEED_TRANSPORT",
       "experiment_stage":"post-2005-pilot temporal descriptive sensitivity",
       "source_dataset":"NOAA ERDDAP yearly_336c_0b32_9cd3 LEGACY OSCAR",
       "not_NASA_OSCAR_v2":True,"original_2005_run":43896071940,
       "geometry_source_run":43895414730,"EuDiS_source_run":43892007427,
       "time_index_range":[TIME_INDEX_START,TIME_INDEX_END],
       "calendar_years":[YEARS[0],YEARS[-1]],"queries":len(queries),
       "fixed_positions":positions,"verified_anchors":chosen,
       "all_source_cell_time_rows":len(cells),
       "missing_source_cell_time_rows":sum(x["status"]=="NOAA_SOURCE_MISSING" for x in cells),
       "routes":routes,
       "limitations":[
        "This is a fixed 3-point great-circle directional current diagnostic, NOT an integrated Lagrangian trajectory.",
        "ERA/year block bootstrap intervals describe source current variability, NOT plant-arrival effects or independent replicated islands.",
        "Missing/coastal source flow remains missing rather than zero.",
        "Previous 2005 pilot was inspected before the 1993–2010 robustness window and is not an independent pre-registered confirmation.",
        "No historic plant source/sink/taxon introduction, floating seed viability, vessel origin/destination traffic, or butterfly establishment outcome were tested.",
        "Reverse-direction projection is exactly the negative of forward projection and is not an independent falsification null.",
       ],
    }
    (a.out_dir/"legacy_1993_2010_temporal_alignment_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("OSCAR_18YR_HISTORICAL_CORRIDOR_ALIGNMENT",json.dumps({
      "status":audit["status"],"queries":len(queries),"source_rows":len(cells),
      "missing_rows":audit["missing_source_cell_time_rows"],"routes":routes
    },ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--strict-source-ports",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    run(ap.parse_args())
