#!/usr/bin/env python3
"""NASA POWER 1981-2020 thermal-baseline feasibility for fixed island-region cases.

This is a SOURCE/DATE AUDIT, NOT butterfly niche suitability or timing of
resource-gate opening. Region points are within the official WGSRPD3 shape,
not necessarily the actual physical island first-record site.
"""
from __future__ import annotations
import argparse,csv,json,calendar,statistics,time,urllib.parse,urllib.request
from pathlib import Path
from shapely.geometry import shape
from collections import defaultdict

API="https://power.larc.nasa.gov/api/temporal/monthly/point"
CODES=("JAM","PUE","NZN","CUB","SEY")
PERIODS=((1981,2000),(2001,2020))
AGENT="chocho-three-gates-climate-coverage/0.1"

def get_json(params):
    url=API+"?"+urllib.parse.urlencode(params)
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"}),timeout=90) as response:
        if response.status!=200:raise RuntimeError(f"NASA POWER status {response.status}")
        return json.load(response),url

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--wgsrpd3-geojson",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    features=json.loads(a.wgsrpd3_geojson.read_text())["features"]
    sites={}
    for f in features:
        code=str((f.get("properties") or {}).get("LEVEL3_COD") or "")
        if code in CODES:
            pt=shape(f["geometry"]).representative_point()
            sites[code]=(float(pt.x),float(pt.y))
    if set(sites)!=set(CODES):raise RuntimeError("Missing frozen island-region climate coordinates")
    annual_rows=[];audit=[]
    for idx,code in enumerate(CODES):
        lon,lat=sites[code]
        params={"parameters":"T2M","community":"AG","longitude":f"{lon:.5f}",
                "latitude":f"{lat:.5f}","format":"JSON","start":"1981","end":"2020"}
        try:
            reply,url=get_json(params)
            if reply.get("messages") and any("error" in str(x).lower() for x in reply["messages"]):
                raise RuntimeError(f"NASA reports error: {reply['messages']}")
            values=reply.get("properties",{}).get("parameter",{}).get("T2M",{})
            if not isinstance(values,dict):raise RuntimeError("No NASA T2M monthly series")
            valid=defaultdict(list)
            for yyyymm,value in values.items():
                s=str(yyyymm)
                if len(s)!=6 or not s.isdigit():continue
                y,m=int(s[:4]),int(s[4:])
                if y<1981 or y>2020 or not 1<=m<=12:continue
                try:v=float(value)
                except (TypeError,ValueError):continue
                if v<=-990 or not -90<=v<=70:continue
                days=calendar.monthrange(y,m)[1]
                valid[y].append((m,v,days))
            if any(len(valid[y])<12 for y in range(1981,2021)):
                miss=[y for y in range(1981,2021) if len(valid[y])<12]
                raise RuntimeError(f"NASA incomplete monthly T2M: {miss[:20]}")
            year_temp={}
            for y in range(1981,2021):
                annual=sum(x[1]*x[2] for x in valid[y])/sum(x[2] for x in valid[y])
                year_temp[y]=annual
                annual_rows.append({"area_code_l3":code,"year":y,
                                    "annual_mean_T2M_C":annual,
                                    "n_months":len(valid[y]),"longitude":lon,"latitude":lat})
            a1981=statistics.mean(year_temp[y] for y in range(1981,2001))
            b2001=statistics.mean(year_temp[y] for y in range(2001,2021))
            audit.append({"code":code,"query_status":"COMPLETE_1981_2020",
                          "longitude":lon,"latitude":lat,
                          "years_complete":len(year_temp),
                          "T2M_1981_2000_C":a1981,
                          "T2M_2001_2020_C":b2001,
                          "delta_recent_minus_early_C":b2001-a1981,
                          "api_url":url})
        except Exception as e:
            audit.append({"code":code,"query_status":"FAILED_SOURCE",
                          "longitude":lon,"latitude":lat,"error":str(e)})
        print("NASA_POWER",code,audit[-1]["query_status"],flush=True)
        time.sleep(1)
    if annual_rows:
        with (a.out_dir/"power_region_annual_T2M_1981_2020.csv").open("w",newline="",encoding="utf8") as h:
            wr=csv.DictWriter(h,fieldnames=list(annual_rows[0]));wr.writeheader();wr.writerows(annual_rows)
    result={"schema":"chocho_three_gate_nasa_power_monthly_temp_1981_2020_pilot_v0.1",
            "source":"NASA POWER monthly T2M at official WGSRPD3 representative points",
            "source_docs":"https://power.larc.nasa.gov/docs/services/api/temporal/monthly/",
            "regions":audit,
            "complete_region_count":sum(r["query_status"]=="COMPLETE_1981_2020" for r in audit),
            "limitations":[
             "This does not estimate species thermal suitability or the date when a butterfly was physiologically able to establish.",
             "Static representative climate points are a coarse sample of large WGSRPD3 islands, not individual colony or port localities.",
             "NASA POWER meteorology starts in 1981 and cannot retrospectively date 19th-century climatic suitability.",
             "Two-period temperature change is an observational trend, not a threshold or butterfly colonization mechanism.",
             "Species-specific seasonal minimum temperature, larvae degree days, overwinter survival, climatic niche and spatial microclimate remain unspecified."
            ]}
    (a.out_dir/"power_thermal_gate_feasibility_audit.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("NASA_POWER_AUDIT",json.dumps(result,ensure_ascii=False))
if __name__=="__main__":
    main()
