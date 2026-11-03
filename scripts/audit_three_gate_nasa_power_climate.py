#!/usr/bin/env python3
"""NASA POWER 1981–2025 annual T2M climate-context audit for dated island butterfly literature leads.

Not climate suitability, no species-native niche threshold, no "climate gate"
crossing inferred. WGSRPD3 region-representative points DO NOT reconstruct
individual physical island climate for archipelagos. 1981–2025 means 2026
first records do not have complete-event-year POWER data.
"""
from __future__ import annotations
import argparse,csv,json,re,time,statistics
from pathlib import Path
from collections import defaultdict,Counter
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError,URLError
from shapely.geometry import shape
from shapely.ops import unary_union

BASE="https://power.larc.nasa.gov/api/temporal/monthly/point"
YEARS=(1981,2025)
AGENT="chocho-three-gates-NASA-POWER-screen/0.1 (https://github.com/zuizui0223/chocho)"
def records(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def writecsv(p,rr):
    if not rr:return
    with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]))
        w.writeheader();w.writerows(rr)
def coord(geo,pinned):
    j=json.loads(geo.read_text(encoding="utf-8"))
    shapes=defaultdict(list)
    for f in j["features"]:
        code=str((f.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if code in pinned and f.get("geometry"):
            shapes[code].append(shape(f["geometry"]))
    out={}
    for code in pinned:
        if not shapes[code]:raise RuntimeError("No WGSRPD3 geometry "+code)
        allshape=unary_union(shapes[code])
        p=allshape.representative_point()
        out[code]={"lat":float(p.y),"lon":float(p.x),"bounds":list(allshape.bounds)}
    return out

def nasa(lat,lon):
    q={"parameters":"T2M","community":"AG",
       "latitude":f"{lat:.5f}","longitude":f"{lon:.5f}",
       "start":str(YEARS[0]),"end":str(YEARS[1]),"format":"JSON"}
    url=BASE+"?"+urlencode(q)
    error=None
    for attempt in range(3):
        try:
            with urlopen(Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"}),timeout=65) as r:
                response=json.loads(r.read(3_000_000))
            if "messages" in response and response["messages"]:
                raise RuntimeError("NASA API messages: "+str(response["messages"])[:350])
            t=response["properties"]["parameter"]["T2M"]
            if not isinstance(t,dict):raise RuntimeError("NASA T2M missing")
            month=defaultdict(dict)
            for raw,value in t.items():
                if not re.fullmatch(r"[0-9]{6}",str(raw)):continue
                year,m=int(raw[:4]),int(raw[4:])
                if m<1 or m>12 or year<1981 or year>2025:continue
                try:v=float(value)
                except(TypeError,ValueError):continue
                if v<=-900 or v>80:continue
                month[year][m]=v
            # NASA monthly T2M is mean monthly temperature in Celsius,
            # multiply by actual days to obtain annual daily-weighted mean.
            import calendar
            annual={}
            for year,values in month.items():
                if len(values)!=12:continue
                weights=[calendar.monthrange(year,m)[1] for m in range(1,13)]
                annual[year]=sum(values[m]*d for m,d in zip(range(1,13),weights))/sum(weights)
            if len(annual)<35:
                raise RuntimeError(f"NASA only {len(annual)} full years for {lat},{lon}")
            return annual, {"request":q,"service_header":response.get("header",{}),
                           "source_parameter":"T2M monthly","complete_annual_years":len(annual),
                           "first_complete":min(annual),"last_complete":max(annual)}
        except (HTTPError,URLError,TimeoutError,KeyError,TypeError,ValueError,RuntimeError,json.JSONDecodeError) as e:
            error=str(e)
            time.sleep(2*(attempt+1))
    return {},{"status":"NASA_QUERY_FAILED","detail":error,"request":q}

def baseline(annual):
    dates=[annual[y] for y in range(1981,2001) if y in annual]
    return statistics.mean(dates) if len(dates)>=18 else None
def trend(rows,year):
    if year<1992:return None
    xs=[y for y in range(max(1981,year-10),min(year,2025)+1) if y in rows]
    if len(xs)<8:return None
    avgy=sum(xs)/len(xs)
    ay=sum(rows[y] for y in xs)/len(xs)
    ss=sum((x-avgy)**2 for x in xs)
    return (sum((x-avgy)*(rows[x]-ay) for x in xs)/ss)*10 if ss>0 else None

def main():
    ap=argparse.ArgumentParser()
    for x in ("literature-cases","wgsrpd3-geojson","out-dir"):
        ap.add_argument("--"+x,required=True,type=Path)
    a=ap.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    cases=records(a.literature_cases)
    if len(cases)!=47:raise RuntimeError("Original literature reference panel changed")
    wanted={r["wgsrpd3_code"] for r in cases if r["wgsrpd3_code"]}
    climate_points=coord(a.wgsrpd3_geojson,wanted)
    result={}
    audits={}
    for i,(code,point) in enumerate(sorted(climate_points.items()),1):
        yearly,metadata=nasa(point["lat"],point["lon"])
        result[code]=yearly
        audits[code]={**point,**metadata}
        print("NASA_POWER_T2M",i,len(climate_points),code,len(yearly),flush=True)
        time.sleep(.25)
    yearrows=[]
    for code in sorted(result):
        p=climate_points[code]
        for y,t in sorted(result[code].items()):
            yearrows.append({"wgsrpd3_code":code,"year":y,"T2M_mean_C":round(t,5),
                             "point_latitude":p["lat"],"point_longitude":p["lon"]})
    caserows=[]
    for r in cases:
        code=r["wgsrpd3_code"];first=int(r["first_evidence_earliest_year"])
        last=int(r["first_evidence_latest_year"])
        series=result.get(code,{})
        base=baseline(series)
        prefirst=[series[y] for y in range(first-10,first) if y in series]
        recent=statistics.mean(prefirst) if len(prefirst)>=8 else None
        caserows.append({
           "species":r["species"],"physical_island":r["physical_island"],
           "wgsrpd3_code":code,
           "earliest_butterfly_first_documentation":first,
           "latest_butterfly_first_documentation":last,
           "first_record_after_POWER_start":int(first>=1981),
           "climate_coverage":(
             "PRE_1981_OUTSIDE_POWER" if first<1981 else
             "WGSRPD3_ISLAND_UNMAPPED" if not code else
             "NO_POWER_DATA" if not series else
             "TOO_EARLY_FOR_10YR_WINDOW" if recent is None else
             "TEN_PRE_FIRST_YEARS_AVAILABLE"),
           "region_representative_latitude":climate_points.get(code,{}).get("lat",""),
           "region_representative_longitude":climate_points.get(code,{}).get("lon",""),
           "T2M_baseline_1981_2000_C":round(base,4) if base is not None else "",
           "T2M_mean_decade_before_first_record_C":round(recent,4) if recent is not None else "",
           "T2M_decade_before_minus_1981_2000_C":
             round(recent-base,4) if recent is not None and base is not None else "",
           "T2M_decadal_trend_at_first_record_C_per_decade":
             round(trend(series,first),4) if trend(series,first) is not None else "",
           "original_evidence_grade":r["evidence_grade"],
           "proves_climate_suitability_threshold":0,
        })
    writecsv(a.out_dir/"nasa_power_annual_1981_2025_by_region.csv",yearrows)
    writecsv(a.out_dir/"literature_case_power_temperature_context.csv",caserows)
    info={
       "schema":"chocho_nasa_power_1981_2025_island_firstrecord_climate_context_v0.1",
       "status":"REGION_TEMPERATURE_CONTEXT_NOT_SPECIES_CLIMATE_SUITABILITY",
       "source":"NASA POWER monthly point / T2M / AG; published monthly API documentation",
       "source_url":"https://power.larc.nasa.gov/docs/services/api/temporal/monthly/",
       "original_case_events":47,
       "unique_WGSRPD3_codes":len(wanted),
       "region_points_with_complete_T2M":sum(bool(v) for v in result.values()),
       "region_points_with_api_failure":sum(not bool(v) for v in result.values()),
       "cases_first_documentation_1981plus":sum(int(r["first_evidence_earliest_year"])>=1981 for r in cases),
       "cases_with_pre_first_decade_temperature":sum(r["climate_coverage"]=="TEN_PRE_FIRST_YEARS_AVAILABLE" for r in caserows),
       "cases_pre_POWER_1981":sum(r["climate_coverage"]=="PRE_1981_OUTSIDE_POWER" for r in caserows),
       "cases_missing_region":sum(not r["wgsrpd3_code"] for r in cases),
       "climate_source_requests":audits,
       "limitations":[
         "Point T2M values describe a single representative location in an often multi-island WGSRPD3 polygon; not each physical island or microclimate.",
         "Do not use the 1981-2000 reference to claim event time climate suitability or to explain pre-2000 first reports.",
         "A warming trend without a independently calibrated species-specific baseline threshold cannot be called the last limiting gate.",
         "POWER point data are gridded reanalysis, not original local island temperatures or butterfly microhabitat thermal limits.",
         "The butterfly first-record year does not mean arrival or colonization date.",
         "This exploratory climate screen includes non-independent reports from several physical islands in the same broad WGSRPD3 region."
       ]
    }
    (a.out_dir/"nasa_power_case_climate_audit.json").write_text(
        json.dumps(info,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("POWER_CLIMATE_AUDIT",json.dumps({k:v for k,v in info.items() if k!="climate_source_requests"},ensure_ascii=False))
if __name__=="__main__":
    main()
