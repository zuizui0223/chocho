#!/usr/bin/env python3
"""Source-audited Florida cultivated Senna occurrence pilot from public iNaturalist API.

A cultivated plant observation is NOT a confirmed butterfly-host interaction or fitness.
Sampling and taxa were frozen before live observation outcomes (see protocol JSON).
"""
from __future__ import annotations
import argparse,csv,json,math,time,urllib.error,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path

ROOT="https://api.inaturalist.org/v1"
SPECIES=["Senna polyphylla","Senna surattensis","Senna ligustrina","Senna alata"]
PLANT_STATUS={
    "Senna polyphylla":{"id":"2490545","WCVP_FLA":"not_listed"},
    "Senna surattensis":{"id":"2897969","WCVP_FLA":"introduced"},
    "Senna ligustrina":{"id":"2490490","WCVP_FLA":"native"},
    "Senna alata":{"id":"2475546","WCVP_FLA":"introduced"}
}
PLACE_ID=21

def request(endpoint,**params):
    u=ROOT+"/"+endpoint+"?"+urllib.parse.urlencode(params)
    for attempt in range(6):
        try:
            req=urllib.request.Request(u,headers={
                "User-Agent":"chocho-cultivated-senna-data-audit/0.1 (reproducible ecology; contact via github.com/zuizui0223/chocho)",
                "Accept":"application/json"
            })
            with urllib.request.urlopen(req,timeout=50) as response:
                return json.load(response)
        except urllib.error.HTTPError as e:
            if e.code not in (429,500,502,503,504) or attempt==5:
                raise RuntimeError(f"iNaturalist HTTP {e.code} at {endpoint} after {attempt+1} attempts") from e
            time.sleep(min(30,2**attempt+1))
        except (TimeoutError,urllib.error.URLError) as e:
            if attempt==5:raise RuntimeError(f"iNaturalist network error at {endpoint}") from e
            time.sleep(min(30,2**attempt+1))
    raise RuntimeError("inaccessible API")

def one_result(x,label):
    rows=x.get("results") if isinstance(x,dict) else None
    if not isinstance(rows,list):raise RuntimeError(f"missing iNat {label} results")
    return rows

def validate_place():
    res=one_result(request(f"places/{PLACE_ID}"),"place")
    matches=[r for r in res if str(r.get("name") or "").strip().lower()=="florida"]
    if len(matches)!=1:
        raise RuntimeError(f"iNaturalist place {PLACE_ID} is not unambiguously Florida: {str(res)[:250]}")
    r=matches[0]
    disp=str(r.get("display_name") or r.get("name") or "")
    if "florida" not in disp.lower():raise RuntimeError("Florida place label mismatch")
    return {"id":PLACE_ID,"name":r.get("name"),"display_name":disp,"place_type_name":r.get("place_type_name")}

def exact_species(name):
    r=one_result(request("taxa/autocomplete",q=name,per_page=20),"taxon")
    matches=[x for x in r if (x.get("name") or "").strip().lower()==name.lower() and x.get("rank")=="species"]
    if len(matches)!=1:
        raise RuntimeError(f"taxon {name!r} did not resolve to exactly one accepted iNaturalist species: {[x.get('name') for x in r][:15]}")
    x=matches[0]
    return {"name":name,"id":int(x["id"]),"rank":x.get("rank"),"observations_count":x.get("observations_count")}

def valid_observation(o,name,captive):
    if type(o.get("captive")) is not bool or o["captive"]!=captive:
        return None
    taxon=o.get("taxon") or {}
    actual=str(taxon.get("name") or "").strip()
    if not (actual==name or actual.startswith(name+" ")):
        return None
    if not o.get("photos"):return None
    if not o.get("observed_on"):return None
    x=o.get("geojson") or {}
    coords=x.get("coordinates") if isinstance(x,dict) else None
    if not coords or len(coords)!=2:return None
    try:
        lon,lat=[float(v) for v in coords]
        accuracy=o.get("positional_accuracy")
        accuracy=float(accuracy) if accuracy is not None else None
    except (ValueError,TypeError):
        return None
    if not (-90<=lat<=90 and -180<=lon<=180):return None
    if accuracy is None or not (0<=accuracy<=1000):return None
    if o.get("geoprivacy") in ("obscured","private"):return None
    ident=o.get("id")
    if not isinstance(ident,int):return None
    year=str(o["observed_on"])[:4]
    if not (year.isdigit() and 2010<=int(year)<=2025):return None
    # Approximate 5km squares are only an independence screening proxy.
    cell=f"{math.floor(lat/.05)}_{math.floor(lon/.05)}"
    return {
        "observation_id":ident,"source_url":f"https://www.inaturalist.org/observations/{ident}",
        "taxon_requested":name,"taxon_observed":actual,"captive_observation":captive,
        "date":o["observed_on"],"accuracy_m":accuracy,"approx_5km_cell":cell,
        "photo_count":len(o.get("photos") or []),
        "quality_grade":o.get("quality_grade"),
        "cultivated_WCVP_status":PLANT_STATUS[name]["WCVP_FLA"]
    }

def get_observations(taxon_id,captive):
    args={
        "place_id":PLACE_ID,"taxon_id":taxon_id,"captive":str(captive).lower(),
        "d1":"2010-01-01","d2":"2025-12-31","photos":"true",
        "page":1,"per_page":100,"order_by":"observed_on","order":"desc"
    }
    res=request("observations",**args)
    total=res.get("total_results")
    raw=one_result(res,"observations")
    if not isinstance(total,int):raise RuntimeError("missing total_results")
    return total,raw

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--protocol-json",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    p.add_argument("--output-csv",type=Path,required=True)
    a=p.parse_args()
    protocol=json.loads(a.protocol_json.read_text(encoding="utf8"))
    if [r["taxon"] for r in protocol["panel"]["taxa"]]!=SPECIES:
        raise RuntimeError("panel drift from preregistered taxon identities")
    for x in protocol["panel"]["taxa"]:
        name=x["taxon"]
        if x["wcvp_accepted_id"]!=PLANT_STATUS[name]["id"]:
            raise RuntimeError("WCVP identifier drift")
    place=validate_place()
    records=[];audit=[];taxon_lookup={}
    for name in SPECIES:
        tax=exact_species(name)
        taxon_lookup[name]=tax
        for captive in (True,False):
            total,rows=get_observations(tax["id"],captive)
            valid=[]
            for o in rows:
                x=valid_observation(o,name,captive)
                if x is not None:valid.append(x)
            by_id={r["observation_id"]:r for r in valid}
            valid=list(by_id.values())
            cells={r["approx_5km_cell"] for r in valid}
            audit.append({
                "plant_name":name,"frozen_WCVP_Florida_status":PLANT_STATUS[name]["WCVP_FLA"],
                "inat_taxon_id":tax["id"],"captive":captive,
                "api_total_results":total,"sample_returned":len(rows),
                "strict_source_confirmed_in_sample":len(valid),
                "distinct_approx_5km_cells":len(cells),
                "sample_is_exhaustive":total<=len(rows),
                "query_completed":True
            })
            records.extend(valid)
            print(json.dumps({"name":name,"captive":captive,"total":total,"valid_sample":len(valid),"cells":len(cells)}),flush=True)
            time.sleep(.7)
    for name in SPECIES:
        a_total=next(x for x in audit if x["plant_name"]==name and x["captive"])
        b_total=next(x for x in audit if x["plant_name"]==name and not x["captive"])
        print(json.dumps({"taxon":name,"cultivated_sample":a_total["strict_source_confirmed_in_sample"],"wild_sample":b_total["strict_source_confirmed_in_sample"]}),flush=True)
    primary=next(x for x in audit if x["plant_name"]=="Senna polyphylla" and x["captive"])
    pass_gate=primary["strict_source_confirmed_in_sample"]>=3 and primary["distinct_approx_5km_cells"]>=2
    out={
       "schema":"chocho_cultivated_senna_exposure_pilot_v0.1",
       "status":"EXPLORATORY_INDEPENDENT_CULTIVATED_PLANT_OCCURRENCE_AUDIT",
       "protocol":"docs/exploratory/CULTIVATED_SENNA_EXPOSURE_PROTOCOL_V01.json",
       "place":place,"species_taxonomy":taxon_lookup,"observations_audited":audit,
       "primary_gate":{"target":"Senna polyphylla",
           "WCVP_Florida_status":"not_listed",
           "minimum_exact_cultivated_photo_geotag_observations":3,
           "minimum_approx_5km_cells":2,
           "passed":pass_gate},
       "method_limitations":[
         "iNaturalist has volunteer-biased botanical reports and incomplete planted/wild labels. Captive count is NOT horticultural abundance.",
         "Only recent 100 photo observations for each species/status are source-sampled, with public accurate GPS. Coverage cannot establish zeros or proportions for real plants.",
         "Self-annotated captive flag, photos and IDs are not externally adjudicated as planted.",
         "Individual plant observations do not identify butterflies feeding or successfully reproducing at those sites.",
         "Koptur et al. 2024 observed adult emergence on planted S.polyphylla in Miami, but this source sample is independent and need not match those precise sites.",
         "This does not prove S.polyphylla is or is not naturalized in Florida; WCVP absence is only absence from fixed regional inventory.",
         "No claim about butterfly habitat rescue, predation, colonization or geographic generality."
       ],
       "sampled_source_observations":len(records)
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["observation_id","source_url","taxon_requested","taxon_observed","captive_observation","date","accuracy_m","approx_5km_cell","photo_count","quality_grade","cultivated_WCVP_status"]
    with a.output_csv.open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows(sorted(records,key=lambda x:(x["taxon_requested"],str(x["captive_observation"]),x["date"],x["observation_id"])))
    print(json.dumps({"primary_gate_passed":pass_gate,"total_strict_sample_records":len(records)}),flush=True)

if __name__=="__main__":main()
