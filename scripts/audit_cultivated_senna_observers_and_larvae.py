#!/usr/bin/env python3
"""iNaturalist source audit: cultivated Senna independence vs actual Phoebis larval feeding."""
from __future__ import annotations
import argparse,csv,json,math,re,time,urllib.parse,urllib.request,urllib.error
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path
BASE="https://api.inaturalist.org/v1"
TARGET="Senna polyphylla"
BUTTERFLIES=("Phoebis sennae","Phoebis philea")
SENNA=("Senna polyphylla","Senna surattensis","Senna ligustrina","Senna alata")
def request(path,**params):
    url=BASE+"/"+path+("?"+urllib.parse.urlencode(params) if params else "")
    for attempt in range(6):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"chocho-florida-senna-source-audit/0.2 (github.com/zuizui0223/chocho)"})
            with urllib.request.urlopen(req,timeout=45) as f:return json.load(f)
        except urllib.error.HTTPError as e:
            if e.code not in (429,500,502,503,504) or attempt==5:raise
        except (TimeoutError,urllib.error.URLError):
            if attempt==5:raise
        time.sleep(min(25,2**attempt+1))
    raise RuntimeError("API retry exhausted")
def payload_rows(response):
    rows=response.get("results")
    if not isinstance(rows,list):raise RuntimeError("missing results")
    return rows
def chunked(items,n):
    for i in range(0,len(items),n):yield items[i:i+n]
def matching_taxon(observation,name):
    a=str((observation.get("taxon") or {}).get("name") or "")
    return a==name or a.startswith(name+" ")
def geo_point(obs):
    g=obs.get("geojson") or {}
    q=g.get("coordinates") if isinstance(g,dict) else None
    if q is None or len(q)!=2:return None
    try:return float(q[1]),float(q[0])
    except (ValueError,TypeError):return None
def source_consistent(obs,row):
    if not matching_taxon(obs,row["taxon_requested"]):return "taxon_changed"
    if str(obs.get("observed_on") or "")!=row["date"]:return "date_changed"
    if obs.get("captive") is not (row["captive_observation"]=="True"):return "cultivation_flag_changed"
    if not obs.get("photos"):return "photo_removed"
    if geo_point(obs) is None:return "no_public_coordinates"
    try:a=float(obs.get("positional_accuracy"))
    except (TypeError,ValueError):return "accuracy_not_available"
    if not 0<=a<=1000:return "accuracy_too_low"
    return "consistent"
def get_fields(o):
    unique=set()
    for collection in ("ofvs","observation_fields"):
        for field in o.get(collection) or []:
            if not isinstance(field,dict):continue
            name=str(field.get("name") or (field.get("observation_field") or {}).get("name") or "")
            value=str(field.get("value") or "")
            if name or value:unique.add((name,value))
    return sorted(unique)
def source_trophic_candidates(o,plant_ids):
    hits=[];other=[]
    for name,value in get_fields(o):
        label=name.casefold()
        if not any(x in label for x in ("feeding on","eating","herbivore of","associated observation","partner observation")):
            continue
        plant=TARGET.casefold() in value.casefold() or bool(re.search(r"(?<![0-9])335152(?![0-9])",value))
        partner=any(re.search(r"(?:observations/|observation\s*#)"+str(pid)+r"(?![0-9])",value,re.I) for pid in plant_ids)
        if plant or partner:hits.append((name,value[:180]))
        elif any(sp.casefold() in value.casefold() for sp in SENNA):other.append((name,value[:180]))
    return hits,other
def haversine_km(a,b):
    lat1,lon1=map(math.radians,a);lat2,lon2=map(math.radians,b)
    h=math.sin((lat2-lat1)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 12742.0176*math.asin(min(1,math.sqrt(h)))
def describe(rows):
    users=Counter(x["user"] for x in rows if x["user"] is not None)
    years={x["year"] for x in rows}
    cells={x["cell"] for x in rows}
    pairs={(x["user"],x["cell"]) for x in rows if x["user"] is not None}
    n=len(rows);max_n=max(users.values(),default=0)
    effective=1/sum((v/n)**2 for v in users.values()) if n and users else 0.0
    return {"source_consistent":n,"observers":len(users),"largest_observer_records":max_n,
            "largest_observer_fraction":max_n/n if n else None,
            "effective_observers_inverse_simpson":effective,"approx_5km_cells":len(cells),
            "observer_cell_pairs":len(pairs),"years":len(years),
            "year_range":[min(years),max(years)] if years else [],
            "gate_pass":len(users)>=10 and len(years)>=4 and len(pairs)>=10 and max_n<=n/2}
def exact_taxon(name):
    r=payload_rows(request("taxa/autocomplete",q=name,per_page=30))
    matches=[x for x in r if str(x.get("name") or "").casefold()==name.casefold() and x.get("rank")=="species"]
    if len(matches)!=1:raise RuntimeError(f"ambiguous iNat taxon {name}")
    return int(matches[0]["id"])
def main():
    ap=argparse.ArgumentParser()
    for p in ("source_csv","protocol_json","output_json","output_candidates_csv"):
        ap.add_argument("--"+p.replace("_","-"),type=Path,required=True)
    a=ap.parse_args()
    protocol=json.loads(a.protocol_json.read_text())
    if tuple(protocol["consumer_specificity"]["butterfly_species"])!=BUTTERFLIES:raise RuntimeError("butterfly panel drift")
    with a.source_csv.open(newline="",encoding="utf8") as f:frozen=list(csv.DictReader(f))
    if len(frozen)!=573:raise RuntimeError(f"expected 573 original observations, got {len(frozen)}")
    original={int(r["observation_id"]):r for r in frozen}
    if len(original)!=573:raise RuntimeError("duplicate frozen IDs")
    if sum(r["taxon_requested"]==TARGET and r["captive_observation"]=="True" for r in frozen)!=59:
        raise RuntimeError("fixed 59 cultivated Senna changed")
    observed={}
    ids=sorted(original)
    for batch in chunked(ids,30):
        payload=request("observations/"+",".join(map(str,batch)))
        for obs in payload_rows(payload):
            if obs.get("id") in original:observed[obs["id"]]=obs
        time.sleep(.65)
    consistent=defaultdict(list);exclusions=Counter()
    for oid,row in original.items():
        obs=observed.get(oid)
        if obs is None:exclusions["no_longer_accessible"]+=1;continue
        reason=source_consistent(obs,row)
        if reason!="consistent":exclusions[reason]+=1;continue
        point=geo_point(obs);usr=(obs.get("user") or {}).get("id")
        consistent[(row["taxon_requested"],row["captive_observation"])].append({
            "oid":oid,"user":usr if isinstance(usr,int) else None,"cell":row["approx_5km_cell"],
            "year":int(row["date"][:4]),"day":row["date"],"point":point})
    groups=[{"plant":sp,"cultivated":c=="True",**describe(consistent[(sp,c)])}
        for sp in SENNA for c in ("True","False")]
    planted=consistent[(TARGET,"True")]
    if not planted:raise RuntimeError("no target plant records survive source audit")
    plant_ids={x["oid"] for x in planted}
    planted_dates=[(x,date.fromisoformat(x["day"])) for x in planted]
    print(json.dumps({"plant_source_verified":len(observed),"target":describe(planted),"exclusions":dict(exclusions)}),flush=True)
    larvae=[];candidates=[];matched_nearness=set()
    for butterfly in BUTTERFLIES:
        taxon_id=exact_taxon(butterfly)
        total=None;scanned=0;has_fields=0;direct=set();other_host=set()
        for page in range(1,6):
            params={"place_id":21,"taxon_id":taxon_id,"d1":"2010-01-01","d2":"2025-12-31",
                    "photos":"true","captive":"false","term_id":1,"term_value_id":6,
                    "per_page":200,"page":page,"order_by":"observed_on","order":"desc"}
            result=request("observations",**params)
            n=result.get("total_results")
            if not isinstance(n,int):raise RuntimeError("larval total missing")
            if total is None:total=n
            rows=payload_rows(result)
            if not rows:break
            for obs in rows:
                if not matching_taxon(obs,butterfly) or not obs.get("photos"):continue
                ident=obs.get("id")
                if not isinstance(ident,int):continue
                scanned+=1
                fs=get_fields(obs)
                if fs:has_fields+=1
                good,other=source_trophic_candidates(obs,plant_ids)
                if good:
                    direct.add(ident)
                    for name,value in good:
                        candidates.append({"butterfly":butterfly,"larval_observation_id":ident,
                            "source_url":f"https://www.inaturalist.org/observations/{ident}",
                            "field_name":name,"field_value":value,
                            "status":"OBSERVER_ASSERTION_REQUIRES_PHOTO_ADJUDICATION"})
                if other:other_host.add(ident)
                pt=geo_point(obs);dstr=obs.get("observed_on")
                if pt is None or not dstr:continue
                try:day=date.fromisoformat(dstr)
                except ValueError:continue
                if any(abs((day,plantdate)[0].toordinal()-plantdate.toordinal())<=365 and
                       haversine_km(pt,p["point"])<=2 for p,plantdate in planted_dates):
                    matched_nearness.add((butterfly,ident))
            if len(rows)<200 or page*200>=total:break
            time.sleep(.8)
        larvae.append({"species":butterfly,"iNat_taxon_id":taxon_id,
            "API_total_annotated_larvae":total,"photographed_annotated_exact_species_scanned":scanned,
            "fully_scanned":total<=1000,"any_observation_fields":has_fields,
            "original_field_link_to_S_polyphylla":len(direct),
            "original_field_link_to_other_Senna":len(other_host),
            "warning":"A zero direct-link count is not evidence of zero larval feeding. Photos and observer annotations require source audit."})
        print(json.dumps(larvae[-1]),flush=True)
    summary={"schema":"chocho_cultivated_senna_source_and_larval_audit_v0.2",
        "protocol":"docs/exploratory/CULTIVATED_SENNA_OBSERVER_AND_DIRECT_LARVAL_USE_PROTOCOL_V02.json",
        "original_plant_ids":len(original),"original_plant_records_retrieved":len(observed),
        "source_consistency_exclusions":dict(exclusions),
        "cultivated_S_polyphylla_independent_reporting":describe(planted),
        "all_frozen_plant_groups":groups,
        "consumer_larvae":larvae,
        "direct_field_link_candidates":len(candidates),
        "larval_records_near_any_planted_record_2km_365d":len(matched_nearness),
        "interpretation":"Cultivated observations document reported host plant presence, not local butterfly use. Larval co-proximity is contextual and is not a trophic record. Original observation fields are not independent proof of consumption or survival. The native/exotic host rearing comparison in Koptur 2024 sampled cultivated plants in both groups."}
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(summary,indent=2)+"\n",encoding="utf8")
    a.output_candidates_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["butterfly","larval_observation_id","source_url","field_name","field_value","status"]
    with a.output_candidates_csv.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(candidates)
    print(json.dumps({"observer_gate":describe(planted)["gate_pass"],
                      "larval_source_candidates":len(candidates),
                      "nearby_larval_records_not_host_use":len(matched_nearness)}),flush=True)
if __name__=="__main__":main()
