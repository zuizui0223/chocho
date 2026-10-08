#!/usr/bin/env python3
"""Assess PUBLIC CRG butterfly rearing source completeness, not butterfly survival.

Reference: CRG iNaturalist project 2025 field instructions and 2026 report.
The ~9,519 project-wide rearings are NOT presumed to be available through
the iNaturalist project API. Restrict ecological inferences to source records.
"""
from __future__ import annotations
import argparse,csv,json,random,re,time,urllib.error,urllib.parse,urllib.request
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path

API="https://api.inaturalist.org/v1"
BUTTERFLY_FAMILIES={"Papilionidae","Pieridae","Nymphalidae","Lycaenidae","Hesperiidae","Riodinidae"}
UNKNOWN={"","NR","N/A","NONE","UNKNOWN","NOT RECORDED","NOT APPLICABLE","?","NA","-"}
PER_PAGE=200
MAX_PAGES=6

def get(endpoint,**params):
    url=API+"/"+endpoint+("?"+urllib.parse.urlencode(params) if params else "")
    for attempt in range(6):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"chocho-CRG-rearing-source-feasibility/0.1 (github.com/zuizui0223/chocho)",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=60) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as exc:
            if exc.code not in (429,500,502,503,504) or attempt==5:
                raise RuntimeError(f"CRG API error HTTP{exc.code}: {endpoint}") from exc
        except (TimeoutError,urllib.error.URLError) as exc:
            if attempt==5:raise RuntimeError(f"CRG API network timeout: {endpoint}") from exc
        time.sleep(min(30,2**attempt+1))
    raise RuntimeError("CRG API retry exhausted")

def rows(payload,where):
    result=payload.get("results")
    if not isinstance(result,list):raise RuntimeError(f"API missing results list at {where}")
    return result

def unique_project():
    xs=rows(get("projects/autocomplete",q="Caterpillar Rearing Group",per_page=30),"projects/autocomplete")
    exact=[v for v in xs if str(v.get("slug") or "")=="caterpillar-rearing-group"]
    if len(exact)!=1:raise RuntimeError(f"CRG project resolution ambiguous: {[x.get('slug') for x in xs]}")
    pid=int(exact[0]["id"])
    checked=rows(get(f"projects/{pid}"),"projects/id")
    if len(checked)!=1 or str(checked[0].get("slug") or "")!="caterpillar-rearing-group":
        raise RuntimeError("project ID/slug drift")
    return pid,exact[0].get("title") or exact[0].get("name")

def fields(observation):
    accum=defaultdict(list)
    for key in ("ofvs","observation_fields"):
        for item in observation.get(key) or []:
            if not isinstance(item,dict):continue
            name=str(item.get("name") or (item.get("observation_field") or {}).get("name") or "").strip()
            val=str(item.get("value") or "").strip()
            if name and val and val not in accum[name.lower()]:accum[name.lower()].append(val)
    return dict(accum)

def v(fieldmap,key):
    values=fieldmap.get(key.lower(),[])
    return values[0].strip() if values else ""

def present(s):
    return bool(s and s.strip().upper() not in UNKNOWN)

def valid_binomial(name):
    text=re.sub(r"\s+"," ",name.strip())
    return bool(re.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+",text))

def valid_date(x):
    x=x.strip()
    if not present(x):return None
    for sep in ("/","-","."):
        if sep in x:
            parts=x.split(sep)
            if len(parts)==3:
                try:
                    y,m,d=map(int,parts)
                    if 1900<=y<=2026:return date(y,m,d).isoformat()
                except ValueError:pass
    return None

def list_focal_species(path):
    res=set()
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                if float(r["host_family_count"])>0 and int(r["host_wgsrpd3_unit_count"])>0:
                    res.add(r["species"].strip())
            except (ValueError,KeyError):continue
    if len(res)!=239:raise RuntimeError(f"expected 239 focal species, got {len(res)}")
    return res

def butterfly_taxonomy(path):
    out={}
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            name=str(r.get("Species") or "").strip()
            family=str(r.get("Family") or "").strip()
            if name and family in BUTTERFLY_FAMILIES:out[name]=family
    return out

def collect_project(pid,protocol):
    sampled={}
    queries=[]
    for direction in protocol["prespecified_sample"]["ordering"]:
        order="asc" if direction.startswith("oldest") else "desc"
        for page in range(1,MAX_PAGES+1):
            params={"project_id":pid,"d1":"2012-01-01","d2":"2025-12-31",
                "per_page":PER_PAGE,"page":page,"order_by":"observed_on","order":order}
            payload=get("observations",**params)
            total=payload.get("total_results")
            if type(total) is not int:raise RuntimeError("CRG API total_results unavailable")
            result=rows(payload,f"observations {order} p{page}")
            queries.append({"order":order,"page":page,"total_results":total,"returned":len(result)})
            for item in result:
                obsid=item.get("id")
                if not isinstance(obsid,int):continue
                if obsid not in sampled:sampled[obsid]=item
            print(json.dumps({"order":order,"page":page,"total_results":total,"returned":len(result),
                "unique_sample":len(sampled)}),flush=True)
            if len(result)<PER_PAGE or page*PER_PAGE>=total:break
            time.sleep(.7)
    totals={x["total_results"] for x in queries}
    if len(totals)!=1:
        raise RuntimeError(f"project total changed during fixed sample: {sorted(totals)}")
    return sampled,queries,next(iter(totals))

def rearing_row(o,tax,focal):
    f=fields(o)
    scientific=str((o.get("taxon") or {}).get("name") or "").strip()
    rank=str((o.get("taxon") or {}).get("rank") or "")
    crg_species=v(f,"CRG Species")
    sp=crg_species if valid_binomial(crg_species) else scientific
    butterfly_family=v(f,"CRG Family").strip()
    if butterfly_family not in BUTTERFLY_FAMILIES:
        butterfly_family=tax.get(sp,"")
    butterfly=butterfly_family in BUTTERFLY_FAMILIES
    food=v(f,"CRG Food Plant Species")
    food_family=v(f,"CRG Food Plant Family")
    collection=valid_date(v(f,"CRG Date of Collection"))
    emergence=valid_date(v(f,"CRG Date of Emergence"))
    adult_url=v(f,"CRG Moth/Butterfly URL")
    caterpillar_url=v(f,"CRG Caterpillar URL")
    parasitism=v(f,"CRG Parasitoid details")
    parasitoid_url=v(f,"CRG Parasitoid URL")
    ref=v(f,"CRG Ref No.")
    url=f"https://www.inaturalist.org/observations/{o['id']}"
    uid=(o.get("user") or {}).get("id")
    adult_link=(adult_url.startswith("https://www.inaturalist.org/observations/")
                or adult_url.startswith("http://www.inaturalist.org/observations/"))
    larval_link=(caterpillar_url.startswith("https://www.inaturalist.org/observations/")
                 or caterpillar_url.startswith("http://www.inaturalist.org/observations/"))
    adult_emergence=bool(emergence and adult_link)
    parasitoid_emergence=bool(present(parasitism))
    has_larval_source=bool(larval_link or v(f,"Eating: (Interaction)") or v(f,"Eating"))
    return {
        "observation_id":o["id"],"observation_url":url,
        "original_observation_taxon":scientific,"rearing_taxon_name":sp,"taxon_rank":rank,
        "butterfly_family":butterfly_family,"is_butterfly":butterfly,
        "in_frozen_239_butterflies":sp in focal,
        "host_plant_verbatim":food,"host_plant_family_verbatim":food_family,
        "host_plant_strict_binomial":valid_binomial(food),
        "observer_id":uid if isinstance(uid,int) else "",
        "rearing_reference":ref,
        "crg_field_names_count":len(f),
        "collection_date":collection or "","adult_emergence_date":emergence or "",
        "adult_original_url":adult_url,"caterpillar_original_url":caterpillar_url,
        "parasitoid_original_url":parasitoid_url,
        "adult_emergence_with_original_link":adult_emergence,
        "parasitoid_details_present":parasitoid_emergence,
        "larval_source_link_present":has_larval_source,
        "crg_locality":v(f,"CRG Locality")[:200],
        "observation_year":str(o.get("observed_on") or "")[:4],
        "observation_photos_present":bool(o.get("photos")),
        "all_observation_fields":sorted(f)
    }

def main():
    ap=argparse.ArgumentParser()
    for arg in ("protocol_json","descriptors_csv","leptraits_csv","output_json","output_records_csv"):
        ap.add_argument("--"+arg.replace("_","-"),type=Path,required=True)
    a=ap.parse_args()
    protocol=json.loads(a.protocol_json.read_text(encoding="utf-8"))
    if protocol["status"]!="OUTCOME_BLIND_FEASIBILITY_PROTOCOL_FROZEN_BEFORE_PUBLIC_CRG_API_EXTRACTION":
        raise RuntimeError("CRG protocol changed")
    focal=list_focal_species(a.descriptors_csv)
    taxonomy=butterfly_taxonomy(a.leptraits_csv)
    pid,project_name=unique_project()
    raw,queries,total=collect_project(pid,protocol)
    rr=[rearing_row(o,taxonomy,focal) for o in raw.values()]
    rr.sort(key=lambda x:x["observation_id"])
    dedup={};collisions=Counter()
    for row in rr:
        key=(row["observer_id"],row["rearing_reference"]) if row["observer_id"] and present(row["rearing_reference"]) else ("obs",row["observation_id"])
        if key not in dedup:dedup[key]=row
        else:
            collisions["repeated_observer_rearing_reference"]+=1
            # prefer the record with stronger outcome fields and host species name, never fabricate links
            old=dedup[key]
            old_score=sum(bool(old[k]) for k in ("host_plant_strict_binomial","adult_emergence_with_original_link","parasitoid_details_present"))
            new_score=sum(bool(row[k]) for k in ("host_plant_strict_binomial","adult_emergence_with_original_link","parasitoid_details_present"))
            if new_score>old_score:dedup[key]=row
    unique=list(dedup.values())
    butterflies=[r for r in unique if r["is_butterfly"]]
    viable=[r for r in butterflies if r["host_plant_strict_binomial"] and (r["adult_emergence_with_original_link"] or r["parasitoid_details_present"])]
    distinct_species={r["rearing_taxon_name"] for r in viable}
    hosts={r["host_plant_family_verbatim"] for r in viable if present(r["host_plant_family_verbatim"])}
    panel_overlap={r["rearing_taxon_name"] for r in viable if r["in_frozen_239_butterflies"]}
    parasitoid=[r for r in butterflies if r["host_plant_strict_binomial"] and r["parasitoid_details_present"]]
    gate=len(viable)>=30 and len(distinct_species)>=5 and len(hosts)>=2
    source={
        "schema":"chocho_crg_butterfly_rearing_feasibility_v0.1",
        "protocol":"docs/exploratory/CRG_BUTTERFLY_REARING_FEASIBILITY_PROTOCOL_V01.json",
        "project":{"id":pid,"label":project_name,"source_url":"https://www.inaturalist.org/projects/caterpillar-rearing-group"},
        "public_observation_API_total_results":total,
        "observations_sampled_unique":len(rr),"deduplicated_rearing_keys":len(unique),
        "sampling_may_be_incomplete":total>len(rr),
        "query_coverage":queries,"collisions":dict(collisions),
        "all_sampled_CRg_field_observation_count":sum(bool(r["crg_field_names_count"]) for r in rr),
        "butterfly":{
            "unique_rearing_keys":len(butterflies),
            "butterfly_family_distribution":dict(Counter(r["butterfly_family"] for r in butterflies)),
            "with_strict_host_binomial":sum(r["host_plant_strict_binomial"] for r in butterflies),
            "with_dated_and_linked_adult_emergence":sum(r["adult_emergence_with_original_link"] for r in butterflies),
            "with_parasitoid_notes":sum(r["parasitoid_details_present"] for r in butterflies),
            "with_original_larval_source_link":sum(r["larval_source_link_present"] for r in butterflies),
            "with_strict_host_and_either_confirmed_emergence":len(viable),
            "species_with_valid_endpoint":len(distinct_species),
            "host_families_with_valid_endpoint":len(hosts),
            "frozen_239_focal_butterfly_species_with_valid_endpoint":len(panel_overlap),
            "strict_host_plus_parasitoid_details_events":len(parasitoid),
            "field_coverage_gate_passed":gate
        },
        "stop_if_not_feasible":not gate,
        "qualifications":[
            "This is a project-specific source availability pilot, not a representative sample of all project rearings.",
            "Records without food plants or adult URLs are treated as insufficiently documented, NOT biological rearing failures.",
            "Only confirmed positive adult-emergence or parasitoid-emergence reports are counted; there is NO denominator of attempted rearings.",
            "With unknown failure reporting, adult/parasitoid success probabilities and relative survival are non-identifiable.",
            "Host plant scientific names are verbatim from CRG fields, without external taxonomic or geographic validation.",
            "Parasitism details field is observer reporting, not independently confirmed parasitoid taxon identification.",
            "No claim of native/exotic host origin, causal benefit, local populations or multi-site experiments."
        ]
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(source,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    a.output_records_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["observation_id","observation_url","original_observation_taxon","rearing_taxon_name","taxon_rank",
      "butterfly_family","is_butterfly","in_frozen_239_butterflies","host_plant_verbatim",
      "host_plant_family_verbatim","host_plant_strict_binomial","observer_id","rearing_reference",
      "collection_date","adult_emergence_date","adult_original_url","caterpillar_original_url",
      "parasitoid_original_url","adult_emergence_with_original_link","parasitoid_details_present",
      "larval_source_link_present","crg_locality","observation_year","observation_photos_present"]
    with a.output_records_csv.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for record in rr:writer.writerow({k:record.get(k) for k in fields})
    print(json.dumps({"project_id":pid,"public_api_total":total,"sample":len(rr),
          "butterfly_rearing_keys":len(butterflies),"valid_events":len(viable),
          "species_with_endpoint":len(distinct_species),"frozen_239_species":len(panel_overlap),
          "feasibility_pass":gate}),flush=True)
if __name__=="__main__":main()
