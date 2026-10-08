#!/usr/bin/env python3
"""Nonconfirmatory feasibility audit of independently observed butterfly larval host use."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,re,time
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from shapely.geometry import Point,shape
from shapely.strtree import STRtree

SPECIES=[
    "Anteos clorinde",
    "Anteos maerula",
    "Catopsilia pomona",
    "Catopsilia pyranthe",
    "Colias alexandra",
    "Colias christina",
    "Colias croceus",
    "Colias erate",
    "Colias eurytheme",
    "Colias hecla",
    "Colias hyale",
    "Colias nastes",
    "Colias philodice",
    "Delias hyparete",
    "Eurema arbela",
    "Eurema hecabe",
    "Eurema mexicana",
    "Leptidea sinapis",
    "Phoebis agarithe",
    "Phoebis neocypris",
    "Phoebis philea",
    "Phoebis sennae",
    "Pieris brassicae",
    "Zerene cesonia"
]
URL="https://api.globalbioticinteractions.org/interaction.csv"
FIELDS=[
    "source_taxon_name","target_taxon_name","target_taxon_path",
    "interaction_type","event_date","latitude","longitude",
    "source_specimen_life_stage","study_source_citation","study_citation",
    "study_external_id","source_specimen_occurrence_id","locality"
]
HERBIVORE_TERMS={"eats","feeds on","feedsOn","herbivory","consumes","eats plants"}
BAD_PROVENANCE=["globalbioticinteractions/hosts","nhm hosts database","nhm-hosts"]

def load_sources(host_path,native_path,contemp_path):
    host=defaultdict(set);all_host=defaultdict(set)
    with host_path.open(newline="",encoding="utf-8") as h:
        for r in csv.DictReader(h):
            sp=str(r["insect_species"]).strip()
            if sp not in SPECIES:continue
            pid=str(r["accepted_plant_name_id"]).strip()
            all_host[sp].add(pid)
            if str(r.get("family") or "").strip()=="Fabaceae":
                host[(sp,str(r["accepted_name"]).strip())].add(pid)
    def read_distributions(path):
        out=defaultdict(set)
        with path.open(newline="",encoding="utf-8") as f:
            for r in csv.DictReader(f):
                out[str(r["accepted_plant_name_id"]).strip()].add(str(r["area_code_l3"]).strip())
        return out
    native=read_distributions(native_path);contemp=read_distributions(contemp_path)
    native_union={sp:set().union(*(native.get(h,set()) for h in ids)) if ids else set() for sp,ids in all_host.items()}
    return host,native,contemp,native_union

def geometry_index(path):
    data=json.loads(path.read_text(encoding="utf-8"))
    geoms=[];codes=[]
    for feat in data["features"]:
        code=str((feat.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if code and feat.get("geometry"):
            geoms.append(shape(feat["geometry"]));codes.append(code)
    if not geoms:raise RuntimeError("empty WGSRPD3 geometry")
    return geoms,codes,STRtree(geoms)

def map_code(lat,lon,geoms,codes,tree):
    pt=Point(float(lon),float(lat))
    hits=[]
    for idx in tree.query(pt):
        idx=int(idx)
        if geoms[idx].covers(pt):hits.append(codes[idx])
    return "" if not hits else sorted(set(hits))[0]

def fetch_page(species,offset,limit,attempts=4):
    query=urlencode({
        "sourceTaxon":species,
        "interactionType":"eats",
        "includeObservations":"true",
        "limit":limit,
        "offset":offset,
        "fields":",".join(FIELDS)
    })
    url=f"{URL}?{query}";error=None
    for attempt in range(attempts):
        try:
            req=Request(url,headers={"User-Agent":"chocho-independent-larval-evidence-probe/0.1","Accept":"text/csv"})
            with urlopen(req,timeout=75) as response:
                raw=response.read(16*1024*1024+1)
                if len(raw)>16*1024*1024:raise RuntimeError("API page exceeds 16 MiB cap")
                return raw,url,None
        except (HTTPError,URLError,TimeoutError,RuntimeError) as exc:
            error=str(exc)
            if attempt+1<attempts:
                time.sleep(min(20,3*(attempt+1)))
    return None,url,error

def read_year(text):
    match=re.search(r"(?<!\d)(1[8-9]\d{2}|20[0-2]\d)(?!\d)",str(text or ""))
    if not match:return None
    year=int(match.group(1))
    return year if 1800<=year<=2025 else None

def canonical_match(actual,wanted):
    actual=" ".join(str(actual or "").strip().split())
    return actual==wanted or actual.startswith(wanted+" ")

def plausible_plant(path):
    low=str(path or "").lower()
    return "plantae" in low or "tracheophyta" in low or "viridiplantae" in low

def stage_is_larva(row):
    stage=str(row.get("source_specimen_life_stage") or "").lower()
    return any(token in stage for token in ("larva","caterpillar"))

def main():
    ap=argparse.ArgumentParser()
    for opt in ("insect_host_csv","native_distribution_csv","contemporary_distribution_csv","level3_geojson","output_json","output_csv","raw_dir"):
        ap.add_argument("--"+opt.replace("_","-"),required=True,type=Path)
    ap.add_argument("--page-limit",type=int,default=512)
    ap.add_argument("--pages-per-species",type=int,default=2)
    a=ap.parse_args()
    hosts,native,contemp,native_union=load_sources(a.insect_host_csv,a.native_distribution_csv,a.contemporary_distribution_csv)
    geoms,codes,tree=geometry_index(a.level3_geojson)
    a.raw_dir.mkdir(parents=True,exist_ok=True)
    found=[];query_audit=[];stages=Counter();types=Counter();studies=Counter();query_completion=[]
    for species in SPECIES:
        n_species=0;any_failed=False;maybe_capped=False
        for page in range(a.pages_per_species):
            offset=page*a.page_limit
            raw,url,error=fetch_page(species,offset,a.page_limit)
            if raw is None:
                query_audit.append({"species":species,"offset":offset,"records":0,"error":error,"url":url})
                any_failed=True
                break
            hashval=hashlib.sha256(raw).hexdigest()
            filename=f"{species.replace(' ','_')}_{offset:05d}.csv.gz"
            with gzip.open(a.raw_dir/filename,"wb") as target:target.write(raw)
            reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
            rows=list(reader)
            query_audit.append({"species":species,"offset":offset,"records":len(rows),"sha256":hashval,"file":filename,"url":url,"error":None})
            if page==a.pages_per_species-1 and len(rows)>=a.page_limit:maybe_capped=True
            n_species+=len(rows)
            for r in rows:
                if not canonical_match(r.get("source_taxon_name"),species):continue
                stages["exact_butterfly_taxon"]+=1
                relation=str(r.get("interaction_type") or "").strip()
                types[relation]+=1
                if relation.lower() not in {x.lower() for x in HERBIVORE_TERMS}:continue
                stages["herbivory_relation"]+=1
                if not plausible_plant(r.get("target_taxon_path")):continue
                stages["plant_target"]+=1
                year=read_year(r.get("event_date"))
                try:
                    lat=float(r.get("latitude") or "nan");lon=float(r.get("longitude") or "nan")
                    valid=(-90<=lat<=90) and (-180<=lon<=180) and not (lat==0 and lon==0)
                except (ValueError,TypeError):valid=False
                if year is None or not valid:continue
                stages["geo_and_dated"]+=1
                if not stage_is_larva(r):continue
                stages["explicit_larva"]+=1
                provenance=" | ".join(str(r.get(k) or "") for k in ["study_source_citation","study_citation","study_external_id"])
                if not provenance.strip(" |"):continue
                if any(token in provenance.lower() for token in BAD_PROVENANCE):continue
                stages["non_hosts_provenance_larval"]+=1
                code=map_code(lat,lon,geoms,codes,tree)
                if not code:continue
                stages["mapped_to_wgsrpd3"]+=1
                studies[provenance[:140]]+=1
                plant=str(r.get("target_taxon_name") or "").strip()
                matched=hosts.get((species,plant),set())
                introduced_ids=sorted(h for h in matched if code in contemp.get(h,set()) and code not in native.get(h,set()))
                introduced_only=bool(introduced_ids and code not in native_union.get(species,set()))
                record={
                    "butterfly_species":species,
                    "plant_name":plant,
                    "interaction_type":relation,
                    "larval_stage":str(r.get("source_specimen_life_stage") or ""),
                    "year":year,"lat":lat,"lon":lon,"wgsrpd3_code":code,
                    "source_occurrence_id":str(r.get("source_specimen_occurrence_id") or ""),
                    "study_source_citation":str(r.get("study_source_citation") or ""),
                    "study_citation":str(r.get("study_citation") or ""),
                    "study_external_id":str(r.get("study_external_id") or ""),
                    "wcvp_documented_known_host_exact_name":int(bool(matched)),
                    "introduced_known_host_in_record_region":int(bool(introduced_ids)),
                    "introduced_only_resource_cell":int(introduced_only)
                }
                found.append(record)
            if len(rows)<a.page_limit:break
        query_completion.append({"species":species,"rows":n_species,"query_failed":any_failed,"possibly_capped":maybe_capped})
        print(json.dumps({"species":species,"response_rows":n_species,"failed":any_failed,"possibly_capped":maybe_capped}),flush=True)
        stages["species_query_finished"]+=1
    dedup={}
    for r in found:
        key=(r["butterfly_species"],r["plant_name"],r["year"],round(r["lat"],4),round(r["lon"],4),r["source_occurrence_id"],r["study_external_id"])
        dedup[key]=r
    records=sorted(dedup.values(),key=lambda r:(r["butterfly_species"],r["wgsrpd3_code"],r["plant_name"],r["year"]))
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["butterfly_species","plant_name","interaction_type","larval_stage","year","lat","lon","wgsrpd3_code","source_occurrence_id","study_source_citation","study_citation","study_external_id","wcvp_documented_known_host_exact_name","introduced_known_host_in_record_region","introduced_only_resource_cell"]
    with a.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
    # Independent evidence of Faba ceae feeding is a provenance-filtered
    # interaction observation. It is NOT survival or local colonization.
    positive=[r for r in records if r["introduced_known_host_in_record_region"]]
    introduced_only=[r for r in records if r["introduced_only_resource_cell"]]
    score={
        "strict_geo_dated_non_hosts_larval_observations":len(records),
        "exact_known_host_record_observations":sum(r["wcvp_documented_known_host_exact_name"] for r in records),
        "introduced_known_host_region_records":sum(r["introduced_known_host_in_record_region"] for r in records),
        "introduced_fabaceae_larval_record_count":len(positive),
        "introduced_fabaceae_larval_species":len({r["butterfly_species"] for r in positive}),
        "introduced_fabaceae_larval_regions":len({r["wgsrpd3_code"] for r in positive}),
        "strict_introduced_only_resource_region_records":len(introduced_only),
        "strict_introduced_only_species":len(set(r["butterfly_species"] for r in introduced_only)),
        "strict_introduced_only_regions":len(set(r["wgsrpd3_code"] for r in introduced_only))
    }
    result={
        "schema":"chocho_pieridae_fabaceae_independent_larval_use_v0.1",
        "status":"EXPLORATORY_TAXON_BLIND_FABACEAE_FEEDING_FEASIBILITY",
        "protocol":"docs/exploratory/PIERIDAE_FABACEAE_LARVAL_USE_PILOT_PROTOCOL_V01.json",
        "species_panel":SPECIES,
        "page_limit":a.page_limit,
        "pages_per_species":a.pages_per_species,
        "api_restriction":"interactionType=eats, not all GloBI interaction types",
        "frozen_selection":"all 24 Pieridae with documented HOSTS Fabaceae use, source panel frozen before querying GloBI",
        "Fabaceae_host_join":"only exact accepted Fabaceae names; all native host alternatives retained when flagging introduced-only cells",
        "query_completion":query_completion,
        "possibly_capped_species":[x["species"] for x in query_completion if x["possibly_capped"]],
        "query_audit":query_audit,
        "stages":dict(stages),
        "interaction_types":dict(types.most_common(25)),
        "source_citations":dict(studies.most_common(20)),
        "coverage":score,
        "gate_pass":bool(score["introduced_fabaceae_larval_record_count"]>=10 and score["introduced_fabaceae_larval_species"]>=3 and score["introduced_fabaceae_larval_regions"]>=3),
        "claim_boundary":[
            "This is NOT a random sample and GloBI API pagination can be incomplete.",
            "GloBI may re-index HOSTS-like data; source citations require manual provenance verification.",
            "Exact WCVP accepted names and source life stage and ingestion are required; ambiguous aliases are not assumed.",
            "Coordinates may be centroids or imprecise; geographic assignment requires source review.",
            "Presence of a documented host relationship is NOT evidence of local use without an independent larval observation.",
            "Even original-source larval feeding annotations do not prove survival, population growth or colonization.",
            "Only a feasibility gate; no butterfly population inference, causal effect, or formal test."
        ]
    }
    a.output_json.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"coverage":score,"gate_pass":result["gate_pass"],"query_failures":sum(bool(x.get("error")) for x in query_audit)},indent=2))
if __name__=="__main__":main()
