#!/usr/bin/env python3
"""Nonselected original230 butterfly HOSTS x dated GBIF port plants x butterfly recurrence.

Data are source-pinned and author-side. DO NOT republish original butterfly-host
edges, restricted WCVP host concepts or individual occurrence coordinates.
GBIF botanical records are observations, NOT plant absence/larval feeding,
introduction date, true butterfly colonisation, or evidence of ship travel.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,random,time
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

R=Path("provenance/reviewer_defenses/results")
EVENTS=R/"island_host_continuity_20261030/island_host_continuity_events.csv"
CANDIDATES=R/"island_host_network_20261030/island_host_network_all225.csv"
FLORA=R/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
PINS={
 EVENTS:"efbc770469a720e2dea82478db2f6196ba03b353",
 CANDIDATES:"39fe84981cab2b0e47063399f0a5c0c9731e800d",
 FLORA:"54c4515947d2a3ee7545e4007985ef82a2af6868",
}
HOSTS_SHA="e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22"
PLANT_CANDIDATES_SHA="a1cc87b20fbd19415fa81ca4269b8414b0c32167bb6033ca0425354c89c1906f"
ORIGINAL23_SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
NEW15_SHA=None  # loaded from source-pinned FLORA
MAPPING_API="https://api.gbif.org/v1/species/match"
SEED=20261031
DRAWS=299
BURN_PER_EDGE=25
BETWEEN_PER_EDGE=4
MAX_WORKERS=6

def sha(data):return hashlib.sha256(data).hexdigest()
def readcsv(data):return list(csv.DictReader(io.StringIO(data.decode("utf-8") if isinstance(data,bytes) else data)))
def source(root,p):
    data=(root/p).read_bytes()
    s=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if s!=PINS[p]:raise ValueError(f"Tracked source blob SHA mismatch {p}: {s}")
    return data
def read_source_plant_gzip(path,expected,old):
    data=path.read_bytes()
    if sha(data)!=expected:raise ValueError("Original source GBIF speciesKey×port×year GZIP SHA changed")
    a=readcsv(gzip.decompress(data))
    exp=9788 if old else 327
    if len(a)!=exp:raise ValueError(f"Unexpected archived speciesKey×port×year count {len(a)} != {exp}")
    out=[]
    for x in a:
        k=x["gbif_speciesKey"] if old else x["speciesKey"]
        n=x["indexed_record_count"] if old else x["indexed_records"]
        y=int(x["year"])
        if y<2001 or y>2010 or int(n)<1:raise ValueError("Historical plant facet year/count invalid")
        out.append((str(x["island_id"]),str(k),y,str(x["port_id"]),int(n)))
    return out

def single_name_match(name):
    params={"name":name,"kingdom":"Plantae","rank":"SPECIES","strict":"true"}
    url=MAPPING_API+"?"+urlencode(params)
    for attempt in range(5):
        try:
            with urlopen(Request(url,headers={"User-Agent":"chocho-original230-historical-host-source-audit/1.0",
                                              "Accept":"application/json"}),timeout=45) as response:
                data=json.load(response)
            if not isinstance(data,dict):raise ValueError("GBIF taxonomy returned non-object")
            if data.get("matchType")=="NONE":
                return name,None,"NONE"
            exact=(str(data.get("matchType"))=="EXACT" and str(data.get("rank"))=="SPECIES"
                 and str(data.get("canonicalName","")).casefold()==name.casefold()
                 and int(data.get("confidence",0))>=90)
            if not exact:return name,None,"NOT_STRICT_EXACT_SPECIES"
            if data.get("status")!="ACCEPTED":
                return name,None,"EXACT_BUT_NOT_ACCEPTED"
            key=data.get("speciesKey") or data.get("usageKey")
            if not str(key).isdigit():return name,None,"ACCEPTED_NONNUMERIC_KEY"
            return name,str(key),"EXACT_ACCEPTED"
        except (HTTPError,URLError,TimeoutError) as ex:
            if isinstance(ex,HTTPError) and ex.code not in (408,425,429,500,502,503,504):
                raise
            if attempt==4:raise
            time.sleep(min(20,2**attempt)+random.Random(SEED+len(name)+attempt).random()/2)
    raise AssertionError("Unreachable")

def resolve_all_source_plant_names(rows):
    names=sorted({x["query_name"] for x in rows if x["query_eligible"]=="True" and x["query_name"]})
    if len(names)!=1642:raise ValueError("Original 1642 source-eligible plant binomials changed")
    output={}
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures={pool.submit(single_name_match,name):name for name in names}
        for i,future in enumerate(as_completed(futures),1):
            name,key,status=future.result()
            output[name]={"gbif_speciesKey":key,"gate":status}
            if i%200==0:print("GBIF_HOST_TAXON_GATE",i,"/",len(names),flush=True)
    if set(output)!=set(names):raise ValueError("GBIF taxonomy response coverage mismatch")
    reverse=defaultdict(set)
    for name,entry in output.items():
        if entry["gbif_speciesKey"] is not None:reverse[entry["gbif_speciesKey"]].add(name)
    ambiguous={k for k,v in reverse.items() if len(v)>1}
    for name,val in output.items():
        if val["gbif_speciesKey"] in ambiguous:
            val["gbif_speciesKey"]=None
            val["gate"]="DUPLICATE_SOURCE_NAMES_TO_ONE_GBIF_KEY"
    # Strong source identity controls for known matched hosts.
    for known,expected in [("Medicago sativa","9151957"),("Alnus incana","9148577")]:
        if output.get(known,{}).get("gbif_speciesKey")!=expected:
            raise ValueError("Known focal plant taxonomy changed: "+known)
    return output

def classify(events,edges,plants,sites):
    by=defaultdict(set)
    for butterfly,plant,family in edges:
        by[butterfly].add(plant)
    samples=[]
    for ev in events:
        keys=by.get(ev["original_butterfly_name"],set())
        if not keys:
            continue
        matches=keys&plants[ev["landmass_id"]]
        samples.append({**ev,"source_resolved_host_count":len(keys),
            "has_prior_recorded_compatible_host":bool(matches),
            "n_prior_recorded_compatible_hosts":len(matches),
            "site_early_plant_species_keys":sites[ev["landmass_id"]]["n_indexed_species_keys"],
            "site_early_botanical_years":sites[ev["landmass_id"]]["years_with_indexed_botanical_records"]})
    return samples

def summary(samples):
    t=len(samples)
    yes=[x for x in samples if x["has_prior_recorded_compatible_host"]]
    no=[x for x in samples if not x["has_prior_recorded_compatible_host"]]
    def rate(rows,col):
        return sum(x[col] for x in rows)/len(rows) if rows else None
    a=rate(yes,"any_re_record");b=rate(no,"any_re_record")
    fivea=rate(yes,"record_year_five");fiveb=rate(no,"record_year_five")
    return {
      "eligible_first_record_events":t,
      "species":len({x["original_butterfly_name"] for x in samples}),
      "islands":len({x["landmass_id"] for x in samples}),
      "historically_recorded_compatible_host_events":len(yes),
      "no_indexed_compatible_host_events_NOT_TRUE_ABSENCE":len(no),
      "any_later_re_record": {"matched_n":len(yes),"matched_yes":sum(x["any_re_record"] for x in yes),
          "no_index_n":len(no),"no_index_yes":sum(x["any_re_record"] for x in no),
          "matched_rate":a,"no_index_rate":b,
          "difference":a-b if a is not None and b is not None else None},
      "re_record_in_fifth_year":{"matched_rate":fivea,"no_index_rate":fiveb,
         "difference":fivea-fiveb if fivea is not None and fiveb is not None else None},
    }

def restricted_cohorts(samples):
    all_=summary(samples)
    def filt(c,y):
        return summary([x for x in samples if x["site_early_plant_species_keys"]>=c
            and x["site_early_botanical_years"]>=y])
    return {"all_source_eligible":all_,
        "botanical_20keys_2years":filt(20,2),
        "botanical_50keys_3years_exploratory":filt(50,3)}

def conditional_strata(samples):
    contrasts=[]
    for var in ["landmass_id","original_butterfly_name"]:
        group=defaultdict(list)
        for ev in samples:group[ev[var]].append(ev)
        scores=[]
        for case in group.values():
            a=summary(case)["any_later_re_record"]["difference"]
            if a is not None:scores.append(a)
        contrasts.append({"conditioning":var,"groups_with_both_detected_and_undetected":len(scores),
            "unweighted_group_mean_re_record_rate_difference":sum(scores)/len(scores) if scores else None,
            "positive_group_contrasts":sum(x>0 for x in scores),
            "negative_group_contrasts":sum(x<0 for x in scores),
            "zero_group_contrasts":sum(x==0 for x in scores)})
    return contrasts

def island_bootstrap(samples):
    islands=sorted({x["landmass_id"] for x in samples},key=int)
    by=defaultdict(list)
    for x in samples:by[x["landmass_id"]].append(x)
    rng=random.Random(SEED)
    draws=[]
    for _ in range(2000):
        chosen=[rng.choice(islands) for _ in islands]
        s=summary([v for iid in chosen for v in by[iid]])
        delta=s["any_later_re_record"]["difference"]
        if delta is not None:draws.append(delta)
    draws.sort()
    return {"island_count":len(islands),"repetitions_attempted":2000,"valid_repetitions":len(draws),
        "fixed_event_bootstrap_without_model_refit":True,
        "interval95":[draws[int(.025*(len(draws)-1))],draws[int(.975*(len(draws)-1))]] if len(draws)>50 else None}

def source_network_null(samples,edges):
    """Family and per-butterfly/per-plant degree preserving rewires, not causal."""
    rng=random.Random(SEED)
    # Sort source associations to eliminate Python set/hash-seed dependent MCMC ordering.
    ed=[[s,p,f] for s,p,f in sorted(edges)]
    used={(s,p) for s,p,f in ed}
    if len(used)!=len(ed):raise ValueError("Duplicate source GBIF resolved host edges")
    family=defaultdict(list)
    for i,(s,p,f) in enumerate(ed):family[f].append(i)
    family_names=sorted(f for f,x in family.items() if len(x)>1)
    original_bfly=Counter(s for s,p,f in ed)
    original_plant=Counter(p for s,p,f in ed)
    original_bflyfamily=Counter((s,f) for s,p,f in ed)
    def shuffle(n):
        accepted=0
        for _ in range(n):
            if not family_names:break
            f=rng.choice(family_names);pool=family[f]
            i=rng.choice(pool);j=rng.choice(pool)
            a,x,_=ed[i];b,y,_=ed[j]
            if a==b or x==y or (a,y) in used or (b,x) in used:continue
            used.remove((a,x));used.remove((b,y))
            ed[i][1],ed[j][1]=y,x
            used.add((a,y));used.add((b,x));accepted+=1
        return accepted
    observed=summary(samples)
    target_fraction=observed["historically_recorded_compatible_host_events"]/len(samples)
    target_delta=observed["any_later_re_record"]["difference"]
    accepted=shuffle(BURN_PER_EDGE*len(ed))
    fractions=[];deltas=[]
    for _ in range(DRAWS):
        accepted+=shuffle(BETWEEN_PER_EDGE*len(ed))
        if len(used)!=len(ed):raise ValueError("Null lost unique butterfly-host edges")
        by=defaultdict(set)
        for s,p,_ in ed:by[s].add(p)
        flags=[]
        for s in samples:
            new=dict(s)
            new["has_prior_recorded_compatible_host"]=bool(by[s["original_butterfly_name"]]&s["_island_plants"])
            flags.append(new)
        res=summary(flags)
        fractions.append(res["historically_recorded_compatible_host_events"]/len(flags))
        if res["any_later_re_record"]["difference"] is not None:
            deltas.append(res["any_later_re_record"]["difference"])
    if Counter(s for s,p,f in ed)!=original_bfly or Counter(p for s,p,f in ed)!=original_plant or Counter((s,f) for s,p,f in ed)!=original_bflyfamily:
        raise ValueError("Structural host-null degree or family count was NOT preserved")
    def mean(vals):return sum(vals)/len(vals) if vals else None
    def quantile(vals,q):
        return sorted(vals)[int(q*(len(vals)-1))] if vals else None
    return {"algorithm":"fixed source family-constrained bipartite host edge swaps",
        "seed":SEED,"n_draws":DRAWS,"edges_in_resolved_network":len(ed),
        "accepted_edge_swaps":accepted,
        "butterfly_degree_host_degree_and_species_x_family_breadth_preserved":True,
        "observed_prior_host_index_event_fraction":target_fraction,
        "null_fraction_mean":mean(fractions),
        "null_fraction_95pct":[quantile(fractions,.025),quantile(fractions,.975)],
        "null_draws_ge_observed_fraction":sum(v>=target_fraction for v in fractions),
        "observed_followup_rate_delta":target_delta,
        "null_followup_delta_mean":mean(deltas),
        "null_followup_delta_95pct":[quantile(deltas,.025),quantile(deltas,.975)],
        "null_draws_ge_observed_delta":sum(v>=target_delta for v in deltas) if target_delta is not None else None,
        "valid_null_delta_draws":len(deltas),
        "CANNOT_DISTINGUISH_BIOLOGICAL_HOST_COMMONNESS_FROM_REPORTING_PROMINENCE":True,
        "NOT_TREATING_NO_HOST_INDEXED_AS_TRUE_ABSENCE":True}

def run(args):
    root=args.root
    original_edges=(args.author_roster_dir/"original230_host_edges.csv").read_bytes()
    original_hosts=(args.author_roster_dir/"plant_taxon_query_candidates.csv").read_bytes()
    if sha(original_edges)!=HOSTS_SHA or sha(original_hosts)!=PLANT_CANDIDATES_SHA:
        raise ValueError("Original author 230-source host catalogue bytes differ")
    rawedges=readcsv(original_edges)
    queries=readcsv(original_hosts)
    if len(rawedges)!=2524 or len(queries)!=1677:
        raise ValueError("Author original host roster size changed")
    records=read_source_plant_gzip(args.original23,ORIGINAL23_SHA,True)
    src=json.loads(source(root,FLORA))
    records+=read_source_plant_gzip(args.new15,src["new_15port_gzip_sha256"],False)
    site={x["island_id"]:x for x in src["islands"]}
    ids=set(site)
    if (len(ids),len(records))!=(17,10115):
        raise ValueError("Wrong 17-node historical plant index population")
    plant_set=defaultdict(set)
    plant_year=defaultdict(lambda:defaultdict(set))
    for i,k,y,pid,n in records:
        if i not in ids:raise ValueError("Plant port assigned to unknown target island")
        plant_set[i].add(k);plant_year[i][k].add(y)
    for iid,info in site.items():
        if len(plant_set[iid])!=info["n_indexed_species_keys"]:
            raise ValueError("Plant source species-key set does not match original17 source report")
    events=readcsv(source(root,EVENTS))
    study=[x for x in events if x["landmass_id"] in ids and 2011<=int(x["first_record_year"])<=2017]
    if len(events)!=5665 or len(study)!=121 or len({x["species"] for x in study})!=64:
        raise ValueError("Original event target population changed")
    candidate=readcsv(source(root,CANDIDATES))
    if len(candidate)!=225:raise ValueError("Butterfly reference species changed")
    candidates={x["species"] for x in candidate}
    if any(x["species"] not in candidates for x in study):
        raise ValueError("New butterfly records not in fixed 225-species source roster")
    if args.taxonomy_from_archive:
        taxonomy_bytes=args.taxonomy_from_archive.read_bytes()
        payload=json.loads(gzip.decompress(taxonomy_bytes).decode("utf8"))
        if payload.get("schema")!="chocho_accepted_original_host_species_gbif_crosswalk_v1":
            raise ValueError("Historic GBIF botanical name crosswalk archive has incorrect schema")
        if payload.get("source_author_plant_candidates_sha256")!=PLANT_CANDIDATES_SHA:
            raise ValueError("Historic plant source identity changed")
        names=payload["name_resolution"]
    else:
        names=resolve_all_source_plant_names(queries)
        payload={"schema":"chocho_accepted_original_host_species_gbif_crosswalk_v1",
            "source_author_plant_candidates_sha256":PLANT_CANDIDATES_SHA,
            "gbif_api":MAPPING_API,
            "strict_taxonomic_gate":"EXACT ACCEPTED SPECIES >=90 confidence, canonical binomial match",
            "name_resolution":names}
        taxonomy_bytes=gzip.compress(json.dumps(payload,sort_keys=True,ensure_ascii=False,separators=(",",":")).encode("utf8"),mtime=0)
    if len(names)!=1642:
        raise ValueError("GBIF historical host name crosswalk coverage not complete")
    if args.taxonomy_archive:
        args.taxonomy_archive.parent.mkdir(parents=True,exist_ok=True)
        args.taxonomy_archive.write_bytes(taxonomy_bytes)
    taxonomy_sha=sha(taxonomy_bytes)
    # Mark original taxon collisions as unresolvable, NOT evidence for absence.
    host_edges=set()
    catalogue_species=defaultdict(set)
    source_family=defaultdict(set)
    for e in rawedges:
        sp=e["original_name"];source_family[sp].add(e["family"])
        if e["query_eligible"]!="True" or not e["query_name"]:
            continue
        key=names[e["query_name"]]["gbif_speciesKey"]
        if key is not None:host_edges.add((sp,key,e["family"]))
        catalogue_species[sp].add(e["accepted_name"])
    if len({s for s,k,f in host_edges})<200:
        raise ValueError("Unexpectedly poor original butterfly host taxonomic resolution")
    # Original butterfly synonyms not adopted without a unique name identity.
    bysource={s for s,k,f in host_edges}
    accepted_alias=defaultdict(set)
    for e in rawedges:
        accepted_alias[e["current_accepted_name"]].add(e["original_name"])
    normalized=[]
    for x in study:
        sp=x["species"]
        if sp in bysource:
            srcsp=sp
            match="ORIGINAL_EXACT"
        elif len(accepted_alias.get(sp,set())&bysource)==1:
            srcsp=next(iter(accepted_alias[sp]&bysource))
            match="ACCEPTED_UNIQUE"
        else:
            srcsp=None
            match="NO_UNAMBIGUOUS_SOURCE_CATALOGUE_TAXON"
        normalized.append({
            "original_butterfly_name":srcsp,
            "source_name_match_type":match,
            "landmass_id":x["landmass_id"],
            "first_indexed_butterfly_year":int(x["first_record_year"]),
            "any_re_record":x["any_followup_record"]=="True",
            "record_year_five":x["record_in_year5"]=="True",
        })
    viable=[x for x in normalized if x["original_butterfly_name"] is not None]
    scored=classify(viable,host_edges,plant_set,site)
    if len(scored)!=len(viable):raise ValueError("Score frame changed after source butterfly name gate")
    for e in scored:e["_island_plants"]=plant_set[e["landmass_id"]]
    # Derived summaries must not serialize raw host links or butterfly-plant keyed rows.
    rates=restricted_cohorts(scored)
    condition=conditional_strata(scored)
    boots=island_bootstrap(scored)
    null=source_network_null(scored,host_edges)
    taxonomy=Counter(x["gate"] for x in names.values())
    reported={"schema":"chocho_source_original230_dated_plant_host_fiveyear_butterfly_recording_v1",
        "status":"EXPLORATORY_RECORDED_COMPATIBLE_HOST_NOT_CAUSAL_RESCUE",
        "audit_date_utc":datetime.now(timezone.utc).isoformat(),
        "author_original_host_edges_sha256":HOSTS_SHA,
        "author_plant_candidates_sha256":PLANT_CANDIDATES_SHA,
        "public_repo_git_blobs":{str(p):v for p,v in PINS.items()},
        "original_23_port_gzip_sha256":ORIGINAL23_SHA,
        "new_15_port_gzip_sha256":src["new_15port_gzip_sha256"],
        "taxonomic_match":{"eligible_plant_names":len(names),
            "resolved_name_to_gbifkey_frozen_archive_sha256":taxonomy_sha,
            "resolution_status_counts":dict(taxonomy),
            "source_original_butterflies":len({x["original_name"] for x in rawedges}),
            "source_original_edges":len(rawedges),
            "resolved_butterfly_host_edges":len(host_edges),
            "resolved_original_butterfly_species":len({x[0] for x in host_edges}),
            "match_algorithm":"GBIF v1 exact accepted binomial SPECIES >=90 confidence; remove GBIF accepted-key collisions"},
        "butterfly_event_cohort":{"event_rows_original_all_islands":len(events),
            "preselected_near_islands":len(ids),"events_in_17_islands_2011_2017":len(study),
            "species_with_2011_2017_new_index":len({x["species"] for x in study}),
            "target_with_resolved_source_catalogue":len(scored),
            "source_lost_or_unresolved":len(study)-len(scored),
            "outcome_any_followup_is_indexed_2011_2017_first_year_plus_5":True,
            "retrospective_host_group_column_NOT_USED_TO_ASSIGN_HISTORIC_HOST":True},
        "rates":rates,
        "within_group_contrasts":condition,
        "fixed_event_island_cluster_bootstrap":boots,
        "host_network_structural_null":null,
        "original_focal_gbif_taxa_match_expected":True,
        "no_unobserved_plant_treated_as_absent":True,
        "no_island_specific_historical_native_introduced_status":True,
        "cannot_identify_local_larval_use_or_ship_introduction":True,
        "does_not_identify_true_colonization_or_island_source":True,
        "not_an_independent_predeclared_confirmatory_test":True,
        "raw_original_host_edges_or_individual_plant_point_occurrences_published":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(reported,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("ORIGINAL230_PRIOR_HOST_FIVEYEAR",json.dumps({
        "taxonomy":reported["taxonomic_match"],
        "event":reported["butterfly_event_cohort"],
        "rates":rates,
        "null":null,
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--author-roster-dir",type=Path,required=True)
    p.add_argument("--original23",type=Path,required=True)
    p.add_argument("--new15",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--taxonomy-archive",type=Path,default=None,
        help="Store just name-to-GBIF-taxonomy crosswalk in action artifact; not butterfly-host edges")
    p.add_argument("--taxonomy-from-archive",type=Path,default=None,
        help="Replay source-frozen GBIF name mapping offline with exact source checks")
    run(p.parse_args())
