#!/usr/bin/env python3
"""Paired same-island WPI-port vs nonport-coast GBIF index pilot.

Plant data: GBIF 2001-2010, same 5-km radius at WPI2019 port and
a *newly reconstructed*, distant coastal control on the same original GSHHG
sibling-ID landmass. Butterfly data: same circles 2001-2010, 2011-2017,
2018-2022; predefined source 225 butterfly species. Reports INDEXES, not
absences, plant introductions, host use, colonization or ship transport.

New coast controls are NOT the unreleased original 2021 WorldCover points.
Do not label this an exact re-use of earlier 273 site landcover pairs.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,sys,time,zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import geopandas as gpd
import shapely
from sklearn.neighbors import BallTree
from scipy.stats import rankdata,spearmanr,binomtest
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import build_public51_butterfly_risk_20261109 as public51
BASE=ROOT/"provenance/reviewer_defenses/results"
JOINT=BASE/"island_port_first_paper_20261031/joint_79port_WPI_plant_landcover_support.json"
BUTTERFLIES=BASE/"island_host_network_20261030/island_host_network_all225.csv"
PLANTS=BASE/"port_alien_broad_inventory_20261030/port_alien_host_broad_species.csv"
FOCAL=BASE/"dated_focal_host_taxonomy_20261031/focal_host_prior_to_butterfly_20261031.json"
PINS={
 JOINT:"", # Source file SHA256 recorded dynamically in receipts; pinned upstream inputs verified by original gate
 BUTTERFLIES:"39fe84981cab2b0e47063399f0a5c0c9731e800d",
}
GSHHG_SHA=public51.GSHHG_SHA256
GBIF="https://api.gbif.org/v1"
PERIODS={"pre":(2001,2010),"early":(2011,2017),"late":(2018,2022)}
RADIUS_KM=5
MIN_CONTROL_PORT_KM=15
MAX_PAIR_KM=125
MAX_PORTS_PER_ISLAND=3
MIN_PORT_PREPLANT_RECORDS=30

def rows(path):
    return list(csv.DictReader(path.read_text(encoding="utf-8-sig").splitlines()))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def source_check():
    raw=BUTTERFLIES.read_bytes()
    gitsha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if gitsha!=PINS[BUTTERFLIES]:
        raise ValueError("Original source 225 butterfly key catalogue changed "+gitsha)
    v=json.loads(JOINT.read_text())
    assert v["identical_WPI_port_and_landmass_ID_intersection"]==79
    assert v["identical_WPI_ports_in_distinct_islands"]==36
    assert len(rows(PLANTS))==687
    assert len(rows(BUTTERFLIES))==225
    return v
def haversine_km(lat,lon,lat2,lon2):
    x=np.radians(np.array(lat2,dtype=float)-lat)
    y=np.radians(np.array(lon2,dtype=float)-lon)
    z=np.sin(x/2)**2+np.cos(np.radians(lat))*np.cos(np.radians(lat2))*np.sin(y/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.maximum(0,np.minimum(1,z))))

def choose_pairs(gshhg_zip,wpi_geojson,all_joint):
    if sha(gshhg_zip)!=GSHHG_SHA:raise ValueError("GSHHG 2.3.7 archive SHA mismatch")
    with zipfile.ZipFile(gshhg_zip) as z:
        members=[p for p in z.namelist() if "/h/GSHHS_h_L1." in "/"+p
                 and p.endswith((".shp",".shx",".dbf",".prj"))]
        if len(members)<3:raise ValueError("Missing GSHHG high-resolution original coast")
        folder=gshhg_zip.parent/"shore"
        for p in members:z.extract(p,folder)
    land=gpd.read_file(next((gshhg_zip.parent/"shore").rglob("GSHHS_h_L1.shp")))
    land["landmass_id"]=land.sibling_id.astype(int).astype(str)
    features=json.loads(wpi_geojson.read_text())["features"]
    assert len(features)==3630
    coords={int(f["properties"]["OBJECTID"]):tuple(map(float,f["geometry"]["coordinates"][:2])) for f in features}
    assert len(coords)==3630
    all_xy=np.array(list(coords.values()))
    nearest=BallTree(np.radians(all_xy[:,::-1]),metric="haversine")
    accepted=defaultdict(list)
    for r in all_joint["exact_linked_source_2021_landcover_and_2001_2024_botanical_port_IDs"]:
        if int(r["original_port_2001_2010_GBIF_indexed_plants_query_count"])>=MIN_PORT_PREPLANT_RECORDS:
            accepted[r["landmass_id"]].append(r)
    selected=[]
    for lid,rr in sorted(accepted.items(),key=lambda q:int(q[0])):
        selected.extend(sorted(rr,key=lambda x:(-int(x["original_port_2001_2010_GBIF_indexed_plants_query_count"]),int(x["port_id"])))[:MAX_PORTS_PER_ISLAND])
    selected_lids={r["landmass_id"] for r in selected}
    polys={lid:shapely.union_all(land.loc[land.landmass_id==lid,"geometry"].to_numpy()) for lid in selected_lids}
    if any(poly.is_empty or not poly.is_valid for poly in polys.values()):
        raise ValueError("Invalid original GSHHG sibling-ID landmass geometry")
    # Controls were fixed BEFORE GBIF response checks. Uniform coastline
    # interpolation is an index-matched spatial comparator, NOT urban matching.
    available={}
    for lid,poly in polys.items():
        boundary=poly.boundary
        count=1400
        candidates=[boundary.interpolate((k+.5)/count,normalized=True) for k in range(count)]
        xy=np.array([(p.x,p.y) for p in candidates])
        distance=nearest.query(np.radians(xy[:,::-1]),k=1)[0][:,0]*6371.0088
        available[lid]=[tuple(z) for z,d in zip(xy,distance) if d>=MIN_CONTROL_PORT_KM]
    pairs=[];excluded=[];used=defaultdict(list)
    for r in sorted(selected,key=lambda x:(int(x["landmass_id"]),-int(x["original_port_2001_2010_GBIF_indexed_plants_query_count"]))):
        lid=r["landmass_id"];pid=int(r["port_id"])
        if pid not in coords:raise ValueError("Missing NOAA source WPI port "+str(pid))
        lon,lat=coords[pid]
        options=[]
        for ctrl_lon,ctrl_lat in available[lid]:
            km=float(haversine_km(lat,lon,ctrl_lat,ctrl_lon))
            if not (MIN_CONTROL_PORT_KM<=km<=MAX_PAIR_KM):continue
            if any(float(haversine_km(ctrl_lat,ctrl_lon,ex[1],ex[0]))<12 for ex in used[lid]):continue
            options.append((km,ctrl_lon,ctrl_lat))
        if not options:
            excluded.append({"landmass_id":lid,"port_id":pid,"reason":"NO_REMOTE_SAME_ISLAND_COAST_15_125KM"})
            continue
        km,ctrl_lon,ctrl_lat=min(options)
        used[lid].append((ctrl_lon,ctrl_lat))
        pairs.append({"landmass_id":lid,"port_id":pid,
          "port_lat":lat,"port_lon":lon,"control_lat":ctrl_lat,"control_lon":ctrl_lon,
          "same_island_port_control_distance_km":round(km,5),
          "old_port_2001_2010_plant_indexed_records":int(r["original_port_2001_2010_GBIF_indexed_plants_query_count"]),
          "prior_2021_landcover_port_minus_own_old_nonport_built_fraction":float(r["landcover_2021_built_difference_1km_near_minus_own_island_nonport"]),
          "not_exact_old_2021_control_point":True})
    return pairs,excluded,{"source_79_WPI_ports":79,"original_islands_with_79":36,
      "early_plant_min30_candidate_ports":len(selected),
      "early_plant_min30_candidate_islands":len(selected_lids),
      "new_geometric_controls":len(pairs),"new_geometric_islands":len(set(x["landmass_id"] for x in pairs)),
      "controls_distinct_from_2021_worldcover_coastal_samples":True}

def match_plant(row):
    name=row["host_name"]
    try:
        x=public51.get_json(GBIF+"/species/match",{"name":name,"rank":"SPECIES","kingdom":"Plantae","strict":"true"})
        canon=str(x.get("canonicalName") or "").split(" ")[0:2]
        orig=name.split(" ")[0:2]
        good=(str(x.get("matchType","")).upper()=="EXACT"
          and str(x.get("rank","")).upper()=="SPECIES"
          and len(canon)==2 and [s.lower() for s in canon]==[s.lower() for s in orig])
        key=x.get("acceptedUsageKey") or x.get("usageKey")
        return name,str(key) if good and key else "",str(x.get("matchType",""))
    except Exception as e:
        return name,"","ERROR:"+str(e)[:110]
def plants_crosswalk(limit):
    species=rows(PLANTS)[:limit]
    result=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        todo=[ex.submit(match_plant,s) for s in species]
        for f in as_completed(todo):result.append(f.result())
    matches={name:key for name,key,st in result if key}
    return matches,{"source_candidate_plants":len(rows(PLANTS)),"queried_candidate_plants":len(species),
      "strict_matched_accepted_gbif_keys":len(matches),
      "unmatched_plants":[{"name":a,"status":c} for a,b,c in sorted(result) if not b]}

def site_facets(lat,lon,period,taxonfilter):
    params={"geoDistance":f"{lat:.7f},{lon:.7f},{RADIUS_KM}km",
      "year":f"{PERIODS[period][0]},{PERIODS[period][1]}",
      "hasCoordinate":"true","hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
      "facet":"speciesKey","facetLimit":"1000","limit":"0"}
    if taxonfilter=="plants":params["taxonKey"]="6"
    else:
        if taxonfilter["method"]=="GBIF_SUPERFAMILY_TAXONKEY":
            params["taxonKey"]=taxonfilter["GBIF_taxonKey"]
        else:
            params["taxonKey"]=list(taxonfilter["family_taxonkeys"].values())
    k={};offset=0;total=None
    for i in range(40):
        data=public51.get_json(GBIF+"/occurrence/search",{**params,"facetOffset":offset})
        n=int(data["count"])
        if total is not None and n!=total:raise ValueError("GBIF changed inside facet retrieval")
        total=n
        facets=next((x.get("counts",[]) for x in data.get("facets",[])
                   if str(x.get("field","")).upper()=="SPECIES_KEY"),[])
        if n>0 and not facets:raise ValueError("Nonzero source index returned no speciesKey facet")
        for v in facets:
            key=str(v.get("name",""))
            if key in k:raise ValueError("GBIF duplicate speciesKey facets")
            if not key.isdecimal():raise ValueError("GBIF facet not speciesKey: "+key)
            k[key]=int(v["count"])
        if len(facets)<1000:break
        offset+=1000
    else:raise ValueError("More than 40 facet pages")
    if sum(k.values())>total:raise ValueError("Facet counts exceed GBIF records")
    return k,total

def summarize(pairs,flora,known):
    bflies={v["species"]:v["taxon_key"] for v in rows(BUTTERFLIES)}
    focal=json.loads(FOCAL.read_text())["focal_links"]
    site=[];pairrows=[];focalrows=[]
    def measure(siteid):
        facts=flora[siteid]
        pre,pre_total=facts["plant_pre"]
        bs,bs_total=facts["butterfly_pre"]
        early,early_total=facts["butterfly_early"]
        late,late_total=facts["butterfly_late"]
        kk=set(known.values())
        found=kk&set(pre)
        pre_target=sum(1 for k in bflies.values() if k in bs)
        newly=[key for key in bflies.values() if key not in bs and key in early]
        last=[key for key in bflies.values() if key not in bs and key not in early and key in late]
        return {"gbif_plant_pre_total_records":pre_total,
         "gbif_plant_pre_specieskeys":len(pre),
         "gbif_host_candidate_pre_specieskeys":len(found),
         "gbif_host_candidate_pre_fraction":len(found)/len(pre) if pre else None,
         "gbif_butterfly_pre_total_records":bs_total,
         "gbif_butterfly_pre_225_species":pre_target,
         "gbif_butterfly_early_total_records":early_total,
         "gbif_butterfly_new_2011_2017_225_species":len(newly),
         "gbif_butterfly_late_total_records":late_total,
         "gbif_butterfly_new_2018_2022_225_species":len(last)}
    for p in pairs:
        pid=str(p["port_id"])
        l=measure(pid+":port");r=measure(pid+":nonport")
        site.append({"landmass_id":p["landmass_id"],"port_id":p["port_id"],"side":"port",**l})
        site.append({"landmass_id":p["landmass_id"],"port_id":p["port_id"],"side":"nonport",**r})
        def delta(col):return float(np.log1p(l[col])-np.log1p(r[col]))
        pairrows.append({"landmass_id":p["landmass_id"],"port_id":p["port_id"],
          "distance_port_nonport_km":p["same_island_port_control_distance_km"],
          "delta_log_pre_host_candidate_species":delta("gbif_host_candidate_pre_specieskeys"),
          "delta_log_pre_all_plant_species":delta("gbif_plant_pre_specieskeys"),
          "delta_log_pre_plant_total_records":delta("gbif_plant_pre_total_records"),
          "delta_log_pre_butterfly_225_species":delta("gbif_butterfly_pre_225_species"),
          "delta_log_new_butterfly_2011_2017_225_species":delta("gbif_butterfly_new_2011_2017_225_species"),
          "delta_log_new_butterfly_2018_2022_225_species":delta("gbif_butterfly_new_2018_2022_225_species"),
          "port_pre_host_species":l["gbif_host_candidate_pre_specieskeys"],
          "nonport_pre_host_species":r["gbif_host_candidate_pre_specieskeys"],
          "port_later_new_butterflies":l["gbif_butterfly_new_2011_2017_225_species"],
          "nonport_later_new_butterflies":r["gbif_butterfly_new_2011_2017_225_species"]})
        for v in focal:
            butterfly=v["butterfly"];host=v["gbif_accepted_plant_speciesKey"];bk=bflies.get(butterfly)
            if not bk:continue
            for side in ("port","nonport"):
                f=flora[pid+":"+side]
                focalrows.append({"landmass_id":p["landmass_id"],"port_id":p["port_id"],"side":side,
                  "butterfly":butterfly,"compatible_plant":v["source_catalogue_compatible_host"],
                  "compatible_plant_recorded_2001_2010":host in f["plant_pre"][0],
                  "butterfly_indexed_pre2001_2010":bk in f["butterfly_pre"][0],
                  "butterfly_indexed_2011_2017":bk in f["butterfly_early"][0],
                  "butterfly_indexed_2018_2022":bk in f["butterfly_late"][0]})
    by_island=defaultdict(list)
    for p in pairrows:by_island[p["landmass_id"]].append(p)
    fields=["delta_log_pre_host_candidate_species","delta_log_pre_all_plant_species",
        "delta_log_pre_plant_total_records","delta_log_pre_butterfly_225_species",
        "delta_log_new_butterfly_2011_2017_225_species",
        "delta_log_new_butterfly_2018_2022_225_species"]
    island=[{"landmass_id":i,"matched_ports":len(v),
             **{f:float(np.mean([z[f] for z in v])) for f in fields}}
             for i,v in sorted(by_island.items(),key=lambda z:int(z[0]))]
    def corr(xfield,yfield,permutation=5000):
        xx=np.array([q[xfield] for q in island]);yy=np.array([q[yfield] for q in island])
        if len(xx)<5 or len(set(xx))<2 or len(set(yy))<2:return {"n_islands":len(xx),"rho":None,"p_island_shuffle":None}
        obs=float(spearmanr(xx,yy).statistic)
        rng=np.random.default_rng(20261109);null=[]
        for i in range(permutation):
            v=spearmanr(xx,rng.permutation(yy)).statistic
            null.append(abs(float(v))>=abs(obs)-1e-12)
        return {"n_islands":len(xx),"rho":obs,"p_island_shuffle":(1+sum(null))/(1+permutation)}
    signs=[np.sign(v["delta_log_pre_host_candidate_species"]) for v in island]
    pos=sum(s>0 for s in signs);neg=sum(s<0 for s in signs)
    return site,pairrows,focalrows,island,{
      "matched_islands":len(island),"matched_ports":len(pairrows),
      "host_advantage_islands_pos":pos,"host_advantage_islands_neg":neg,
      "host_advantage_islands_tied":len(island)-pos-neg,
      "host_advantage_sign_test_p_two_sided":float(binomtest(pos,pos+neg,.5).pvalue) if pos+neg else None,
      "host_delta_vs_later_new_butterfly_delta":corr(fields[0],fields[4]),
      "host_delta_vs_2018_2022_new_butterfly_delta":corr(fields[0],fields[5]),
      "plant_total_delta_vs_later_new_butterfly_delta_negative_control":corr(fields[2],fields[4]),
      "prior_butterfly_delta_vs_later_new_butterfly_delta_negative_control":corr(fields[3],fields[4]),
      "focal_compatible_host_and_subsequent_butterfly_rows":len(focalrows),
      "focal_host_pre_and_bfly_new_2011_2017":sum(x["compatible_plant_recorded_2001_2010"] and
                  not x["butterfly_indexed_pre2001_2010"] and x["butterfly_indexed_2011_2017"] for x in focalrows)}

def savecsv(path,records):
    if not records:return
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gshhg",type=Path,required=True)
    p.add_argument("--wpi",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--candidate-limit",type=int,default=687)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    src=source_check()
    pairs,excluded,support=choose_pairs(args.gshhg,args.wpi,src)
    savecsv(args.out/"fixed_same_island_port_nonport_sites.csv",pairs)
    if len(set(z["landmass_id"] for z in pairs))<5:
        raise ValueError("Less than 5 island pairs; do not run an under-supported lag test")
    print("MATCHED_SITE_GEOMETRY",json.dumps(support),flush=True)
    matched,crosswalk=plants_crosswalk(args.candidate_limit)
    # The two previously identified focal taxa remain strictly independent
    # of the broad global candidate crosswalk; their original GBIF keys are
    # retained only for the transparent 3-butterfly focal index audit.
    filter_=public51.resolve_butterfly_taxon_filter()
    flora={};done=0
    for x in pairs:
        for side in ("port","nonport"):
            siteid=str(x["port_id"])+":"+side
            lat=x["port_lat"] if side=="port" else x["control_lat"]
            lon=x["port_lon"] if side=="port" else x["control_lon"]
            flora[siteid]={}
            for kind,period in [("plant","pre"),("butterfly","pre"),
                                ("butterfly","early"),("butterfly","late")]:
                v,total=site_facets(lat,lon,period,"plants" if kind=="plant" else filter_)
                flora[siteid][kind+"_"+period]=(v,total)
            done+=1
            print("GBIF_SITE",done,"/",2*len(pairs),siteid,
                  "plants",flora[siteid]["plant_pre"][1],
                  "butterflies",flora[siteid]["butterfly_early"][1],flush=True)
    site,pair_rows,focal,islands,model=summarize(pairs,flora,matched)
    savecsv(args.out/"port_nonport_plant_butterfly_site_indices.csv",site)
    savecsv(args.out/"matched_pair_host_lag_differences.csv",pair_rows)
    savecsv(args.out/"island_equal_weight_host_lag_differences.csv",islands)
    savecsv(args.out/"focal_named_butterfly_host_time_order.csv",focal)
    audit={"schema":"chocho_same_island_port_vs_nonport_coast_historical_host_butterfly_index_pilot_v1",
       "status":"EMPIRICAL_OBSERVATION_INDEX_NOT_HOST_USE_OR_COLONIZATION",
       "query_time_utc":datetime.now(timezone.utc).isoformat(),
       "data_source_sha256":{"joint79":sha(JOINT),"source225":sha(BUTTERFLIES),
                             "candidate687":sha(PLANTS),"focal":sha(FOCAL),
                             "GSHHG237":sha(args.gshhg),"NOAA_WPI_2019":sha(args.wpi)},
       "design":support,"excluded_no_geographic_control":excluded,
       "plant_catalogue_taxonomy":crosswalk,
       "gbif_butterfly_family_filter":filter_,
       "radius_km_port_and_nonport":5,"plant_baseline_period":"2001-2010",
       "butterfly_periods":["2001-2010","2011-2017","2018-2022"],
       "analysis":model,
       "limitations":["Candidate introduced host catalogue is worldwide source label; neither local alien status nor historical introductions are independently identified",
         "Broad 687 host catalogue combines all Lepidoptera; host identities for all 225 butterfly species cannot be recovered from public aggregate: broad host-richness correlation is NOT a species-by-species compatible-host test",
         "Three focal butterfly-host links are species-specific but observational and sparse",
         "New controls are newly reconstructed remote coastline points, NOT the earlier unreleased WorldCover point-matched controls",
         "GBIF species name index changes, positional errors and GBIF survey bias remain",
         "No observed larval feeding, butterfly breeding, ship cargo, boat itinerary, true absence or causal colonisation is measured",
         "Port vs remote-coast habitats differ; land cover, effort, climate and human disturbance not fully matched",
         "Island is the independent bootstrap/permutation sampling unit; multiple ports within same island are not independent"]}
    (args.out/"audit.json").write_text(json.dumps(audit,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print("ISLAND_PORT_NONPORT_LAG_PILOT",json.dumps(model,sort_keys=True),flush=True)
if __name__=="__main__":main()
