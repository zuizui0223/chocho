#!/usr/bin/env python3
"""Source-locked *spatial accounting*, not causal attribution, of LGM climate land masks.

Exact identity for every plant / insect-hostset, with native WCVP modern-trained
climate-suitable AND PMIP4 icefree land:
A21 - A0 = (A21 on land-both - A0 on land-both)
          + A21 on land-only-21ka - A0 on land-only-0BP.
This DOES NOT reconstruct historical plant occupancy, species dispersal,
genetic breaks, or causal climate-vs-geography interaction terms.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path
import numpy as np

SHA={
 "plant_primary":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
 "plant_extrapolated":"67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a",
 "PMIP4_terrain":"d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db",
 "insect_host_manifest":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
}
GRID=(360,720)
CELL_N=360*720

def sha256_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def accepted_id(value):
    d=Decimal(str(value).strip())
    if not d.is_finite() or d<=0 or d!=d.to_integral_value():
        raise ValueError("invalid WCVP accepted host taxon ID")
    return str(int(d))

def area_by_flat_cell():
    centers=np.deg2rad(89.75-.5*np.arange(GRID[0]))
    dlat=np.deg2rad(.25)
    per_row=6371.0088**2*np.deg2rad(.5)*(np.sin(centers+dlat)-np.sin(centers-dlat))
    return np.repeat(per_row,GRID[1])

def partition_area(mask0,mask21,physical0,physical21,area):
    h0=np.asarray(mask0,bool).reshape(-1)
    h21=np.asarray(mask21,bool).reshape(-1)
    land0=np.asarray(physical0,bool).reshape(-1)
    land21=np.asarray(physical21,bool).reshape(-1)
    w=np.asarray(area,float).reshape(-1)
    if not (h0.shape==h21.shape==land0.shape==land21.shape==w.shape==(CELL_N,)):
        raise ValueError("source 0.5degree grids not aligned")
    if np.any(h0&~land0) or np.any(h21&~land21):
        raise ValueError("source climatic mask includes terrain-blocked cells")
    common=land0&land21
    former_only=land0&~land21
    glacial_only=land21&~land0
    measure=lambda x:float(w[x].sum())
    common0=measure(h0&common)
    common21=measure(h21&common)
    former_lost=measure(h0&former_only)
    glacial_gained=measure(h21&glacial_only)
    total0=measure(h0)
    total21=measure(h21)
    diff=total21-total0
    common_climate=common21-common0
    domain_net=glacial_gained-former_lost
    remainder=diff-common_climate-domain_net
    if abs(remainder)>2e-5:
        raise ValueError("exact spatial area accounting identity failed")
    return {
       "potential_area_0BP_km2":total0,"potential_area_21ka_km2":total21,
       "paired_total_change_km2":diff,
       "unchanged_physically_passable_domain_climate_change_km2":common_climate,
       "21ka_only_passable_land_suitable_area_km2":glacial_gained,
       "0BP_only_passable_land_suitable_area_km2":former_lost,
       "net_passable_terrain_replacement_km2":domain_net,
       "area_identity_remainder_km2":remainder
    }

def describe(rows):
    if not rows:raise ValueError("no qualifying plant or complete known hostset")
    get=lambda k:np.array([float(r[k]) for r in rows],float)
    whole=get("paired_total_change_km2")
    common=get("unchanged_physically_passable_domain_climate_change_km2")
    domain=get("net_passable_terrain_replacement_km2")
    def negative(k):return int((get(k)<0).sum())
    return {"n":len(rows),
        "full_suitable_area_decrease":int((whole<0).sum()),
        "full_suitable_area_increase":int((whole>0).sum()),
        "unchanged_land_climate_decrease":int((common<0).sum()),
        "unchanged_land_climate_increase":int((common>0).sum()),
        "terrain_replacement_net_negative":int((domain<0).sum()),
        "terrain_replacement_net_positive":int((domain>0).sum()),
        "both_climate_and_terrain_negative":int(((common<0)&(domain<0)).sum()),
        "climate_negative_terrain_positive":int(((common<0)&(domain>0)).sum()),
        "climate_positive_terrain_negative":int(((common>0)&(domain<0)).sum()),
        "climate_positive_terrain_positive":int(((common>0)&(domain>0)).sum()),
        "climate_term_abs_greater_than_domain_term":int((np.abs(common)>np.abs(domain)).sum()),
        "median_full_paired_area_delta_km2":float(np.median(whole)),
        "median_common_land_climate_delta_km2":float(np.median(common)),
        "median_net_terrain_replacement_delta_km2":float(np.median(domain)),
        "median_21ka_only_passable_land_suitable_km2":float(np.median(get("21ka_only_passable_land_suitable_area_km2"))),
        "median_0BP_only_passable_land_suitable_km2":float(np.median(get("0BP_only_passable_land_suitable_area_km2"))),
        "max_abs_identity_remainder_km2":float(np.max(np.abs(get("area_identity_remainder_km2"))))}

def rows_to_csv(path,rows):
    if not rows:raise ValueError("empty source response")
    fields=sorted(set().union(*(x.keys() for x in rows)))
    with Path(path).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader();w.writerows(rows)

def source_plan(args):
    for label,path in (("plant_primary",args.plant_primary),("plant_extrapolated",args.plant_extrapolated),
                       ("PMIP4_terrain",args.terrain),("insect_host_manifest",args.manifest)):
        if sha256_file(path)!=SHA[label]:
            raise ValueError("original SHA-locked source changed: "+label)
    with np.load(args.plant_primary,allow_pickle=False) as q:
        plant_ids=[accepted_id(x) for x in q["plant_wcvp_ids"]]
        m0=np.asarray(q["mask_0BP"],dtype=np.uint8)
        m21=np.asarray(q["mask_21ka"],dtype=np.uint8)
        assert bool(q["nonanalog_excluded"])
        if tuple(q["grid_shape"].tolist())!=GRID or int(q["land_cut"])!=5:
            raise ValueError("source terrain and plant threshold changed")
    with np.load(args.plant_extrapolated,allow_pickle=False) as q:
        ex_ids=[accepted_id(x) for x in q["plant_wcvp_ids"]]
        e0=np.asarray(q["mask_0BP"],dtype=np.uint8)
        e21=np.asarray(q["mask_21ka"],dtype=np.uint8)
        assert not bool(q["nonanalog_excluded"])
    if (not(len(plant_ids)==len(set(plant_ids))==416) or plant_ids!=ex_ids
        or m0.shape!=m21.shape!=e0.shape!=e21.shape!=(416,32400)):
        raise ValueError("modelled 416 accepted host sources and scenarios not aligned")
    if np.any(m0&~e0) or np.any(m21&~e21):
        raise ValueError("extrapolation scenario is not a superset of original source habitat")
    with np.load(args.terrain,allow_pickle=False) as q:
        p0=np.asarray(q["pmip0_icefree_subcells"]).reshape(-1)>=5
        p21=np.asarray(q["pmip21_icefree_subcells"]).reshape(-1)>=5
    if not (p0.shape==p21.shape==(CELL_N,) and (int(p0.sum()),int(p21.sum()))==(58765,57186)):
        raise ValueError("PMIP4 exact source area mask changed")
    links=defaultdict(set)
    with args.manifest.open("r",newline="",encoding="utf8") as stream:
        for r in csv.DictReader(stream):
            links[r["insect_species"]].add(accepted_id(r["wcvp_id_as_frozen"]))
    if len(links)!=641:raise ValueError("source 641 insect species incomplete")
    index={x:i for i,x in enumerate(plant_ids)}
    by_hostset=defaultdict(list)
    for insect,hostset in links.items():
        if hostset<=index.keys():
            by_hostset[tuple(sorted(hostset,key=int))].append(insect)
    if not(sum(len(v) for v in by_hostset.values())==513 and len(by_hostset)==349):
        raise ValueError("fixed 513 insects or 349 unique accepted known host combinations drift")
    return plant_ids,(m0,m21),(e0,e21),(p0,p21),links,index,by_hostset

def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    plants,primary,extrap,physical,links,idx,groups=source_plan(args)
    w=area_by_flat_cell()
    all_plants=[];all_insects=[];unique=[]
    for scenario,matrix in (("analog",primary),("extrapolated",extrap)):
        x0=np.unpackbits(matrix[0],axis=1,count=CELL_N).astype(bool)
        x21=np.unpackbits(matrix[1],axis=1,count=CELL_N).astype(bool)
        for pi,pid in enumerate(plants):
            all_plants.append({"wcvp_id":pid,"scenario":scenario,
                **partition_area(x0[pi],x21[pi],*physical,w)})
        complete=set()
        for hosts,species in sorted(groups.items()):
            idlist=[idx[h] for h in hosts]
            metric=partition_area(np.any(x0[idlist],axis=0),np.any(x21[idlist],axis=0),*physical,w)
            key=";".join(hosts)
            unique.append({"hostset":key,"scenario":scenario,
                           "insects_reusing_hostset":len(species),**metric})
            for sp in species:
                complete.add(sp)
                all_insects.append({"insect_species":sp,"status":"COMPLETE_ALL_RECORDED_HOST_MODELS",
                                    "hostset":key,"scenario":scenario,
                                    "host_taxa_n":len(hosts),"insects_reusing_hostset":len(species),
                                    **metric})
        for sp in sorted(set(links)-complete):
            all_insects.append({"insect_species":sp,"status":"INCOMPLETE_RECORDED_HOST_MODELS_UNSCORED",
                                "hostset":";".join(sorted(links[sp],key=int)),"scenario":scenario,
                                "host_taxa_n":len(links[sp]),"insects_reusing_hostset":None})
    if not(len(all_plants)==832 and len(all_insects)==1282 and len(unique)==698):
        raise ValueError("original WCVP/TTF source qualification not preserved")
    report={"schema":"chocho_513_insects_349_hostsets_416_plants_0BP21ka_common_land_accounting_v07",
       "status":"EXPLORATORY_EXACT_SPATIAL_ACCOUNTING_NOT_CAUSAL_FACTORIAL_OR_HISTORICAL_GENETICS",
       "original_insects":641,"qualifying_insects":513,
       "insects_missing_any_known_host_model":128,
       "qualified_distinct_hostsets":349,"modelled_native_plants":416,
       "grid_cells":CELL_N,
       "icefree_land_0BP_cells":int(physical[0].sum()),
       "icefree_land_21ka_cells":int(physical[1].sum()),
       "icefree_land_common_cells":int((physical[0]&physical[1]).sum()),
       "icefree_land_only_0BP_cells":int((physical[0]&~physical[1]).sum()),
       "icefree_land_only_21ka_cells":int((~physical[0]&physical[1]).sum()),
       "spatial_accounting_identity":"Total area(21ka)-area(0BP) = (21ka-0BP suitability on COMMON 0BP+21ka ice-free cells) + (21ka suitability on only-21ka ice-free cells) - (0BP suitability on only-0BP ice-free cells)",
       "NOT_causal_climate_vs_land_two_by_two":True,
       "source_sha256":dict(SHA),
       "strata":{},
       "limits":["Exact terrain-domain accounting is not a causal climate-vs-land/ice main-effects or interaction decomposition",
          "Modern WCVP native-trained plant climate suitability is hypothetical; PMIP4 icefree cell majority masks are source coarse-grid land estimates",
          "No ancient plant or insect occupation, source genetics, dispersal or host-specific genetic effect read",
          "513 insects depend on only 349 accepted host combinations; all within-insect edge correlations remain",
          "Original TTF actual-host vs modern-climate-matched other-plant negative v0.5 is unchanged and not rescued"]}
    for scenario in ("analog","extrapolated"):
        plants0=[r for r in all_plants if r["scenario"]==scenario]
        insects0=[r for r in all_insects if r["scenario"]==scenario and r["status"].startswith("COMPLETE")]
        groups0=[r for r in unique if r["scenario"]==scenario]
        report["strata"][scenario]={"plants":describe(plants0),
            "insects":describe(insects0),"distinct_hostsets":describe(groups0)}
    assert report["strata"]["analog"]["insects"]["full_suitable_area_decrease"]==399
    assert report["strata"]["extrapolated"]["insects"]["full_suitable_area_decrease"]==322
    assert report["strata"]["analog"]["plants"]["full_suitable_area_decrease"]==289
    assert report["strata"]["extrapolated"]["plants"]["full_suitable_area_decrease"]==243
    rows_to_csv(args.out/"all_416_plants_land_climate_partition_v07.csv",all_plants)
    rows_to_csv(args.out/"all_641_insects_land_climate_partition_v07.csv",all_insects)
    rows_to_csv(args.out/"349_unique_hostsets_land_climate_partition_v07.csv",unique)
    (args.out/"summary.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))

def main():
    p=argparse.ArgumentParser()
    for field in ("plant-primary","plant-extrapolated","terrain","manifest","out"):
        p.add_argument("--"+field,required=True,type=Path)
    run(p.parse_args())
if __name__=="__main__":main()
