#!/usr/bin/env python3
"""Within-plant-family butterfly-host link rewiring; potential-resource geography."""
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict,Counter
from pathlib import Path
from statistics import median
from shapely.geometry import shape

def desc(path):
    with path.open(newline="",encoding="utf8") as f:
        return {r["species"]:(float(r["host_family_count"]),int(r["host_wgsrpd3_unit_count"])) for r in csv.DictReader(f)}
def interactions(path,panel):
    by=defaultdict(set); family={}
    with path.open(newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            sp=r["insect_species"].strip();h=r["accepted_plant_name_id"].strip()
            if sp not in panel or not h:continue
            fam=r["family"].strip()
            if h in family and family[h]!=fam:raise RuntimeError("conflicting botanical family")
            family[h]=fam;by[sp].add(h)
    return by,family
def load_butterfly_taxonomy(path):
    taxa={}
    with path.open(newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            name=str(r.get("Species") or "").strip()
            fam=str(r.get("Family") or "").strip()
            if name and fam:taxa[name]=fam
    return taxa

def dist(path):
    out=defaultdict(set)
    with path.open(newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            h=r["accepted_plant_name_id"].strip();u=r["area_code_l3"].strip()
            if h and u:out[h].add(u)
    return out
def geo(path):
    out={}
    for f in json.load(path.open(encoding="utf8"))["features"]:
        p=f["properties"];name=str(p["LEVEL3_COD"]);point=shape(f["geometry"]).representative_point()
        out[name]=(str(p["LEVEL1_COD"]),point.y,point.x)
    return out
def climate(path):
    vals={}
    with path.open(newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            k=r["wgsrpd3_code"].strip()
            try:
                x=(float(r["bio1"]),float(r["bio7"]),float(r["bio12"]),float(r["bio15"]))
            except (ValueError,TypeError,KeyError):continue
            if x[2]>=65535 or x[2]<0:continue
            x=(x[0],x[1],math.log1p(x[2]),x[3])
            if k in vals and any(abs(a-b)>1e-8 for a,b in zip(x,vals[k])):raise RuntimeError("climate duplicate drift")
            vals[k]=x
    if len(vals)<340:raise RuntimeError("climate coverage drift")
    mu=[sum(v[j] for v in vals.values())/len(vals) for j in range(4)]
    sd=[math.sqrt(sum((v[j]-mu[j])**2 for v in vals.values())/(len(vals)-1)) for j in range(4)]
    return {k:tuple((v[j]-mu[j])/sd[j] for j in range(4)) for k,v in vals.items()}
def hav(a,b):
    lat1=math.radians(a[1]);lat2=math.radians(b[1])
    dlat=lat2-lat1;dlon=math.radians(((b[2]-a[2]+180)%360)-180)
    h=math.sin(dlat/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return 2*6371.0088*math.asin(min(1,math.sqrt(h)))
def q(a,p):
    s=sorted(a);x=(len(s)-1)*p;i=int(math.floor(x));j=int(math.ceil(x))
    return s[i] if i==j else s[i]*(j-x)+s[j]*(x-i)
def region_masks(rows,plant_native,plant_contemp,nreg):
    n=[0]*nreg;c=[0]*nreg
    for s,hosts in enumerate(rows):
        nn=cc=0
        for h in hosts:
            nn|=plant_native[h];cc|=plant_contemp[h]
        bit=1<<s
        while nn:
            x=nn&-nn;n[x.bit_length()-1]|=bit;nn-=x
        while cc:
            x=cc&-cc;c[x.bit_length()-1]|=bit;cc-=x
    return n,c
def gain(n,c,pairs):
    a=b=0.0
    for i,j in pairs:
        x=n[i];y=n[j];den=(x|y).bit_count()
        a+=(x&y).bit_count()/den if den else 0.0
        x=c[i];y=c[j];den=(x|y).bit_count()
        b+=(x&y).bit_count()/den if den else 0.0
    return (b-a)/len(pairs)
def main():
    ap=argparse.ArgumentParser()
    for opt in ["descriptors_csv","insect_host_csv","native_distribution_csv","contemporary_distribution_csv","level3_geojson","climate_crossfit_csv","output_json"]:
        ap.add_argument("--"+opt.replace("_","-"),type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=499)
    ap.add_argument("--seed",type=int,default=20261008)
    ap.add_argument('--leptraits-csv',type=Path,required=True)
    a=ap.parse_args()
    d=desc(a.descriptors_csv);species=[s for s,(family_count,n) in d.items() if family_count>0 and n>0]
    if len(species)!=239:raise RuntimeError("239-species panel drift")
    links,family=interactions(a.insect_host_csv,set(species))
    butterfly_tax=load_butterfly_taxonomy(a.leptraits_csv)
    if any(s not in butterfly_tax for s in species):raise RuntimeError('missing butterfly family')
    plant_ids=sorted(set().union(*(links[s] for s in species)))
    if len(plant_ids)!=1706:raise RuntimeError(f"plant panel drift: {len(plant_ids)}")
    plants={p:i for i,p in enumerate(plant_ids)}
    nr=dist(a.native_distribution_csv);cr=dist(a.contemporary_distribution_csv)
    native_by_sp={}
    for sp in species:
        footprint=set()
        for h in links[sp]:footprint|=nr.get(h,set())
        if len(footprint)!=d[sp][1]:raise RuntimeError("native species footprint drift")
        native_by_sp[sp]=footprint
    regions=sorted(set().union(*native_by_sp.values()))
    if len(regions)!=355:raise RuntimeError(f"region panel drift {len(regions)}")
    reg_ix={r:i for i,r in enumerate(regions)}
    pn=[];pc=[]
    for h in plant_ids:
        native=sum(1<<reg_ix[u] for u in nr.get(h,set()) if u in reg_ix)
        contemp=sum(1<<reg_ix[u] for u in cr.get(h,set()) if u in reg_ix)
        if native&~contemp:raise RuntimeError("nonmonotonic plant geography")
        pn.append(native);pc.append(contemp)
    g=geo(a.level3_geojson);env=climate(a.climate_crossfit_csv)
    far=[]
    for i,r in enumerate(regions):
        for j in range(i+1,len(regions)):
            s=regions[j]
            if r not in env or s not in env:continue
            if g[r][0]==g[s][0] or hav(g[r],g[s])<3000:continue
            distance=math.sqrt(sum((x-y)**2 for x,y in zip(env[r],env[s])))
            far.append((distance,i,j))
    lo=q([x[0] for x in far],.25);hi=q([x[0] for x in far],.75)
    analog=[(i,j) for dist,i,j in far if dist<=lo]
    discord=[(i,j) for dist,i,j in far if dist>=hi]
    if len(analog)!=12944 or len(discord)!=12944:raise RuntimeError("climate group drift")
    rows=[{plants[h] for h in links[s]} for s in species]
    edges=[(i,h) for i,r in enumerate(rows) for h in sorted(r)]
    if len(edges)!=2596:raise RuntimeError("host link count drift")
    fam_groups=defaultdict(list)
    for edge_index,(i,h) in enumerate(edges):fam_groups[(butterfly_tax[species[i]],family[plant_ids[h]])].append(edge_index)
    eligible=[ids for ids in fam_groups.values() if len(ids)>=2]
    rng=random.Random(a.seed)
    def jmean(masks,pairs):
        total=0.0
        for i,j in pairs:
            aa=masks[i];bb=masks[j];den=(aa|bb).bit_count()
            total+=(aa&bb).bit_count()/den if den else 0
        return total/len(pairs)
    def components(rowsets):
        n,c=region_masks(rowsets,pn,pc,len(regions))
        na=jmean(n,analog);nd=jmean(n,discord)
        ca=jmean(c,analog);cd=jmean(c,discord)
        return {
            "native_analog":na,"native_discordant":nd,
            "contemporary_analog":ca,"contemporary_discordant":cd,
            "native_selectivity":na-nd,
            "contemporary_selectivity":ca-cd,
            "introduction_delta_selectivity":(ca-cd)-(na-nd)
        }
    observed_components=components(rows)
    observed=observed_components["introduction_delta_selectivity"]
    if abs(observed-0.09843700070646616)>1e-9:raise RuntimeError(f"original network contrast drift {observed}")
    original_family_degree=[Counter(family[plant_ids[h]] for h in r) for r in rows]
    original_host_degree=Counter(h for _,h in edges)
    accepted=0
    def swaps(n):
        nonlocal accepted
        for _ in range(n):
            ids=eligible[rng.randrange(len(eligible))]
            k1=ids[rng.randrange(len(ids))];k2=ids[rng.randrange(len(ids))]
            if k1==k2:continue
            i,h1=edges[k1];j,h2=edges[k2]
            if i==j or h1==h2 or h2 in rows[i] or h1 in rows[j]:continue
            rows[i].remove(h1);rows[i].add(h2)
            rows[j].remove(h2);rows[j].add(h1)
            edges[k1]=(i,h2);edges[k2]=(j,h1)
            accepted+=1
    burn=max(10000,15*len(edges));spacing=max(2000,3*len(edges))
    swaps(burn)
    sims=[];null_components=defaultdict(list)
    for step in range(a.permutations):
        swaps(spacing)
        comp=components(rows)
        sims.append(comp["introduction_delta_selectivity"])
        for key,val in comp.items():null_components[key].append(val)
        if (step+1)%100==0:print(json.dumps({"completed":step+1,"accepted_swaps":accepted}),flush=True)
    if accepted<1000:raise RuntimeError(f'insufficient taxonomically constrained edge swaps: {accepted}')
    if [Counter(family[plant_ids[h]] for h in r) for r in rows]!=original_family_degree:
        raise RuntimeError("butterfly x botanical family host degree drift")
    if Counter(h for _,h in edges)!=original_host_degree:
        raise RuntimeError("plant host-user degree drift")
    if any(len(r)!=len(set(r)) for r in rows):raise RuntimeError("duplicate host")
    decomposition={}
    for key,val in observed_components.items():
        sims_for_key=null_components[key]
        med=median(sims_for_key)
        decomposition[key]={
            "observed":val,
            "null_median":med,
            "observed_minus_null_median":val-med,
            "null_ci95":[q(sims_for_key,.025),q(sims_for_key,.975)],
            "p_high":(1+sum(v>=val-1e-14 for v in sims_for_key))/(1+len(sims_for_key))
        }
    payload={
        "schema":"chocho_strict_host_identity_temporal_contrast_decomposition_v0.1",
        "status":"POSTHOC_DECOMPOSITION_OF_FIXED_PHYLOGENETIC_HOST_LINK_NULL",
        "protocol":"docs/exploratory/butterfly_within_consumer_family_host_identity_null_protocol_v0.2.json",
        "panel":{"butterflies":len(species),"plants":len(plant_ids),"host_links":len(edges),"botanical_families":len(fam_groups),"regions":len(regions)},
        "climate_pairs":{"far_climate_analogue":len(analog),"far_climate_discordant":len(discord),"q25":lo,"q75":hi},
        "primary":{
            "observed":observed,"null_median":median(sims),"null_ci95":[q(sims,.025),q(sims,.975)],
            "excess":observed-median(sims),
            "p_one_sided":(1+sum(x>=observed-1e-14 for x in sims))/(len(sims)+1)
        },
        "native_vs_contemporary_decomposition":decomposition,
        "interpretation_note":"An introduced-resource-specific residual is the difference between contemporary and native conditional host-link selectivity. Native host identity already shapes climate-selective resource geography; do not attribute that pre-existing signal to recent plant globalization.",
        "null":{"permutations":a.permutations,"seed":a.seed,"burn_attempts":burn,"spacing_attempts":spacing,"accepted_swaps":accepted,
                "preserves":["all host plant native and contemporary regional distributions","each butterfly host-family degree","each plant number of consumers within each butterfly Family"]},
        "interpretation_boundary":"Does not preserve each butterfly's regional opportunity count; therefore it complements rather than replaces the fixed-butterfly-margin null. Network membership is based on documented potential hosts, not realized regional host use."
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))
if __name__=="__main__":main()
