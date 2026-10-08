#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict
from pathlib import Path
from statistics import median
from shapely.geometry import shape

def load_desc(p):
    out={}
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            out[r['species']]={'b':float(r['host_family_count']),'n':int(r['host_wgsrpd3_unit_count'])}
    return out

def load_pairs(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            s=r['insect_species'].strip(); h=r['accepted_plant_name_id'].strip()
            if s and h: out[s].add(h)
    return out

def load_units(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            h=r['accepted_plant_name_id'].strip(); u=r['area_code_l3'].strip()
            if h and u: out[h].add(u)
    return out

def load_level1(p):
    data=json.load(open(p,encoding='utf-8'))
    out={}
    for feat in data.get('features',[]):
        props={str(k).upper():v for k,v in (feat.get('properties') or {}).items()}
        l3=props.get('LEVEL3_COD') or props.get('LEVEL3_CODE') or props.get('LEVEL3')
        l1=props.get('LEVEL1_COD') or props.get('LEVEL1_CODE') or props.get('LEVEL1')
        if l3 is None or l1 is None: continue
        l3=str(l3).strip()
        # Some GIS exports may render numeric codes as floats; keep level1 as a string label only.
        l1=str(l1).strip()
        if l3:
            if l3 in out and out[l3]!=l1: raise RuntimeError(f'conflicting level1 for {l3}')
            out[l3]=l1
    if not out: raise RuntimeError('no WGSRPD level3-to-level1 mapping recovered')
    return out

def masks(species,regions,by_species):
    idx={u:i for i,u in enumerate(regions)}
    m=[0]*len(regions)
    for i,s in enumerate(species):
        bit=1<<i
        for u in by_species[s]:
            if u in idx:m[idx[u]]|=bit
    return m

def mean_jaccard_masks(ms):
    vals=[]
    for i,a in enumerate(ms):
        for b in ms[i+1:]:
            z=(a|b).bit_count()
            if z: vals.append((a&b).bit_count()/z)
    return sum(vals)/len(vals)

def mean_jaccard_sets(ss):
    vals=[]
    for i,a in enumerate(ss):
        for b in ss[i+1:]:
            z=len(a|b)
            if z: vals.append(len(a&b)/z)
    return sum(vals)/len(vals)

def q(xs,p):
    xs=sorted(xs); z=p*(len(xs)-1); i=int(math.floor(z)); j=int(math.ceil(z))
    if i==j:return xs[i]
    a=z-i; return xs[i]*(1-a)+xs[j]*a


def load_geography(path):
    data=json.load(open(path,encoding="utf-8"))
    out={}
    for f in data["features"]:
        p=f["properties"]; code=str(p["LEVEL3_COD"]).strip()
        pos=shape(f["geometry"]).representative_point()
        out[code]={"level1":str(p["LEVEL1_COD"]),"lat":pos.y,"lon":pos.x}
    return out

def haversine(a,b):
    p1=math.radians(a["lat"]);p2=math.radians(b["lat"])
    dl=p2-p1;dn=math.radians(((b["lon"]-a["lon"]+180)%360)-180)
    h=math.sin(dl/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dn/2)**2
    return 2*6371.0088*math.asin(min(1,math.sqrt(h)))

def grouped_jaccard(masks,groups):
    means={}
    for name,arr in groups.items():
        total=0.0
        for i,j in arr:
            a=masks[i];b=masks[j];u=(a|b).bit_count()
            total+=(a&b).bit_count()/u if u else 0.0
        means[name]=total/len(arr)
    return means


def load_climate_zscores(path):
    # Cross-fit climate records are identical within each WGSRPD3 region.
    raw={}
    with open(path,newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            code=r["wgsrpd3_code"].strip()
            if not code:continue
            try:
                bio1=float(r["bio1"]);bio7=float(r["bio7"])
                bio12=float(r["bio12"]);bio15=float(r["bio15"])
            except (KeyError,TypeError,ValueError):continue
            if bio12>=65535 or bio12<0:continue
            vals=(bio1,bio7,math.log1p(bio12),bio15)
            if code in raw and any(abs(x-y)>1e-8 for x,y in zip(raw[code],vals)):
                raise RuntimeError("inconsistent duplicate climate: "+code)
            raw[code]=vals
    if len(raw)<340:raise RuntimeError(f"insufficient climate regions: {len(raw)}")
    mu=[sum(x[j] for x in raw.values())/len(raw) for j in range(4)]
    sd=[math.sqrt(sum((x[j]-mu[j])**2 for x in raw.values())/(len(raw)-1)) for j in range(4)]
    if any(x<=0 for x in sd):raise RuntimeError("zero climate spread")
    return {code:tuple((v[j]-mu[j])/sd[j] for j in range(4)) for code,v in raw.items()}

def climate_distance(a,b):
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv','level3_geojson']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--permutations',type=int,default=499)
    ap.add_argument('--seed',type=int,default=20261007)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--climate-crossfit-csv',type=Path,required=True)
    args=ap.parse_args()

    desc=load_desc(args.descriptors_csv); pairs=load_pairs(args.insect_host_csv)
    nat=load_units(args.native_distribution_csv); con=load_units(args.contemporary_distribution_csv)
    l1map=load_level1(args.level3_geojson)
    spp=[s for s,d in desc.items() if d['b']>0 and d['n']>0]
    assert len(spp)==239
    nb={}; cb={}; ab={}
    for s in spp:
        n=set(); c=set()
        for h in pairs[s]:
            n |= nat.get(h,set()); c |= con.get(h,set())
        assert len(n)==desc[s]['n']
        nb[s]=n; cb[s]=c; ab[s]=c-n
    assert sum(len(x) for x in ab.values())==14553

    allr=sorted(set().union(*nb.values(),*cb.values()))
    active=sorted(u for u in allr if any(u in nb[s] for s in spp))
    ridx={u:i for i,u in enumerate(active)}
    missing=[u for u in active if u not in l1map]
    if missing: raise RuntimeError(f'missing level1 mapping for {missing[:10]} (n={len(missing)})')
    rl1=[l1map[u] for u in active]

    nsets=[{ridx[u] for u in nb[s] if u in ridx} for s in spp]
    asets=[{ridx[u] for u in ab[s] if u in ridx} for s in spp]
    native_masks=masks(spp,active,nb)
    obs_masks=masks(spp,active,cb)
    obs_region=mean_jaccard_masks(obs_masks)
    obs_species=mean_jaccard_sets([nsets[i]|asets[i] for i in range(len(spp))])
    geo=load_geography(args.level3_geojson)
    geographic_groups=defaultdict(list)
    for i,a_name in enumerate(active):
        a=geo[a_name]
        for j in range(i+1,len(active)):
            b=geo[active[j]]
            dist=haversine(a,b); latdiff=abs(abs(a["lat"])-abs(b["lat"]))
            cross=a["level1"]!=b["level1"]
            if cross and dist>=3000:
                geographic_groups["far_cross_all"].append((i,j))
                if latdiff<=10:
                    geographic_groups["far_cross_lat_analog"].append((i,j))
                    if a["lat"]*b["lat"]>0: geographic_groups["far_cross_analog_same_hemisphere"].append((i,j))
                if latdiff>=20: geographic_groups["far_cross_lat_discordant"].append((i,j))
            elif not cross and dist<=1000:
                geographic_groups["near_within_all"].append((i,j))
    # Climate group definitions use environmental covariates ONLY.
    climate_z=load_climate_zscores(args.climate_crossfit_csv)
    far_climate=[]
    for i,j in geographic_groups["far_cross_all"]:
        a=climate_z.get(active[i]);b=climate_z.get(active[j])
        if a is None or b is None:continue
        far_climate.append((climate_distance(a,b),i,j))
    if len(far_climate)<1000:raise RuntimeError("too few climate-comparable pairs")
    threshold_low=q([d for d,_,_ in far_climate],.25)
    threshold_high=q([d for d,_,_ in far_climate],.75)
    geographic_groups["far_cross_climate_analog"]=[(i,j) for d,i,j in far_climate if d<=threshold_low]
    geographic_groups["far_cross_climate_discordant"]=[(i,j) for d,i,j in far_climate if d>=threshold_high]
    group_names=["far_cross_all","far_cross_lat_analog","far_cross_lat_discordant",
                 "far_cross_analog_same_hemisphere","near_within_all",
                 "far_cross_climate_analog","far_cross_climate_discordant"]
    for g in group_names:
        if len(geographic_groups[g])<100: raise RuntimeError("insufficient fixed group "+g)
    geographic_groups={g:geographic_groups[g] for g in group_names}
    native_geo=grouped_jaccard(native_masks,geographic_groups)
    observed_geo=grouped_jaccard(obs_masks,geographic_groups)
    observed_gains={g:observed_geo[g]-native_geo[g] for g in group_names}
    observed_primary=observed_gains["far_cross_lat_analog"]-observed_gains["far_cross_lat_discordant"]
    observed_climate_primary=observed_gains["far_cross_climate_analog"]-observed_gains["far_cross_climate_discordant"]

    rows=[set(x) for x in asets]
    edges=[(i,c) for i,cs in enumerate(rows) for c in cs]
    by_l1=defaultdict(list)
    for k,(i,c) in enumerate(edges): by_l1[rl1[c]].append(k)
    eligible_groups=[g for g,ids in by_l1.items() if len(ids)>=2]
    if not eligible_groups: raise RuntimeError('no level1 groups eligible for swaps')

    added_masks=[0]*len(active)
    for i,cs in enumerate(rows):
        bit=1<<i
        for c in cs: added_masks[c]|=bit

    rng=random.Random(args.seed)
    def swaps(attempts):
        acc=0
        for _ in range(attempts):
            g=rng.choice(eligible_groups); ids=by_l1[g]
            k1=ids[rng.randrange(len(ids))]; k2=ids[rng.randrange(len(ids))]
            if k1==k2: continue
            i,a=edges[k1]; j,b=edges[k2]
            if i==j or a==b: continue
            if b in rows[i] or a in rows[j]: continue
            if b in nsets[i] or a in nsets[j]: continue
            rows[i].remove(a); rows[j].remove(b); rows[i].add(b); rows[j].add(a)
            edges[k1]=(i,b); edges[k2]=(j,a)
            bi=1<<i; bj=1<<j
            added_masks[a]^=bi; added_masks[a]^=bj; added_masks[b]^=bj; added_masks[b]^=bi
            acc+=1
        return acc

    burn=max(10000,15*len(edges)); accepted=swaps(burn); spacing=max(2000,3*len(edges))
    reg=[]; spe=[]; geo_null=defaultdict(list); primary_null=[]; climate_primary_null=[]
    for _ in range(args.permutations):
        accepted+=swaps(spacing)
        pseudo_masks=[native_masks[j]|added_masks[j] for j in range(len(active))]
        reg.append(mean_jaccard_masks(pseudo_masks))
        spe.append(mean_jaccard_sets([nsets[i]|rows[i] for i in range(len(spp))]))
        pseudo_geo=grouped_jaccard(pseudo_masks,geographic_groups)
        pseudo_gains={g:pseudo_geo[g]-native_geo[g] for g in group_names}
        for g in group_names:geo_null[g].append(pseudo_gains[g])
        primary_null.append(pseudo_gains["far_cross_lat_analog"]-pseudo_gains["far_cross_lat_discordant"])
        climate_primary_null.append(pseudo_gains["far_cross_climate_analog"]-pseudo_gains["far_cross_climate_discordant"])

    # Audit preservation of row x level1 counts.
    def row_l1_counts(rowsets):
        out=[]
        for cs in rowsets:
            d=defaultdict(int)
            for c in cs:d[rl1[c]]+=1
            out.append(dict(d))
        return out
    original_counts=row_l1_counts(asets); final_counts=row_l1_counts(rows)
    if original_counts!=final_counts: raise RuntimeError('species x level1 margins not preserved')

    geo_result={
        "schema":"chocho_selective_geographic_bridging_v0.1",
        "protocol":"docs/exploratory/butterfly_selective_geographic_bridging_protocol_v0.1.json",
        "status":"EXPLORATORY_NOT_FOR_MAIN_GEB",
        "primary":{
            "observed":observed_primary,
            "null_median":median(primary_null),
            "null_ci95":[q(primary_null,.025),q(primary_null,.975)],
            "observed_minus_null_median":observed_primary-median(primary_null),
            "p_one_sided":(1+sum(v>=observed_primary-1e-14 for v in primary_null))/(len(primary_null)+1)
        },
        "groups":{g:{
            "pair_count":len(geographic_groups[g]),
            "native_mean_jaccard":native_geo[g],
            "contemporary_mean_jaccard":observed_geo[g],
            "observed_gain":observed_gains[g],
            "null_median_gain":median(geo_null[g]),
            "excess":observed_gains[g]-median(geo_null[g]),
            "p_one_sided":(1+sum(v>=observed_gains[g]-1e-14 for v in geo_null[g]))/(len(geo_null[g])+1)
        } for g in group_names},
        "geographic_definition":{
            "cross_level1":"different continent codes, not formal biogeographic realms",
            "far_km_min":3000,"near_km_max":1000,
            "latitude_analogue_absolute_latitude_difference_max_degrees":10,
            "latitude_discordant_absolute_latitude_difference_min_degrees":20,
            "boundary":"Absolute latitude is not a climate measurement; region representative points approximate distance."
        },
        "claim_boundary":"This is conditional potential larval-resource geography, not realized use or colonization; positive claims require real climate and plant-only comparators."
    }
    geo_result["climate_validation"]={
        "protocol":"docs/exploratory/butterfly_climate_analog_bridging_protocol_v0.2.json",
        "status":"CLIMATE_VALIDATION_FROZEN_BEFORE_OUTCOME",
        "observed":observed_climate_primary,
        "null_median":median(climate_primary_null),
        "null_ci95":[q(climate_primary_null,.025),q(climate_primary_null,.975)],
        "observed_minus_null_median":observed_climate_primary-median(climate_primary_null),
        "p_one_sided":(1+sum(v>=observed_climate_primary-1e-14 for v in climate_primary_null))/(len(climate_primary_null)+1),
        "environmental_source":"TTF independent crossfit BIO1 BIO7 log1p(BIO12) BIO15",
        "climate_regions":len(climate_z),
        "eligible_far_pairs":len(far_climate),
        "climate_distance_q25":threshold_low,
        "climate_distance_q75":threshold_high,
        "interpretation_boundary":"Coarse four-variable regional climate summary; needs a plant-flora comparison before a consumer-specific homogenization claim."
    }
    payload={
      'schema':'chocho_butterfly_resource_homogenization_level1_null_v0.1',
      'status':'SUCCESS_LEVEL1_CONSTRAINED_FIXED_MARGIN_NULL',
      'observed':{
        'regional_mean_jaccard_contemporary':obs_region,
        'butterfly_resource_geography_mean_jaccard_contemporary':obs_species
      },
      'null':{
        'permutations':args.permutations,'seed':args.seed,'added_edges_in_active_regions':len(edges),
        'accepted_swaps_total':accepted,'eligible_level1_groups':len(eligible_groups),
        'regional_mean_jaccard_median':median(reg),'regional_mean_jaccard_ci95':[q(reg,.025),q(reg,.975)],
        'regional_observed_minus_null_median':obs_region-median(reg),
        'regional_p_greater_or_equal':(1+sum(x>=obs_region for x in reg))/(len(reg)+1),
        'butterfly_mean_jaccard_median':median(spe),'butterfly_mean_jaccard_ci95':[q(spe,.025),q(spe,.975)],
        'butterfly_observed_minus_null_median':obs_species-median(spe),
        'butterfly_p_greater_or_equal':(1+sum(x>=obs_species for x in spe))/(len(spe)+1),
        'preserves':[
          'each WGSRPD3 region column total of introduced-added butterfly resource incidences',
          'each butterfly total number of added WGSRPD3 regions',
          'each butterfly number of added WGSRPD3 regions within each WGSRPD level-1 continent',
          'native resource incidences as forbidden placements'
        ]
      },
      'selective_geographic_bridging':geo_result,
      'claim_boundary':'This post-hoc null tests whether observed resource homogenization exceeds expectations after controlling broad continental allocation as well as row and column margins. It does not establish realized competition, colonization, or fitness effects.'
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
