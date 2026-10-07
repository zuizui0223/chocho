#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict
from pathlib import Path
from statistics import median

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

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv','level3_geojson']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--permutations',type=int,default=499)
    ap.add_argument('--seed',type=int,default=20261007)
    ap.add_argument('--output-json',type=Path,required=True)
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
    reg=[]; spe=[]
    for _ in range(args.permutations):
        accepted+=swaps(spacing)
        reg.append(mean_jaccard_masks([native_masks[j]|added_masks[j] for j in range(len(active))]))
        spe.append(mean_jaccard_sets([nsets[i]|rows[i] for i in range(len(spp))]))

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
      'claim_boundary':'This post-hoc null tests whether observed resource homogenization exceeds expectations after controlling broad continental allocation as well as row and column margins. It does not establish realized competition, colonization, or fitness effects.'
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
