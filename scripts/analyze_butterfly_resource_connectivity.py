#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,random
from collections import defaultdict,deque
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
            s=r['insect_species'].strip();h=r['accepted_plant_name_id'].strip()
            if s and h:out[s].add(h)
    return out

def load_units(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            h=r['accepted_plant_name_id'].strip();u=r['area_code_l3'].strip()
            if h and u:out[h].add(u)
    return out

def load_graph(p):
    d=json.load(open(p))
    geoms={};level1={}
    for f in d['features']:
        pr=f.get('properties') or {}; code=str(pr.get('LEVEL3_COD') or '').strip()
        if not code or not f.get('geometry'):continue
        geoms[code]=shape(f['geometry'])
        level1[code]=str(pr.get('LEVEL1_COD') or code[:2]).strip()
    codes=sorted(geoms); adj={c:set() for c in codes}
    for i,a in enumerate(codes):
        ga=geoms[a]
        for b in codes[i+1:]:
            gb=geoms[b]
            if ga.touches(gb):
                adj[a].add(b);adj[b].add(a)
    return codes,level1,adj

def comps(nodes,adj):
    nodes=set(nodes);seen=set();out=[]
    for s in nodes:
        if s in seen:continue
        q=[s];seen.add(s);cc=set()
        while q:
            x=q.pop();cc.add(x)
            for y in adj.get(x,()):
                if y in nodes and y not in seen:seen.add(y);q.append(y)
        out.append(cc)
    return out

def bridge_metrics(native,contemporary,adj):
    nc=comps(native,adj); cc=comps(contemporary,adj)
    if not native:return {'native_components':0,'contemporary_components':len(cc),'bridged_native_pairs':0,'possible_cross_component_native_pairs':0,'bridged_fraction':None,'native_components_merged':0}
    nlab={}
    for i,c in enumerate(nc):
        for u in c:nlab[u]=i
    clab={}
    for i,c in enumerate(cc):
        for u in c:clab[u]=i
    sizes=[len(c) for c in nc]
    possible=0
    for i in range(len(sizes)):
        for j in range(i+1,len(sizes)):possible+=sizes[i]*sizes[j]
    bridged=0
    comp_ids=range(len(nc))
    for i in comp_ids:
        for j in range(i+1,len(nc)):
            # all pairs between native components i,j are connected if any nodes share same contemporary component
            # components cannot partially merge in an undirected graph: if one pair connects, the full components merge.
            ci=next(iter(nc[i])); cj=next(iter(nc[j]))
            if clab.get(ci)==clab.get(cj):bridged+=len(nc[i])*len(nc[j])
    # how many native components participate in a merged contemporary component
    merged_groups=defaultdict(set)
    for u in native: merged_groups[clab[u]].add(nlab[u])
    merged=sum(max(0,len(v)-1) for v in merged_groups.values())
    return {'native_components':len(nc),'contemporary_components':len(cc),'bridged_native_pairs':bridged,
            'possible_cross_component_native_pairs':possible,'bridged_fraction':bridged/possible if possible else 0.0,
            'native_components_merged':merged,'largest_native_component':max(map(len,nc)) if nc else 0,
            'largest_contemporary_component':max(map(len,cc)) if cc else 0}

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv','level3_geojson']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--permutations',type=int,default=499)
    ap.add_argument('--seed',type=int,default=20261007)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-species-csv',type=Path,required=True)
    args=ap.parse_args()
    d=load_desc(args.descriptors_csv);pairs=load_pairs(args.insect_host_csv)
    nat=load_units(args.native_distribution_csv);con=load_units(args.contemporary_distribution_csv)
    codes,l1,adj=load_graph(args.level3_geojson);code_set=set(codes)
    spp=[s for s,x in d.items() if x['b']>0 and x['n']>0];assert len(spp)==239
    rows=[];native_sets={};added_sets={}
    for s in spp:
        n=set();c=set()
        for h in pairs[s]:n|=nat.get(h,set());c|=con.get(h,set())
        n&=code_set;c&=code_set
        assert len(n)==d[s]['n']
        native_sets[s]=n;added_sets[s]=c-n
        m=bridge_metrics(n,c,adj)
        rows.append({'species':s,'host_family_count':d[s]['b'],'native_regions':len(n),'added_regions':len(c-n),**m})
    eligible=[r for r in rows if r['possible_cross_component_native_pairs']>0 and r['added_regions']>0]
    obs_total=sum(r['bridged_native_pairs'] for r in eligible)
    obs_possible=sum(r['possible_cross_component_native_pairs'] for r in eligible)
    obs_frac=obs_total/obs_possible if obs_possible else 0
    species_any=sum(r['bridged_native_pairs']>0 for r in eligible)
    rng=random.Random(args.seed)
    null_tot=[];null_species=[]
    pools=defaultdict(list)
    for u in codes:pools[l1[u]].append(u)
    for rep in range(args.permutations):
        total=0;ns=0
        for r in eligible:
            s=r['species']; n=native_sets[s]; a=added_sets[s]
            counts=defaultdict(int)
            for u in a:counts[l1[u]]+=1
            sampled=set()
            for lev,k in counts.items():
                avail=[u for u in pools[lev] if u not in n]
                if len(avail)<k: raise RuntimeError(f'insufficient null pool {s} {lev}')
                sampled.update(rng.sample(avail,k))
            m=bridge_metrics(n,n|sampled,adj)
            total+=m['bridged_native_pairs'];ns+=m['bridged_native_pairs']>0
        null_tot.append(total);null_species.append(ns)
    null_tot.sort();null_species.sort()
    def q(xs,p):return xs[int(p*(len(xs)-1))]
    payload={'schema':'chocho_butterfly_resource_connectivity_v0.1','status':'SUCCESS_POSTHOC_RESOURCE_CONNECTIVITY',
      'definition':'WGSRPD3 polygons sharing a boundary are adjacent. A native resource-region pair is bridged when its two originally disconnected native components become connected after introduced host-resource regions are added.',
      'panel':{'butterflies':239,'species_with_multiple_native_components_and_added_resources':len(eligible)},
      'observed':{'species_with_any_native_component_bridging':species_any,
                  'total_bridged_native_region_pairs':obs_total,'total_possible_cross_component_native_region_pairs':obs_possible,
                  'bridged_fraction':obs_frac,
                  'median_species_bridged_fraction':median([r['bridged_fraction'] for r in eligible]),
                  'species_with_fewer_components_contemporary':sum(r['contemporary_components']<r['native_components'] for r in rows),
                  'species_with_larger_largest_component':sum(r['largest_contemporary_component']>r['largest_native_component'] for r in rows)},
      'level1_constrained_null':{'permutations':args.permutations,'seed':args.seed,
             'preserves':'For each butterfly, the exact number of introduced-added WGSRPD3 regions within each WGSRPD level-1 region; native resource regions remain fixed.',
             'bridged_pairs_median':median(null_tot),'bridged_pairs_ci95':[q(null_tot,.025),q(null_tot,.975)],
             'p_null_ge_observed':(1+sum(x>=obs_total for x in null_tot))/(args.permutations+1),
             'species_with_any_bridge_median':median(null_species),'species_with_any_bridge_ci95':[q(null_species,.025),q(null_species,.975)]},
      'claim_boundary':'Resource connectivity is a coarse geographic opportunity metric, not dispersal or colonization itself. Boundary adjacency is conservative for over-water dispersal and WGSRPD3 is too coarse for local corridor management.'}
    args.output_json.parent.mkdir(parents=True,exist_ok=True);args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    with open(args.output_species_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
