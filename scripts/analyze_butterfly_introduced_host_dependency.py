#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def load_desc(p):
    out={}
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            out[r['species']]={'b':float(r['host_family_count']),'n':int(r['host_wgsrpd3_unit_count'])}
    return out

def load_pairs(p):
    by=defaultdict(set);meta={}
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            s=r['insect_species'].strip();h=r['accepted_plant_name_id'].strip()
            if s and h:
                by[s].add(h);meta.setdefault(h,{'accepted_name':r['accepted_name'].strip(),'family':r['family'].strip()})
    return by,meta

def load_units(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            h=r['accepted_plant_name_id'].strip();u=r['area_code_l3'].strip()
            if h and u:out[h].add(u)
    return out

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-csv',type=Path,required=True)
    args=ap.parse_args()
    d=load_desc(args.descriptors_csv);pairs,meta=load_pairs(args.insect_host_csv)
    nat=load_units(args.native_distribution_csv);con=load_units(args.contemporary_distribution_csv)
    spp=[s for s,x in d.items() if x['b']>0 and x['n']>0];assert len(spp)==239
    rows=[];expanded=0
    for s in spp:
        hs=sorted(pairs[s]);n=set();c=set();added_by={}
        for h in hs:
            hn=set(nat.get(h,set()));hc=set(con.get(h,set()));n|=hn;c|=hc;added_by[h]=hc-hn
        assert len(n)==d[s]['n'];add=c-n
        if not add:continue
        expanded+=1
        unique=defaultdict(int);frac=defaultdict(float);regions_by=defaultdict(list)
        for u in add:
            cs=[h for h in hs if u in added_by[h]]
            sh=1/len(cs)
            for h in cs:frac[h]+=sh
            if len(cs)==1:
                unique[cs[0]]+=1;regions_by[cs[0]].append(u)
        for h in sorted(set(frac)|set(unique)):
            rows.append({'species':s,'host_family_count':d[s]['b'],'added_resource_regions':len(add),
                         'accepted_plant_name_id':h,'accepted_name':meta[h]['accepted_name'],'plant_family':meta[h]['family'],
                         'fractional_contribution_units':frac[h],'fractional_contribution_share':frac[h]/len(add),
                         'uniquely_supported_added_regions':unique[h],
                         'unique_dependency_fraction':unique[h]/len(add),
                         'unique_region_codes':';'.join(sorted(regions_by[h]))})
    assert expanded==206
    critical25=[r for r in rows if r['unique_dependency_fraction']>=.25]
    critical50=[r for r in rows if r['unique_dependency_fraction']>=.5]
    critical100=[r for r in rows if r['unique_dependency_fraction']>=1]
    sp25={r['species'] for r in critical25};sp50={r['species'] for r in critical50};sp100={r['species'] for r in critical100}
    ranked=sorted(rows,key=lambda r:(-r['unique_dependency_fraction'],-r['uniquely_supported_added_regions'],r['species'],r['accepted_name']))
    payload={'schema':'chocho_butterfly_introduced_host_dependency_ledger_v0.1','status':'SUCCESS_POSTHOC_DEPENDENCY_LEDGER',
      'expanded_butterflies':206,'species_with_any_host_uniquely_supporting_ge_25pct_added_regions':len(sp25),
      'species_with_any_host_uniquely_supporting_ge_50pct_added_regions':len(sp50),
      'species_with_any_host_uniquely_supporting_100pct_added_regions':len(sp100),
      'critical_species_host_links_ge_25pct':len(critical25),'critical_species_host_links_ge_50pct':len(critical50),
      'top_dependency_links':ranked[:50],
      'management_interpretation':'Use this ledger only to prioritize locality-level checks before removal or restoration decisions. A coarse-region dependency means no other recorded host in the reconstruction supplies that butterfly in that WGSRPD3 unit; it does not prove local larval use, demographic dependence, host quality, or that the focal introduced plant should be retained.',
      'claim_boundary':'HOSTS and WCVP are incomplete and WGSRPD3 is coarse. Unique support is reconstruction-level resource redundancy, not realized ecological dependence.'}
    args.output_json.parent.mkdir(parents=True,exist_ok=True);args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    with open(args.output_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(ranked)
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
