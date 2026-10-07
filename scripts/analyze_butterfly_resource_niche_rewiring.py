#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict,Counter
from pathlib import Path
from statistics import median

def load_descriptors(path):
    out={}
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            out[row['species']]={'host_family_count':float(row['host_family_count']),'native_resource_units':int(row['host_wgsrpd3_unit_count'])}
    return out

def load_taxonomy(path):
    out={}
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            s=row['Species'].strip(); fam=row['Family'].strip()
            if s and fam: out[s]=fam
    return out

def load_pairs(path):
    by=defaultdict(set); meta={}
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            s=row['insect_species'].strip(); h=row['accepted_plant_name_id'].strip()
            if s and h:
                by[s].add(h)
                meta.setdefault(h,{'accepted_name':row['accepted_name'].strip(),'family':row['family'].strip()})
    return by,meta

def load_units(path):
    out=defaultdict(set)
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            h=row['accepted_plant_name_id'].strip(); u=row['area_code_l3'].strip()
            if h and u: out[h].add(u)
    return out

def q(xs,p):
    xs=sorted(xs); pos=p*(len(xs)-1); lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    if lo==hi:return xs[lo]
    a=pos-lo; return xs[lo]*(1-a)+xs[hi]*a

def pairwise_summary(native_cells,cont_cells,taxonomy):
    spp=list(native_cells)
    dn=[]; dc=[]; delta=[]; same=[]; cross=[]; inc=dec=eq=0
    pairs_with_new_shared=0; cross_new=0
    for i,a in enumerate(spp):
        A=native_cells[a]; C=cont_cells[a]
        for b in spp[i+1:]:
            B=native_cells[b]; D=cont_cells[b]
            un=A|B; uc=C|D
            jn=len(A&B)/len(un) if un else 0.0
            jc=len(C&D)/len(uc) if uc else 0.0
            d=jc-jn
            dn.append(jn);dc.append(jc);delta.append(d)
            (same if taxonomy[a]==taxonomy[b] else cross).append(d)
            if d>1e-15: inc+=1
            elif d<-1e-15: dec+=1
            else:eq+=1
            if (C&D)-(A&B):
                pairs_with_new_shared+=1
                if taxonomy[a]!=taxonomy[b]:cross_new+=1
    return {
        'butterfly_pairs':len(delta),
        'mean_jaccard_native':sum(dn)/len(dn),'mean_jaccard_contemporary':sum(dc)/len(dc),
        'delta_mean_jaccard':sum(dc)/len(dc)-sum(dn)/len(dn),
        'median_jaccard_native':median(dn),'median_jaccard_contemporary':median(dc),
        'pairs_increased':inc,'pairs_decreased':dec,'pairs_unchanged':eq,
        'fraction_pairs_increased':inc/len(delta),
        'median_delta_same_family':median(same),'median_delta_cross_family':median(cross),
        'mean_delta_same_family':sum(same)/len(same),'mean_delta_cross_family':sum(cross)/len(cross),
        'pairs_gaining_new_shared_host_region_cells':pairs_with_new_shared,
        'cross_family_pairs_gaining_new_shared_cells':cross_new,
        'fraction_new_shared_pairs_cross_family':cross_new/pairs_with_new_shared if pairs_with_new_shared else None,
    }

def random_removal(units,hosts,ks,reps,seed):
    rng=random.Random(seed); out={}; units=list(units)
    for k in ks:
        losses=[]
        for _ in range(reps):
            rem=set(rng.sample(hosts,k))
            losses.append(sum(contrib<=rem for _,_,contrib in units))
        out[str(k)]={'reps':reps,'median_lost_units':median(losses),'ci95_lost_units':[q(losses,.025),q(losses,.975)]}
    return out

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','leptraits_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--random-reps',type=int,default=499)
    ap.add_argument('--seed',type=int,default=20261007)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-host-csv',type=Path,required=True)
    args=ap.parse_args()
    desc=load_descriptors(args.descriptors_csv); tax=load_taxonomy(args.leptraits_csv); pairs,meta=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv); cont=load_units(args.contemporary_distribution_csv)
    spp=[s for s,d in desc.items() if d['host_family_count']>0 and d['native_resource_units']>0]
    assert len(spp)==239
    for s in spp: assert s in tax
    host_idx={h:i for i,h in enumerate(sorted(meta))}
    all_regions=sorted(set().union(*native.values(),*cont.values())); reg_idx={u:i for i,u in enumerate(all_regions)}
    stride=len(all_regions)
    ncells={}; ccells={}; unit_contrib=[]; host_credit=defaultdict(float); host_consumers=defaultdict(set)
    native_paircells=0; intro_paircells=0; intro_cross_paircells=0; pair_gain=set(); cross_pair_gain=set()
    eligible=set(spp)
    consumers_by_host=defaultdict(set)
    for s in spp:
        for h in pairs[s]:
            consumers_by_host[h].add(s)
    for h,consset in consumers_by_host.items():
        consumers=sorted(consset)
        n=len(consumers)
        if n<2: continue
        choose=n*(n-1)//2
        native_paircells += choose*len(native.get(h,set()))
        intr=set(cont.get(h,set()))-set(native.get(h,set()))
        intro_paircells += choose*len(intr)
        fam_counts=Counter(tax[s] for s in consumers)
        same_pairs=sum(v*(v-1)//2 for v in fam_counts.values())
        cross_pairs=choose-same_pairs
        intro_cross_paircells += cross_pairs*len(intr)
        if intr:
            for i,a in enumerate(consumers):
                for b in consumers[i+1:]:
                    key=(a,b); pair_gain.add(key)
                    if tax[a]!=tax[b]:cross_pair_gain.add(key)
    for s in spp:
        hs=sorted(pairs[s]); n_union=set(); c_union=set(); added_by={}; ncell=set(); ccell=set()
        for h in hs:
            hn=set(native.get(h,set())); hc=set(cont.get(h,set()))
            n_union |= hn; c_union |= hc; added_by[h]=hc-hn
            hi=host_idx[h]
            ncell.update(hi*stride+reg_idx[u] for u in hn)
            ccell.update(hi*stride+reg_idx[u] for u in hc)
            host_consumers[h].add(s)
        assert len(n_union)==desc[s]['native_resource_units']
        ncells[s]=ncell;ccells[s]=ccell
        for u in c_union-n_union:
            cs=frozenset(h for h in hs if u in added_by[h]); assert cs
            unit_contrib.append((s,u,cs))
            share=1/len(cs)
            for h in cs:host_credit[h]+=share
    assert len(unit_contrib)==14553
    ps=pairwise_summary(ncells,ccells,tax)

    # Fixed-demand null allocation model.
    # Each butterfly has total larval demand = 1 and allocates it uniformly
    # across all exact host×region resource cells available to it. Under this
    # deliberately simple scenario, sum(p^2) is within-species concentration,
    # while sum(p_i * p_j) over shared cells is between-species co-allocation.
    self_native=[]
    self_contemporary=[]
    inter_native=[]
    inter_contemporary=[]
    inter_delta=[]
    for sp in spp:
        self_native.append(1.0/len(ncells[sp]))
        self_contemporary.append(1.0/len(ccells[sp]))
    for i,a in enumerate(spp):
        A=ncells[a]; C=ccells[a]
        for b in spp[i+1:]:
            B=ncells[b]; D=ccells[b]
            on=len(A&B)/(len(A)*len(B))
            oc=len(C&D)/(len(C)*len(D))
            inter_native.append(on)
            inter_contemporary.append(oc)
            inter_delta.append(oc-on)

    competition_null={
        'assumption':'Each butterfly has fixed total larval demand = 1 and distributes it uniformly across all available exact host×WGSRPD3 resource cells. All resource cells are equal quality and no density response is fitted.',
        'within_species_resource_concentration':{
            'mean_native':sum(self_native)/len(self_native),
            'mean_contemporary':sum(self_contemporary)/len(self_contemporary),
            'median_native':median(self_native),
            'median_contemporary':median(self_contemporary),
            'relative_change_mean':(sum(self_contemporary)/sum(self_native))-1.0,
        },
        'between_species_exact_resource_coallocation':{
            'mean_pair_overlap_native':sum(inter_native)/len(inter_native),
            'mean_pair_overlap_contemporary':sum(inter_contemporary)/len(inter_contemporary),
            'sum_pair_overlap_native':sum(inter_native),
            'sum_pair_overlap_contemporary':sum(inter_contemporary),
            'relative_change_sum':(sum(inter_contemporary)/sum(inter_native)-1.0) if sum(inter_native)>0 else None,
            'pairs_increased':sum(d>1e-15 for d in inter_delta),
            'pairs_decreased':sum(d<-1e-15 for d in inter_delta),
            'pairs_unchanged':sum(abs(d)<=1e-15 for d in inter_delta),
        },
        'inter_to_intra_pressure_ratio':{
            'native':sum(inter_native)/sum(self_native),
            'contemporary':sum(inter_contemporary)/sum(self_contemporary),
            'relative_change':(sum(inter_contemporary)/sum(self_contemporary))/(sum(inter_native)/sum(self_native))-1.0,
        },
        'boundary':'This is a null allocation scenario, not a demographic competition model. It assumes equal resource quality, fixed abundance and uniform use, and omits phenology, host biomass, preference, natural enemies and local density dependence.'
    }

    ps['native_shared_host_region_pair_cells']=native_paircells
    ps['introduced_added_shared_host_region_pair_cells']=intro_paircells
    ps['relative_increase_shared_resource_cells']=intro_paircells/native_paircells if native_paircells else None
    ps['introduced_added_cross_family_pair_cells']=intro_cross_paircells
    ps['fraction_introduced_pair_cells_cross_family']=intro_cross_paircells/intro_paircells if intro_paircells else None
    ps['distinct_pairs_with_any_introduced_shared_host']=len(pair_gain)
    ps['distinct_cross_family_pairs_with_any_introduced_shared_host']=len(cross_pair_gain)
    redundancy=Counter(len(c) for _,_,c in unit_contrib)
    single=redundancy.get(1,0)
    contributing_hosts=sorted(host_credit,key=lambda h:(-host_credit[h],h))
    unique_loss=Counter()
    for s,u,c in unit_contrib:
        if len(c)==1:unique_loss[next(iter(c))]+=1
    ks=[1,5,10,20,38,50,100]
    targeted={}
    for k in ks:
        rem=set(contributing_hosts[:k]); lost=[(s,u,c) for s,u,c in unit_contrib if c<=rem]
        lost_by_sp=Counter(s for s,_,_ in lost); total_by_sp=Counter(s for s,_,_ in unit_contrib)
        complete=sum(total_by_sp[s]>0 and lost_by_sp[s]==total_by_sp[s] for s in total_by_sp)
        targeted[str(k)]={'lost_units':len(lost),'fraction_added_opportunity_lost':len(lost)/len(unit_contrib),
                          'butterflies_losing_any_added_region':len(lost_by_sp),'butterflies_losing_all_added_regions':complete}
    rnd=random_removal(unit_contrib,contributing_hosts,ks,args.random_reps,args.seed)
    host_rows=[]
    for rank,h in enumerate(contributing_hosts,1):
        host_rows.append({'rank':rank,'accepted_plant_name_id':h,'accepted_name':meta[h]['accepted_name'],'plant_family':meta[h]['family'],
                          'fractional_credit':host_credit[h],'unique_support_units':unique_loss[h],
                          'consumer_species':len(host_consumers[h]),'consumer_families':len({tax[s] for s in host_consumers[h]})})
    species_total=Counter(s for s,_,_ in unit_contrib); species_single=Counter(s for s,_,c in unit_contrib if len(c)==1)
    vulnerable=[{'species':s,'butterfly_family':tax[s],'added_units':species_total[s],'single_host_supported_units':species_single[s],
                 'single_host_dependency_fraction':species_single[s]/species_total[s]} for s in species_total if species_total[s]>=10]
    vulnerable.sort(key=lambda x:(-x['single_host_dependency_fraction'],-x['added_units'],x['species']))
    payload={'schema':'chocho_butterfly_resource_niche_rewiring_v0.1','status':'SUCCESS_POSTHOC_RESOURCE_NICHE_REWIRING_AND_REMOVAL_STRESS_TEST',
             'panel':{'butterflies':239,'added_butterfly_x_region_units':len(unit_contrib),'contributing_introduced_hosts':len(contributing_hosts)},
             'resource_niche_overlap':ps,
             'fixed_demand_competition_null':competition_null,
             'added_opportunity_redundancy':{'contributor_count_distribution':dict(sorted(redundancy.items())),
                 'single_host_supported_units':single,'fraction_single_host_supported':single/len(unit_contrib),
                 'median_contributing_hosts_per_added_unit':median([len(c) for _,_,c in unit_contrib])},
             'targeted_global_removal_stress_test':{'removal_order':'descending fractional contribution to added butterfly×region opportunity','scenarios':targeted},
             'random_host_removal_null':rnd,
             'most_single_host_dependent_butterflies':vulnerable[:30],
             'top_hosts':host_rows[:30],
             'claim_boundary':'Resource-niche overlap and shared host×region cells quantify potential resource co-use and interaction exposure, not realized competition. The removal analysis is a global counterfactual stress test, not a site-specific management prescription. Competition intensity, host quality, abundance, phenology, natural enemies and local replacement resources are not observed here.'}
    args.output_json.parent.mkdir(parents=True,exist_ok=True); args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    with open(args.output_host_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(host_rows[0]));w.writeheader();w.writerows(host_rows)
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
