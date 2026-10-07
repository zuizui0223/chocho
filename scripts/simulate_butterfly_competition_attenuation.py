#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
from statistics import median

ALPHAS=(0.0,0.05,0.1,0.25,0.5,1.0,2.0)

def load_descriptors(path):
    out={}
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            out[row['species']]={'host_family_count':float(row['host_family_count']),
                                 'native_resource_units':int(row['host_wgsrpd3_unit_count'])}
    return out

def load_pairs(path):
    by=defaultdict(set)
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            s=row['insect_species'].strip(); h=row['accepted_plant_name_id'].strip()
            if s and h: by[s].add(h)
    return by

def load_units(path):
    out=defaultdict(set)
    with open(path,newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        for row in r:
            h=row['accepted_plant_name_id'].strip(); u=row['area_code_l3'].strip()
            if h and u: out[h].add(u)
    return out

def summarize(values):
    xs=sorted(values)
    n=len(xs)
    return {'n':n,'mean':sum(xs)/n if n else None,
            'median':median(xs) if n else None,
            'q25':xs[int(.25*(n-1))] if n else None,
            'q75':xs[int(.75*(n-1))] if n else None}

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-species-csv',type=Path,required=True)
    args=ap.parse_args()

    desc=load_descriptors(args.descriptors_csv)
    pairs=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv)
    cont=load_units(args.contemporary_distribution_csv)
    spp=[s for s,d in desc.items() if d['host_family_count']>0 and d['native_resource_units']>0]
    assert len(spp)==239
    focal=set(spp)

    consumers=defaultdict(set)
    for s in spp:
        for h in pairs[s]: consumers[h].add(s)

    # For each species and region, build competitor exposure under native and contemporary
    native_regions={}; cont_regions={}; native_comp={}; cont_comp={}
    for s in spp:
        nr=set(); cr=set()
        ncomp=defaultdict(set); ccomp=defaultdict(set)
        for h in pairs[s]:
            others=consumers[h]-{s}
            for u in native.get(h,set()):
                nr.add(u); ncomp[u].update(others)
            for u in cont.get(h,set()):
                cr.add(u); ccomp[u].update(others)
        assert len(nr)==desc[s]['native_resource_units']
        native_regions[s]=nr; cont_regions[s]=cr
        native_comp[s]={u:len(ncomp[u]) for u in nr}
        cont_comp[s]={u:len(ccomp[u]) for u in cr}

    rows=[]
    grid={}
    for alpha in ALPHAS:
        gains=[]; ratios=[]; positive=zero=negative=0
        for s in spp:
            nscore=sum(1/(1+alpha*native_comp[s][u]) for u in native_regions[s])
            cscore=sum(1/(1+alpha*cont_comp[s][u]) for u in cont_regions[s])
            gain=cscore-nscore
            ratio=cscore/nscore if nscore else None
            gains.append(gain); ratios.append(ratio)
            if gain>1e-12: positive+=1
            elif gain<-1e-12: negative+=1
            else: zero+=1
            rows.append({'species':s,'host_family_count':desc[s]['host_family_count'],
                         'alpha':alpha,'native_effective_opportunity':nscore,
                         'contemporary_effective_opportunity':cscore,
                         'effective_gain':gain,'effective_ratio':ratio,
                         'native_regions':len(native_regions[s]),
                         'contemporary_regions':len(cont_regions[s]),
                         'mean_native_competitor_exposure':sum(native_comp[s].values())/len(native_regions[s]),
                         'mean_contemporary_competitor_exposure':sum(cont_comp[s].values())/len(cont_regions[s])})
        grid[str(alpha)]={'species_positive_gain':positive,'species_zero_gain':zero,'species_negative_gain':negative,
                          'median_effective_ratio':median(ratios),'median_effective_gain':median(gains),
                          'effective_ratio_summary':summarize(ratios)}

    # new regions are more or less crowded?
    native_load=[]; added_load=[]
    for s in spp:
        native_load.extend(native_comp[s][u] for u in native_regions[s])
        added=cont_regions[s]-native_regions[s]
        added_load.extend(cont_comp[s][u] for u in added)

    # Identify robust winners and crowding-sensitive species
    by_species=defaultdict(dict)
    for row in rows: by_species[row['species']][row['alpha']]=row
    sensitive=[]
    for s in spp:
        base=by_species[s][0.0]
        first_flip=None
        for a in ALPHAS[1:]:
            if by_species[s][a]['effective_gain']<=0:
                first_flip=a; break
        sensitive.append({'species':s,'host_family_count':desc[s]['host_family_count'],
                          'raw_added_regions':base['contemporary_regions']-base['native_regions'],
                          'mean_native_competitor_exposure':base['mean_native_competitor_exposure'],
                          'mean_contemporary_competitor_exposure':base['mean_contemporary_competitor_exposure'],
                          'first_alpha_with_nonpositive_effective_gain':first_flip,
                          'ratio_alpha_0.25':by_species[s][0.25]['effective_ratio'],
                          'ratio_alpha_1':by_species[s][1.0]['effective_ratio'],
                          'ratio_alpha_2':by_species[s][2.0]['effective_ratio']})
    sensitive.sort(key=lambda x:(x['first_alpha_with_nonpositive_effective_gain'] is None,
                                 x['first_alpha_with_nonpositive_effective_gain'] if x['first_alpha_with_nonpositive_effective_gain'] is not None else 99,
                                 x['ratio_alpha_1']))

    payload={'schema':'chocho_butterfly_competition_attenuation_scenarios_v0.1',
             'status':'SUCCESS_SCENARIO_ANALYSIS_NOT_REALIZED_COMPETITION',
             'model':{'effective_region_value':'1 / (1 + alpha * C_sr)',
                      'C_sr':'number of other focal butterflies documented to use at least one exact same host species present in region r',
                      'alphas':ALPHAS,
                      'interpretation':'alpha is an unknown scenario parameter, not an estimated competition coefficient.'},
             'panel':{'butterflies':239},
             'competitor_exposure':{'native_resource_units':summarize(native_load),
                                    'introduced_added_resource_units':summarize(added_load),
                                    'mean_added_minus_native':sum(added_load)/len(added_load)-sum(native_load)/len(native_load)},
             'competition_attenuation_grid':grid,
             'species_flipping_to_nonpositive_gain_by_alpha':{
                 str(a):sum(x['first_alpha_with_nonpositive_effective_gain'] is not None and x['first_alpha_with_nonpositive_effective_gain']<=a for x in sensitive)
                 for a in ALPHAS[1:]},
             'most_crowding_sensitive_species':sensitive[:30],
             'claim_boundary':'This simulation asks whether the geographic opportunity created by host redistribution remains positive after imposing hypothetical interspecific crowding penalties. It does not estimate realized competition, abundance, carrying capacity, demographic growth, or colonization probability. Intraspecific competition is not modelled because abundance is unavailable.'}
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    with open(args.output_species_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
