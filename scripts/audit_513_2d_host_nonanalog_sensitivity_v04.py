#!/usr/bin/env python3
"""Source-verified 513-insect 349-host-combination extrapolation sensitivity.

Compares the unchanged, source-approved native plant 0BP/21ka climate landscape
with source 0BP logistic-function-extrapolated novel climates. No occupied
historical plants, dispersal barriers or genetic responses are inferred.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from audit_513_2d_resource_topology_v03 import (
    accepted,compute_masks,group_metrics,summarize
)

SHA={
"original_mask":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
"extrapolated_mask":"67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a",
"original513":"838b73df7f9ba7cd23cc33b981d3a434931ab1a6ed61f4bd27621f009d129020",
"original_insect_host":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044",
"source_logit_recovery_summary":"3ac3d4cd6020e22377315359a2d50ddd461cbcd892ecb49b8f34f0a90a7f32d6"
}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as file:
        for chunk in iter(lambda:file.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def read_rows(path):
    with Path(path).open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))

def count_sign_transitions(a,b,field):
    if set(a)!=set(b):raise ValueError('source species mismatch for comparison')
    count=Counter()
    for species in a:
        p=float(a[species][field]);q=float(b[species][field])
        if not np.isfinite(p) or not np.isfinite(q):raise ValueError('nonfinite original/extended score')
        before='decline' if p<0 else 'nondecline'
        after='decline' if q<0 else 'nondecline'
        count[before+'_to_'+after]+=1
    return dict(count)

def source_comparison(primary,alternate,original513):
    # Inputs equal-species source-verified, each has 513 records.
    if not (len(primary)==len(alternate)==len(original513)==513):
        raise ValueError('not 513 prequalified insect species')
    orig={r['insect_species']:r for r in original513}
    if set(orig)!={r['insect_species'] for r in primary}:
        raise ValueError('original frozen insect source identities not preserved')
    for row in primary:
        earlier=orig[row['insect_species']]
        for key in ('area_delta_km2','component_delta_4','component_delta_8',
                    'largest_fraction_delta_4','largest_fraction_delta_8'):
            if abs(float(row[key])-float(earlier[key]))>1e-7:
                raise ValueError('v0.3 original 4- and 8-neighbor source results not replayed')
    return {'primary_source_replayed':True,
            'area_direction_transitions':count_sign_transitions(orig,
               {r['insect_species']:r for r in alternate},'area_delta_km2'),
            'largest_share4_direction_transitions':count_sign_transitions(orig,
               {r['insect_species']:r for r in alternate},'largest_fraction_delta_4'),
            'largest_share8_direction_transitions':count_sign_transitions(orig,
               {r['insect_species']:r for r in alternate},'largest_fraction_delta_8')}

def run(plant_main,plant_extended,prior513,links,model_summary):
    for key,path in (('original_mask',plant_main),('extrapolated_mask',plant_extended),
                     ('original513',prior513),('original_insect_host',links),
                     ('source_logit_recovery_summary',model_summary)):
        if sha(path)!=SHA[key]:raise ValueError('source SHA mismatch '+key)
    original=read_rows(prior513)
    if len(original)!=513 or len({r['insect_species'] for r in original})!=513:
        raise ValueError('source original 513-species model panel drift')
    q=json.loads(Path(model_summary).read_text())
    if (q['modeled_plant_taxa']!=416 or
        q['source_holdout_max_absolute_error_21ka']>1e-5 or
        not q['original_416_main_masks_inherited_byte_for_byte'] or
        not q['source_21ka_predictions_used_only_for_validation']):
        raise ValueError('416 source model functions failed original 21ka verification')
    bysp=defaultdict(set)
    for r in read_rows(links):
        bysp[r['insect_species']].add(accepted(r['wcvp_id_as_frozen']))
    if len(bysp)!=641:raise ValueError('original all-641 host identities changed')
    original_species=set(r['insect_species'] for r in original)
    groups=defaultdict(list)
    for sp in original_species:
        groups[tuple(sorted(bysp[sp],key=int))].append(sp)
    if len(groups)!=349:raise ValueError('original 349 recorded host combinations changed')
    main_rows=[];extra_rows=[];per_species=[]
    with np.load(plant_main,allow_pickle=False) as first, np.load(plant_extended,allow_pickle=False) as second:
        ids=[accepted(x) for x in first['plant_wcvp_ids']]
        extra_ids=[accepted(x) for x in second['plant_wcvp_ids']]
        if (ids!=extra_ids or len(ids)!=416 or
            first['mask_0BP'].shape!=first['mask_21ka'].shape!=(416,32400) or
            second['mask_0BP'].shape!=second['mask_21ka'].shape!=(416,32400) or
            not bool(first['nonanalog_excluded']) or bool(second['nonanalog_excluded'])):
            raise ValueError('plant-model identity/count/original analog exclusion drift')
        for era in ('0BP','21ka'):
            a=np.asarray(first['mask_'+era]);b=np.asarray(second['mask_'+era])
            if np.any((a & (~b))!=0):
                raise ValueError('extended source has eliminated originally supported plant habitat')
        idx={plant:i for i,plant in enumerate(ids)}
        if len(idx)!=416:raise ValueError('duplicate plant model IDs')
        lat=np.asarray(first['lat_centers'],float)
        for j,(hostset,sp_list) in enumerate(sorted(groups.items())):
            if not set(hostset)<=set(idx):
                raise ValueError('qualified insect has an unmodeled known host')
            order=[idx[x] for x in hostset]
            original_grids=compute_masks(order,first['mask_0BP'],first['mask_21ka'])
            extrap_grids=compute_masks(order,second['mask_0BP'],second['mask_21ka'])
            if any(np.any(a&~b) for a,b in zip(original_grids,extrap_grids)):
                raise ValueError('extrapolated insect climate habitat is not a superset')
            main=group_metrics(*original_grids,lat)
            ext=group_metrics(*extrap_grids,lat)
            for sp in sp_list:
                base={'insect_species':sp,
                      'hostset':'/'.join(hostset),
                      'hostset_shared_insects':len(sp_list),
                      'recorded_host_n':len(hostset)}
                main_rows.append({**base,**main})
                extra_rows.append({**base,**ext})
                per_species.append({**base,
                    **{'main_'+k:v for k,v in main.items()},
                    **{'extended_'+k:v for k,v in ext.items()}})
            if (j+1)%75==0:
                print(json.dumps({'unique_host_combinations_processed':j+1,'total':349}),flush=True)
    main_rows.sort(key=lambda r:r['insect_species'])
    extra_rows.sort(key=lambda r:r['insect_species'])
    per_species.sort(key=lambda r:r['insect_species'])
    paired=source_comparison(main_rows,extra_rows,original)
    representative=lambda rows:list({r['hostset']:r for r in rows}.values())
    out={'schema':'chocho_513insect_349hostset_frozen0BP_native_logit_nonanalog_v04',
        'status':'EXPLORATORY_CLIMATE_EXTRAPOLATION_SENSITIVITY_NOT_ANCIENT_OCCUPANCY',
        'species_qualified':513,'original_species':641,'unique_known_hostsets':349,
        'source_validation':paired,
        'main_species_equal_weight':summarize(main_rows),
        'extended_species_equal_weight':summarize(extra_rows),
        'main_distinct_hostset_equal_weight':summarize(representative(main_rows)),
        'extended_distinct_hostset_equal_weight':summarize(representative(extra_rows)),
        'source_sha256':dict(SHA),
        'limits':['Source 0BP logistic functions reconstructed from saved 0BP probabilities, heldout 21ka scores verified against originals',
            'Including environmentally nonanalog climate predictions measures extrapolation sensitivity, NOT recovered 21ka plant occupancy',
            'Four and eight neighbors refer to hypothetical coarse icefree climate cells, not dispersal routes or observed genetic breaks',
            'Insect host records not documented as complete, and many of the 513 insects reuse one of 349 identical host combination climate landscapes',
            'The older negative TTF v0.5 recorded-host versus matched-other-plant null stays unchanged']}
    return per_species,out

def main():
    p=argparse.ArgumentParser()
    for key in ('plant-main','plant-extended','prior513','insect-host','model-summary','out'):
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows,q=run(a.plant_main,a.plant_extended,a.prior513,a.insect_host,a.model_summary)
    with (a.out/'all_513_main_vs_extrapolated_climatic_host_graphs_v04.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    (a.out/'summary.json').write_text(json.dumps(q,indent=2,sort_keys=True)+'\n')
    print(json.dumps(q,sort_keys=True))
if __name__=='__main__':main()
