#!/usr/bin/env python3
"""Post-hoc climate-conditioned host-introduction null; fixed submission untouched.

Reuse the tested plant-family/continent edge-switch engine. Add within-continent
quantile bins of region-level raw CHELSA BIO1 and BIO12. Missing-climate regions
are each frozen separately, not pooled or dropped. Region climate is deduplicated
across species; no butterfly occurrence, fitted niche score, or response enters
stratification. This is coarse climate conditioning, not climate equivalence.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median, mean


def load_climate(path: Path):
    climate, counts, invalid = {}, {}, set()
    nrows = 0
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'wgsrpd3_code', 'bio1', 'bio12', 'climate_sample_points'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('climate CSV schema drift')
        for row in reader:
            nrows += 1
            code = row['wgsrpd3_code'].strip()
            if not code:
                raise ValueError('empty region code')
            try:
                n = int(row['climate_sample_points'])
                vector = (float(row['bio1']), float(row['bio12']))
            except (ValueError, TypeError):
                invalid.add(code)
                continue
            if n < 1 or not all(math.isfinite(v) for v in vector):
                invalid.add(code)
                continue
            if code in climate and any(abs(a-b)>1e-8 for a,b in zip(climate[code], vector)):
                raise ValueError('conflicting raw regional climate: '+code)
            if code in counts and counts[code] != n:
                raise ValueError('conflicting climate point count: '+code)
            climate[code], counts[code] = vector, n
    if invalid & set(climate):
        raise ValueError('conflicting valid/missing regional climate')
    if not climate:
        raise ValueError('no finite regional climate')
    return climate, {
        'csv_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_rows': nrows, 'unique_valid_regions': len(climate),
        'invalid_regions_in_csv': sorted(invalid),
        'min_valid_points': min(counts.values()),
        'median_valid_points': median(counts.values()),
        'raw_variables': ['bio1','bio12'],
    }


def quantile(values, probability):
    vals=sorted(values)
    z=(len(vals)-1)*probability
    lo=int(z); hi=min(lo+1,len(vals)-1)
    return vals[lo]+(vals[hi]-vals[lo])*(z-lo)


def climate_strata(continents, climate, bins):
    if bins not in (0,2,3):
        raise ValueError('bins must be 0, 2 or 3')
    by_continent=defaultdict(list)
    for region in sorted(continents):
        if region in climate:
            by_continent[str(continents[region])].append(region)
    labels, cutpoints = {}, {}
    for continent, regions in sorted(by_continent.items()):
        cuts=[]
        for axis in range(2):
            values=[climate[r][axis] for r in regions]
            cuts.append([quantile(values,j/bins) for j in range(1,bins)] if bins else [])
        cutpoints[continent]={'regions':len(regions),'bio1':cuts[0],'bio12':cuts[1]}
        for r in regions:
            t=bisect.bisect_right(cuts[0],climate[r][0])
            p=bisect.bisect_right(cuts[1],climate[r][1])
            labels[r]=f'{continent}|T{t}|P{p}'
    missing=sorted(set(continents)-set(climate))
    for r in missing:
        labels[r]=f'{continents[r]}|missing:{r}'
    return labels, {
        'bins_per_climate_axis':bins,
        'cutpoints_within_continent_raw_raster_units':cutpoints,
        'missing_climate_regions':missing,
        'missing_handling':'Each missing-climate region is a unique fixed block, in every model.',
        'block_region_counts':dict(sorted(Counter(labels.values()).items())),
    }


def trace_diagnostics(values):
    x=[float(v) for v in values]
    m=mean(x); den=sum((v-m)**2 for v in x)
    rho=sum((a-m)*(b-m) for a,b in zip(x[:-1],x[1:]))/den if den else None
    # Descriptive first-order ESS approximation, not a convergence certificate.
    ess=min(len(x),max(1.,len(x)*(1-rho)/(1+rho))) if rho is not None else None
    return {'lag1_autocorrelation':rho,'ar1_effective_samples_approx':ess,
            'first_half_mean':mean(x[:len(x)//2]),'second_half_mean':mean(x[len(x)//2:])}


def compare(observed, values):
    n=len(values)
    return {'observed':observed, 'null_mean':mean(values), 'null_median':median(values),
            'null_q025':quantile(values,.025),'null_q975':quantile(values,.975),
            'difference_from_null_mean':observed-mean(values),
            'p_upper':(1+sum(v>=observed for v in values))/(n+1),
            'p_lower':(1+sum(v<=observed for v in values))/(n+1),
            'trace_diagnostics':trace_diagnostics(values)}


def execute_chain(engine, species_hosts, family, native, contemporary, panel, selected,
                  strata, observed, climate, seed, replicas, burn, gap):
    introduced, blocks=engine.init_edges(native,contemporary,family,strata)
    movable=engine.build_swappable_blocks(blocks)
    signature=engine.degree_signature(blocks)
    original={edge for edges in blocks.values() for edge in edges}
    focal=set().union(*(species_hosts[sp] for sp in panel))
    focal_original={e for e in original if e[0] in focal}
    rng=random.Random(seed)
    attempts,accepted=engine.swaps(movable,native,introduced,rng,burn)
    collected={tag:[] for tag in selected}
    changed=[]
    for j in range(replicas):
        at,ac=engine.swaps(movable,native,introduced,rng,gap)
        attempts+=at; accepted+=ac
        summary,_=engine.geographic_resources(panel,species_hosts,native,introduced,selected)
        for tag,row in summary.items():
            total=row['introduced_host_region_incidences']
            if total != observed[tag]['introduced_host_region_incidences']:
                raise RuntimeError('fixed incidence total drift')
            if total != (row['added_units']+row['introduced_incidences_already_covered_by_native_hosts']+
                         row['duplicated_introduced_incidences_outside_native_resource']):
                raise RuntimeError('overlap decomposition drift')
            collected[tag].append(row)
        if j in (0,replicas//2,replicas-1):
            now={e for edges in blocks.values() for e in edges}
            changed.append({'replicate':j+1,'changed_all_edges':len(original-now),
                            'changed_focal_host_edges':len(focal_original-now)})
    if signature!=engine.degree_signature(blocks):
        raise RuntimeError('fixed climate/family margin drift')
    if any(introduced[h] & native.get(h,set()) for h in introduced):
        raise RuntimeError('native range exclusion drift')
    keys=['added_units','single_host_share','mean_species_single_share',
          'introduced_incidences_already_covered_by_native_hosts',
          'duplicated_introduced_incidences_outside_native_resource']
    comparisons={tag:{k:compare(observed[tag][k],[r[k] for r in rows]) for k in keys}
                 for tag,rows in collected.items()}
    return {'seed':seed,'replicas':replicas,'burn_sweeps':burn,'gap_sweeps':gap,
            'comparisons':comparisons,
            'mixing':{'introduced_edges':len(original),'focal_host_edges':len(focal_original),
                      'edges_in_eligible_blocks':sum(len(x) for x in movable.values()),
                      'eligible_blocks':len(movable), 'attempted_swaps':attempts,
                      'accepted_swaps':accepted,'acceptance_rate':accepted/attempts if attempts else 0.,
                      'changed_edges':changed,'invariants_verified':True},
            'trace':collected}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--climate-root',type=Path,required=True)
    ap.add_argument('--insect-host-csv',type=Path,required=True)
    ap.add_argument('--native-csv',type=Path,required=True)
    ap.add_argument('--contemporary-csv',type=Path,required=True)
    ap.add_argument('--expansion-csv',type=Path,default=Path('data/frozen/figure_sources/anthropogenic_species_metrics.csv'))
    ap.add_argument('--lookup-json',type=Path,default=Path('data/frozen/wgsrpd3_level1_lookup_v0.1.json'))
    ap.add_argument('--previous-json',type=Path,default=Path('provenance/reviewer_defenses/results/butterfly_resource_redundancy_geonull_decomposition_20261011.json'))
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--replicas',type=int,default=199)
    args=ap.parse_args()
    if args.replicas<19:
        raise ValueError('at least 19 replicas needed')
    import analyze_butterfly_resource_redundancy_geonull as engine
    climate_files=list(args.climate_root.rglob('unit_climate_crossfit.csv'))
    if len(climate_files)!=1:
        raise ValueError(f'exactly one frozen regional climate CSV required, found {len(climate_files)}')
    climate,clim_source=load_climate(climate_files[0])
    continents=json.loads(args.lookup_json.read_text())['units']
    hosts,family=engine.read_hosts(args.insect_host_csv)
    native=engine.read_units(args.native_csv)
    contemporary=engine.read_units(args.contemporary_csv)
    with args.expansion_csv.open(newline='') as f:
        frozen={r['species']:r for r in csv.DictReader(f)}
    panel=sorted(frozen)
    adequate=[s for s in panel if int(frozen[s]['introduced_added_units'])>0 and
              frozen[s]['host_taxonomy_lower_bound_adequate']=='True']
    rich=[s for s in adequate if int(frozen[s]['resolved_host_species'])>=10]
    if (len(panel),len(adequate),len(rich))!=(239,191,80):
        raise RuntimeError('focal panel drift')
    selected={'all_239':panel,'adequate_191':adequate,'rich_hosts_80':rich}
    introductions,_=engine.init_edges(native,contemporary,family,continents)
    observed,per_species=engine.geographic_resources(panel,hosts,native,introductions,selected)
    previous=json.loads(args.previous_json.read_text())['baseline']
    for tag in selected:
        for key,expected in previous[tag].items():
            if abs(observed[tag][key]-expected)>1e-10:
                raise RuntimeError(f'frozen baseline drift {tag}/{key}')
    for s in panel:
        if per_species[s][0]!=int(frozen[s]['native_resource_units']):
            raise RuntimeError('native union drift '+s)
    output=args.output_dir; output.mkdir(parents=True,exist_ok=True)
    with (output/'regional_climate_unique.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['region','level1','bio1_raw','bio12_raw'])
        for r in sorted(climate):
            w.writerow([r,continents.get(r,''),*climate[r]])
    # Same missing-climate support in all three null definitions.
    configurations=[('matched_support_no_climate',0,2026101101,15,3),
                    ('climate_2x2',2,2026101102,15,3),
                    ('climate_3x3_primary',3,2026101103,15,3),
                    ('climate_3x3_longer_mixing',3,2026101104,30,6)]
    report={'status':'EXPLORATORY_POSTHOC_CLIMATE_CONDITIONED_SPATIAL_NULL',
            'climate_source':clim_source,'baseline':observed,
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [args.insect_host_csv,args.native_csv,args.contemporary_csv,args.expansion_csv,args.lookup_json,args.previous_json]},
            'models':{},'limitation':'BIO1/BIO12 regional medians and quantile bins are coarse climate controls, not exact climatic matching. No causal/demographic inference. Missing regions frozen; shared-data post-hoc tests, not independent biological replication.'}
    all_edges=[(h,r) for h,regs in introductions.items() for r in regs]
    report['climate_source']['introduced_edges_without_climate']=sum(r not in climate for h,r in all_edges)
    report['climate_source']['introduced_edges_total']=len(all_edges)
    for name,bins,seed,burn,gap in configurations:
        strata,receipt=climate_strata(continents,climate,bins)
        print('START_MODEL '+name,flush=True)
        result=execute_chain(engine,hosts,family,native,contemporary,panel,selected,strata,observed,
                             climate,seed,args.replicas,burn,gap)
        trace=result.pop('trace')
        (output/(name+'_trace.json')).write_text(json.dumps(trace,indent=2)+'\n')
        result['stratification']=receipt
        report['models'][name]=result
        (output/'climate_geonull_result.json').write_text(json.dumps(report,indent=2)+'\n')
        compact={tag:result['comparisons'][tag] for tag in selected}
        print('MODEL_RESULT '+json.dumps({'model':name,'comparisons':compact,'mixing':result['mixing']}),flush=True)
    print('CLIMATE_NULL_SUMMARY_BEGIN',flush=True)
    print(json.dumps(report,indent=2),flush=True)
    print('CLIMATE_NULL_SUMMARY_END',flush=True)


if __name__=='__main__':
    main()
