#!/usr/bin/env python3
"""Post-hoc seasonality/subcontinental refinements of the fixed-introduction null.

Runs the existing tested swap engine. No host identities, native ranges, or
observed response values change. Seasonal proxies are BIO7 and BIO15, not
phenology. WGSRPD level2 is coarse regional conditioning, not distance matching.
"""
from __future__ import annotations
import argparse
import bisect
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
import analyze_butterfly_resource_climate_geonull as previous
import analyze_butterfly_resource_redundancy_geonull as engine

VARIABLES = ('bio1', 'bio12', 'bio7', 'bio15')
CONFIGURATIONS = {
    'baseline': (False, False, 2026101121, 15, 3),
    'seasonal': (True, False, 2026101122, 15, 3),
    'geographic': (False, True, 2026101123, 15, 3),
    'combined': (True, True, 2026101124, 15, 3),
    'combined_long': (True, True, 2026101125, 30, 6),
}


def load_four_climate(path: Path):
    """Deduplicate regional vectors; keep flagged raw values for audit."""
    values, points, invalid, flagged = {}, {}, set(), set()
    nrows = 0
    with path.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        required = {'wgsrpd3_code', 'climate_sample_points', *VARIABLES}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError('climate CSV schema drift')
        for row in reader:
            nrows += 1
            region = row['wgsrpd3_code'].strip()
            if not region:
                raise ValueError('empty climate region')
            try:
                vector = tuple(float(row[k]) for k in VARIABLES)
                n = int(row['climate_sample_points'])
            except (ValueError, TypeError):
                invalid.add(region)
                continue
            if n < 1 or not all(math.isfinite(v) for v in vector):
                invalid.add(region)
                continue
            if region in values and (values[region] != vector or points[region] != n):
                raise ValueError('conflicting regional climate ' + region)
            values[region], points[region] = vector, n
            if any(v == 65535 for v in vector):
                flagged.add(region)
    if invalid & set(values):
        raise ValueError('conflicting valid/missing climate')
    if not values:
        raise ValueError('no regional climate')
    return values, {
        'source_rows': nrows, 'unique_regions': len(values),
        'raw_variables': list(VARIABLES),
        'invalid_regions_in_csv': sorted(invalid),
        'ceiling_flagged_regions': sorted(flagged),
        'min_points': min(points.values()), 'median_points': median(points.values()),
        'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'ceiling_policy': 'Possible storage ceiling, not confirmed nodata; retained but individually frozen in all models, excluded from climate cutpoint calibration.',
    }


def load_subregions(path: Path, continents: dict[str, str]) -> dict[str, str]:
    """Read official geometry properties without approximating region names."""
    result = {}
    data = json.loads(path.read_text(encoding='utf-8'))
    for feature in data.get('features', []):
        props = feature.get('properties', {})
        r = str(props.get('LEVEL3_COD', '')).strip()
        if r not in continents:
            continue
        c = str(int(props['LEVEL1_COD']))
        sub = str(int(props['LEVEL2_COD']))
        if c != str(continents[r]) or not sub.startswith(c):
            raise ValueError('continent mismatch for ' + r)
        if r in result and result[r] != sub:
            raise ValueError('conflicting geography ' + r)
        result[r] = sub
    if set(result) != set(continents):
        raise ValueError('geography coverage differs from frozen lookup')
    return result


def refine_strata(continents, subregions, climate, *, seasonal, geographic, quarantine=()):
    if geographic and set(continents) - set(subregions):
        raise ValueError('missing geography')
    quarantine = set(quarantine)
    usable = set(climate) - quarantine
    labels, cuts = {}, {}
    grouped = defaultdict(list)
    for r in sorted(continents):
        if r in usable:
            grouped[str(continents[r])].append(r)
    for continent, regions in sorted(grouped.items()):
        bounds = []
        for axis, bins in enumerate((3, 3, 2, 2)):
            vals = [climate[r][axis] for r in regions]
            bounds.append([previous.quantile(vals, i / bins) for i in range(1, bins)])
        cuts[continent] = dict(zip(VARIABLES, bounds))
        for r in regions:
            indices = [bisect.bisect_right(b, climate[r][j]) for j,b in enumerate(bounds)]
            label = f'{continent}|T{indices[0]}|P{indices[1]}'
            if seasonal:
                label += f'|A{indices[2]}|S{indices[3]}'
            if geographic:
                label += '|L2:' + str(subregions[r])
            labels[r] = label
    frozen = sorted(set(continents) - usable)
    for r in frozen:
        labels[r] = f'{continents[r]}|frozen:{r}'
    counts = Counter(labels.values())
    return labels, {
        'mean_climate_variables': ['bio1', 'bio12'],
        'mean_climate_bins': [3,3],
        'seasonal_variables': ['bio7','bio15'] if seasonal else [],
        'seasonal_bins': [2,2] if seasonal else [],
        'geography_level': 2 if geographic else 1,
        'frozen_regions': frozen,
        'cutpoints_within_continent': cuts,
        'region_block_counts': dict(sorted(counts.items())),
        'region_blocks': len(counts),
        'single_region_blocks': sum(n == 1 for n in counts.values()),
    }


def run(args):
    climate, audit = load_four_climate(args.climate_csv)
    continents = json.loads(args.lookup_json.read_text())['units']
    subregions = load_subregions(args.level3_geojson, continents)
    hosts, family = engine.read_hosts(args.insect_host_csv)
    native = engine.read_units(args.native_csv)
    contemporary = engine.read_units(args.contemporary_csv)
    with args.expansion_csv.open(newline='') as f:
        frozen = {r['species']: r for r in csv.DictReader(f)}
    panel = sorted(frozen)
    adequate = [s for s in panel if int(frozen[s]['introduced_added_units']) > 0
                and frozen[s]['host_taxonomy_lower_bound_adequate'] == 'True']
    rich = [s for s in adequate if int(frozen[s]['resolved_host_species']) >= 10]
    if (len(panel), len(adequate), len(rich)) != (239,191,80):
        raise ValueError('frozen panel drift')
    selected = {'all_239':panel, 'adequate_191':adequate, 'rich_hosts_80':rich}
    introductions,_ = engine.init_edges(native,contemporary,family,continents)
    observed, per_species = engine.geographic_resources(panel,hosts,native,introductions,selected)
    receipt = json.loads(args.previous_result.read_text())
    for tag, expected in receipt['baseline'].items():
        for key, val in expected.items():
            if abs(observed[tag][key] - val) > 1e-9:
                raise ValueError(f'baseline drift {tag}/{key}')
    quarantine = set(audit['ceiling_flagged_regions'])
    seasonal, geographic, seed, burn, gap = CONFIGURATIONS[args.model]
    labels, stratification = refine_strata(continents,subregions,climate,
        seasonal=seasonal,geographic=geographic,quarantine=quarantine)
    all_edges = {(h,r) for h,regs in introductions.items() for r in regs}
    focal_hosts = set().union(*(hosts[s] for s in panel))
    audit['quarantined_introduction_edges'] = sum(r in quarantine for h,r in all_edges)
    audit['quarantined_focal_host_introduction_edges'] = sum(r in quarantine and h in focal_hosts for h,r in all_edges)
    result = previous.execute_chain(engine,hosts,family,native,contemporary,panel,selected,
        labels,observed,climate,seed,args.replicas,burn,gap)
    trace = result.pop('trace')
    result.update({
        'schema': 'chocho_seasonal_subcontinental_geonull_20261011',
        'status': 'EXPLORATORY_POSTHOC_CONDITIONAL_NULL',
        'model': args.model, 'baseline': observed, 'climate_audit': audit,
        'stratification': stratification,
        'source_sha256': {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (args.climate_csv,args.lookup_json,args.level3_geojson,args.insect_host_csv,
             args.native_csv,args.contemporary_csv,args.expansion_csv)},
        'limitation': 'Coarse seasonal-climate and level2 geographic constraints, not phenological synchrony, exact distance, host performance, causality, or population persistence. More constraints reduce the randomization space. Report movement diagnostics, not just p values.',
    })
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/(args.model+'_trace.json')).write_text(json.dumps(trace,indent=2)+'\n')
    (args.output_dir/(args.model+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')
    key='introduced_incidences_already_covered_by_native_hosts'
    print(json.dumps({'model':args.model, 'primary':result['comparisons']['adequate_191'][key],
                     'mixing':result['mixing']},indent=2),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',choices=CONFIGURATIONS,required=True)
    ap.add_argument('--climate-csv',type=Path,required=True)
    ap.add_argument('--level3-geojson',type=Path,required=True)
    ap.add_argument('--insect-host-csv',type=Path,required=True)
    ap.add_argument('--native-csv',type=Path,required=True)
    ap.add_argument('--contemporary-csv',type=Path,required=True)
    ap.add_argument('--expansion-csv',type=Path,default=Path('data/frozen/figure_sources/anthropogenic_species_metrics.csv'))
    ap.add_argument('--lookup-json',type=Path,default=Path('data/frozen/wgsrpd3_level1_lookup_v0.1.json'))
    ap.add_argument('--previous-result',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--replicas',type=int,default=199)
    args=ap.parse_args()
    if args.replicas < 19:
        raise ValueError('at least 19 sampled matrices required')
    run(args)

if __name__ == '__main__':
    main()
