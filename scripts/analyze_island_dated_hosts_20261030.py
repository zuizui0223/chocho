#!/usr/bin/env python3
"""Join source-verified host catalogues to complete pre2018 island plant facets.

The unit is the original frozen species–target-island candidate. Confirmed
catalogue-compatible plant recording is distinguished from no matched record,
unresolved catalogue identity, and missing source queries. No feeding,
introduced-host, arrival, or establishment endpoint is inferred.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

import acquire_island_early_plants_20261030 as plants
import analyze_island_source_exact_null_20261030 as exact
import bridge_island_plant_taxonomy_20261030 as taxonomy

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'data/exploratory/island_dated_host_20261030'


def sha(body):
    return hashlib.sha256(body).hexdigest()


def read_csv(path):
    with path.open(newline='', encoding='utf8') as handle:
        return list(csv.DictReader(handle))


def truth(value):
    if value is True or value == 'True':
        return True
    if value is False or value == 'False':
        return False
    raise ValueError('Boolean field must be explicit True or False')


def build_catalogues(edges, bridge):
    by_id = {r['wcvp_accepted_id']: r for r in bridge}
    if len(by_id) != len(bridge):
        raise ValueError('Duplicated source plant taxonomy identity')
    if {r['wcvp_accepted_id'] for r in edges} - set(by_id):
        raise ValueError('Missing source plant taxonomy identity')
    grouped = defaultdict(set)
    for edge in edges:
        grouped[edge['original_name']].add(edge['wcvp_accepted_id'])
    catalogues = {}
    for name, ids in grouped.items():
        resolved, keys, missing = [], set(), []
        for pid in sorted(ids):
            row = by_id[pid]
            status = row['query_status']
            verification = row.get('accepted_verification_query_status', '')
            if status not in ('SUCCESS', 'FAILED', 'UNQUERIED_SOURCE_NAME'):
                raise ValueError('Unknown taxonomy query status')
            if verification not in ('', 'SUCCESS', 'FAILED'):
                raise ValueError('Unknown accepted-verification status')
            source_missing = status == 'FAILED' or verification == 'FAILED'
            if source_missing:
                missing.append(pid)
            if truth(row['conservative_join_eligible']):
                key = row['current_accepted_key']
                if status != 'SUCCESS' or source_missing or not key:
                    raise ValueError('Taxonomy eligible row has incomplete source/key')
                resolved.append(pid)
                keys.add(key)
        catalogues[name] = {
            'catalogue_host_ids': len(ids), 'resolved_source_host_ids': len(resolved),
            'resolved_host_concepts': len(keys), 'unresolved_host_ids': len(ids) - len(resolved),
            'missing_taxonomy_source_host_ids': len(missing), 'taxonomy_source_complete': not missing,
            'current_host_keys': sorted(keys), 'unresolved_source_host_ids': sorted(ids - set(resolved)),
        }
    return catalogues


def measure_hosts(catalogue, island):
    result = {k: v for k, v in catalogue.items() if k not in ('current_host_keys', 'unresolved_source_host_ids')}
    complete = island['status'] == 'COMPLETE'
    result['plant_query_complete'] = complete
    if not complete:
        return {**result, 'recorded_host_concepts': None, 'indexed_host_records': None,
                'confirmed_host_record': None, 'recorded_current_host_keys': '',
                'host_record_status': 'PLANT_QUERY_MISSING'}
    counts = island['facets']['SPECIES_KEY']
    if any(type(n) is not int or n <= 0 for n in counts.values()):
        raise ValueError('Complete species facet must contain positive integer counts')
    keys = sorted(set(catalogue['current_host_keys']) & set(counts))
    if keys:
        status = 'RECORDED_KNOWN_HOST'
    elif catalogue['unresolved_host_ids']:
        status = 'NO_MATCH_UNRESOLVED_HOST_CATALOGUE'
    else:
        status = 'NO_MATCH_COMPLETE_RESOLVED_CATALOGUE'
    return {**result, 'recorded_host_concepts': len(keys),
            'indexed_host_records': sum(counts[k] for k in keys),
            'confirmed_host_record': bool(keys), 'recorded_current_host_keys': '|'.join(keys),
            'host_record_status': status}


def recording_summary(rows):
    complete = [r for r in rows if r['plant_query_complete']]
    positives = [r for r in rows if int(r['late_positive']) == 1]
    positive_complete = [r for r in positives if r['plant_query_complete']]
    confirmed = sum(r['confirmed_host_record'] is True for r in complete)
    positive_confirmed = sum(r['confirmed_host_record'] is True for r in positive_complete)
    return {'all_candidate_cells': len(rows), 'complete_plant_query_cells': len(complete),
            'plant_query_missing_cells': len(rows) - len(complete),
            'confirmed_host_cells': confirmed,
            'confirmed_host_fraction_among_complete': confirmed / len(complete) if complete else None,
            'all_later_positive_cells': len(positives),
            'complete_plant_query_later_positive_cells': len(positive_complete),
            'later_positive_plant_query_missing_cells': len(positives) - len(positive_complete),
            'later_positive_confirmed_host_cells': positive_confirmed,
            'later_positive_confirmed_host_fraction_among_complete': positive_confirmed / len(positive_complete) if positive_complete else None,
            'record_status_counts': dict(Counter(r.get('host_record_status', '') for r in rows)),
            'descriptive_only_no_independent_host_effect_estimated': True}


def load_sources(root, output, plants_manifest_path):
    _, edges, source_pins = taxonomy.load_sources(root)
    targets = plants.load_targets(root)
    taxonomy_manifest = json.loads((output / 'plant_taxonomy_manifest.json').read_bytes())
    if taxonomy_manifest['status'] != 'FINISHED':
        raise ValueError('Taxonomy source acquisition is not finished')
    for name, expected in taxonomy_manifest['output_sha256'].items():
        if sha((output / name).read_bytes()) != expected:
            raise ValueError('Taxonomy output hash differs: ' + name)
    bridge = read_csv(output / 'plant_taxonomy_bridge.csv')
    if len(bridge) != 1062 or {r['wcvp_accepted_id'] for r in bridge} != {r['wcvp_accepted_id'] for r in edges}:
        raise ValueError('Not exactly1062 early source plant identities')
    manifest = json.loads(plants_manifest_path.read_bytes())
    if manifest['status'] != 'FINISHED':
        raise ValueError('Early plant source acquisition is not finished')
    scope = manifest['query_scope']
    if scope != {'kingdom': 'Plantae', 'taxonKey': 'P', 'checklistKey': plants.pilot.CHECKLIST,
                 'year_min': 2010, 'year_max': 2017, 'individual_record_limit': 0}:
        raise ValueError('Plant source scope changed')
    results = manifest['results']
    by_island = {r['island_key']: r for r in results}
    if len(results) != 28 or set(by_island) != {r['island_key'] for r in targets}:
        raise ValueError('Fixed28 complete target status roster required')
    for target in targets:
        result = by_island[target['island_key']]
        if result['geometry_sha256'] != target['geometry_sha256']:
            raise ValueError('Joined botanical geometry differs from butterfly geometry')
        if result['status'] == 'COMPLETE':
            facets = result['facets']
            if any(not 2010 <= int(y) <= 2017 for y in facets['YEAR']):
                raise ValueError('Plant recording dates exceed early window')
            if any(sum(facets[f].values()) != result['count'] for f in ('YEAR', 'DATASET_KEY')):
                raise ValueError('Complete botanical source margins differ')
        elif result.get('count') is not None:
            raise ValueError('Missing botanical query cannot carry complete count')
    _, risks = exact.load_inputs(root)
    return edges, bridge, by_island, risks, source_pins, taxonomy_manifest, manifest


def run(root, output, plants_manifest_path):
    edges, bridge, islands, risks, pins, tax_manifest, plant_manifest = load_sources(root, output, plants_manifest_path)
    catalogues = build_catalogues(edges, bridge)
    if len(catalogues) != 93 or sum(c['catalogue_host_ids'] for c in catalogues.values()) != 1440:
        raise ValueError('Fixed93 butterfly host catalogue totals differ')
    rows = []
    for cell in risks['primary']['cells']:
        name, island = cell['original_name'], cell['target_island']
        support = plants.plant_support(islands[island])
        rows.append({'original_name': name, 'current_butterfly_key': cell['species_key'],
                     'target_island': island, 'longitude_block': cell['longitude_block'],
                     'late_positive': int(cell['late_positive']),
                     'target_late_focal_records': int(cell['target_late_focal_records']),
                     'nearest_source_upper_km': float(cell['nearest_source_upper_km']),
                     'all_plant_indexed_records': support['record_count'],
                     'plant_botanical_support_screen': support['botanical_support_screen'],
                     **measure_hosts(catalogues[name], islands[island])})
    rows.sort(key=lambda r: (r['original_name'], r['target_island']))
    if len(rows) != 2385 or sum(r['late_positive'] for r in rows) != 66:
        raise ValueError('Original primary candidate/positive counts changed')
    risk_path = output / 'dated_host_primary_risk_cells.csv'
    plants.pilot.write_csv(risk_path, rows, list(rows[0]))
    coverage = []
    for name, cat in sorted(catalogues.items()):
        coverage.append({'original_name': name, **{k: ('|'.join(v) if isinstance(v, list) else v) for k,v in cat.items()}})
    coverage_path = output / 'dated_host_butterfly_catalogue_coverage.csv'
    plants.pilot.write_csv(coverage_path, coverage, list(coverage[0]))
    by_target = []
    for island in sorted(islands):
        subset = [r for r in rows if r['target_island'] == island]
        summary = recording_summary(subset)
        by_target.append({'target_island': island, **plants.plant_support(islands[island]),
                          **{k: v for k,v in summary.items() if k != 'record_status_counts'}})
    target_path = output / 'dated_host_target_support.csv'
    plants.pilot.write_csv(target_path, by_target, list(by_target[0]))
    report = {'schema': 'chocho_original93_fixed28_early_host_recording_20261030_v1',
              'status': 'COMPLETE_DERIVED_RECORDING_JOIN_WITH_EXPLICIT_MISSING_STATES',
              'source_pins': pins,
              'inputs': {p.name: sha(p.read_bytes()) for p in [output / 'plant_taxonomy_manifest.json',
                         output / 'plant_taxonomy_bridge.csv', plants_manifest_path]},
              'plan_sha256': sha((root / plants.PLAN_PATH).read_bytes()),
              'analysis_script_sha256': sha(Path(__file__).read_bytes()),
              'plant_acquisition_summary': plant_manifest['summary'],
              'taxonomy_acquisition_summary': tax_manifest['summary'],
              'primary_summary': recording_summary(rows),
              'catalogue_summary': {'butterflies': len(catalogues), 'source_links': len(edges),
                                    'all_source_names_resolved_butterflies': sum(c['unresolved_host_ids'] == 0 for c in catalogues.values()),
                                    'butterflies_with_at_least_one_resolved_host': sum(c['resolved_host_concepts'] > 0 for c in catalogues.values()),
                                    'taxonomy_source_missing_butterflies': sum(not c['taxonomy_source_complete'] for c in catalogues.values())},
              'limitations': [
                  'Species and year facets are separate margins; records are known within2010–2017, not species-specific dates or simultaneous butterfly–plant observations.',
                  'Known host means source-catalogue compatibility plus a current exact accepted-concept plant record; no observed feeding or larval development.',
                  'No matched indexed host record is not biological absence, including when all source names resolve.',
                  'Unresolved source names, uncertain geography, incomplete recording and cultivation remain possible.',
                  'Introduced/native plant status and historical port/shipping exposure were not established.',
                  'Descriptive co-recording fractions are not an independent host effect, especially with geographically structured ranges and recording.'],
              'output_sha256': {p.name: sha(p.read_bytes()) for p in [risk_path, coverage_path, target_path]}}
    plants.write_json(output / 'dated_host_recording_summary.json', report)
    print(json.dumps({'primary_summary': report['primary_summary'], 'catalogue_summary': report['catalogue_summary']}, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    parser.add_argument('--plants-manifest', type=Path)
    args = parser.parse_args()
    run(args.root, args.output_dir, args.plants_manifest or args.output_dir / 'early_plants_resolved_manifest.json')
