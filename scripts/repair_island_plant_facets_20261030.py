#!/usr/bin/env python3
"""Repair tied-count facet pagination without changing scientific predicates."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import acquire_island_early_plants_20261030 as original

ROOT = Path(__file__).resolve().parents[1]
PLAN = 'provenance/archive/host_hindcast/ISLAND_PLANT_FACET_REPAIR_20261030.md'
FIELDS = original.pilot.FIELDS
LIMIT = 10000


def maps(page, max_bins=LIMIT):
    original.validate_count_response(original.pilot.canonical(page))
    facets = page.get('facets')
    if not isinstance(facets, list):
        raise ValueError('Missing facets')
    result = {f: {} for f in FIELDS}
    seen = set()
    for facet in facets:
        field = facet['field']
        if field not in result or field in seen:
            raise ValueError('Duplicate or unknown facet field')
        seen.add(field)
        counts = facet['counts']
        if not isinstance(counts, list) or len(counts) >= max_bins:
            raise ValueError('Facet reached bounded single-page capacity')
        for entry in counts:
            key, value = entry['name'], entry['count']
            if not isinstance(key, str) or not key or type(value) is not int or value <= 0:
                raise ValueError('Invalid positive facet bin')
            if key in result[field]:
                raise ValueError('Duplicate facet key within single page')
            result[field][key] = value
    return result


def collect_single_page(fetch):
    first = fetch(None, 0, LIMIT)
    facets = maps(first)
    count = first['count']
    for field in FIELDS:
        terminal = fetch(field, len(facets[field]), 1)
        original.validate_count_response(original.pilot.canonical(terminal))
        if terminal['count'] != count:
            raise ValueError('Total count changed during terminal probes')
        terminal_fields = terminal.get('facets', [])
        if not isinstance(terminal_fields, list):
            raise ValueError('Malformed terminal facets')
        if any(f.get('field') != field or not isinstance(f.get('counts'), list) or f['counts']
               for f in terminal_fields):
            raise ValueError('Nonempty or malformed terminal probe; possible silent source cap')
        if len(terminal_fields) > 1:
            raise ValueError('Duplicate terminal facet fields')
    again = fetch(None, 0, LIMIT, recheck=True)
    if again['count'] != count or maps(again) != facets:
        raise ValueError('Full source facet key-to-count maps changed on recheck')
    residual = {f: count - sum(values.values()) for f, values in facets.items()}
    if any(value < 0 for value in residual.values()) or residual['YEAR'] or residual['DATASET_KEY']:
        raise ValueError('Complete source facet margins do not reconcile')
    if any(not key.isdigit() or not 2010 <= int(key) <= 2017 for key in facets['YEAR']):
        raise ValueError('Plant year outside fixed early window')
    return {'count': count, 'facets': facets, 'unassigned_facet_records': residual,
            'facet_pages': 1, 'first_page_recheck_passed': True,
            'terminal_probes_passed': True, 'recheck_compares_key_count_maps_not_tie_order': True}


def merge_results(initial, repairs):
    by_key = {r['island_key']: r for r in initial}
    if len(by_key) != len(initial) or len({r['island_key'] for r in repairs}) != len(repairs):
        raise ValueError('Duplicate island identities')
    merged = copy.deepcopy(by_key)
    for result in repairs:
        key = result['island_key']
        if key not in by_key or by_key[key]['status'] == 'COMPLETE':
            raise ValueError('Cannot replace original complete result or add target')
        merged[key] = copy.deepcopy(result)
    return [merged[key] for key in sorted(merged)]


def acquire_repair(row, wkt, cache, client):
    folder = cache / row['island_key']
    folder.mkdir(parents=True, exist_ok=True)
    receipts = []
    def fetch(field, offset, limit, recheck=False):
        payload = original.plant_payload(wkt, 0)
        payload.update(facetLimit=limit, facetOffset=offset)
        if field:
            payload['facets'] = [field]
        prefix = ('terminal_' + field.lower()) if field else ('full_recheck' if recheck else 'full')
        request = original.pilot.canonical(payload)
        if (folder / (prefix + '_receipt.json')).exists():
            page, receipt = original.read_count_cache(folder, prefix, request)
        else:
            page, receipt = original.save_count_response(folder, prefix, request, client.post(request))
        receipts.append({**receipt, 'facet': field, 'facet_offset': offset, 'facet_limit': limit,
                         'is_first_page_recheck': recheck,
                         'cache_prefix': row['island_key'] + '/' + prefix})
        return page
    started = original.pilot.utc()
    result = original.pilot.query_safely(lambda: collect_single_page(fetch))
    return {'island_key': row['island_key'], 'geometry_sha256': row['geometry_sha256'],
            'started_utc': started, 'finished_utc': original.pilot.utc(),
            'source_pages': receipts, **result}


def run(output, cache, source_cache, root=ROOT):
    initial_path = output / 'early_plants_manifest.json'
    initial_bytes = initial_path.read_bytes()
    initial = json.loads(initial_bytes)
    if initial['status'] != 'FINISHED' or len(initial['results']) != 28:
        raise ValueError('Initial fixed28 acquisition must finish first')
    initial_archive = output / initial['source_archive']['file']
    if original.sha(initial_archive.read_bytes()) != initial['source_archive']['sha256']:
        raise ValueError('Initial source archive changed')
    failed = {r['island_key'] for r in initial['results'] if r['status'] != 'COMPLETE'}
    if len(failed) != 3:
        raise ValueError('Repair expected exactly the original3failed queries')
    targets = [r for r in original.load_targets(root) if r['island_key'] in failed]
    lock_path = output / 'early_plants_repair_run_lock.json'
    spec = {'initial_manifest_sha256': original.sha(initial_bytes),
            'repair_plan_sha256': original.sha((root / PLAN).read_bytes()),
            'repair_script_sha256': original.sha(Path(__file__).read_bytes()),
            'initial_script_sha256': original.sha(Path(original.__file__).read_bytes()),
            'eligible_failed_islands': sorted(failed), 'full_facet_limit': LIMIT,
            'terminal_probe_each_facet': True, 'scientific_predicate_unchanged': True}
    if lock_path.exists():
        raise ValueError('Repair first-run lock exists; do not rerun after inspecting results')
    lock = {'locked_before_repair_queries_utc': original.pilot.utc(), 'spec': spec}
    original.write_json(lock_path, lock)
    cache.mkdir(parents=True, exist_ok=True)
    wkts = original.prepare_geometries(source_cache, targets)
    client, repairs = original.SourceClient(), []
    for target in targets:
        result = acquire_repair(target, wkts[target['island_key']], cache, client)
        repairs.append(result)
        print(json.dumps({'stage': 'repair', 'island_key': target['island_key'],
                          'status': result['status'], 'count': result.get('count'),
                          'error': result.get('error')}), flush=True)
    archive = original.archive_count_sources(cache, output / 'early_plants_repair_raw_facets.tar.gz')
    repair_manifest = {'schema': 'chocho_early_plants_tied_facet_repair_20261030_v1',
                       'status': 'FINISHED', 'lock': lock, 'results': repairs,
                       'source_archive': archive, 'http_attempts': client.attempts,
                       'http_retries': client.retries, 'finished_utc': original.pilot.utc()}
    repair_path = output / 'early_plants_repair_manifest.json'
    original.write_json(repair_path, repair_manifest)
    merged = merge_results(initial['results'], repairs)
    complete = [r for r in merged if r['status'] == 'COMPLETE']
    support = [{'island_key': r['island_key'], 'status': r['status'], **original.plant_support(r)} for r in merged]
    summary = {'requested_islands': 28, 'completed_islands': len(complete),
               'missing_islands': 28 - len(complete), 'pending_islands': 0,
               'complete_zero_islands': sum(r['count'] == 0 for r in complete),
               'indexed_plant_records': sum(r['count'] for r in complete),
               'botanical_support_islands': sum(r['botanical_support_screen'] for r in support)}
    merged_manifest = {'schema': 'chocho_early_plants_original_plus_repaired_source_20261030_v1',
                       'status': 'FINISHED', 'query_scope': initial['query_scope'],
                       'summary': summary, 'results': merged,
                       'source_manifests': {initial_path.name: original.sha(initial_bytes),
                                            repair_path.name: original.sha(repair_path.read_bytes())},
                       'source_archives': [initial['source_archive'], archive],
                       'source_route_by_island': {r['island_key']: ('REPAIR_SINGLE_PAGE_TERMINAL_PROBES' if r['island_key'] in failed else 'ORIGINAL_COMPLETE_PAGINATION') for r in merged},
                       'initial_complete_result_rows_preserved': all(next(v for v in merged if v['island_key'] == r['island_key']) == r for r in initial['results'] if r['status'] == 'COMPLETE'),
                       'interpretation': initial['interpretation'], 'assembled_utc': original.pilot.utc()}
    support_path = output / 'early_plants_resolved_island_support.csv'
    original.pilot.write_csv(support_path, support, list(support[0]))
    facet_paths = []
    for field in FIELDS:
        path = output / ('early_plants_resolved_' + field.lower() + '_facets.csv')
        rows = [{'island_key': r['island_key'], 'facet_key': k, 'record_count': n}
                for r in complete for k,n in sorted(r['facets'][field].items())]
        original.pilot.write_csv(path, rows, ['island_key', 'facet_key', 'record_count'])
        facet_paths.append(path)
    merged_manifest['output_sha256'] = {p.name: original.sha(p.read_bytes()) for p in [support_path, *facet_paths]}
    original.write_json(output / 'early_plants_resolved_manifest.json', merged_manifest)
    if original.sha(initial_path.read_bytes()) != spec['initial_manifest_sha256']:
        raise ValueError('Initial manifest changed during repair')
    print(json.dumps({'stage': 'merged_finished', **summary}), flush=True)
    return merged_manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=original.OUTPUT)
    parser.add_argument('--cache-dir', type=Path, default=Path('/tmp/chocho-early-plants-repair-cache'))
    parser.add_argument('--source-cache', type=Path, default=Path('/tmp/chocho-island-sources'))
    args = parser.parse_args()
    run(args.output_dir, args.cache_dir, args.source_cache)
