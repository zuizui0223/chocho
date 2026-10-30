#!/usr/bin/env python3
"""Complete early Plantae facets on the original 28 target islands.

Only count-only public responses are acquired. This records indexed botanical
evidence, not plant absence, wild establishment, feeding or uncertainty-disk
containment. No later plant observations or individual records are requested.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import acquire_island_butterfly_panel_20261029 as pilot

ROOT = Path(__file__).resolve().parents[1]
TARGET_PATH = 'data/exploratory/island_recorded_sources_20261029/recorded_source_target_support.csv'
SELECTION_PATH = 'data/exploratory/island_expansion_20261029/expanded_island_panel_selection.json'
PLAN_PATH = 'provenance/archive/host_hindcast/ISLAND_DATED_HOST_LATITUDE_PLAN_20261030.md'
PINS = {
    TARGET_PATH: '9bbf0d50f1a42b82879f23ef966a6de4a6dc36f5dc03cc92e5c511fc6bad3e74',
    SELECTION_PATH: '5ccc868b89973b8142a7fdb4b3014c4394f134596ae2ef58e1179eaab68f5243',
}
OUTPUT = ROOT / 'data/exploratory/island_dated_host_20261030'
WINDOW = (2010, 2017)


def sha(body):
    return hashlib.sha256(body).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + '\n')


def load_targets(root=ROOT):
    root = Path(root)
    for path, expected in PINS.items():
        if sha((root / path).read_bytes()) != expected:
            raise ValueError('Frozen source hash differs: ' + path)
    with (root / TARGET_PATH).open(newline='') as handle:
        support = list(csv.DictReader(handle))
    selected = json.loads((root / SELECTION_PATH).read_bytes())['selection']['selected']
    keyed = {r['island_key']: r for r in selected}
    if len(keyed) != 96 or len(support) != 96 or len({r['island_key'] for r in support}) != 96:
        raise ValueError('Original96 island identity changed')
    if {r['island_key'] for r in support} != set(keyed):
        raise ValueError('Target support and selection identities differ')
    targets = [keyed[r['island_key']] for r in support if r['fixed_target_recording_support'] == '1']
    if len(targets) != 28:
        raise ValueError('Exactly28 frozen target islands required')
    return sorted(targets, key=lambda r: r['island_key'])


def plant_payload(wkt, offset=0):
    if offset < 0 or offset % 100:
        raise ValueError('Facet offset must be a nonnegative multiple of100')
    payload = pilot.payload(wkt, WINDOW, offset)
    payload['predicate']['predicates'][0]['value'] = 'P'
    return payload


def validate_count_response(raw):
    result = json.loads(raw)
    if (not isinstance(result, dict) or result.get('results') != []
            or type(result.get('count')) is not int or result['count'] < 0):
        raise ValueError('Source violated count-only response contract; response not retained')
    return result


def save_count_response(folder, prefix, request, raw):
    result = validate_count_response(raw)
    folder.mkdir(parents=True, exist_ok=True)
    receipt = {'request_sha256': sha(request), 'response_sha256': sha(raw),
               'retrieved_utc': pilot.utc(), 'response_bytes': len(raw),
               'count_only_validated': True}
    (folder / (prefix + '_request.json')).write_bytes(request)
    (folder / (prefix + '_response.json')).write_bytes(raw)
    write_json(folder / (prefix + '_receipt.json'), receipt)
    return result, receipt


def read_count_cache(folder, prefix, request):
    saved = (folder / (prefix + '_request.json')).read_bytes()
    receipt = json.loads((folder / (prefix + '_receipt.json')).read_bytes())
    if request != saved or receipt['request_sha256'] != sha(request):
        raise ValueError('Cached request identity differs')
    raw = (folder / (prefix + '_response.json')).read_bytes()
    if receipt['response_sha256'] != sha(raw):
        raise ValueError('Cached response identity differs')
    return validate_count_response(raw), receipt


def collect_early_facets(fetch):
    result = pilot.collect_facets(fetch, page_size=100, max_pages=100)
    years = result['facets']['YEAR']
    if any(not y.isdigit() or not 2010 <= int(y) <= 2017 for y in years):
        raise ValueError('Source contains a year outside the fixed early window')
    return result


def plant_support(result):
    if result['status'] != 'COMPLETE':
        return {'query_complete': False, 'record_count': None, 'recorded_species': None,
                'species_assignment_fraction': None, 'positive_years': None,
                'dataset_count': None, 'botanical_support_screen': False}
    count = result['count']
    facets = result['facets']
    fraction = sum(facets['SPECIES_KEY'].values()) / count if count else None
    years = len(facets['YEAR'])
    return {'query_complete': True, 'record_count': count,
            'recorded_species': len(facets['SPECIES_KEY']),
            'species_assignment_fraction': fraction, 'positive_years': years,
            'dataset_count': len(facets['DATASET_KEY']),
            'botanical_support_screen': bool(count >= 20 and fraction >= .9 and years >= 3)}


def checked_wkt(geometry, expected_sha):
    from shapely import orient_polygons, to_wkt
    geometry = orient_polygons(geometry, exterior_cw=False)
    if geometry.is_empty or not geometry.is_valid or geometry.geom_type not in ('Polygon', 'MultiPolygon'):
        raise ValueError('Invalid original island geometry')
    wkt = to_wkt(geometry, rounding_precision=-1, trim=True)
    if sha(wkt.encode()) != expected_sha:
        raise ValueError('Frozen exact geometry hash differs; no replacement permitted')
    return wkt


def prepare_geometries(source_cache, targets):
    from measure_island_port_geography_20261029 import source_landmasses, derived_key
    land = source_landmasses(source_cache)
    expected = {r['island_key']: r for r in targets}
    result = {}
    for _, row in land[(land.area_km2 >= 20) & (land.area_km2 < 20000)].iterrows():
        key = derived_key(row.parent_source_ID)
        if key in expected:
            result[key] = checked_wkt(row.geometry, expected[key]['geometry_sha256'])
    if set(result) != set(expected):
        raise ValueError('An original target geometry is missing')
    return result


class SourceClient:
    def __init__(self):
        self.last_start = 0.
        self.attempts = 0
        self.retries = 0

    def post(self, data):
        for attempt in range(3):
            time.sleep(max(0., 1.5 - (time.monotonic() - self.last_start)))
            self.last_start = time.monotonic()
            self.attempts += 1
            try:
                req = Request(pilot.ENDPOINT, data=data, headers={
                    'Content-Type': 'application/json',
                    'User-Agent': 'chocho-dated-island-plants/1.0 (+https://github.com/zuizui0223/chocho)'})
                with urlopen(req, timeout=45) as response:
                    raw = response.read(5_000_001)
                if len(raw) > 5_000_000:
                    raise ValueError('Count-only source response exceeds5MB bound')
                validate_count_response(raw)
                return raw
            except HTTPError as exc:
                if attempt == 2 or exc.code not in (429, 500, 502, 503, 504):
                    raise RuntimeError(f'GBIF HTTP{exc.code}; query remains missing') from None
                delay = 60 if exc.code == 429 else 5
                retry = exc.headers.get('Retry-After', '')
                if retry.isdigit():
                    delay = max(delay, int(retry))
                self.retries += 1
                print(json.dumps({'stage': 'source_retry', 'http': exc.code, 'delay_seconds': delay}), flush=True)
                deadline = time.monotonic() + delay
                while time.monotonic() < deadline:
                    time.sleep(min(30, max(0, deadline - time.monotonic())))
            except (URLError, TimeoutError, OSError):
                if attempt == 2:
                    raise RuntimeError('GBIF source transport failed3attempts; query remains missing') from None
                self.retries += 1
                time.sleep(5)
        raise RuntimeError('Source retry budget exhausted')


def acquire_target(row, wkt, cache, client):
    folder = cache / row['island_key']
    folder.mkdir(parents=True, exist_ok=True)
    receipts = []
    def fetch(offset, recheck=False):
        prefix = f'page_{offset:05d}' + ('_recheck' if recheck else '')
        request = pilot.canonical(plant_payload(wkt, offset))
        if (folder / (prefix + '_receipt.json')).exists():
            data, receipt = read_count_cache(folder, prefix, request)
        else:
            data, receipt = save_count_response(folder, prefix, request, client.post(request))
        receipts.append({**receipt, 'facet_offset': offset, 'is_first_page_recheck': recheck,
                         'cache_prefix': row['island_key'] + '/' + prefix})
        return data
    started = pilot.utc()
    result = pilot.query_safely(lambda: collect_early_facets(fetch))
    return {'island_key': row['island_key'], 'geometry_sha256': row['geometry_sha256'],
            'started_utc': started, 'finished_utc': pilot.utc(), 'source_pages': receipts, **result}


def archive_count_sources(cache, destination):
    files = sorted(cache.rglob('*.json'))
    for path in files:
        if path.name.endswith('_response.json'):
            validate_count_response(path.read_bytes())
    with destination.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as archive:
            for path in files:
                body = path.read_bytes()
                info = tarfile.TarInfo(path.relative_to(cache).as_posix())
                info.size, info.mtime, info.mode = len(body), 0, 0o644
                archive.addfile(info, io.BytesIO(body))
    return {'file': destination.name, 'bytes': destination.stat().st_size,
            'sha256': sha(destination.read_bytes()), 'files': len(files)}


def publish(output, lock, results, client, final=False):
    results = sorted(results, key=lambda r: r['island_key'])
    complete = [r for r in results if r['status'] == 'COMPLETE']
    summaries = [{'island_key': r['island_key'], 'status': r['status'], **plant_support(r)} for r in results]
    manifest = {'schema': 'chocho_fixed28_early_plants_20261030_v1',
                'status': 'FINISHED' if final else 'RUNNING', 'updated_utc': pilot.utc(),
                'query_scope': {'kingdom': 'Plantae', 'taxonKey': 'P', 'checklistKey': pilot.CHECKLIST,
                                'year_min': 2010, 'year_max': 2017, 'individual_record_limit': 0},
                'lock': lock, 'http_attempts_this_invocation': client.attempts,
                'http_retries_this_invocation': client.retries,
                'summary': {'requested_islands': 28, 'completed_islands': len(complete),
                            'missing_islands': len(results) - len(complete), 'pending_islands': 28 - len(results),
                            'complete_zero_islands': sum(r['count'] == 0 for r in complete),
                            'indexed_plant_records': sum(r['count'] for r in complete),
                            'botanical_support_islands': sum(r['botanical_support_screen'] for r in summaries)},
                'interpretation': 'Retrospective indexed plant recording in2010–2017. Missing matched records are not biological absence; point membership does not certify reported uncertainty-disk containment or wild/cultivated status.',
                'results': results}
    write_json(output / 'early_plants_manifest.json', manifest)
    if summaries:
        pilot.write_csv(output / 'early_plants_island_support.csv', summaries, list(summaries[0]))
    for field in pilot.FIELDS:
        rows = [{'island_key': r['island_key'], 'facet_key': k, 'record_count': n}
                for r in complete for k, n in sorted(r['facets'][field].items())]
        pilot.write_csv(output / ('early_plants_' + field.lower() + '_facets.csv'), rows,
                        ['island_key', 'facet_key', 'record_count'])
    return manifest


def run(source_cache, cache, output, root=ROOT):
    targets = load_targets(root)
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    spec = {'source_pins': PINS, 'plan_sha256': sha((root / PLAN_PATH).read_bytes()),
            'acquisition_script_sha256': sha(Path(__file__).read_bytes()),
            'target_geometry_sha256': {r['island_key']: r['geometry_sha256'] for r in targets},
            'query_contract': plant_payload('<exact frozen unsimplified CCW WKT>'),
            'target_roster_already_selected_using_later_butterfly_richness': True}
    lock_path = output / 'early_plants_run_lock.json'
    if lock_path.exists():
        lock = json.loads(lock_path.read_bytes())
        if lock['spec'] != spec:
            raise ValueError('Source acquisition specification changed since first lock')
    else:
        lock = {'locked_before_plant_occurrence_queries_utc': pilot.utc(), 'spec': spec}
        write_json(lock_path, lock)
    print(json.dumps({'stage': 'restore_exact28_geometries'}), flush=True)
    wkts = prepare_geometries(source_cache, targets)
    client, results = SourceClient(), []
    publish(output, lock, results, client)
    for index, row in enumerate(targets, 1):
        result = acquire_target(row, wkts[row['island_key']], cache, client)
        results.append(result)
        publish(output, lock, results, client)
        print(json.dumps({'stage': 'early_plants', 'completed': index, 'total': 28,
                          'island_key': row['island_key'], 'status': result['status'],
                          'count': result.get('count'), 'http_attempts': client.attempts}), flush=True)
    manifest = publish(output, lock, results, client, final=True)
    archive = archive_count_sources(cache, output / 'early_plants_raw_facets.tar.gz')
    manifest['source_archive'] = archive
    manifest['output_sha256'] = {p.name: sha(p.read_bytes()) for p in sorted(output.glob('early_plants_*.csv'))}
    write_json(output / 'early_plants_manifest.json', manifest)
    print(json.dumps({'stage': 'finished', **manifest['summary'], 'archive': archive}), flush=True)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-cache', type=Path, default=Path('/tmp/chocho-island-sources'))
    parser.add_argument('--cache-dir', type=Path, default=Path('/tmp/chocho-early-plants-cache'))
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.source_cache, args.cache_dir, args.output_dir)
