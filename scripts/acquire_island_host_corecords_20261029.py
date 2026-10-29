#!/usr/bin/env python3
"""Fresh original26 host-source acquisition and exact island co-recording audit.

This is a new transient-memory local pipeline, not a replay of the old
Actions-only point-proximity pipeline. Raw biological coordinates, occurrence
identifiers and response bodies are never written to cache, logs or artifacts.
The fixed96 primary panel and the independent full8260 island inventory are
frozen before source acquisition. Co-recording never establishes host use.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qsl, urlencode, urlparse
from urllib.request import Request, urlopen

from shapely import contains_properly, orient_polygons, to_wkt
from shapely.geometry import Point, box
from shapely.strtree import STRtree

from analyze_butterfly_host_point_proximity_20261026 import (
    build_queries as old_build_queries,
    records_url as old_records_url,
    source_inputs,
    valid_point as old_valid_point,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'provenance/reviewer_defenses/results/gbif_country_20261020'
SELECTION = ROOT / 'data/exploratory/island_expansion_20261029/expanded_island_panel_selection.json'
SELECTION_FILE_SHA = '5ccc868b89973b8142a7fdb4b3014c4394f134596ae2ef58e1179eaab68f5243'
SELECTION_SHA = 'db0f49c244309716b22e05ea5eeddd7d48749e46f46d2757114d98cf0df4df72'
CHECKLIST = 'd7dddbf4-2cf0-4f39-9b2a-bb099caae36c'
COL_CHECKLIST = '7ddf754f-d193-4cc9-b351-99906754a03b'
SCHEMA = 'chocho_original26_host_exact_island_corecords_20261029_v1'
WINDOWS = {'2010_2017': (2010, 2017), '2018_2025': (2018, 2025)}
PAGE_LIMIT, PER_QUERY_CAP, TOTAL_RECORD_CAP, MAX_HTTP_ATTEMPTS = 300, 5000, 12000, 200
MIN_WAIT_SECONDS = 1.5
# WGS84's minimum meridional radius is 6,335,439.327 m; rounding DOWN
# makes a spherical cap an outer bound for the same-radius ellipsoidal disk.
MIN_WGS84_CURVATURE_RADIUS_M = 6_335_439.0
COMPLETE = 'COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT'


def sha(body):
    return hashlib.sha256(body).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def utc():
    return datetime.now(timezone.utc).isoformat()


def load_frozen_sources():
    return source_inputs(BASE / 'gbif_all26_source_host_country_inventory_20261023.json',
                         BASE / 'gbif_positive48_gadm_lifestage_facets_20261020.json',
                         BASE / 'gbif_all26_host_gadm_source_20261023.json')


def build_source_queries(cases, plants):
    queries = old_build_queries(cases, plants)
    for key, q in queries.items():
        q['query_id'] = '|'.join(key)
        q['source_identity'] = key[1]
        q['checklist_key'] = CHECKLIST
        if key[0] == 'plant':
            q['wcvp_accepted_id'] = key[1]
    return queries


def records_url(query, offset=0, limit=PAGE_LIMIT):
    """Validate the original contract before adding its explicit namespace."""
    if not 0 <= limit <= PAGE_LIMIT or offset < 0:
        raise ValueError('Unbounded page size or negative offset')
    validated = old_records_url(query['source_url'], query['gbif_species_key'],
                                query['country'], offset)
    pairs = [(k, str(limit) if k == 'limit' else v)
             for k, v in parse_qsl(urlparse(validated).query)]
    pairs.append(('checklistKey', CHECKLIST))
    return 'https://api.gbif.org/v1/occurrence/search?' + urlencode(pairs)


def uncertainty_envelope(lon, lat, uncertainty_m):
    """Conservative bounding rectangle of a reported WGS84 uncertainty disk.

    A sphere using the minimum WGS84 curvature radius has a length metric
    no larger than WGS84's at every latitude. Its same-radius disk therefore
    encloses the WGS84 disk. The rectangle enclosing that spherical cap is a
    sufficient, deliberately conservative full-disk containment test. Wraps
    across the dateline and polar caps are not certified by this method.
    Source shoreline/coordinate accuracy beyond the reported radius is unknown.
    """
    if uncertainty_m is None or not math.isfinite(uncertainty_m) or uncertainty_m < 0:
        return None
    if uncertainty_m == 0:
        return Point(lon, lat)
    radius = uncertainty_m / MIN_WGS84_CURVATURE_RADIUS_M
    phi = math.radians(lat)
    if radius >= math.pi / 2 or abs(phi) + radius >= math.pi / 2:
        return None
    delta_lat = math.degrees(radius)
    delta_lon = math.degrees(math.asin(min(1.0, math.sin(radius) / math.cos(phi))))
    # Outward floating point guard; it can only reject additional strict cases.
    lo, hi = lon - delta_lon - 1e-10, lon + delta_lon + 1e-10
    if lo <= -180 or hi >= 180:
        return None
    return box(lo, lat - delta_lat - 1e-10, hi, lat + delta_lat + 1e-10)


class IslandJoiner:
    def __init__(self, islands, selected_keys=None):
        self.islands = islands
        self.geometries = [r['geometry'] for r in islands]
        self.tree = STRtree(self.geometries)
        self.selected_keys = (set(selected_keys) if selected_keys is not None
                              else {r['island_key'] for r in islands})

    def join(self, point):
        coordinate = Point(point['lon'], point['lat'])
        matched = [int(i) for i in self.tree.query(coordinate)
                   if self.geometries[int(i)].covers(coordinate)]
        if len(matched) > 1:
            raise ValueError('Point lies in multiple independent source island polygons')
        uncertainty = point['uncertainty']
        result = {'island_key': None, 'in_frozen96': False,
                  'year': point['year'], 'dataset': point.get('dataset'),
                  'on_polygon_boundary': False, 'strict_uncertainty_inside': False,
                  'strict_uncertainty_le1km_inside': False,
                  'uncertainty_status': 'UNKNOWN' if uncertainty is None else 'KNOWN',
                  'uncertainty_reported': uncertainty is not None}
        if not matched:
            return result
        ix = matched[0]
        geometry = self.geometries[ix]
        result['island_key'] = self.islands[ix]['island_key']
        result['in_frozen96'] = result['island_key'] in self.selected_keys
        result['on_polygon_boundary'] = not geometry.contains(coordinate)
        envelope = uncertainty_envelope(point['lon'], point['lat'], uncertainty)
        strict = (not result['on_polygon_boundary'] and envelope is not None
                  and bool(contains_properly(geometry, envelope)))
        result['strict_uncertainty_inside'] = strict
        result['strict_uncertainty_le1km_inside'] = strict and uncertainty <= 1000
        if result['on_polygon_boundary']:
            result['uncertainty_status'] = 'ON_POLYGON_BOUNDARY'
        elif uncertainty is not None:
            result['uncertainty_status'] = ('STRICT_KNOWN_INSIDE' if strict
                                            else 'KNOWN_NOT_CERTIFIED')
        return result


class RecordBudget:
    def __init__(self, limit=TOTAL_RECORD_CAP):
        self.limit = limit
        self.used = 0
        self.hard_stopped = False

    @property
    def remaining(self):
        return self.limit - self.used


class Throttle:
    """Public source calls only; raw response bytes die after JSON/hash use."""
    def __init__(self):
        self.last = None
        self.requests = 0
        self.retries = 0
        self.first_query_utc = None
        self.last_query_utc = None

    def fetch(self, url):
        for attempt in range(4):
            if self.requests >= MAX_HTTP_ATTEMPTS:
                raise RuntimeError('Global HTTP attempt cap reached')
            if self.last is not None:
                time.sleep(max(0, MIN_WAIT_SECONDS - (time.monotonic() - self.last)))
            self.last = time.monotonic()
            started = utc()
            self.first_query_utc = self.first_query_utc or started
            self.last_query_utc = started
            self.requests += 1
            try:
                request = Request(url, headers={'Accept': 'application/json',
                    'User-Agent': 'chocho-original26-island-corecords/1.0 (+https://github.com/zuizui0223/chocho)'})
                with urlopen(request, timeout=40) as response:
                    raw = response.read(20_000_001)
                    if response.status != 200 or len(raw) > 20_000_000:
                        raise ValueError('Source response HTTP status or byte bound failed')
                data = json.loads(raw)
                if not isinstance(data, dict) or 'count' not in data or 'results' not in data:
                    raise ValueError('GBIF occurrence response schema mismatch')
                receipt = {'body_sha256': sha(raw), 'response_bytes': len(raw),
                           'request_started_utc': started, 'retrieved_utc': utc(), 'http_status': 200}
                del raw
                return data, receipt
            except (HTTPError, URLError, TimeoutError) as exc:
                if isinstance(exc, HTTPError) and exc.code not in (429, 500, 502, 503, 504):
                    raise RuntimeError(f'Nonretryable GBIF HTTP {exc.code}') from None
                if attempt == 3:
                    raise RuntimeError('Repeated GBIF transport failure') from None
                self.retries += 1
                time.sleep(min(2 ** (attempt + 1), 12))
        raise RuntimeError('GBIF source fetch exhausted')


def validated_point(record, query):
    key = record.get('speciesKey', record.get('acceptedTaxonKey'))
    try:
        if key is None or int(key) != query['gbif_species_key']:
            return None, 'accepted_species_key_mismatch'
    except (TypeError, ValueError):
        return None, 'accepted_species_key_mismatch'
    return old_valid_point(record, query['gbif_species_key'], query['country'])


def retrieve_query(query, fetcher, joiner, budget):
    """Complete query or unknown; an exceeded cap never yields a partial zero."""
    receipt = {'query_id': query['query_id'], 'kind': query['kind'],
               'source_identity': query['source_identity'], 'taxon_label': query['label'],
               'country': query['country'], 'gbif_species_key': query['gbif_species_key'],
               'checklist_key': CHECKLIST,
               'original_facet_count': query['original_faceted_index_count'],
               'query_url': records_url(query), 'current_count': None,
               'valid_points': None, 'selected_island_points': None,
               'all8260_island_points': None, 'pages': [], 'started_utc': utc()}
    if query['kind'] == 'plant':
        receipt['wcvp_accepted_id'] = query['wcvp_accepted_id']
    if budget.hard_stopped:
        return {**receipt, 'status': 'NOT_STARTED_AFTER_HARD_CAP'}, []
    points, seen, excluded = [], set(), Counter()
    dataset_year = Counter()
    try:
        probe, meta = fetcher.fetch(records_url(query, limit=0))
        count = int(probe['count'])
        receipt['count_probe'] = meta
        receipt['current_count'] = count
        receipt['delta_vs_original_facet_count'] = count - receipt['original_facet_count']
        if count < 0 or probe['results'] != []:
            raise ValueError('Count-only source probe inconsistent')
        if count > PER_QUERY_CAP or count > budget.remaining:
            budget.hard_stopped = True
            return {**receipt, 'status': 'OVER_BOUNDED_ACQUISITION_CAP',
                    'reason': 'Complete query would exceed per-query or total row cap',
                    'finished_utc': utc()}, []
        for offset in range(0, count, PAGE_LIMIT):
            limit = min(PAGE_LIMIT, count - offset, budget.remaining)
            page, meta = fetcher.fetch(records_url(query, offset=offset, limit=limit))
            rows = page['results']
            if not isinstance(rows, list):
                raise ValueError('Occurrence results are not a list')
            budget.used += len(rows)
            receipt['pages'].append({'offset': offset, 'limit': limit, 'rows': len(rows), **meta})
            if budget.used > budget.limit:
                budget.hard_stopped = True
                raise ValueError('Source violated requested page and total row bound')
            if int(page['count']) != count or len(rows) != limit:
                raise ValueError('Count drift or short/overfull page during fresh pagination')
            for record in rows:
                identity = record.get('key')
                if identity is None or str(identity) in seen:
                    raise ValueError('Missing or duplicate occurrence identity during pagination')
                seen.add(str(identity))
                dataset_year[(str(record.get('datasetKey') or 'MISSING'), int(record.get('year') or 0))] += 1
                point, reason = validated_point(record, query)
                excluded[reason] += 1
                if point is not None:
                    points.append(joiner.join(point))
            # Only derived island/year/uncertainty values survive between pages.
            del rows, page
        if len(seen) != count:
            raise ValueError('Exhaustive source count and unique retrieved rows differ')
        receipt.update(status=COMPLETE, valid_points=len(points),
                       selected_island_points=sum(p['in_frozen96'] for p in points),
                       all8260_island_points=sum(p['island_key'] is not None for p in points),
                       strict_uncertainty_inside_points=sum(p['strict_uncertainty_inside'] for p in points),
                       uncertainty_unknown_points=sum(not p['uncertainty_reported'] for p in points),
                       excluded=dict(excluded),
                       dataset_year_counts=[{'dataset_key': d, 'year': y, 'indexed_records': n}
                                            for (d, y), n in sorted(dataset_year.items())],
                       finished_utc=utc())
        return receipt, points
    except (ValueError, TypeError, KeyError, HTTPError, URLError, RuntimeError, TimeoutError) as exc:
        if 'Global HTTP attempt cap' in str(exc):
            budget.hard_stopped = True
        # Error strings originate in validators or transport only, never a record body.
        receipt.update(status='INCOMPLETE_GBIF_SOURCE', reason=str(exc)[:240], finished_utc=utc())
        return receipt, []


def measure_case(butterflies, plants, island, window, source_complete):
    lo, hi = WINDOWS[window]
    fields = ('butterfly_points', 'plant_points', 'strict_butterfly_points',
              'strict_plant_points', 'corecorded_butterfly_indices_within3years',
              'corecorded_butterfly_indices_within3years_same_window',
              'corecorded_butterfly_indices_within3years_cross_window_only',
              'corecorded_butterfly_indices_same_year', 'corecorded_butterfly_indices_any_year',
              'strict_corecorded_butterfly_indices_within3years',
              'strict_le1km_corecorded_butterfly_indices_within3years')
    if not source_complete:
        return {'status': 'INCOMPLETE_REQUIRED_SOURCE', **dict.fromkeys(fields, None)}
    b = [p for p in butterflies if p['island_key'] == island and lo <= p['year'] <= hi]
    p = [p for p in plants if p['island_key'] == island]
    out = {'status': 'COMPLETE_OBSERVATIONAL_INDEX_ONLY', 'butterfly_points': len(b),
           'plant_points': sum(lo <= r['year'] <= hi for r in p),
           'strict_butterfly_points': sum(r['strict_uncertainty_inside'] for r in b),
           'strict_plant_points': sum(r['strict_uncertainty_inside'] and lo <= r['year'] <= hi for r in p)}
    for label, span in (('within3years', 3), ('same_year', 0), ('any_year', 100)):
        out['corecorded_butterfly_indices_' + label] = sum(
            any(abs(r['year'] - host['year']) <= span for host in p) for r in b)
    same_window = sum(any(abs(r['year'] - host['year']) <= 3 and lo <= host['year'] <= hi
                          for host in p) for r in b)
    out['corecorded_butterfly_indices_within3years_same_window'] = same_window
    out['corecorded_butterfly_indices_within3years_cross_window_only'] = (
        out['corecorded_butterfly_indices_within3years'] - same_window)
    for label, field in (('strict', 'strict_uncertainty_inside'),
                         ('strict_le1km', 'strict_uncertainty_le1km_inside')):
        out[label + '_corecorded_butterfly_indices_within3years'] = sum(
            r[field] and any(host[field] and abs(r['year'] - host['year']) <= 3 for host in p)
            for r in b)
    return out


def load_island_geometries(source_cache):
    from measure_island_port_geography_20261029 import source_landmasses, derived_key, SHA, URL
    selection_raw = SELECTION.read_bytes()
    lock = json.loads(selection_raw)
    if sha(selection_raw) != SELECTION_FILE_SHA or lock['selection_sha256'] != SELECTION_SHA \
            or sha(canonical(lock['selection'])) != SELECTION_SHA:
        raise ValueError('Frozen96 island selection identity changed')
    selected = {r['island_key']: r for r in lock['selection']['selected']}
    if len(selected) != 96:
        raise ValueError('Not exactly96 original selected islands')
    land = source_landmasses(Path(source_cache))
    universe = land[(land.area_km2 >= 5) & (land.area_km2 < 7_000_000)]
    if len(universe) != 8260:
        raise ValueError('Original >=5km2 physical-island universe differs from8260')
    islands = []
    for _, row in universe.iterrows():
        key = derived_key(row.parent_source_ID)
        geometry = orient_polygons(row.geometry, exterior_cw=False)
        if geometry.is_empty or not geometry.is_valid or geometry.geom_type not in ('Polygon', 'MultiPolygon'):
            raise ValueError('Invalid exact source island geometry')
        geometry_sha = sha(to_wkt(geometry, rounding_precision=-1, trim=True).encode())
        if key in selected and geometry_sha != selected[key]['geometry_sha256']:
            raise ValueError('Exact96 unsimplified polygon bytes differ from frozen selection')
        islands.append({'island_key': key, 'geometry': geometry, 'geometry_sha256': geometry_sha,
                        'area_km2': float(row.area_km2), 'in_frozen96': key in selected})
    if {r['island_key'] for r in islands if r['in_frozen96']} != set(selected):
        raise ValueError('Selected island missing from original source universe')
    public = [{k: v for k, v in r.items() if k != 'geometry'} for r in islands]
    public.sort(key=lambda r: r['island_key'])
    source = {'selection_file': str(SELECTION.relative_to(ROOT)),
              'selection_file_sha256': SELECTION_FILE_SHA, 'selection_sha256': SELECTION_SHA,
              'gshhg_url': URL, 'gshhg_zip_sha256': SHA, 'resolution': 'h', 'level': 1,
              'primary_frozen96_count': 96, 'diagnostic_island_universe_count': 8260,
              'diagnostic_island_rule': 'All original dissolved GSHHG L1 landmasses with area>=5km2 and <7000000km2; five continental units excluded',
              'diagnostic_universe_sha256': sha(canonical(public)),
              'geometry_policy': 'Exact unsimplified source polygons and frozen selection hashes; no outcome-based island replacement'}
    return islands, selected, public, source


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + '\n')


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def source_edges(cases, plants):
    return [{'butterfly': sp, 'country': country, 'wcvp_accepted_id': pid,
             'candidate_host_binomial': plants[(pid, country)]['plant_name']}
            for (sp, country), case in sorted(cases.items()) for pid, _ in case['hosts']]


def freeze_inputs(output, cases, plants, source, queries, geography, universe, plan):
    contract = {'schema': SCHEMA, 'source': source, 'geography': geography,
                'implementation_script_sha256': sha(Path(__file__).read_bytes()),
                'implementation_tests_sha256': sha((ROOT / 'tests/test_island_host_corecords_20261029.py').read_bytes()),
                'source_queries': list(queries.values()), 'original_source_edges': source_edges(cases, plants),
                'taxonomy': {'original_source_checklist_key': CHECKLIST,
                             'current_COL_bridge_checklist_key_not_queried_here': COL_CHECKLIST,
                             'old_numeric_keys_never_replaced_by_COL_keys': True},
                'new_runtime': 'Transient-memory local subprocess, separately authorized continuation; not old Actions pipeline',
                'raw_biological_response_or_coordinate_persistence': False,
                'window_rule': WINDOWS, 'matching_rule': 'Same exact physical island and absolute source year difference<=3; butterfly year assigns output window; hosts may cross window boundary within2010-2025',
                'strict_uncertainty_rule': 'Reported WGS84 uncertainty disk conservatively enclosed by min-curvature spherical-cap rectangle; entire rectangle strictly inside exact source polygon; missing is unknown; dateline/polar wraps uncertified',
                'caps': {'queries': 70, 'per_query_records': PER_QUERY_CAP, 'total_retrieved_records': TOTAL_RECORD_CAP,
                         'http_attempts_including_retries': MAX_HTTP_ATTEMPTS, 'minimum_call_start_spacing_seconds': MIN_WAIT_SECONDS},
                'plan_file': str(plan.relative_to(ROOT)) if plan.is_relative_to(ROOT) else str(plan),
                'plan_sha256': sha(plan.read_bytes())}
    output.mkdir(parents=True, exist_ok=True)
    path = output / 'source_island_query_lock.json'
    contract_sha = sha(canonical(contract))
    if path.exists():
        lock = json.loads(path.read_text())
        if lock['contract_sha256'] != contract_sha or sha(canonical(lock['contract'])) != contract_sha:
            raise ValueError('Pre-acquisition contract changed; do not overwrite previous lock')
    else:
        lock = {'locked_before_first_source_query_utc': utc(), 'contract_sha256': contract_sha, 'contract': contract}
        write_json(path, lock)
        write_csv(output / 'frozen_original_island_universe.csv', universe)
    return lock


def aggregate_source_island_year(queries, receipts, records):
    rows = []
    for key, q in queries.items():
        if receipts[key]['status'] != COMPLETE:
            continue
        bins = defaultdict(list)
        for point in records[key]:
            if point['island_key'] is not None:
                bins[(point['island_key'], point['year'])].append(point)
        for (island, year), points in sorted(bins.items()):
            rows.append({'query_id': q['query_id'], 'kind': q['kind'], 'taxon_label': q['label'],
                         'source_identity': q['source_identity'], 'country': q['country'],
                         'gbif_species_key': q['gbif_species_key'], 'island_key': island,
                         'in_frozen96': points[0]['in_frozen96'], 'year': year,
                         'window': '2010_2017' if year <= 2017 else '2018_2025',
                         'point_record_count': len(points),
                         'strict_uncertainty_inside_count': sum(p['strict_uncertainty_inside'] for p in points),
                         'strict_uncertainty_le1km_inside_count': sum(p['strict_uncertainty_le1km_inside'] for p in points),
                         'unknown_uncertainty_count': sum(not p['uncertainty_reported'] for p in points),
                         'on_polygon_boundary_count': sum(p['on_polygon_boundary'] for p in points),
                         'distinct_dataset_count': len({p['dataset'] for p in points if p['dataset'] is not None})})
    return rows


def build_case_tables(cases, plants, selected, receipts, records):
    primary, diagnostic, edges = [], [], []
    for (sp, country), case in sorted(cases.items()):
        bkey = ('butterfly', sp, country)
        pkeys = [('plant', pid, country) for pid, _ in case['hosts']]
        complete = all(receipts[key]['status'] == COMPLETE for key in [bkey, *pkeys])
        b = records.get(bkey, [])
        p = [r for key in pkeys for r in records.get(key, [])]
        present = {r['island_key'] for r in b + p if r['island_key'] is not None}
        base = {'butterfly': sp, 'country': country, 'candidate_plant_taxa': len(pkeys),
                'candidate_wcvp_ids': '|'.join(key[1] for key in pkeys)}
        for island in sorted(set(selected) | present):
            for window in WINDOWS:
                row = {**base, 'island_key': island, 'in_frozen96': island in selected,
                       'window': window, **measure_case(b, p, island, window, complete)}
                if island in selected:
                    primary.append(row)
                if island in present:
                    diagnostic.append(row)
        for pkey in pkeys:
            host = records.get(pkey, [])
            edge_islands = {r['island_key'] for r in b + host if r['island_key'] is not None}
            edge_complete = all(receipts[key]['status'] == COMPLETE for key in (bkey, pkey))
            for island in sorted(edge_islands):
                for window in WINDOWS:
                    edges.append({'butterfly': sp, 'country': country, 'wcvp_accepted_id': pkey[1],
                                  'candidate_host_binomial': plants[(pkey[1], country)]['plant_name'],
                                  'island_key': island, 'in_frozen96': island in selected, 'window': window,
                                  **measure_case(b, host, island, window, edge_complete)})
    return primary, diagnostic, edges


def scope_summary(rows, query_receipts, records, selected_only):
    complete_rows = [r for r in rows if r['status'] == 'COMPLETE_OBSERVATIONAL_INDEX_ONLY']
    positive = [r for r in complete_rows if r['corecorded_butterfly_indices_within3years'] > 0]
    complete_all = all(r['status'] == COMPLETE for r in query_receipts.values())
    counts = Counter()
    for key, points in records.items():
        for point in points:
            if point['island_key'] is not None and (not selected_only or point['in_frozen96']):
                counts[key[0]] += 1
    subset = {
            'case_island_windows_with_any3year_corecording': len(positive),
            'source_cases_with_any3year_corecording': len({(r['butterfly'], r['country']) for r in positive}),
            'islands_with_any3year_corecording': len({r['island_key'] for r in positive}),
            'butterfly_indices_with_any3year_source_host': sum(r['corecorded_butterfly_indices_within3years'] for r in complete_rows),
            'strict_butterfly_indices_with_any3year_strict_source_host': sum(r['strict_corecorded_butterfly_indices_within3years'] for r in complete_rows),
            'same_year_corecorded_butterfly_indices': sum(r['corecorded_butterfly_indices_same_year'] for r in complete_rows)}
    return {'all70_sources_complete': complete_all,
            'butterfly_indices_in_scope_complete_sources': counts['butterfly'],
            'unique_plant_country_indices_in_scope_complete_sources': counts['plant'],
            'complete_case_island_window_rows': len(complete_rows),
            **(subset if complete_all else dict.fromkeys(subset, None)),
            'complete_source_subset': subset,
            'unknown_sources_are_never_encoded_as_zero': True}


def run(args):
    cases, plants, original = load_frozen_sources()
    queries = build_source_queries(cases, plants)
    islands, selected, universe, geography = load_island_geometries(args.source_cache)
    lock = freeze_inputs(args.output, cases, plants, original, queries, geography, universe, args.plan)
    print(json.dumps({'stage': 'PREQUERY_SOURCE_AND_GEOMETRY_LOCK_VERIFIED',
                      'contract_sha256': lock['contract_sha256'], 'primary_islands': 96,
                      'diagnostic_islands': 8260, 'source_queries': 70}), flush=True)
    if args.prepare_only:
        return
    if (args.output / 'host_corecord_manifest.json').exists():
        raise ValueError('Acquisition output exists; never silently overwrite a live snapshot')
    joiner = IslandJoiner(islands, selected_keys=selected)
    throttle, budget = Throttle(), RecordBudget()
    receipts, records = {}, {}
    for number, (key, q) in enumerate(queries.items(), 1):
        receipts[key], records[key] = retrieve_query(q, throttle, joiner, budget)
        print(json.dumps({'query_number': number, 'query_id': q['query_id'],
                          'status': receipts[key]['status'], 'current_count': receipts[key]['current_count'],
                          'rows_acquired_total': budget.used}), flush=True)
    source_rows = aggregate_source_island_year(queries, receipts, records)
    primary, diagnostic, edges = build_case_tables(cases, plants, selected, receipts, records)
    all_complete = all(r['status'] == COMPLETE for r in receipts.values())
    write_json(args.output / 'source_query_receipts.json', {'schema': SCHEMA,
        'first_source_query_utc': throttle.first_query_utc, 'last_source_query_utc': throttle.last_query_utc,
        'contract_sha256': lock['contract_sha256'], 'queries': list(receipts.values())})
    write_csv(args.output / 'source_island_year_counts.csv', source_rows,
              list(source_rows[0]) if source_rows else ['query_id', 'island_key', 'year', 'point_record_count'])
    write_csv(args.output / 'frozen96_case_island_window_corecords.csv', primary)
    write_csv(args.output / 'original8260_positive_source_case_inventory.csv', diagnostic,
              list(diagnostic[0]) if diagnostic else list(primary[0]))
    write_csv(args.output / 'original46_edge_island_window_inventory.csv', edges,
              list(edges[0]) if edges else ['butterfly', 'country', 'wcvp_accepted_id', 'island_key', 'window'])
    manifest = {'schema': SCHEMA, 'status': ('COMPLETE_NEW_LIVE_70_SOURCE_ISLAND_JOIN' if all_complete
                                           else 'INCOMPLETE_NEW_LIVE_SOURCE_ISLAND_JOIN'),
                'executed_utc': utc(), 'first_source_query_utc': throttle.first_query_utc,
                'last_source_query_utc': throttle.last_query_utc, 'contract_sha256': lock['contract_sha256'],
                'source': original, 'geography': geography,
                'acquisition': {'total_source_queries': 70, 'complete_source_queries': sum(r['status'] == COMPLETE for r in receipts.values()),
                                'http_attempts_including_retries': throttle.requests, 'rate_limit_retries': throttle.retries,
                                'retrieved_occurrence_rows': budget.used, 'hard_stopped': budget.hard_stopped,
                                'new_retrieval_not_original_snapshot_replay': True,
                                'raw_coordinates_occurrence_ids_or_response_bodies_preserved': False,
                                'runtime': lock['contract']['new_runtime'],
                                'source_queries_with_count_drift': sum(r.get('delta_vs_original_facet_count', 0) != 0 for r in receipts.values()),
                                'complete_butterfly_source_count': sum(r['current_count'] for r in receipts.values() if r['kind'] == 'butterfly' and r['status'] == COMPLETE),
                                'complete_plant_source_count': sum(r['current_count'] for r in receipts.values() if r['kind'] == 'plant' and r['status'] == COMPLETE)},
                'primary_frozen96': scope_summary(primary, receipts, records, True),
                'diagnostic_original8260': scope_summary(diagnostic, receipts, records, False),
                'interpretation': [
                    'Primary96 remains the pre-existing geography-selected panel with no replacement.',
                    'All8260 is a separate descriptive inventory of records from already selected positive26 country cases; it is not a new representative cohort or effect test.',
                    'Any-point containment establishes source spatial record presence only; an uncertain point can lie outside the true island.',
                    'Strict containment is conditional on the reported coordinate uncertainty covering actual error and on supplied GSHHG geometry; it does not cover source shoreline age, map accuracy or omitted lakes.',
                    'Same island and within3calendar years is observational co-recording, not feeding, host use, wild plant availability, arrival, recruitment, climate niche expansion or causality.',
                    'The two output windows follow butterfly year; candidate plant years may cross the2017/2018 boundary within the fixed2010-2025 acquisition range.',
                    'No matched source records is no indexed co-detection, not biological absence. Failed or capped sources remain incomplete.',
                    'The current COL XR island panel is a different taxonomic namespace; this acquisition retains original numeric GBIF Backbone keys explicitly.'],
                'outputs': {}}
    for path in sorted(args.output.iterdir()):
        if path.is_file() and path.name != 'host_corecord_manifest.json':
            manifest['outputs'][path.name] = {'sha256': sha(path.read_bytes()), 'bytes': path.stat().st_size}
    write_json(args.output / 'host_corecord_manifest.json', manifest)
    print(json.dumps({'stage': 'DERIVED_EXPORT_COMPLETE', 'status': manifest['status'],
                      'primary_frozen96': manifest['primary_frozen96'],
                      'diagnostic_original8260': manifest['diagnostic_original8260']}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-cache', type=Path, default=Path('/tmp/chocho-island-sources'))
    parser.add_argument('--output', type=Path, default=ROOT / 'data/exploratory/island_host_corecords_20261029')
    parser.add_argument('--plan', type=Path, required=True, help='Frozen root implementation/source plan, created before any source query')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
