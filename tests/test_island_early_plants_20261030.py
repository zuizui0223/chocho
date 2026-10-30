"""Protect early-only, exact-island and no-point acquisition contracts."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / 'scripts/acquire_island_early_plants_20261030.py'


def module():
    assert PATH.exists(), 'Early plant acquisition implementation is missing'
    sys.path.insert(0, str(PATH.parent))
    spec = importlib.util.spec_from_file_location('early_plants_test', PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def page(count=2, results=None):
    return {'count': count, 'results': [] if results is None else results,
            'facets': [
                {'field': 'SPECIES_KEY', 'counts': [{'name': 'A', 'count': count}] if count else []},
                {'field': 'BASIS_OF_RECORD', 'counts': [{'name': 'HUMAN_OBSERVATION', 'count': count}] if count else []},
                {'field': 'DATASET_KEY', 'counts': [{'name': 'd', 'count': count}] if count else []},
                {'field': 'YEAR', 'counts': [{'name': '2017', 'count': count}] if count else []}]}


def test_payload_cannot_retrieve_later_years_or_points():
    m = module()
    wkt = 'POLYGON ((0 0, 1 0, 1 1, 0 0))'
    p = m.plant_payload(wkt, 100)
    filters = p['predicate']['predicates']
    assert filters[0] == {'type': 'equals', 'key': 'TAXON_KEY', 'value': 'P',
                          'checklistKey': '7ddf754f-d193-4cc9-b351-99906754a03b'}
    assert {'type': 'greaterThanOrEquals', 'key': 'YEAR', 'value': '2010'} in filters
    assert {'type': 'lessThanOrEquals', 'key': 'YEAR', 'value': '2017'} in filters
    assert {'type': 'within', 'geometry': wkt} in filters
    assert p['limit'] == 0 and p['facetOffset'] == 100
    assert p['parameters']['CHECKLIST_KEY'] == ['7ddf754f-d193-4cc9-b351-99906754a03b']


def test_actual_target_roster_is_exactly_frozen28():
    rows = module().load_targets(ROOT)
    assert len(rows) == len({r['island_key'] for r in rows}) == 28
    assert all(len(r['geometry_sha256']) == 64 for r in rows)


def test_changed_source_hash_fails_instead_of_reselecting(tmp_path):
    m = module()
    p = tmp_path / m.TARGET_PATH
    p.parent.mkdir(parents=True)
    p.write_text('changed')
    with pytest.raises(ValueError, match='source hash'):
        m.load_targets(tmp_path)


def test_response_containing_individual_records_is_never_cached(tmp_path):
    m = module()
    raw = json.dumps(page(results=[{'decimalLatitude': 5}])).encode()
    with pytest.raises(ValueError, match='count-only'):
        m.save_count_response(tmp_path, 'p', b'{}', raw)
    assert not list(tmp_path.iterdir())


def test_cached_response_and_request_identity_are_both_checked(tmp_path):
    m = module()
    raw = json.dumps(page()).encode()
    m.save_count_response(tmp_path, 'p', b'{"a":1}', raw)
    assert m.read_count_cache(tmp_path, 'p', b'{"a":1}')[0]['count'] == 2
    with pytest.raises(ValueError, match='request'):
        m.read_count_cache(tmp_path, 'p', b'{"a":2}')
    (tmp_path / 'p_response.json').write_text(json.dumps(page(0)))
    with pytest.raises(ValueError, match='response'):
        m.read_count_cache(tmp_path, 'p', b'{"a":1}')


def test_stable_zero_query_is_complete_not_missing():
    m = module()
    result = m.collect_early_facets(lambda *args, **kwargs: page(0))
    assert result['count'] == 0 and result['first_page_recheck_passed']
    summary = m.plant_support({'status': 'COMPLETE', **result})
    assert summary['query_complete'] and summary['record_count'] == 0
    assert summary['species_assignment_fraction'] is None
    assert summary['botanical_support_screen'] is False


def test_changed_total_and_out_of_window_year_fail():
    m = module()
    with pytest.raises(ValueError, match='changed'):
        m.collect_early_facets(lambda offset, recheck=False: page(3 if recheck else 2))
    bad = page()
    bad['facets'][-1]['counts'][0]['name'] = '2018'
    with pytest.raises(ValueError, match='early'):
        m.collect_early_facets(lambda *args, **kwargs: bad)


def test_missing_query_cannot_be_converted_to_zero():
    s = module().plant_support({'status': 'MISSING_QUERY_FAILED', 'count': None})
    assert s['query_complete'] is False and s['record_count'] is None
    assert s['recorded_species'] is None and s['botanical_support_screen'] is False


def test_support_screen_uses_species_coverage_and_three_years():
    m = module()
    base = {'status': 'COMPLETE', 'count': 20, 'facets': {
        'SPECIES_KEY': {'A': 18}, 'YEAR': {'2010': 6, '2011': 7, '2017': 7},
        'DATASET_KEY': {'d': 20}, 'BASIS_OF_RECORD': {'HUMAN_OBSERVATION': 20}}}
    assert m.plant_support(base)['botanical_support_screen'] is True
    base['facets']['YEAR'] = {'2010': 10, '2017': 10}
    assert m.plant_support(base)['botanical_support_screen'] is False


def test_geometry_hash_mismatch_is_never_replaced_by_rectangle():
    m = module()
    from shapely.geometry import Polygon
    polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 0)])
    with pytest.raises(ValueError, match='geometry'):
        m.checked_wkt(polygon, '0' * 64)
