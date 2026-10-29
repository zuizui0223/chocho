"""Geometric and acquisition controls for the frozen 26-source / 96-island join."""
import copy
import importlib.util
import json
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest
from shapely.geometry import Polygon, box

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


@pytest.fixture
def m():
    path = ROOT / 'scripts/acquire_island_host_corecords_20261029.py'
    assert path.exists(), 'fresh source-locked acquisition implementation is missing'
    spec = importlib.util.spec_from_file_location('host_corecords_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def query(m):
    cases, plants, source = m.load_frozen_sources()
    return next(iter(m.build_source_queries(cases, plants).values()))


def raw_record(q, key=1, lon=10.5, lat=10.5, year=2017, uncertainty=100):
    return {'key': key, 'decimalLongitude': lon, 'decimalLatitude': lat,
            'year': year, 'countryCode': q['country'],
            'speciesKey': q['gbif_species_key'], 'datasetKey': 'synthetic-dataset',
            'coordinateUncertaintyInMeters': uncertainty,
            'basisOfRecord': 'HUMAN_OBSERVATION', 'occurrenceStatus': 'PRESENT'}


class SyntheticPages:
    """Only the unavailable public network is replaced; join and validation are real."""
    def __init__(self, pages):
        self.pages = iter(pages)
        self.calls = 0

    def fetch(self, url):
        parameters = parse_qs(urlparse(url).query)
        assert parameters['checklistKey'] == ['d7dddbf4-2cf0-4f39-9b2a-bb099caae36c']
        data = next(self.pages)
        self.calls += 1
        return data, {'body_sha256': '0' * 64, 'retrieved_utc': '2026-10-29T00:00:00+00:00',
                      'response_bytes': 100, 'http_status': 200}


def test_source_keys_are_not_reinterpreted_as_current_col_keys(m):
    cases, plants, source = m.load_frozen_sources()
    queries = m.build_source_queries(cases, plants)
    assert len(cases) == 26 and len(plants) == 44 and len(queries) == 70
    assert source['frozen_butterfly_source_count'] == 467
    assert source['frozen_plant_source_count'] == 3595
    q = query(m)
    params = parse_qs(urlparse(m.records_url(q, limit=0)).query)
    assert params['taxonKey'] == ['1944574']
    assert params['checklistKey'] == ['d7dddbf4-2cf0-4f39-9b2a-bb099caae36c']
    assert params['year'] == ['2010,2025']
    assert params['country'] == ['CU']
    assert params['limit'] == ['0'] and 'facet' not in params
    changed = copy.deepcopy(q)
    changed['source_url'] = changed['source_url'].replace('2010%2C2025', '2018%2C2025')
    with pytest.raises(ValueError):
        m.records_url(changed)


def test_membership_keeps_unknown_uncertainty_and_boundary_separate(m):
    joiner = m.IslandJoiner([{'island_key': 'A', 'geometry': box(10, 10, 11, 11)}])
    p = {'lat': 10.5, 'lon': 10.5, 'year': 2017, 'dataset': 'd', 'uncertainty': 100}
    inside = joiner.join(p)
    assert inside['island_key'] == 'A' and inside['strict_uncertainty_inside']
    assert inside['strict_uncertainty_le1km_inside']
    unknown = joiner.join({**p, 'uncertainty': None})
    assert unknown['island_key'] == 'A' and not unknown['strict_uncertainty_inside']
    assert unknown['uncertainty_status'] == 'UNKNOWN'
    boundary = joiner.join({**p, 'lon': 10, 'uncertainty': 0})
    assert boundary['island_key'] == 'A' and boundary['on_polygon_boundary']
    assert not boundary['strict_uncertainty_inside']
    edge = joiner.join({**p, 'lon': 10.0001, 'uncertainty': 100})
    assert edge['island_key'] == 'A' and not edge['strict_uncertainty_inside']
    assert joiner.join({**p, 'lon': 12})['island_key'] is None
    assert 'lat' not in inside and 'lon' not in inside


def test_polygon_holes_and_large_uncertainty_do_not_become_strict_presence(m):
    geometry = Polygon([(10,10), (11,10), (11,11), (10,11)],
                       holes=[[(10.4,10.4), (10.6,10.4), (10.6,10.6), (10.4,10.6)]])
    joiner = m.IslandJoiner([{'island_key': 'A', 'geometry': geometry}])
    point = {'lat': 10.5, 'lon': 10.5, 'year': 2017, 'dataset': None, 'uncertainty': 0}
    assert joiner.join(point)['island_key'] is None
    assert not joiner.join({**point, 'lon': 10.2, 'uncertainty': 50_000})['strict_uncertainty_inside']
    assert m.uncertainty_envelope(179.99, 0, 10_000) is None
    with pytest.raises(ValueError, match='multiple'):
        m.IslandJoiner([{'island_key': 'A', 'geometry': box(10,10,11,11)},
                        {'island_key': 'B', 'geometry': box(10,10,11,11)}]).join(point)


def test_corecording_requires_same_island_and_preserves_butterfly_window(m):
    def p(island, year, strict=True):
        return {'island_key': island, 'year': year, 'strict_uncertainty_inside': strict,
                'strict_uncertainty_le1km_inside': strict}
    b = [p('A', 2017), p('A', 2018), p('B', 2017), p(None, 2017)]
    plants = [p('A', 2020), p('B', 2021)]
    out = m.measure_case(b, plants, 'A', '2010_2017', source_complete=True)
    assert out['butterfly_points'] == 1 and out['plant_points'] == 0
    assert out['corecorded_butterfly_indices_within3years'] == 1
    assert out['corecorded_butterfly_indices_within3years_same_window'] == 0
    assert out['corecorded_butterfly_indices_within3years_cross_window_only'] == 1
    assert out['corecorded_butterfly_indices_same_year'] == 0
    assert out['strict_corecorded_butterfly_indices_within3years'] == 1
    assert m.measure_case(b, plants, 'B', '2010_2017', True)['corecorded_butterfly_indices_within3years'] == 0
    incomplete = m.measure_case(b, plants, 'A', '2010_2017', source_complete=False)
    assert incomplete['corecorded_butterfly_indices_within3years'] is None
    assert incomplete['status'] == 'INCOMPLETE_REQUIRED_SOURCE'


def test_true_zero_and_transport_failure_are_distinct(m):
    q = query(m)
    joiner = m.IslandJoiner([])
    zero, pts = m.retrieve_query(q, SyntheticPages([{'count': 0, 'results': []}]), joiner, m.RecordBudget(12000))
    assert zero['status'] == 'COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT'
    assert zero['current_count'] == 0 and pts == []
    bad, pts = m.retrieve_query(q, SyntheticPages([{'count': 2, 'results': []},
                                                 {'count': 2, 'results': [raw_record(q)]}]), joiner, m.RecordBudget(12000))
    assert bad['status'] == 'INCOMPLETE_GBIF_SOURCE'
    assert bad['valid_points'] is None and pts == []


def test_complete_query_emits_receipts_without_raw_record_identity(m):
    q = query(m)
    joiner = m.IslandJoiner([{'island_key': 'A', 'geometry': box(10,10,11,11)}])
    pages = SyntheticPages([{'count': 2, 'results': []},
                            {'count': 2, 'results': [raw_record(q, key=999999999999999),
                                                      raw_record(q, key=888888888888888, lon=12, year=2021)]}])
    budget = m.RecordBudget(12000)
    receipt, points = m.retrieve_query(q, pages, joiner, budget)
    assert receipt['status'] == 'COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT'
    assert receipt['current_count'] == 2 and receipt['valid_points'] == 2
    assert receipt['selected_island_points'] == 1 and budget.used == 2
    assert receipt['dataset_year_counts'] == [
        {'dataset_key': 'synthetic-dataset', 'year': 2017, 'indexed_records': 1},
        {'dataset_key': 'synthetic-dataset', 'year': 2021, 'indexed_records': 1}]
    assert receipt['delta_vs_original_facet_count'] == 0
    text = json.dumps(receipt)
    assert all(x not in text for x in ('decimalLatitude', 'decimalLongitude', '999999999999999', '888888888888888'))
    assert all('lat' not in p and 'lon' not in p for p in points)


def test_caps_prevent_partial_source_acquisition_and_count_drift_fails_closed(m):
    q = query(m)
    pages = SyntheticPages([{'count': 5001, 'results': []}])
    budget = m.RecordBudget(12000)
    receipt, points = m.retrieve_query(q, pages, m.IslandJoiner([]), budget)
    assert receipt['status'] == 'OVER_BOUNDED_ACQUISITION_CAP'
    assert receipt['valid_points'] is None and points == [] and budget.hard_stopped
    pages = SyntheticPages([{'count': 3, 'results': []}])
    budget = m.RecordBudget(2)
    receipt, points = m.retrieve_query(q, pages, m.IslandJoiner([]), budget)
    assert receipt['status'] == 'OVER_BOUNDED_ACQUISITION_CAP' and budget.used == 0
    pages = SyntheticPages([{'count': 1, 'results': []}, {'count': 2, 'results': [raw_record(q)]}])
    receipt, points = m.retrieve_query(q, pages, m.IslandJoiner([]), m.RecordBudget(12000))
    assert receipt['status'] == 'INCOMPLETE_GBIF_SOURCE' and points == []


def test_duplicate_keys_and_taxon_mismatch_do_not_inflate_membership(m):
    q = query(m)
    pages = SyntheticPages([{'count': 2, 'results': []},
                            {'count': 2, 'results': [raw_record(q), raw_record(q)]}])
    receipt, points = m.retrieve_query(q, pages, m.IslandJoiner([]), m.RecordBudget(12000))
    assert receipt['status'] == 'INCOMPLETE_GBIF_SOURCE' and points == []
    wrong = raw_record(q)
    wrong['speciesKey'] = 'COL-WRONG-NAMESPACE'
    pages = SyntheticPages([{'count': 1, 'results': []}, {'count': 1, 'results': [wrong]}])
    receipt, points = m.retrieve_query(q, pages, m.IslandJoiner([]), m.RecordBudget(12000))
    assert receipt['status'] == 'COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT'
    assert receipt['valid_points'] == 0 and receipt['excluded']['accepted_species_key_mismatch'] == 1
    assert points == []


def test_incomplete_scope_summary_does_not_report_missing_sources_as_zero(m):
    receipt = {('butterfly', 'synthetic', 'CU'): {'status': 'INCOMPLETE_GBIF_SOURCE'}}
    result = m.scope_summary([], receipt, {}, selected_only=True)
    assert result['all70_sources_complete'] is False
    assert result['butterfly_indices_with_any3year_source_host'] is None
    assert result['strict_butterfly_indices_with_any3year_strict_source_host'] is None
    assert result['source_cases_with_any3year_corecording'] is None
    assert result['complete_source_subset']['butterfly_indices_with_any3year_source_host'] == 0
