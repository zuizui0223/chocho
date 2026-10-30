"""Catalogue-compatible recording is not host absence or causal establishment."""
import importlib.util
from pathlib import Path
import sys

import pytest

PATH = Path(__file__).resolve().parents[1] / 'scripts/analyze_island_dated_hosts_20261030.py'


def module():
    assert PATH.exists(), 'Dated-host integration implementation is missing'
    sys.path.insert(0, str(PATH.parent))
    spec = importlib.util.spec_from_file_location('dated_hosts_test', PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def bridge(pid, key='', eligible=True, status='SUCCESS', verification=''):
    return {'wcvp_accepted_id': pid, 'current_accepted_key': key,
            'conservative_join_eligible': str(eligible), 'query_status': status,
            'accepted_verification_query_status': verification}


def test_two_source_ids_one_current_concept_are_not_double_counted():
    m = module()
    edges = [{'original_name': 'Butterfly A', 'wcvp_accepted_id': '1'},
             {'original_name': 'Butterfly A', 'wcvp_accepted_id': '2'}]
    catalogue = m.build_catalogues(edges, [bridge('1', 'plant'), bridge('2', 'plant')])['Butterfly A']
    result = m.measure_hosts(catalogue, {'status': 'COMPLETE', 'facets': {'SPECIES_KEY': {'plant': 5}}})
    assert result['catalogue_host_ids'] == 2
    assert result['resolved_host_concepts'] == 1
    assert result['recorded_host_concepts'] == 1
    assert result['indexed_host_records'] == 5
    assert result['host_record_status'] == 'RECORDED_KNOWN_HOST'


def test_query_missing_is_null_and_not_no_match():
    m = module()
    cat = m.build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}], [bridge('1', 'p')])['A']
    result = m.measure_hosts(cat, {'status': 'MISSING_QUERY_FAILED'})
    assert result['recorded_host_concepts'] is None
    assert result['indexed_host_records'] is None
    assert result['confirmed_host_record'] is None
    assert result['host_record_status'] == 'PLANT_QUERY_MISSING'


def test_unresolved_catalogue_no_match_is_not_complete_catalogue_no_match():
    m = module()
    cat = m.build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}],
                             [bridge('1', eligible=False, status='UNQUERIED_SOURCE_NAME')])['A']
    result = m.measure_hosts(cat, {'status': 'COMPLETE', 'facets': {'SPECIES_KEY': {'other': 10}}})
    assert result['host_record_status'] == 'NO_MATCH_UNRESOLVED_HOST_CATALOGUE'
    assert result['unresolved_host_ids'] == 1
    assert result['missing_taxonomy_source_host_ids'] == 0
    assert result['confirmed_host_record'] is False


def test_failed_accepted_verification_is_explicit_missing_taxonomy():
    m = module()
    cat = m.build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}],
                             [bridge('1', 'p', eligible=False, verification='FAILED')])['A']
    assert cat['missing_taxonomy_source_host_ids'] == 1
    assert not cat['taxonomy_source_complete']


def test_matching_confirmed_host_survives_other_unresolved_links():
    m = module()
    edges = [{'original_name': 'A', 'wcvp_accepted_id': '1'}, {'original_name': 'A', 'wcvp_accepted_id': '2'}]
    cat = m.build_catalogues(edges, [bridge('1', 'p'), bridge('2', eligible=False, status='FAILED')])['A']
    result = m.measure_hosts(cat, {'status': 'COMPLETE', 'facets': {'SPECIES_KEY': {'p': 2}}})
    assert result['confirmed_host_record'] is True
    assert result['unresolved_host_ids'] == 1
    assert result['missing_taxonomy_source_host_ids'] == 1


def test_complete_resolved_query_with_no_match_means_record_zero_only():
    m = module()
    cat = m.build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}], [bridge('1', 'p')])['A']
    result = m.measure_hosts(cat, {'status': 'COMPLETE', 'facets': {'SPECIES_KEY': {}}})
    assert result['host_record_status'] == 'NO_MATCH_COMPLETE_RESOLVED_CATALOGUE'
    assert result['recorded_host_concepts'] == 0 and result['indexed_host_records'] == 0
    assert 'absent' not in result


@pytest.mark.parametrize('rows', [[bridge('1', 'p'), bridge('1', 'q')], []])
def test_duplicate_or_missing_taxonomy_source_identity_fails(rows):
    with pytest.raises(ValueError):
        module().build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}], rows)


def test_eligible_row_with_failed_source_cannot_be_used():
    with pytest.raises(ValueError, match='eligible'):
        module().build_catalogues([{'original_name': 'A', 'wcvp_accepted_id': '1'}],
                                  [bridge('1', 'p', eligible=True, status='FAILED')])


def test_missing_record_values_are_excluded_from_complete_summary_denominator():
    rows = [{'late_positive': 1, 'plant_query_complete': True, 'confirmed_host_record': True},
            {'late_positive': 1, 'plant_query_complete': False, 'confirmed_host_record': None},
            {'late_positive': 0, 'plant_query_complete': True, 'confirmed_host_record': False}]
    result = module().recording_summary(rows)
    assert result['all_candidate_cells'] == 3
    assert result['complete_plant_query_cells'] == 2
    assert result['all_later_positive_cells'] == 2
    assert result['complete_plant_query_later_positive_cells'] == 1
    assert result['later_positive_confirmed_host_cells'] == 1
    assert result['later_positive_confirmed_host_fraction_among_complete'] == 1.0
