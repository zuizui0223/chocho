import importlib.util
from pathlib import Path
import sys

import pytest

PATH = Path(__file__).resolve().parents[1] / 'scripts/repair_island_plant_facets_20261030.py'


def module():
    assert PATH.exists(), 'Facet repair implementation is missing'
    sys.path.insert(0, str(PATH.parent))
    spec = importlib.util.spec_from_file_location('repair_plants_test', PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def response(facets, count=4):
    return {'count': count, 'results': [], 'facets': [
        {'field': f, 'counts': [{'name': k, 'count': n} for k,n in counts]}
        for f, counts in facets.items()]}


def base():
    return {'SPECIES_KEY': [('a', 2), ('b', 2)], 'DATASET_KEY': [('d', 4)],
            'BASIS_OF_RECORD': [('HUMAN_OBSERVATION', 4)], 'YEAR': [('2010', 2), ('2017', 2)]}


def test_tied_order_changes_are_accepted_only_with_complete_empty_terminal_probes():
    m = module()
    calls = []
    def fetch(field, offset, limit, recheck=False):
        calls.append((field, offset, limit, recheck))
        if field:
            return response({field: []})
        facets = base()
        if recheck:
            facets['SPECIES_KEY'].reverse()
        return response(facets)
    result = m.collect_single_page(fetch)
    assert result['facets']['SPECIES_KEY'] == {'a': 2, 'b': 2}
    assert result['terminal_probes_passed']
    assert set(c[0] for c in calls if c[0]) == {'SPECIES_KEY', 'BASIS_OF_RECORD', 'DATASET_KEY', 'YEAR'}
    assert ('SPECIES_KEY', 2, 1, False) in calls
    assert ('DATASET_KEY', 1, 1, False) in calls


def test_silent_server_cap_is_detected_at_actual_bin_count():
    def fetch(field, offset, limit, recheck=False):
        return response({field: [('hidden', 1)]}) if field else response(base())
    with pytest.raises(ValueError, match='terminal'):
        module().collect_single_page(fetch)


def test_recheck_changed_map_is_rejected_even_with_same_total():
    def fetch(field, offset, limit, recheck=False):
        if field:
            return response({field: []})
        f = base()
        if recheck:
            f['SPECIES_KEY'] = [('a', 3), ('b', 1)]
        return response(f)
    with pytest.raises(ValueError, match='changed'):
        module().collect_single_page(fetch)


def test_terminal_count_change_is_rejected():
    def fetch(field, offset, limit, recheck=False):
        return response({field: []}, count=5) if field else response(base())
    with pytest.raises(ValueError, match='count'):
        module().collect_single_page(fetch)


def test_duplicate_bin_within_single_page_is_rejected():
    f = base()
    f['SPECIES_KEY'] = [('a', 2), ('a', 2)]
    with pytest.raises(ValueError, match='Duplicate'):
        module().collect_single_page(lambda *args, **kwargs: response(f))


def test_original_complete_results_are_not_replaced():
    m = module()
    original = [{'island_key': 'a', 'status': 'COMPLETE', 'count': 4},
                {'island_key': 'b', 'status': 'MISSING_QUERY_FAILED', 'count': None}]
    replacement = [{'island_key': 'b', 'status': 'COMPLETE', 'count': 8}]
    merged = m.merge_results(original, replacement)
    assert merged[0] == original[0]
    assert merged[1]['count'] == 8
    with pytest.raises(ValueError, match='complete'):
        m.merge_results(original, [{'island_key': 'a', 'status': 'COMPLETE', 'count': 99}])
