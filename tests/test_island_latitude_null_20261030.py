"""Geographic refinement must preserve every frozen recording outcome and score."""
import csv
from functools import lru_cache
import gzip
import importlib.util
import itertools
import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / 'scripts/analyze_island_latitude_null_20261030.py'


@lru_cache(None)
def module():
    assert PATH.exists(), 'Latitude-conditioned comparison is not implemented'
    spec = importlib.util.spec_from_file_location('island_latitude_null', PATH)
    result = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = result
    spec.loader.exec_module(result)
    return result


@pytest.mark.parametrize('latitude,zone', [(-90, 'south_extratropical'), (-23.500001, 'south_extratropical'),
    (-23.5, 'tropical'), (0, 'tropical'), (23.5, 'tropical'), (23.500001, 'north_extratropical'),
    (90, 'north_extratropical')])
def test_rounded_zone_boundaries_and_hemispheres(latitude, zone):
    assert module().latitude_zone(latitude) == zone


@pytest.mark.parametrize('latitude', [float('nan'), float('inf'), -90.1, 90.1])
def test_invalid_latitude_is_not_silently_classified(latitude):
    with pytest.raises(ValueError, match='latitude'):
        module().latitude_zone(latitude)


def toy():
    mask = np.ones((3, 4), dtype=bool)
    mask[2, :2] = False
    state = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1]], dtype=bool)
    values = np.arange(1, 13, dtype=float).reshape(3, 4)
    cells = [{'analysis': 'toy', 'row_index': i, 'column_index': j, 'species_key': f's{i}',
              'original_name': f'Species {i}', 'target_island': f'i{j}', 'longitude_block': 'longitude_0',
              'late_positive': int(state[i, j])} for i, j in zip(*np.where(mask))]
    return {'analysis': 'toy', 'mask': mask, 'state': state, 'values': values,
            'row_margins': state.sum(1), 'column_margins': state.sum(0),
            'row_labels': [(f's{i}', 'longitude_0') for i in range(3)],
            'column_labels': [f'i{j}' for j in range(4)], 'cells': cells,
            'observed_T': math.fsum(values[state]) / 4}, dict(zip([f'i{j}' for j in range(4)], [0, 10, -40, 40]))


def test_refinement_preserves_mask_outcomes_scores_and_original_row_and_island_margins():
    m = module()
    old, latitudes = toy()
    new = m.refine_risk(old, latitudes)
    assert m.assert_refinement(old, new)['all_frozen_cell_scores_identical']
    assert new['mask'].sum() == old['mask'].sum()
    assert new['state'].sum() == old['state'].sum()
    assert np.array_equal(new['column_margins'], old['column_margins'])
    for old_i in range(3):
        selected = [i for i, origin in enumerate(new['original_row_indices']) if origin == old_i]
        assert sum(new['row_margins'][selected]) == old['row_margins'][old_i]
        assert np.array_equal(new['mask'][selected].sum(0), old['mask'][old_i])
        assert np.array_equal(new['state'][selected].sum(0), old['state'][old_i])
    # Species 2 has separate forced south/north positives, never cross-hemisphere swaps.
    positive_labels = [new['row_labels'][i][2] for i in np.where(new['state'][:, 2:].any(1))[0]]
    assert set(positive_labels) == {'south_extratropical', 'north_extratropical'}
    corrupted = dict(new, values=new['values'].copy())
    corrupted['values'][new['mask']] += .01
    with pytest.raises(ValueError, match='score'):
        m.assert_refinement(old, corrupted)


def test_tiny_enumeration_and_forced_statistic_denominator_independently_reconcile(tmp_path):
    m = module()
    old, latitudes = toy()
    new = m.refine_risk(old, latitudes)
    blocks = m.build_blocks(new)
    total = math.prod(b['fiber'].count for b in blocks)
    assert total == 2
    expected = []
    for bits in itertools.product((0, 1), repeat=int(new['mask'].sum())):
        state = np.zeros_like(new['mask'])
        state[new['mask']] = bits
        if np.array_equal(state.sum(1), new['row_margins']) and np.array_equal(state.sum(0), new['column_margins']):
            expected.append(state)
    actual = [m.exact.table_from_global_rank(blocks, k, new['mask'].shape) for k in range(total)]
    assert {s.tobytes() for s in actual} == {s.tobytes() for s in expected}
    for state in actual:
        assert not np.any(state & ~new['mask'])
        assert list(state.sum(1)) == list(new['row_margins'])
        assert list(state.sum(0)) == list(old['column_margins'])
        m.validate_refined_table(state, new, old)
    summary, occupancy = m.summarize_fiber(new, blocks)
    scores = [math.fsum(new['values'][s]) / 4 for s in expected]
    assert summary['null_mean_T'] == pytest.approx(sum(scores) / 2)
    assert summary['null_variance_T'] == pytest.approx(np.var(scores))
    assert summary['support']['forced_positive_cells'] == 2
    assert summary['support']['variable_species'] == 2
    assert summary['support']['variable_islands'] == 2
    assert summary['support']['variable_positive_cells'] == 2
    numerators = np.sum(expected, axis=0)
    assert all(int(c['probability_numerator']) == numerators[c['row_index'], c['column_index']] for c in occupancy)
    sampled = m.exact.sample_scenario(new, blocks, 18, 25, tmp_path / 'samples.csv.gz')
    assert sampled['stored_margin_checks'] == 25
    with gzip.open(tmp_path / 'samples.csv.gz', 'rt') as handle:
        assert all(float(row['mean_log1p_source_gap_km']) in scores for row in csv.DictReader(handle))
    invalid = actual[0].copy()
    i, j = np.argwhere(invalid)[0]
    invalid[i, j] = False
    with pytest.raises(ValueError, match='margin'):
        m.validate_refined_table(invalid, new, old)


def test_real_frozen_inputs_independently_reconstruct_exactly_and_fixed_settings():
    m = module()
    _, risks, assertions = m.load_inputs()
    assert assertions['panel_islands'] == 96
    assert assertions['eligible_target_islands'] == 28
    assert all(a['all_frozen_cell_scores_identical'] for a in assertions['scenarios'].values())
    for name, risk in risks.items():
        expected = {'primary': (2385, 66), 'source_records_ge3': (2197, 62), 'late_records_ge2': (2385, 35)}[name]
        assert (int(risk['mask'].sum()), int(risk['state'].sum())) == expected
        assert sum(risk['row_margins']) == sum(risk['column_margins'])
        assert len(risk['column_labels']) == 28
    assert m.SAMPLES == 100_000
    assert [s[1] for s in m.SCENARIOS] == [2026103011, 2026103012, 2026103013]
    assert m.exact.MAX_STATES == 100_000
    assert m.exact.MAX_TRANSITIONS == 2_000_000


def test_forced_fiber_is_explicitly_no_randomization_contrast():
    m = module()
    old, latitudes = toy()
    new = m.refine_risk(old, {key: [-40, 0, 40, 41][i] for i, key in enumerate(latitudes)})
    blocks = m.build_blocks(new)
    summary, _ = m.summarize_fiber(new, blocks)
    assert summary['fiber_cardinality_decimal'] == '1'
    assert summary['randomization_contrast'] == 'NONE_FULL_FIBER_FORCED'
    assert summary['support']['variable_positive_cells'] == 0
