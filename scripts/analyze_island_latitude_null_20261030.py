#!/usr/bin/env python3
"""Exploratory fixed latitude-zone refinement of frozen recorded-source outcomes.

Exact integer finite-table law reuses the archived exact engine. No source,
recording threshold, structural cell, outcome or floating distance score changes.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ENGINE = 'scripts/analyze_island_source_exact_null_20261030.py'
spec = importlib.util.spec_from_file_location('latitude_reused_exact_engine', ROOT / ENGINE)
exact = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = exact
spec.loader.exec_module(exact)
ZONES = ('south_extratropical', 'tropical', 'north_extratropical')
SCENARIOS = (('primary', 2026103011, 1, 1), ('source_records_ge3', 2026103012, 3, 1),
             ('late_records_ge2', 2026103013, 1, 2))
SAMPLES = 100_000
PANEL = 'data/exploratory/island_expansion_20261029/expanded_island_panel_selection.csv'
SUPPORT = exact.OLD + 'recorded_source_target_support.csv'
OLD_RESULT = 'data/exploratory/island_source_exact_null_20261030/exact_null_diagnostic.json'
PINS = {**exact.PINS,
        PANEL: '05c370d6fb63c623c48714ce14b377639c446d523bc742e4f93dfbb96d562538',
        SUPPORT: '9bbf0d50f1a42b82879f23ef966a6de4a6dc36f5dc03cc92e5c511fc6bad3e74',
        ENGINE: 'ae37f7902910c86d3ec9173824920f110636bcd0989dc9c41249d40c37f349f0',
        OLD_RESULT: '9b2b21478c6a73dce547328995628f58e7d3c0f351649a037c7fee9fa74c939c'}
PLAN_TEXT = '''# Fixed latitude-zone comparison plan — 2026-10-30

This is exploratory reanalysis after the old outcomes and ratios were known.
Use exactly the frozen 96-island panel, 28 supported target islands, 93 early
butterfly concepts and all three old scenario masks, outcomes and upper-gap
scores. Only geography refines the old species-by-longitude row margins.

Zones use the archived representative latitude: south extratropical < -23.5;
tropical -23.5 through +23.5 inclusive; north extratropical > +23.5. Preserve
hemisphere distinction. These rounded tropical boundaries are not measured
climate. No alternative zones, cut points or grouping selection after results.
Each species-longitude row is split into these three zones. Preserve its sum
and the original target-island margins. Exact integer suffix DP counts use the
unchanged archived engine, capped at 100000 memo states and 2000000 attempted
transitions per longitude-zone block. Incomplete runs yield no qualified tail.
The global statistic remains mean log(1 + upper coastline-gap km), including
all forced positives and the full observed positive denominator. Save exact
cell occupancies, block counts, centered floating moments and contributions.

Draw exactly 100000 independent uniform complete-table ranks with replacement
per scenario: primary 2026103011; source_records_ge3 2026103012;
late_records_ge2 2026103013. No rerun or budget extension after results. Save
every sampled rank and statistic. Inclusive closer tail is T_null <= T_observed
+ 1e-12. Report existing individual one-sided 95% and Bonferroni simultaneous
three-scenario 95% binomial upper limits. Fully forced fibers have no
randomization contrast and cannot count as support. Validate every sample
against the full refined mask and both margins. Save plan/script/input hashes,
pre-run lock, actual times and independent assertions of frozen input identity.

Compare old and new descriptive geometric-distance ratios without claiming a
causal proportion explained. Narrower conditional support is not ecological
robustness. Known target selection, indexed-recording processes, omitted
mainland/most island sources and climate/habitat confounding remain. Earlier
nondetection is not absence; later recording is not arrival or establishment.
Static WPI2025 proximity is not historical shipping exposure. No port effect is
estimated. No old files are changed.
'''


def latitude_zone(latitude):
    latitude = float(latitude)
    if not math.isfinite(latitude) or not -90 <= latitude <= 90:
        raise ValueError('A finite latitude between -90 and 90 is required')
    return ZONES[0] if latitude < -23.5 else ZONES[2] if latitude > 23.5 else ZONES[1]


def refine_risk(old, latitudes):
    """Partition every old labeled row; each old allowed cell appears once."""
    if set(old['column_labels']) - set(latitudes):
        raise ValueError('A frozen target lacks an archived latitude')
    column_zones = [latitude_zone(latitudes[key]) for key in old['column_labels']]
    row_labels = [(r[0], r[1], zone) for r in old['row_labels'] for zone in ZONES]
    origins = [i for i in range(len(old['row_labels'])) for _ in ZONES]
    shape = (len(row_labels), len(column_zones))
    mask, state, values = np.zeros(shape, bool), np.zeros(shape, bool), np.full(shape, np.nan)
    cells = []
    for cell in old['cells']:
        i, j = int(cell['row_index']), int(cell['column_index'])
        zone = column_zones[j]
        refined_i = i * 3 + ZONES.index(zone)
        mask[refined_i, j] = old['mask'][i, j]
        state[refined_i, j] = old['state'][i, j]
        values[refined_i, j] = old['values'][i, j]
        cells.append({**cell, 'row_index': refined_i, 'original_row_index': i,
                      'latitude_zone': zone, 'representative_latitude': float(latitudes[old['column_labels'][j]])})
    risk = {**old, 'mask': mask, 'state': state, 'values': values,
            'row_margins': state.sum(1), 'column_margins': old['column_margins'].copy(),
            'row_labels': row_labels, 'original_row_indices': origins,
            'column_zones': column_zones, 'cells': cells}
    assert_refinement(old, risk)
    return risk


def assert_refinement(old, new):
    """Collapse independently; assert all original masks/outcomes/score bits."""
    mask = np.zeros_like(old['mask'], dtype=int)
    state = np.zeros_like(old['state'], dtype=int)
    values = np.full(old['values'].shape, np.nan)
    for i, original_i in enumerate(new['original_row_indices']):
        if tuple(new['row_labels'][i][:2]) != tuple(old['row_labels'][original_i]):
            raise ValueError('Refined original row identity changed')
        zone = new['row_labels'][i][2]
        for j in np.where(new['mask'][i])[0]:
            if new['column_zones'][j] != zone:
                raise ValueError('A refined cell crosses a latitude zone')
            mask[original_i, j] += 1
            state[original_i, j] += int(new['state'][i, j])
            values[original_i, j] = new['values'][i, j]
    if not np.array_equal(mask, old['mask']):
        raise ValueError('Refinement changed a frozen mask cell')
    if not np.array_equal(state, old['state']):
        raise ValueError('Refinement changed a frozen outcome')
    if not np.array_equal(values[old['mask']], old['values'][old['mask']]):
        raise ValueError('Refinement changed a frozen score')
    collapsed_margins = np.zeros(len(old['row_labels']), dtype=int)
    for i, original_i in enumerate(new['original_row_indices']):
        collapsed_margins[original_i] += new['row_margins'][i]
    if not np.array_equal(collapsed_margins, old['row_margins']):
        raise ValueError('Refined row margins fail to sum to the frozen margins')
    if not np.array_equal(new['column_margins'], old['column_margins']):
        raise ValueError('Frozen target-island margins changed')
    exact.validate_table(new['state'], new['mask'], new['row_margins'], new['column_margins'])
    observed = math.fsum(new['values'][new['state']]) / int(new['state'].sum())
    if abs(observed - old['observed_T']) > exact.TIE_TOLERANCE:
        raise ValueError('Refined observed statistic differs')
    return {'all_frozen_cell_scores_identical': True, 'all_frozen_mask_cells_identical': True,
            'all_frozen_outcomes_identical': True, 'original_row_margin_sums_identical': True,
            'target_island_margins_identical': True, 'observed_T_fsum_recomputed': observed,
            'original_mask_sha256': exact.sha(old['mask'].tobytes()),
            'original_outcome_sha256': exact.sha(old['state'].tobytes()),
            'original_allowed_score_bytes_sha256': exact.sha(old['values'][old['mask']].tobytes()),
            'refined_mask_sha256': exact.sha(new['mask'].tobytes()),
            'refined_outcome_sha256': exact.sha(new['state'].tobytes()),
            'allowed_cells': int(mask.sum()), 'positive_cells': int(state.sum())}


def validate_refined_table(state, risk, old):
    exact.validate_table(state, risk['mask'], risk['row_margins'], risk['column_margins'])
    collapsed = np.zeros_like(old['mask'], dtype=int)
    for i, original_i in enumerate(risk['original_row_indices']):
        collapsed[original_i] += state[i]
    if np.any(collapsed > 1):
        raise ValueError('Collapsed table is not binary')
    exact.validate_table(collapsed.astype(bool), old['mask'], old['row_margins'], old['column_margins'])


def load_inputs(root=ROOT):
    root = Path(root)
    exact.verify_pinned_files(root, PINS)
    old_report, originals = exact.load_inputs(root)
    with (root / PANEL).open(newline='', encoding='utf8') as handle:
        panel = list(csv.DictReader(handle))
    latitudes = {r['island_key']: float(r['representative_lat']) for r in panel}
    if len(panel) != 96 or len(latitudes) != 96:
        raise ValueError('Frozen 96-island panel differs')
    with (root / SUPPORT).open(newline='', encoding='utf8') as handle:
        support = list(csv.DictReader(handle))
    targets = {r['island_key'] for r in support if int(r['fixed_target_recording_support']) == 1}
    if len(targets) != 28:
        raise ValueError('Frozen eligible target roster differs')
    # This second read directly compares archived labeled CSV cells to the new
    # matrices, independently of the exact engine loader and collapse loop.
    with (root / exact.OLD / 'recorded_source_risk_cells.csv').open(newline='', encoding='utf8') as handle:
        direct_cells = list(csv.DictReader(handle))
    risks, assertions = {}, {}
    for name, _, _, _ in SCENARIOS:
        old = originals[name]
        if set(old['column_labels']) != targets:
            raise ValueError('Scenario target roster differs')
        new = refine_risk(old, latitudes)
        direct = [r for r in direct_cells if r['analysis'] == name]
        if len(direct) != int(new['mask'].sum()):
            raise ValueError('Direct original cell count differs')
        for cell in direct:
            j = old['column_labels'].index(cell['target_island'])
            zone_index = ZONES.index(latitude_zone(latitudes[cell['target_island']]))
            i = int(cell['row_index']) * 3 + zone_index
            if (not new['mask'][i, j] or int(new['state'][i, j]) != int(cell['late_positive'])
                    or new['values'][i, j] != float(cell['log1p_source_upper_km'])):
                raise ValueError('Independent archived-cell assertion failed')
        assertions[name] = {**assert_refinement(old, new), 'independent_direct_csv_assertions': len(direct)}
        risks[name] = new
    return old_report, risks, {'panel_islands': 96, 'eligible_target_islands': len(targets),
                              'scenarios': assertions,
                              'target_zones': [{'island_key': key, 'representative_latitude': latitudes[key],
                                                'latitude_zone': latitude_zone(latitudes[key])} for key in sorted(targets)]}


def build_blocks(risk):
    blocks, covered = [], np.zeros_like(risk['mask'])
    for longitude in sorted({r[1] for r in risk['row_labels']}):
        for zone in ZONES:
            rows = [i for i, r in enumerate(risk['row_labels']) if r[1:] == (longitude, zone)
                    or tuple(r[1:]) == (longitude, zone)]
            columns = [j for j, z in enumerate(risk['column_zones']) if z == zone
                       and any(risk['mask'][i, j] for i in rows)]
            selector = np.ix_(rows, columns)
            mask = risk['mask'][selector]
            covered[selector] = mask
            fiber = exact.ExactBlockFiber(mask, risk['row_margins'][rows], risk['column_margins'][columns],
                                         max_states=100_000, max_transitions=2_000_000)
            if not fiber.count:
                raise ValueError('Observed feasible refined table has zero completion count')
            blocks.append({'longitude_block': longitude + ':' + zone, 'latitude_zone': zone,
                           'original_longitude_block': longitude, 'rows': rows, 'columns': columns, 'fiber': fiber})
    if not np.array_equal(covered, risk['mask']):
        raise ValueError('Longitude-zone blocks do not partition the frozen mask')
    return blocks


def summarize_fiber(risk, blocks):
    summary, occupancy = exact.summarize_fiber(risk, blocks)
    for record, block in zip(summary['blocks'], blocks):
        record.update(latitude_zone=block['latitude_zone'], original_longitude_block=block['original_longitude_block'])
    for record, cell in zip(occupancy, risk['cells']):
        record.update(latitude_zone=cell['latitude_zone'], original_row_index=cell['original_row_index'],
                      representative_latitude=cell['representative_latitude'])
    variable = [r for r in occupancy if r['occupancy_status'] == 'variable']
    forced = [r for r in occupancy if r['occupancy_status'] == 'forced_positive']
    summary['support'].update(variable_species=len({r['species_key'] for r in variable}),
                              variable_islands=len({r['target_island'] for r in variable}),
                              forced_positive_species=len({r['species_key'] for r in forced}),
                              forced_positive_islands=len({r['target_island'] for r in forced}),
                              refined_rows=len(risk['row_labels']))
    summary['randomization_contrast'] = 'NONE_FULL_FIBER_FORCED' if summary['fiber_cardinality_decimal'] == '1' else 'VARIABLE_CONDITIONAL_FIBER'
    return summary, occupancy


def run(output, root=ROOT):
    output, root = Path(output), Path(root)
    exact.verify_pinned_files(root, PINS)
    output.mkdir(parents=True, exist_ok=False)
    plan_path = output / 'latitude_null_fixed_plan.md'
    plan_path.write_text(PLAN_TEXT, encoding='utf8')
    plan_sha = exact.sha(plan_path.read_bytes())
    script_sha = exact.sha(Path(__file__).read_bytes())
    lock = {'schema': 'chocho_latitude_null_first_run_lock_20261030_v1', 'started_utc': exact.utc(),
            'plan_sha256': plan_sha, 'script_sha256': script_sha, 'input_sha256': PINS,
            'samples_per_scenario': SAMPLES, 'sampling_with_replacement': True,
            'scenario_order_and_seeds': [{'analysis': n, 'seed': seed} for n, seed, _, _ in SCENARIOS],
            'zones': {'south_extratropical': 'latitude < -23.5', 'tropical': '-23.5 <= latitude <= 23.5',
                      'north_extratropical': 'latitude > 23.5'},
            'max_memo_states_per_block': 100_000, 'max_attempted_pattern_transitions_per_block': 2_000_000,
            'python_version': platform.python_version(), 'numpy_version': np.__version__,
            'tail_tie_tolerance': exact.TIE_TOLERANCE, 'no_reruns_or_budget_extensions': True}
    lock_path = output / 'latitude_null_run_lock.json'
    lock_path.write_text(json.dumps(lock, indent=2, sort_keys=True, allow_nan=False) + '\n')
    lock_sha, started = exact.sha(lock_path.read_bytes()), time.perf_counter()
    analyses, occupancy_rows, block_rows, stage = {}, [], [], 'load_inputs'
    try:
        _, risks, assertions = load_inputs(root)
        (output / 'independent_input_assertions.json').write_text(json.dumps(assertions, indent=2, sort_keys=True) + '\n')
        old_exact = json.loads((root / OLD_RESULT).read_text())
        for name, seed, source_min, late_min in SCENARIOS:
            stage = name + ':exact_counts_and_moments'
            counted_at, count_started = exact.utc(), time.perf_counter()
            risk = risks[name]
            blocks = build_blocks(risk)
            summary, cells = summarize_fiber(risk, blocks)
            summary.update(source_min_records=source_min, late_min_records=late_min,
                           counting_started_utc=counted_at, counting_elapsed_seconds=time.perf_counter() - count_started,
                           frozen_input_assertions=assertions['scenarios'][name],
                           old_longitude_only_ratio=old_exact['analyses'][name]['observed_to_null_geometric_mean_one_plus_gap_ratio'])
            analyses[name] = summary
            occupancy_rows.extend(cells)
            block_rows.extend({'analysis': name, **{k: v for k, v in b.items() if not isinstance(v, list)}} for b in summary['blocks'])
            print(json.dumps({'stage': 'exact_counts_complete', 'analysis': name,
                              'fiber_cardinality_decimal': summary['fiber_cardinality_decimal'],
                              'block_counts': [b['fiber_cardinality_decimal'] for b in summary['blocks']],
                              'ratio': summary['observed_to_null_geometric_mean_one_plus_gap_ratio']}), flush=True)
            stage = name + ':fixed_independent_draws'
            summary['sampling_started_utc'] = exact.utc()
            summary['sampling'] = exact.sample_scenario(risk, blocks, seed, SAMPLES, output / f'{name}_rank_T_samples.csv.gz')
            summary['sampling_finished_utc'] = exact.utc()
            print(json.dumps({'stage': 'sampling_complete', 'analysis': name, **summary['sampling']}), flush=True)
        stage = 'final_checks_and_outputs'
        exact.verify_pinned_files(root, PINS)
        if (exact.sha(Path(__file__).read_bytes()) != script_sha or exact.sha(lock_path.read_bytes()) != lock_sha
                or exact.sha(plan_path.read_bytes()) != plan_sha):
            raise ValueError('Script, plan or first-run lock changed during execution')
        if any(a['sampling']['stored_samples'] != SAMPLES or a['sampling']['stored_margin_checks'] != SAMPLES for a in analyses.values()):
            raise ValueError('A fixed sample budget or validation was incomplete')
        exact.write_gzip_csv(output / 'latitude_exact_cell_occupancy_counts.csv.gz', occupancy_rows, list(occupancy_rows[0]))
        exact.write_gzip_csv(output / 'latitude_exact_block_counts_and_moments.csv.gz', block_rows, list(block_rows[0]))
        for summary in analyses.values():
            forced = summary['randomization_contrast'] == 'NONE_FULL_FIBER_FORCED'
            summary['tail'] = exact.tail_summary(summary['sampling']['lower_tail_count_including_ties'], SAMPLES)
            summary['tail']['qualified_bound'] = not forced
            if forced:
                summary['tail']['interpretation'] = 'No randomization contrast; full fiber forced. This cannot count as support.'
            summary['sampling']['qualified_bound'] = not forced
            summary['sampling']['sample_mean_minus_DP_mean_T'] = summary['sampling']['sample_mean_T'] - summary['null_mean_T']
            summary['sampling']['sample_mean_MC_standard_error_from_DP'] = math.sqrt(summary['null_variance_T'] / SAMPLES)
        report = {'schema': 'chocho_latitude_exact_null_20261030_v1',
                  'status': 'COMPLETE_EXACT_COUNT_IID_MONTE_CARLO_BOUNDS',
                  'started_utc': lock['started_utc'], 'finished_utc': exact.utc(),
                  'elapsed_seconds': time.perf_counter() - started,
                  'plan_sha256': plan_sha, 'run_lock_sha256': lock_sha, 'script_sha256': script_sha, 'input_sha256': PINS,
                  'method': {'grouping': 'Species by original longitude block by fixed archived latitude zone; unchanged target-island margins',
                             'counts': 'Reused exact Python-integer suffix DP and product-rank bijection',
                             'scores': 'Unchanged floating log1p upper-gap scores; all forced positives and original positive denominator',
                             'moments': 'Centered floating mixture moments independently checked against exact integer occupancies',
                             'sampling': '100000 independent uniform arbitrary-integer product ranks with replacement per scenario; every complete table validated',
                             'tail': 'Inclusive T_null <= frozen observed T + 1e-12; one-sided exact binomial and three-scenario Bonferroni upper limits',
                             'caps_per_block': {'memo_states': 100_000, 'attempted_transitions': 2_000_000}},
                  'interpretation': ['Exploratory ecological reanalysis after the old outcomes, ratios and concentration were known; not preregistered or independent replication.',
                                     'Rounded tropical boundaries use representative geography and are not measured climate. Hemisphere is retained.',
                                     'A smaller conditional fiber is less randomization support, not evidence of ecological robustness.',
                                     'A ratio change is not a causal proportion explained by latitude.',
                                     'Target selection depends on known later butterfly recording richness. Recording effort, habitat/climate, mainland and most island sources remain omitted.',
                                     'Early nondetection is not absence; later recording is not first arrival, breeding or establishment. No historical port effect is estimated.'],
                  'independent_input_assertions': assertions, 'analyses': analyses,
                  'output_sha256': {p.name: exact.sha(p.read_bytes()) for p in sorted(output.iterdir()) if p.is_file()}}
        (output / 'latitude_null_diagnostic.json').write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')
        return report
    except BaseException as error:
        failure = {'schema': 'chocho_latitude_exact_null_20261030_v1', 'status': 'INCOMPLETE_NO_QUALIFIED_BOUND',
                   'qualified_bounds': False, 'failed_stage': stage, 'error_type': type(error).__name__, 'error': str(error),
                   'started_utc': lock['started_utc'], 'finished_utc': exact.utc(), 'elapsed_seconds': time.perf_counter() - started,
                   'plan_sha256': plan_sha, 'run_lock_sha256': lock_sha, 'script_sha256': script_sha, 'input_sha256': PINS,
                   'partial_output_sha256': {p.name: exact.sha(p.read_bytes()) for p in sorted(output.iterdir()) if p.is_file()}}
        (output / 'latitude_null_diagnostic.json').write_text(json.dumps(failure, indent=2, sort_keys=True, allow_nan=False) + '\n')
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'data/exploratory/island_latitude_null_20261030')
    run(parser.parse_args().output_dir)
