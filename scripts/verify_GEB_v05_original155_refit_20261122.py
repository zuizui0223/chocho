#!/usr/bin/env python3
"""Compare the original fixed-155 fit to a fresh execution. Export aggregate evidence only."""
import argparse
from datetime import datetime
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


COMMIT = '42f8fcc3a388e96456c76180d7e56d75dd2ed48c'
LOCK = '0c7fc34d7858d6f62fe31a914188aee63e56ab1d41872abe40c314f6c5ed8b1e'
PANEL = '5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1'
FIT = 'a34ad4cd1e4f40460c4f6bf6ed0b22ba3885b17d20ac17f74bd6564d0618ce6d'
FROZEN = {
    'temporal': 'f5303e00c52a9d5dea82f990a8e86bff16fea041d2a1c827729cf102f9c83ac2',
    'spatial': '76500f3fd1e059750fd56546718bbd36aaaffe90786c56d76047b532b6ca7ad4',
    'species': '7896720594ab99bbc47437f3222fedce1b2a02a0098f1b1de2a3fd2ef23c1792',
    'spatial_buffer100': '8a157515a85303435ac5c1c86229035f345db67c4ad22614ea49de58da8b7f1e',
}
MODELS = ('persistence', 'M0_background', 'M1_climate', 'M2_landcover',
          'M3_both', 'M4_host_resources', 'M5_host_composition')
SETTINGS = dict(n_estimators=160, max_depth=12, min_samples_leaf=40,
                max_features=.7, random_state=20261029, n_jobs=4)
RTOL, ATOL = 1e-6, 1e-9
KEYS = ['speciesKey', 'cell']
COUNTS = {'rows', 'species', 'initially_unrecorded_rows'}
BUFFER_METRICS = {'rows', 'species', 'species_equal_mse', 'species_equal_mae'}
METRICS = BUFFER_METRICS | {'pooled_mse', 'initially_unrecorded_rows', 'initially_unrecorded_mse'}


class Invalid(Exception):
    """Only pass constant labels or labels formed entirely from trusted constants."""


def require(ok, label):
    if not ok:
        raise Invalid(label)


class Audit:
    def __init__(self):
        self.comparisons = {}

    def compare(self, label, expected, actual, exact=False):
        left = np.asarray(expected, dtype=np.float64)
        right = np.asarray(actual, dtype=np.float64)
        require(left.shape == right.shape and left.size > 0, 'numeric_shape_mismatch')
        finite = np.isfinite(left) & np.isfinite(right)
        with np.errstate(over='ignore'):
            delta = np.abs(right[finite] - left[finite])
        maximum = float(delta.max()) if delta.size else None
        tolerance = 0 if exact else ATOL + RTOL * np.abs(left[finite])
        self.comparisons[label] = {
            'compared': int(left.size),
            'outside_tolerance': int((~finite).sum() + (delta > tolerance).sum()),
            'nonfinite': int((~finite).sum()),
            'max_absolute_difference': maximum if maximum is not None and math.isfinite(maximum) else None,
        }

    @property
    def passed(self):
        return bool(self.comparisons) and all(not c['outside_tolerance'] for c in self.comparisons.values())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def header(result, mode):
    for key, required in dict(status='COMPLETE', mode=mode, train_rows=193703, test_rows=208482,
                              train_species=155, test_species=155, model_lock_sha256=LOCK,
                              settings=SETTINGS).items():
        require(result.get(key) == required, 'inconsistent_' + key)
    require(result.get('sensitivity') is None and result.get('rarefaction_n', 10) == 10,
            'incorrect_primary_fit')


def keyed(frame):
    require(set(KEYS + ['species']).issubset(frame.columns), 'missing_keys')
    require(not frame[KEYS + ['species']].isna().any().any(), 'null_key')
    require(not frame.duplicated(KEYS).any(), 'duplicate_keys')
    require(not frame.duplicated(['species','cell']).any(), 'duplicate_species_cells')
    return frame.set_index(KEYS).sort_index()


def compare_predictions(audit, original, new, mode):
    old, fresh = keyed(original), keyed(new)
    require(old.index.equals(fresh.index), 'prediction_keys_differ')
    require((old.species.to_numpy() == fresh.species.to_numpy()).all(), 'species_labels_differ')
    columns = ['target', 'baseline_probability', *MODELS]
    if mode == 'spatial':
        columns.append('distance_to_training_cell_km')
    for column in columns:
        require(column in old and column in fresh, 'missing_prediction_column_' + column)
        audit.compare('row.' + column, old[column], fresh[column])
    audit.compare('baseline_persistence', fresh.baseline_probability, fresh.persistence)
    for column in ['species_fold', 'spatial_fold', 'region', 'spatial_block', 'climate_id', 'days', 'future_days']:
        if column in old:
            require(column in fresh and not old[column].isna().any() and
                    not fresh[column].isna().any() and
                    (old[column].to_numpy() == fresh[column].to_numpy()).all(),
                    'row_metadata_differ_' + column)
    if mode == 'spatial':
        audit.compare('row.buffer100_membership', old.distance_to_training_cell_km >= 100,
                      fresh.distance_to_training_cell_km >= 100, exact=True)


def scores(frame, model, buffered=False):
    err = frame[model] - frame.target
    errors = pd.DataFrame({'mse': err ** 2, 'mae': abs(err), 'species': frame.species})
    groups = errors.groupby('species')[['mse', 'mae']].mean()
    groups['rows'] = errors.groupby('species').size()
    result = dict(rows=len(frame), species=len(groups), species_equal_mse=groups.mse.mean(),
                  species_equal_mae=groups.mae.mean())
    if not buffered:
        unseen = errors[frame.baseline_probability == 0]
        result.update(pooled_mse=errors.mse.mean(), initially_unrecorded_rows=len(unseen),
                      initially_unrecorded_mse=unseen.groupby('species').mse.mean().mean())
    return result, groups


def compare_species(audit, label, expected, actual):
    require(expected.index.is_unique and actual.index.is_unique and
            not expected.index.hasnans and not actual.index.hasnans, 'duplicate_or_null_species')
    left, right = expected.sort_index(), actual.sort_index()
    require(left.index.equals(right.index), 'species_metric_keys_differ')
    require(set(right.columns) == set(left.columns) == {'mse','mae','rows'}, 'species_metric_columns')
    for key in ['mse','mae','rows']:
        audit.compare(label+'.'+key, left[key], right[key], exact=(key == 'rows'))


def metrics_compare(audit, label, expected, actual, buffered=False):
    require(set(expected) == set(actual) == (BUFFER_METRICS if buffered else METRICS),
            'score_fields_differ')
    for key in sorted(expected):
        require(np.isscalar(expected[key]) and np.isscalar(actual[key]), 'nonscalar_score')
        name = label + '.' + key
        audit.compare(name, expected[key], actual[key], exact=(key in COUNTS))
        if not audit.comparisons[name]['nonfinite']:
            audit.comparisons[name].update(frozen_value=float(expected[key]), new_value=float(actual[key]))


def load_frozen(snapshot, mode):
    file = snapshot / f'forecast_{mode}_results.json'
    require(digest(file) == FROZEN[mode], 'snapshot_hash_mismatch')
    result = read_json(file)
    require(set(result['results']) == set(MODELS), 'original_score_groups_differ')
    if mode == 'spatial_buffer100':
        require(result['status']=='COMPLETE' and result['mode']==mode and result['cells']==1788,
                'buffer_snapshot_header')
    else:
        header(result, mode)
    return result['results']


def timestamp(value):
    if type(value) in (int, float):
        require(math.isfinite(value), 'invalid_execution_timestamp')
        return float(value)
    require(isinstance(value, str), 'invalid_execution_timestamp')
    try:
        parsed = datetime.fromisoformat(value.replace('Z', '+00:00'))
        require(parsed.utcoffset() is not None, 'ambiguous_execution_timestamp')
        return parsed.timestamp()
    except (ValueError, OverflowError):
        raise Invalid('invalid_execution_timestamp') from None


def execution(receipt, mode):
    run_id = receipt.get('run_id')
    require(type(run_id) in (str, int) and re.fullmatch(r'[1-9][0-9]{0,19}', str(run_id)),
            'invalid_workflow_run')
    for key, expected in dict(mode=mode, source_script_sha256=FIT, panel_sha256=PANEL,
                              returncode=0, preexisting_outputs=False).items():
        require(type(receipt.get(key)) is type(expected) and receipt[key] == expected,
                'invalid_execution_' + key)
    start, end = [timestamp(receipt.get(key)) for key in ('started_at', 'completed_at')]
    require(0 < start <= end < time.time() + 2, 'inconsistent_execution_time')
    return str(run_id), start, end


def fresh_file(path, start, end):
    require(path.is_file() and start <= path.stat().st_mtime < end + 1,
            'missing_or_stale_fit_file')
    return path


def validate(args, audit, report):
    mode, run = args.mode, args.fit_root
    outputs = run / 'outputs'
    require(not args.original_predictions.resolve().is_relative_to(run.resolve()) and
            not args.snapshot_dir.resolve().is_relative_to(run.resolve()), 'fresh_outputs_must_be_isolated')
    run_id, start, end = execution(read_json(args.execution_receipt), mode)
    report['workflow_run_id'] = run_id
    require(digest(run/'scripts/temporal/fit_forecast.py') == FIT and
            digest(outputs/'forecast_model_lock.json') == LOCK and
            digest(run/'work/temporal/forecast/species_cell_period_panel.parquet') == PANEL,
            'fresh_fit_script_panel_or_lock_hash_mismatch')
    input_hash = digest(args.original_predictions)
    require(input_hash == args.original_predictions_sha256, 'archived_prediction_hash_mismatch')
    oldscores = load_frozen(args.snapshot_dir, mode)
    newscore_path = fresh_file(outputs/f'forecast_{mode}_results.json', start, end)
    newscore = read_json(newscore_path)
    header(newscore, mode)
    require(set(newscore['results']) == set(MODELS), 'new_score_groups_differ')
    newpath = fresh_file(run/f'work/temporal/forecast/predictions_{mode}.parquet',start,end)
    old, new = pd.read_parquet(args.original_predictions), pd.read_parquet(newpath)
    require(len(new) == 208482 and new.species.nunique() == 155, 'new_test_support_differ')
    compare_predictions(audit, old, new, mode)
    require(all(c['nonfinite']==0 for c in audit.comparisons.values()), 'nonfinite_predictions')
    report['sha256'] = dict(original_predictions=input_hash, new_predictions=digest(newpath),
                            new_scores=digest(newscore_path), frozen_scores=FROZEN[mode],
                            model_lock=LOCK, fit_script=FIT, panel=PANEL)
    report['support'] = dict(train_rows=193703, train_species=155, test_rows=len(new), test_species=155)
    report['settings'] = SETTINGS
    report['direct_archived_species_csv_count'] = 0
    for model in MODELS:
        newmetrics, newgroups = scores(new, model)
        metrics_compare(audit, 'fit_output.'+model, oldscores[model], newscore['results'][model])
        metrics_compare(audit, 'fresh_rows.'+model, oldscores[model], newmetrics)
        if model == 'persistence':
            continue
        name = f'forecast_{mode}_{model}_species_metrics.csv'
        _, oldgroups = scores(old, model)
        freshcsv = pd.read_csv(fresh_file(outputs/name, start, end), index_col='species')
        compare_species(audit, 'by_species.rows.'+model, oldgroups, newgroups)
        compare_species(audit, 'by_species.csv.'+model, oldgroups, freshcsv)
        original_csv_dir = args.original_species_csv_dir or args.snapshot_dir
        if (original_csv_dir/name).is_file():
            require(not (original_csv_dir/name).resolve().is_relative_to(run.resolve()),
                    'original_species_csv_must_be_separate')
            archived_csv = pd.read_csv(original_csv_dir/name,index_col='species')
            compare_species(audit, 'by_species.archived_csv.'+model, archived_csv, freshcsv)
            report['direct_archived_species_csv_count'] += 1
    if mode == 'spatial':
        subset = new[new.distance_to_training_cell_km >= 100]
        require((len(subset), subset.cell.nunique(), subset.species.nunique()) == (102728,1788,155),
                'buffer_support_differ')
        report['buffer100_support'] = dict(rows=102728, cells=1788, species=155)
        buffer = load_frozen(args.snapshot_dir,'spatial_buffer100')
        report['sha256']['frozen_buffer100_scores'] = FROZEN['spatial_buffer100']
        for model in MODELS:
            observed, _ = scores(subset,model,buffered=True)
            metrics_compare(audit,'buffer100.'+model,buffer[model],observed,buffered=True)
    report['comparison_scope'] = ('Fresh fit execution attested by calling workflow. Validator checks output files, '
                                  'including original keyed predictions and all frozen score groups. '
                                  'The 100-km subset reuses the spatial fit.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['temporal','species','spatial'], required=True)
    for key in ['snapshot-dir','fit-root','original-predictions','output']:
        parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--original-predictions-sha256',required=True,
                        help='Hash from the independently verified original recovery inventory')
    parser.add_argument('--original-species-csv-dir', type=Path)
    parser.add_argument('--execution-receipt', required=True, type=Path,
                        help='Receipt from the calling workflow for its freshly executed fit')
    args = parser.parse_args()
    audit=Audit()
    report = dict(schema='GEB_original155_refit_validation_v1', mode=args.mode,
                  status='FAIL', original_commit=COMMIT, model_refit_completed=False,
                  absolute_tolerance=ATOL, relative_tolerance=RTOL)
    try:
        validate(args,audit,report)
        require(audit.passed,'numerical_comparison_failed')
        report.update(status='PASS', model_refit_completed=True)
    except Invalid as error:
        report['failure']=str(error)
    except Exception as error:
        # Exceptions from pandas/Arrow/filesystem parsing may contain private paths, labels or records.
        report['failure']='input_read_or_validation_error_'+type(error).__name__
    report['runtime_versions']={'python':platform.python_version()}
    for package in ['pandas','numpy','scikit-learn','pyarrow','joblib','scipy']:
        try:
            report['runtime_versions'][package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            report['runtime_versions'][package]='unavailable'
    report['comparisons']=audit.comparisons
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('GEB_ORIGINAL155_VALIDATION_'+report['status'])
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':
    sys.exit(main())
