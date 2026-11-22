#!/usr/bin/env python3
"""Controlled M0 species-fold probe, exporting hashes, counts and fixed metadata."""
import argparse
import contextlib
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

import diagnose_GEB_v05_original155_species_20261122 as lineage
import recover_GEB_v05_original155_private_inputs_20261122 as recovery
import verify_GEB_v05_original155_refit_20261122 as check

ARMS = {'baseline': ('2.2.3','23.0.1'), 'pandas_only': ('2.3.3','23.0.1'),
        'arrow_only': ('2.2.3','24.0.0'), 'joint': ('2.3.3','24.0.0')}
COMMON = {'numpy':'2.3.5','scikit-learn':'1.7.2','scipy':'1.18.1','joblib':'1.5.3',
          'threadpoolctl':'3.6.0','python-dateutil':'2.9.0.post0','six':'1.17.0',
          'tzdata':'2025.2','pytz':'2025.2'}
PYTHON = '3.12.16'
ARCHIVED_PREDICTIONS = '85f319512d638f491656c3c0aca2289aca1298b18397254b071a67b9e65ef038'
REFERENCE_COUNT, REFERENCE_MAX = 28575, 0.007789458622820877
THREADS = {'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OMP_NUM_THREADS':'1',
           'PYTHONHASHSEED':'0'}


class ProbeError(Exception):
    """Fixed diagnostic messages suitable for a public runner."""


class FiveFoldsComplete(Exception):
    pass


def require(ok, label):
    if not ok:
        raise ProbeError(label)


def array_signature(array):
    a = np.asarray(array)
    require(a.dtype.kind in 'fibu' and np.isfinite(a).all(),'invalid_numeric_array')
    return dict(sha256=hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest(),
                dtype=a.dtype.str, shape=list(a.shape), strides=list(a.strides))


def key_hash(frame):
    h = hashlib.sha256()
    for key in check.KEYS:
        for value in frame[key].to_numpy():
            if isinstance(value, (str, np.str_)):
                kind = b's'
            elif isinstance(value, (int, np.integer)):
                kind = b'i'
            else:
                raise ProbeError('unsupported_key_type')
            data = str(value).encode('utf-8')
            h.update(kind + len(data).to_bytes(4,'big') + data)
    return h.hexdigest()


def instrument(raw):
    require(hashlib.sha256(raw).hexdigest()==recovery.SCRIPT_SHA256,'incorrect_source_digest')
    text = raw.decode('utf-8')
    replacements = {
        'weights*=len(weights)/weights.sum()':
            'weights*=probe.normalization(len(weights)/weights.sum(), weights)',
        'reg=RandomForestRegressor(**settings)': 'reg=probe.factory(name, fold, tr, te, settings, features)',
    }
    for source, replacement in replacements.items():
        require(text.count(source)==1,'ambiguous_probe_boundary')
        text = text.replace(source,replacement)
    return text


def forest_hash(model):
    h = hashlib.sha256()
    for est in model.estimators_:
        h.update(int(est.random_state).to_bytes(8,'big'))
        for field in ('children_left','children_right','feature','threshold','value',
                      'n_node_samples','weighted_n_node_samples','impurity','missing_go_to_left'):
            array = getattr(est.tree_,field)
            header = json.dumps([field,array.dtype.str,list(array.shape),array.nbytes],
                                separators=(',',':')).encode('ascii')
            h.update(len(header).to_bytes(4,'big'))
            h.update(header)
            h.update(np.ascontiguousarray(array).tobytes())
    return h.hexdigest()


class Observer:
    def __init__(self, archived, fit):
        self.archived = archived
        self.fit_models = fit
        self.folds = []
        self.positions = None
        self.coverage = np.zeros(len(archived),dtype=bool)
        self.predictions = np.full(len(archived),np.nan) if fit else None

    def normalization(self, scalar, raw=None):
        self.weights = dict(normalization_scalar_hex=float(scalar).hex())
        if raw is not None:
            self.weights.update(raw=array_signature(raw), raw_sum_hex=float(raw.sum()).hex())
        return scalar

    def factory(self, name, fold, train, test, settings, features):
        require(name=='M0_background' and settings==recovery.MODEL_SETTINGS,'unexpected_estimator')
        positions = self.archived.index.get_indexer(pd.MultiIndex.from_frame(test[check.KEYS]))
        require((positions>=0).all() and len(np.unique(positions))==len(positions) and
                not self.coverage[positions].any(), 'test_key_mismatch')
        require(np.all(self.archived.species_fold.to_numpy()[positions]==fold),'test_fold_mismatch')
        self.positions = positions
        self.coverage[positions] = True
        self.record = dict(fold=int(fold), feature_order_sha256=hashlib.sha256(
            json.dumps(features,separators=(',',':')).encode()).hexdigest(), train_key_order=key_hash(train),
                           test_key_order=key_hash(test), arrays={})
        if self.fit_models:
            from sklearn.ensemble import RandomForestRegressor
            self.model = RandomForestRegressor(**settings)
        return self

    def fit(self, x, y, sample_weight):
        self.record['weights'] = self.weights
        self.record['arrays'].update(train_features=array_signature(x), target=array_signature(y),
                                     normalized_weights=array_signature(sample_weight))
        if self.fit_models:
            self.model.fit(x,y,sample_weight=sample_weight)
        return self

    def predict(self, x):
        self.record['arrays']['test_features'] = array_signature(x)
        if self.fit_models:
            pred = self.model.predict(x)
            self.predictions[self.positions] = pred
            audit = check.Audit()
            audit.compare('m0',self.archived.M0_background.to_numpy()[self.positions],pred)
            self.record.update(comparison=audit.comparisons['m0'], forest_sha256=forest_hash(self.model),
                               predictions=array_signature(pred))
            del self.model
        else:
            # A placeholder satisfies the source's fold assignment. No score/output is produced.
            pred = np.zeros(len(x))
        self.folds.append(self.record)
        if len(self.folds)==5:
            require([r['fold'] for r in self.folds]==list(range(5)) and self.coverage.all(),
                    'incomplete_fold_coverage')
            raise FiveFoldsComplete()
        return pred


def equivalent_arms(signatures):
    require(set(signatures)==set(ARMS),'missing_preflight_arm')
    equivalent, previous = {}, {}
    for arm in ARMS:
        arrays = [(f['fold'],f['arrays']) for f in signatures[arm]]
        match = next((candidate for candidate, prior in previous.items() if prior==arrays),None)
        equivalent[arm] = match or arm
        if match is None:
            previous[arm] = arrays
    return equivalent


def select_arms(signatures):
    return [arm for arm,representative in equivalent_arms(signatures).items()
            if arm==representative]


def packages(arm):
    pandas_version,arrow_version=ARMS[arm]
    return dict(COMMON,pandas=pandas_version,pyarrow=arrow_version)


def check_environment(arm):
    versions = {name:importlib.metadata.version(name) for name in packages(arm)}
    require(versions==packages(arm) and platform.python_version()==PYTHON,'runtime_pin_mismatch')
    require(all(os.environ.get(k)==v for k,v in THREADS.items()),'thread_setting_mismatch')
    require(pd.get_option('io.parquet.engine')=='auto','unexpected_parquet_engine_option')
    from pandas.io.parquet import get_engine
    require(type(get_engine('auto')).__name__=='PyArrowImpl','unexpected_parquet_reader')
    return versions


def run_arm(root, snapshot, arm, fit):
    versions = check_environment(arm)
    panel, preds = lineage.inputs(root)
    require(recovery.sha256_file(preds)==ARCHIVED_PREDICTIONS,'species_prediction_identity')
    source = snapshot/'analysis_source/fit_forecast.py'
    code = instrument(source.read_bytes())
    lock = snapshot/'forecast_model_lock.json'
    require(recovery.sha256_file(lock)==recovery.LOCK_SHA256,'incorrect_model_lock')
    stage = 'fit' if fit else 'preflight'
    run = root/'probe'/arm/stage
    panel_dest = run/recovery.PANEL
    panel_dest.parent.mkdir(parents=True,exist_ok=False)
    os.link(panel,panel_dest)
    (run/'outputs').mkdir()
    shutil.copyfile(lock,run/recovery.LOCK)
    script = run/'scripts/temporal/fit_forecast.py'
    script.parent.mkdir(parents=True)
    shutil.copyfile(source,script)
    old = pd.read_parquet(preds,columns=[*check.KEYS,'species_fold','M0_background',
                                         'target','baseline_probability']).reset_index(drop=True)
    require(len(old)==208482 and not old.duplicated(check.KEYS).any() and
            not old.isna().any().any(),'incorrect_archived_support')
    observer = Observer(old.set_index(check.KEYS),fit)
    sys.argv = [str(script),'species']
    with (run/'source.log').open('w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            exec(compile(code,str(script),'exec'),dict(__file__=str(script),__name__='__main__',probe=observer))
        except FiveFoldsComplete:
            pass
        else:
            raise ProbeError('model_probe_did_not_stop_after_m0')
    ordered = old.sort_values(check.KEYS).reset_index(drop=True)
    oracle = dict(keys=key_hash(ordered), values={k:array_signature(ordered[k].to_numpy())
                  for k in ['species_fold','M0_background','target','baseline_probability']})
    result = dict(arm=arm, stage=stage, versions=versions, folds=observer.folds,
                  ordered_archived_oracle=oracle, instrumented_source_sha256=hashlib.sha256(
                      code.encode('utf-8')).hexdigest())
    if fit:
        require(np.isfinite(observer.predictions).all(),'incomplete_prediction_array')
        # Each environment reads the same archive. Use the same canonical key order for cross-arm comparisons.
        order = old.sort_values(check.KEYS).index.to_numpy()
        pred = observer.predictions[order]
        truth = old.M0_background.to_numpy()[order]
        audit = check.Audit()
        audit.compare('vs_archived',truth,pred)
        if arm=='baseline':
            np.save(run/'baseline.npy',pred,allow_pickle=False)
        else:
            baseline = np.load(root/'probe/baseline/fit/baseline.npy',allow_pickle=False)
            audit.compare('vs_baseline',baseline,pred)
        result['comparisons'] = audit.comparisons
        result['canonical_test_keys_sha256'] = key_hash(old.iloc[order])
    recovery.write_json(run/'aggregate.json',result)


def orchestrate(root, snapshot, output):
    work = root/'probe'
    work.mkdir()
    report = dict(schema='geb_species_controlled_m0_runtime_probe_v1',status='IN_PROGRESS',
                  source_sha256=recovery.SCRIPT_SHA256, panel_sha256=recovery.PANEL_SHA256,
                  archived_predictions_sha256=ARCHIVED_PREDICTIONS, lock_sha256=recovery.LOCK_SHA256,
                  python=PYTHON, architecture=platform.machine(), os=platform.system(), threads=THREADS,
                  settings=recovery.MODEL_SETTINGS, absolute_tolerance=check.ATOL,relative_tolerance=check.RTOL,
                  previous_replay=dict(run_id='44970180129',outside_tolerance=REFERENCE_COUNT,
                                       max_absolute_difference=REFERENCE_MAX), preflight={}, fits={})
    recovery.write_json(output,report)
    def worker(arm, stage):
        executable = work/'env'/arm/'bin/python'
        result_path = work/arm/stage/'aggregate.json'
        subprocess.run([str(executable),str(Path(__file__).resolve()),stage,'--arm',arm,'--root',str(root),
                        '--snapshot',str(snapshot)],check=True)
        return json.loads(result_path.read_text())
    try:
        for arm in ARMS:
            env = work/'env'/arm
            subprocess.run([sys.executable,'-m','venv',str(env)],check=True)
            deps = [f'{k}=={v}' for k,v in packages(arm).items()]
            subprocess.run([str(env/'bin/python'),'-m','pip','install','--disable-pip-version-check',
                            '--only-binary=:all:','--no-deps',*deps],check=True)
            report['preflight'][arm] = worker(arm,'preflight')
            recovery.write_json(output,report)
        oracle = report['preflight']['baseline']['ordered_archived_oracle']
        for r in report['preflight'].values():
            current = r['ordered_archived_oracle']
            require(oracle['keys']==current['keys'] and
                    all(oracle['values'][k]['sha256']==v['sha256'] and
                        oracle['values'][k]['dtype']==v['dtype'] for k,v in current['values'].items()),
                    'archived_oracle_differs_by_reader')
        signatures = {arm:r['folds'] for arm,r in report['preflight'].items()}
        selected = select_arms(signatures)
        report['equivalent_estimator_inputs'] = equivalent_arms(signatures)
        report['selected_fit_arms'] = selected
        for arm in selected:
            r = worker(arm,'fit')
            require(r['ordered_archived_oracle']==report['preflight'][arm]['ordered_archived_oracle'],
                    'real_fit_archived_oracle_differs_from_preflight')
            for pre,real in zip(signatures[arm],r['folds'],strict=True):
                require(all(real[k]==pre[k] for k in pre),'real_fit_inputs_differ_from_preflight')
            if arm=='baseline':
                actual = r['comparisons']['vs_archived']
                r['previous_replay_outside_count_difference'] = actual['outside_tolerance']-REFERENCE_COUNT
                r['previous_replay_max_error_difference'] = actual['max_absolute_difference']-REFERENCE_MAX
            else:
                require(r['canonical_test_keys_sha256']==report['fits']['baseline']['canonical_test_keys_sha256'],
                        'baseline_comparison_key_order')
            report['fits'][arm] = r
            recovery.write_json(output,report)
        report['status'] = 'DIAGNOSTIC_COMPLETED'
    except Exception:
        report['status']='DIAGNOSTIC_INCOMPLETE'
        raise
    finally:
        recovery.write_json(output,report)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('stage',choices=['run','preflight','fit'])
    p.add_argument('--root',required=True,type=Path)
    p.add_argument('--snapshot',required=True,type=Path)
    p.add_argument('--arm',choices=ARMS)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.stage=='run':
        require(a.output is not None,'missing_output')
        orchestrate(a.root.resolve(),a.snapshot.resolve(),a.output)
    else:
        require(a.arm is not None,'missing_arm')
        run_arm(a.root.resolve(),a.snapshot.resolve(),a.arm,a.stage=='fit')


if __name__=='__main__':
    try:
        main()
    except ProbeError as exc:
        sys.exit(str(exc))
    except Exception as exc:
        sys.exit('SPECIES_PROBE_ERROR_'+type(exc).__name__)
