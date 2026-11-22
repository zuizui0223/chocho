import hashlib
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import probe_GEB_v05_species_runtime_20261122 as p


class RuntimeProbe(unittest.TestCase):
    def test_array_signature_preserves_values_dtype_and_original_layout(self):
        a = np.arange(12, dtype=np.float32).reshape((3,4))
        value = p.array_signature(a)
        self.assertEqual(value['sha256'], hashlib.sha256(a.tobytes()).hexdigest())
        self.assertNotEqual(value, p.array_signature(np.asfortranarray(a)))
        self.assertNotEqual(value, p.array_signature(a[::-1]))
        self.assertNotEqual(value, p.array_signature(a.astype('float64')))

    def test_component_arm_selection_uses_estimator_arrays_only(self):
        def fingerprint(values, key='original'):
            return [{'fold':i,'arrays':{'x':n},'train_key_order':key} for i,n in enumerate(values)]
        base = fingerprint([1]*5)
        pandas = fingerprint([2]*5)
        arms = dict(baseline=base,pandas_only=pandas,arrow_only=fingerprint([1]*5,'other'),
                    joint=fingerprint([2]*5))
        self.assertEqual(p.select_arms(arms), ['baseline','pandas_only'])
        self.assertEqual(p.equivalent_arms(arms)['joint'],'pandas_only')
        arms['joint'] = fingerprint([3]*5)
        self.assertEqual(p.select_arms(arms), ['baseline','pandas_only','joint'])

    def test_keys_hash_order_and_boundaries_without_exporting_values(self):
        f = pd.DataFrame({'speciesKey':[12,1],'cell':['private3','private23']})
        self.assertNotEqual(p.key_hash(f),p.key_hash(f.iloc[::-1]))
        self.assertNotIn('private',p.key_hash(f))
        self.assertEqual(len(p.key_hash(f)),64)

    def test_source_binding_and_model_boundary(self):
        public = Path(__file__).resolve().parents[1]/'provenance/reviewer_defenses/results/host_forecast_validation_20261029/analysis_source/fit_forecast.py'
        code = p.instrument(public.read_bytes())
        self.assertIn('probe.factory(name, fold, tr, te, settings, features)',code)
        self.assertIn('weights*=probe.normalization(len(weights)/weights.sum(), weights)',code)
        self.assertIn("reg.fit(tr[features].to_numpy(dtype='float32'),tr.target.to_numpy(),sample_weight=weights)",code)
        with self.assertRaises(p.ProbeError):
            p.instrument(public.read_bytes()+b'#')

    def test_aggregate_checked_predictions_keep_values_private(self):
        te = pd.DataFrame({'speciesKey':[2,1],'cell':['secret-a','secret-b']})
        rows = te.assign(species_fold=1,M0_background=[.3,.4]).set_index(['speciesKey','cell']).iloc[::-1]
        observer = p.Observer(rows, False)
        reg = observer.factory('M0_background',1,te,te,dict(p.recovery.MODEL_SETTINGS),['latitude'])
        x = np.array([[2.],[1.]], dtype='float32')
        observer.normalization(1.)
        reg.fit(x,np.array([.3,.4]),sample_weight=np.ones(2))
        reg.predict(x)
        self.assertNotIn('secret',json.dumps(observer.folds))
        self.assertNotIn('comparison',observer.folds[0])
        self.assertEqual(observer.folds[0]['arrays']['test_features'],p.array_signature(x))

    def test_source_preflight_executes_five_m0_boundaries_without_estimators(self):
        source = Path(__file__).resolve().parents[1]/'provenance/reviewer_defenses/results/host_forecast_validation_20261029/analysis_source/fit_forecast.py'
        labels = {}
        for i in range(100):
            name = 'synthetic_'+str(i)
            fold = int(hashlib.sha256(name.encode()).hexdigest()[:8],16)%5
            labels[fold] = name
            if len(labels) == 5:
                break
        panel = pd.DataFrame([dict(species=name,cell=fold,speciesKey=fold,period=period,
                                   days=12,distance_km=10,recording_probability10=.25,
                                   spatial_block='synthetic',region='synthetic',climate_id=3)
                              for period in range(3) for fold,name in sorted(labels.items())])
        climate=['temperature','precipitation','temperature_range','warmest_month','coldest_month']
        cover=['cropland','crop_mosaic','natural_crop_mosaic','tree','tree_mosaic','herbaceous_mosaic',
               'shrub','grass','wetland','built','bare','sparse','lichen','water','snow']
        for name in (climate+cover+[k+'_trend_per_year' for k in climate+cover]+
                     ['latitude','longitude','provider_dominance','available_hosts','available_families',
                      'global_hosts','available_share','native_host_count','introduced_host_count',
                      'introduced_share']+['season_q'+str(i) for i in range(4)]):
            panel[name] = 0.
        ref = panel[panel.period==1][p.check.KEYS].copy()
        ref['species_fold']=ref.speciesKey
        obs=p.Observer(ref.set_index(p.check.KEYS),False)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'outputs').mkdir()
            (root/p.recovery.LOCK).write_text(json.dumps({'random_forest':p.recovery.MODEL_SETTINGS}))
            script=root/'scripts/temporal/fit_forecast.py'
            with patch.object(pd,'read_parquet',return_value=panel),patch.object(sys,'argv',[str(script),'species']), \
                 patch('sklearn.ensemble.RandomForestRegressor.fit',side_effect=AssertionError('Unexpected fit')), \
                 self.assertRaises(p.FiveFoldsComplete):
                exec(compile(p.instrument(source.read_bytes()),str(script),'exec'),
                     dict(__file__=str(script),__name__='__main__',probe=obs))
            self.assertEqual(list((root/'outputs').iterdir()),[root/p.recovery.LOCK])
        self.assertEqual([r['fold'] for r in obs.folds],list(range(5)))
        self.assertTrue(obs.coverage.all())
        self.assertNotIn('synthetic',json.dumps(obs.folds))


if __name__=='__main__':
    unittest.main()
