"""Temporary synthetic inputs exercise the real aggregator, never production outputs."""
from pathlib import Path
import contextlib
import hashlib
import io
import json
import tempfile
import unittest
import numpy as np
import pandas as pd
import aggregate_route_terrain as target


class AggregationSupportTest(unittest.TestCase):
    def test_missing_support_and_integrity(self):
        original = target.R, target.W, target.O
        with tempfile.TemporaryDirectory() as folder:
            try:
                target.R = Path(folder)
                target.W = w = target.R / 'work'
                target.O = o = target.R / 'outputs'
                w.mkdir(); o.mkdir(); (w/'route_terrain_values').mkdir()
                def save(path, value):
                    path.write_text(json.dumps(value))
                save(o/'route_terrain_design_lock.json', {'synthetic': True})
                lock = target.digest(o/'route_terrain_design_lock.json')
                mapping = pd.DataFrame({'speciesKey':[1,2,1,1,1], 'period':[0]*5,
                    'cell':['a','a','b','c','d'], 'route_key':['r0','r0','r1','r2','r3']})
                mapping.to_parquet(w/'sea_route_mapping.parquet', index=False)
                save(o/'route_terrain_design_audit.json', {'lock_sha256':lock,
                    'source_mapping_sha256':target.digest(w/'sea_route_mapping.parquet'),
                    'planned_land_samples':11})
                pd.DataFrame({'tile':['test'], 'route_land_sample_points':[11]}).to_csv(o/'route_terrain_required_tiles.csv', index=False)
                pd.DataFrame({'route_key':['r0','r1','r2','r3'], 'route_index':range(4)}).to_parquet(w/'route_terrain_routes.parquet', index=False)
                pd.DataFrame({'route_index':range(4),'planned_samples':[5,4,2,6],
                    'land_samples':[5,4,2,0],'sample_step_km':[2]*4}).to_parquet(w/'route_terrain_sampling_support.parquet', index=False)
                q=pd.DataFrame({'route_index':[0]*5+[1]*4+[2]*2,
                    'longitude':np.arange(11,dtype=float), 'latitude':np.zeros(11),
                    'elevation_m':[10,20,30,40,np.nan,10,20,30,np.nan,100,200]})
                data=w/'route_terrain_values/test.parquet'; q.to_parquet(data,index=False)
                receipt={'lock_sha256':lock,'points':11,'status':'COMPLETE','valid_points':9,
                    'data_sha256':target.digest(data),
                    'point_hash':hashlib.sha256(pd.util.hash_pandas_object(q[['route_index','longitude','latitude']],index=False).to_numpy().tobytes()).hexdigest()}
                rp=w/'route_terrain_values/test.json'; save(rp,receipt)
                with self.assertRaisesRegex(RuntimeError,'not finished'):
                    target.main()
                self.assertFalse((w/'route_terrain_features.parquet').exists())
                save(o/'route_terrain_acquisition.json',{'requested_tiles':1,'partial_limit':0})
                save(rp,{**receipt,'status':'ERROR'})
                with self.assertRaisesRegex(RuntimeError,'failed retrieval'):
                    target.main()
                save(rp,{**receipt,'data_sha256':'corrupted'})
                with self.assertRaises(AssertionError):
                    target.main()
                save(rp,receipt)
                with contextlib.redirect_stdout(io.StringIO()):
                    target.main()
                m=pd.read_parquet(w/'route_terrain_metrics.parquet').set_index('route_index')
                self.assertEqual(m.route_terrain_unknown.tolist(),[0,1,1,1])
                self.assertAlmostEqual(m.loc[0,'route_elevation_p90_m'],37)
                self.assertAlmostEqual(m.loc[0,'route_elevation_spread_m'],24)
                self.assertAlmostEqual(m.loc[0,'log_route_elevation_p90'],np.log(38))
                self.assertAlmostEqual(m.loc[1,'valid_land_fraction'],.75)
                self.assertTrue(pd.isna(m.loc[3,'route_elevation_p90_m']))
                self.assertTrue(pd.isna(m.loc[3,'valid_land_fraction']))
                self.assertTrue(m.loc[[1,2,3],'log_route_elevation_p90'].eq(0).all())
                f=pd.read_parquet(w/'route_terrain_features.parquet')
                pd.testing.assert_frame_equal(f[mapping.columns],mapping)
                self.assertEqual(len(f),5)
            finally:
                target.R,target.W,target.O=original


if __name__=='__main__':
    unittest.main()
