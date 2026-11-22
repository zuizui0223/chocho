import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import refit_GEB_v05_original155_sensitivities_20261122 as x
class TestBoundaries(unittest.TestCase):
 def test_exact_sensitivity_cohorts(self):
  self.assertEqual(set(x.S),{'effort20','radius250','exclude_top_provider','precision100'})
  self.assertEqual(x.S['effort20'][:3],(122547,159084,20))
  self.assertEqual(x.S['precision100'][:3],(98838,135532,10))
  for k,v in x.S.items():
   self.assertEqual(len(v[3]),64)
   self.assertTrue(x.source_names(k)[0].endswith(f'panel_{k}.parquet'))
   self.assertTrue(x.source_names(k)[1].endswith(f'predictions_temporal_{k}.parquet'))
 def test_guard_failure(self):
  with self.assertRaises(x.SafeError):x.ok(False,'test')
if __name__=='__main__':unittest.main()
