import sys,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import replay_GEB_v05_original155_all_columns_20261124 as m
class TestComparison(unittest.TestCase):
 def test_critical_count_source_pinned(self):
  self.assertEqual(len(m.EXPECTED_SHA),3)
  self.assertEqual(len(m.EXPECTED_SHA['work/remote_sensing/host_region_counts.csv']),64)
  self.assertEqual(m.SCRIPTHASH,'99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6')
 def test_comparison_preserves_mismatch(self):
  a=pd.DataFrame({'speciesKey':[1],'cell':['1:1'],'period':[0],
   'species':['a'],'native_host_count':[2.], 'introduced_host_count':[3.],
   'introduced_share':[.6]})
  with self.assertRaises(m.GateError):
   m.check_parquets(a,a)
 def test_not_clear_if_count_not_same(self):
  self.assertTrue(m.EXPECTED_SHA.get('work/remote_sensing/host_region_counts.csv'))
if __name__=='__main__':unittest.main()
