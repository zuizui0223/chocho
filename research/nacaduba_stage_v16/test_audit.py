"""Checks source-locked v16 arithmetic and forbids false biological inference."""
import unittest
from audit import audit

class V16Test(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.r=audit()
 def test_source_only_flowering_calendar(self):
  self.assertEqual(self.r["full_calendar_binyears"],57)
  self.assertEqual(self.r["season_unique_g_binyears"],[9,6,8])
 def test_minimal_clones_unique(self):
  self.assertEqual(self.r["three_clone_full_calendar_subsets"],[[146,148,152]])
  self.assertEqual(self.r["three_clones_ramets"],5)
 def test_source_pair_unique(self):
  self.assertEqual(self.r["two_provenance_pairs_complete"],[["g","h"]])
  self.assertEqual(self.r["g_plus_h_ramets"],12)
 def test_density_gate_changes_conclusion(self):
  a=self.r["concurrent_flowering_proxy"]
  self.assertEqual([a[x]["at_least_1"] for x in ("all","g_plus_h","three_clones")],[57,57,57])
  self.assertEqual([a[x]["at_least_3"] for x in ("all","g_plus_h","three_clones")],[48,42,32])
  self.assertEqual([a[x]["at_least_5"] for x in ("all","g_plus_h","three_clones")],[40,29,0])
  self.assertEqual([a[x]["flowering_ramet_decads"] for x in ("all","g_plus_h","three_clones")],[397,249,139])
 def test_size_matched_enumeration(self):
  self.assertEqual(self.r["g_matched_size_subsets_as_extreme"],1)
  self.assertEqual(self.r["seven_clone_matched_portfolios_complete"],2184)
  self.assertEqual(self.r["seven_clone_matched_portfolios_total"],91728)
 def test_no_butterfly_causality_inferred(self):
  for k in ("butterfly_larval_success_effect","flower_bud_mass","biological_colonization_lag","causal_p_value"):
   self.assertIsNone(self.r[k])

if __name__=="__main__":unittest.main()
