"""Research-result regression checks for original flower-data transcription."""
import unittest
from audit import analyze

class KatoProvenance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=analyze()
    def test_table_sample_and_event_accounting(self):
        self.assertEqual((self.d['n_clones'],self.d['n_ramets'],self.d['n_flowering_ramet_year_events']),(23,39,79))
    def test_original_clone_flowering_means(self):
        expected=[65.21739130434783,73.91304347826087,60.869565217391305]
        for a,b in zip([v['mean_clone_flowering_percent'] for v in self.d['years']],expected):
            self.assertAlmostEqual(a,b)
    def test_source_g_sole_coverage(self):
        self.assertEqual([r['thai_only_10day_bins'] for r in self.d['years']],[9,6,8])
        self.assertEqual([r['combined_10day_bins'] for r in self.d['years']],[19,17,21])
    def test_source_g_all_bloom_each_year(self):
        self.assertEqual([r['source_g_flowering_ramets'] for r in self.d['years']],[6,6,6])
    def test_conflicting_reported_duration_left_visible(self):
        x=self.d['source_discrepancy']
        self.assertEqual(x['2009_max_in_paper_prose_10day_bins'],8)
        self.assertEqual(x['2009_max_in_table7_10day_bins'],13)
    def test_no_invented_local_ecology(self):
        for k in ('singapore_daily_buds','butterfly_flowering_mismatch_effect','biological_colonization_lag'):
            self.assertIsNone(self.d[k])

if __name__=='__main__':unittest.main()
