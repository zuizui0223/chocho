"""Regression checks for separate Malaysian 1987-1989 *crown* flowering index."""
import unittest
import analyze
class ACIARPublishedSeries(unittest.TestCase):
 def test_n_months_and_missing_month(self):
  d=analyze.analyze()
  self.assertEqual(d['reported_months'],23)
  self.assertEqual(d['omitted_month'],['1989-02'])
 def test_table_extrema(self):
  d=analyze.analyze()
  self.assertEqual(d['min']['intensity_percent'],.6)
  self.assertEqual(d['max']['intensity_percent'],85.6)
  self.assertEqual(d['same_month_january_1988_vs_1989'],[65.6,.6])
 def test_monthly_score_is_not_daily_larval_food(self):
  d=analyze.analyze()
  self.assertIsNone(d['joint_singapore_butterfly_plant_pheno_effect'])
  self.assertEqual(d['month_categories'],{'little':10,'moderate':11,'heavy':2})
 def test_constructed_equal_monthly_mean_different_feeding_window(self):
  self.assertEqual(analyze.viable_start_days([1]*15+[0]*15,10),6)
  self.assertEqual(analyze.viable_start_days([1,0]*15,10),0)
if __name__=='__main__':unittest.main()
