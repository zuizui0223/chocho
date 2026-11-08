"""NC exploratory public 17-island index-risk/source-access prediction gate."""
import importlib.util
import math
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("ncomms_p17",ROOT/"scripts/run_ncomms_public17_ocean_access_holdout.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class SourceLockedAccessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources,cls.pins=m.load()
        cls.early,cls.early_info=m.prepare(cls.sources,"2011_2017")
        cls.late,cls.late_info=m.prepare(cls.sources,"2018_2022")

    def test_original_source_has_seventeen_islands_225_butterflies(self):
        self.assertEqual(len(self.sources["butterflies"]),3825)
        self.assertEqual(len(self.pins),8)
        self.assertEqual(len({r["landmass_id"] for r in self.sources["butterflies"]}),17)
        self.assertEqual(len({r["speciesKey"] for r in self.sources["butterflies"]}),225)
        self.assertTrue(all(z["git_blob_sha1"]==m.SOURCES[k][1] for k,z in self.pins.items()))

    def test_future_record_periods_observationally_heldout(self):
        self.assertEqual((len(self.early),self.early_info["first_index_cases"]),(3584,80))
        self.assertEqual((len(self.late),self.late_info["first_index_cases"]),(3504,18))
        self.assertEqual(self.early_info["prior_source_rows"],864)
        self.assertEqual(self.late_info["prior_source_rows"],1039)
        self.assertEqual(self.early_info["prior_source_case_count"],57)
        self.assertEqual(self.late_info["prior_source_case_count"],17)

    def test_no_outcome_based_host_or_shipping_features(self):
        for row in self.early+self.late:
            self.assertEqual(set(row["features"]),set(m.FEATURES["B"]+m.FEATURES["D"]+m.FEATURES["C"]))
            self.assertTrue(all(math.isfinite(x) for x in row["features"].values()))
            if not row["prior_source_islands"]:
                self.assertEqual(row["features"]["min_log_geodesic_km"],0)
                self.assertEqual(row["features"]["min_log_uniform_effective_resistance"],0)
        self.assertNotIn("host",str(m.FEATURES).lower())
        self.assertNotIn("ais",str(m.FEATURES).lower())

    def test_stratified_island_region_and_species_folds_exact_once(self):
        self.assertEqual(len({r["landmass_id"] for r in self.early}),17)
        self.assertEqual(len({r["component"] for r in self.early}),6)
        self.assertEqual(len({str(int(r["speciesKey"])%5) for r in self.early}),5)

    def test_original_public_temporal_risk_validity_margins(self):
        for r in self.sources["butterflies"]:
            n0=int(r["pre2001_2010_indexed_records"])
            n1=int(r["future2011_2017_indexed_records"])
            n2=int(r["future2018_2022_indexed_records"])
            self.assertEqual(int(r["eligible_2011_2017"]),int(n0==0))
            self.assertEqual(int(r["first_indexed_in_2011_2017"]),int(n0==0 and n1>0))
            self.assertEqual(int(r["eligible_2018_2022"]),int(n0==0 and n1==0))
            self.assertEqual(int(r["first_indexed_in_2018_2022"]),int(n0==0 and n1==0 and n2>0))

    def test_finite_and_positive_risk_metrics(self):
        y=[0,0,1,1]
        p=[0.1,0.4,0.8,0.9]
        a=m.loss(__import__("numpy").array(y),__import__("numpy").array(p))
        self.assertTrue(0<a["brier"]<1)
        self.assertTrue(0<a["logloss"]<2)


if __name__=="__main__":
    unittest.main()
