import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_GEB_v05_species_manifest_lineage_20261122 as a

class ManifestSafety(unittest.TestCase):
    def test_groups_omit_paths(self):
        records=[
            {"path":"/private/user/fit_forecast.py","sha256":"a"*64,"bytes":20},
            {"path":"C:\\hidden\\fit_forecast_original.py","sha256":"b"*64,"bytes":30},
            {"path":"/private/user/predictions_species.parquet","sha256":"c"*64,"bytes":900},
            {"path":"/private/user/secret.csv","sha256":"d"*64,"bytes":50},
        ]
        got=a.summarize(records)
        s=str(got)
        for private in ["private","hidden","user","secret","original.py"]:
            self.assertNotIn(private,s)
        self.assertEqual(got["fit_forecast_script_candidates"]["distinct_content_hashes"],2)
        self.assertEqual(got["species_predictions"]["entries"],1)

    def test_bad_matching_digest_rejected(self):
        with self.assertRaises(a.ScanError):
            a.summarize([{"path":"/sensitive/fit_forecast.py","sha256":"bad","bytes":1}])

    def test_substring_not_classified(self):
        self.assertIsNone(a.category("/private/fit_forecast.py.backup"))
        self.assertEqual(a.category("/sensitive/model_M0_background.joblib"),"saved_model_candidates")
        self.assertEqual(a.category("/sensitive/fit_forecast_v2.py"),"fit_forecast_script_candidates")

if __name__=="__main__":unittest.main()
