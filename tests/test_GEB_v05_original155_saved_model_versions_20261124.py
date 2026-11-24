import sys
import unittest
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import audit_GEB_v05_original155_saved_model_versions_20261124 as audit


class StaticVersionMetadataTest(unittest.TestCase):
    def test_expected_30_model_names_no_traversal(self):
        paths = [audit.model_relative_path(c, m) for c in audit.CONDITIONS for m in audit.MODELS]
        self.assertEqual(len(paths), 30)
        self.assertEqual(len(set(paths)), 30)
        self.assertTrue(all(x.startswith('work/temporal/forecast/model_') for x in paths))
        self.assertRaises(audit.VersionAuditError, audit.model_relative_path, '../private', 'M0_background')

    def test_bounded_string_parser_on_compressed_synthetic_data(self):
        raw = b'prefix._sklearn_version\x94\x8c\x051.7.2\x94'
        # The scanner is static, it NEVER calls joblib.load or pickle.load.
        self.assertEqual(audit.read_version_safely(zlib.compress(raw, 3)), '1.7.2')
        self.assertEqual(audit.read_version_safely(raw), '1.7.2')

    def test_missing_string_is_unknown_not_compliance(self):
        self.assertIsNone(audit.read_version_safely(b'not a model'))


if __name__ == '__main__':
    unittest.main()
