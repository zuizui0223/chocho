"""Synthetic behavior tests for the independent atlas verifier; no model fits."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("verify_atlas.py")


def fold(island, reference, augmented, *, comparison="C_to_HC", mae_reference=.2, mae_augmented=.1):
    ref, aug = {"C_to_HC": ("C", "HC"), "H_to_HC": ("H", "HC")}[comparison]
    return dict(target="new_butterfly_jaccard", comparison=comparison,
                reference_model=ref, augmented_model=aug, heldout_island_id=str(island),
                train_pairs=1225, test_pairs=50, reference_mse=reference, augmented_mse=augmented,
                delta_mse_augmented_minus_reference=augmented-reference,
                reference_mae=mae_reference, augmented_mae=mae_augmented,
                improved=augmented<reference, worse=augmented>reference, tied=augmented==reference)


def geography_fixture():
    records = [(10, "Port Ten", "AA", 179., 0.), (3, "Port Three", "CC", 0., 10.),
               (2, "Port Two", "BB", -179., 0.)]
    features = [{"properties": {"OBJECTID": key, "PORT_NAME": name, "COUNTRY": country,
                                 "HARBORSIZE": "S", "HARBORTYPE": "CB"},
                 "geometry": {"type": "Point", "coordinates": [lon, lat]}}
                for key, name, country, lon, lat in records]
    ledger = []
    for key, name, _, lon, lat in records:
        included = key != 2
        row = dict(landmass_id="2" if key == 3 else "1", port_id=str(key),
                   longitude=lon, latitude=lat, source_harbor_type="CB", source_port_name=name,
                   included=included, status="INCLUDED" if included else "EXCLUDED",
                   reason="" if included else "fixture exclusion")
        row.update(dict(node=key if included else "", raster_row=1 if included else "",
                        raster_col=2 if included else "", marine_longitude=lon if included else "",
                        marine_latitude=lat if included else "", displacement_km=0. if included else ""))
        ledger.append(row)
    recovery = {
        "source_51_port_groups": {"1": ["10", "2"], "2": ["3"]},
        "verified_source_51_mean_port_points": {
            "1": dict(longitude=180., latitude=0., n_complete_single_landmass_ports=2),
            "2": dict(longitude=0., latitude=10., n_complete_single_landmass_ports=1)},
    }
    roster = [dict(landmass_id="1", eligible=True, ports=2), dict(landmass_id="2", eligible=True, ports=1)]
    prep = {"by_island": {"1": dict(original_ports=2, marine_ports=1, unique_terminal_nodes=1),
                           "2": dict(original_ports=1, marine_ports=1, unique_terminal_nodes=1)}}
    return recovery, {"features": features}, ledger, roster, prep


class IndependentAtlasVerifierTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.is_file(), "Independent atlas verifier is not implemented")
        spec = importlib.util.spec_from_file_location("independent_atlas_verifier", SCRIPT)
        self.audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.audit)

    def test_global_percentage_uses_reference_sum_not_mean_island_percentages(self):
        gains, summary = self.audit.summarize_group([fold("1", .01, .005), fold("2", .09, .081)])
        self.assertAlmostEqual(gains[0]["mse_relative_gain_percent"], 50)
        self.assertAlmostEqual(gains[1]["mse_relative_gain_percent"], 10)
        self.assertAlmostEqual(summary["mse_net_global_pp"], 14)
        self.assertAlmostEqual(sum(row["mse_global_gain_contribution_pp"] for row in gains), 14)
        self.assertAlmostEqual(summary["mse_net_mean_gain"], .007)
        self.assertAlmostEqual(gains[0]["mse_gain"], -gains[0]["delta_mse_augmented_minus_reference"])

    def test_zero_reference_retains_adverse_gain_and_undefined_relative_percentage(self):
        gains, summary = self.audit.summarize_group([fold("1", 0, .1), fold("2", .1, .1)])
        self.assertIsNone(gains[0]["mse_relative_gain_percent"])
        self.assertAlmostEqual(gains[0]["mse_gain"], -.1)
        self.assertAlmostEqual(summary["mse_net_global_pp"], -100)
        self.assertEqual(summary["positive_islands"], 0)
        self.assertIsNone(summary["positive_top1_share"])
        self.assertIsNone(summary["positive_half_mass_count"])

    def test_all_zero_mass_and_reference_keep_undefined_values(self):
        gains, summary = self.audit.summarize_group([fold("1", 0, 0), fold("2", 0, 0)])
        self.assertTrue(all(row["mse_relative_gain_percent"] is None for row in gains))
        self.assertIsNone(summary["mse_net_global_pp"])
        self.assertIsNone(summary["positive_top5_share"])
        self.assertIsNone(summary["adverse_half_mass_count"])
        self.assertEqual(summary["tied_islands"], 2)

    def test_concentration_uses_gross_mass_despite_positive_negative_cancellation(self):
        rows = [fold("1", .8, .2), fold("2", .6, .2), fold("3", .1, 1.), fold("4", .2, .2)]
        _, summary = self.audit.summarize_group(rows)
        self.assertAlmostEqual(summary["mse_positive_mass"], 1.)
        self.assertAlmostEqual(summary["mse_adverse_mass"], .9)
        self.assertAlmostEqual(summary["mse_net_mass"], .1)
        self.assertAlmostEqual(summary["positive_top1_share"], .6)
        self.assertEqual(summary["positive_half_mass_count"], 1)
        self.assertEqual(summary["adverse_top5_share"], 1.)

    def test_exact_gain_ties_use_numeric_ids_not_input_or_lexical_order(self):
        rows = [fold("10", .2, .1), fold("2", .2, .1)]
        _, summary = self.audit.summarize_group(rows)
        self.assertEqual(summary["ranked_positive_island_ids"], ["2", "10"])
        self.assertEqual(summary["positive_top1_share"], .5)

    def test_wrong_saved_delta_is_rejected_before_geographic_accounting(self):
        row = fold("1", .2, .1)
        row["delta_mse_augmented_minus_reference"] = .1
        with self.assertRaisesRegex(ValueError, "delta"):
            self.audit.summarize_group([row])

    def test_joint_sign_categories_preserve_ties_and_mae_disagreement(self):
        gains = []
        signs = [(1, 1), (1, -1), (-1, 1), (-1, -1), (0, 1), (1, 0), (0, 0), (-1, 0), (0, -1)]
        for i, (host, circuit) in enumerate(signs, 1):
            for comparison, sign in [("C_to_HC", host), ("H_to_HC", circuit)]:
                part, _ = self.audit.summarize_group([fold(str(i), .2+.1*sign, .2,
                    comparison=comparison, mae_reference=.1, mae_augmented=.2)])
                gains.extend(part)
        rows = self.audit.conditional_status(gains)
        labels = [row["mse_category"] for row in rows]
        self.assertEqual(labels[:4], ["both_positive", "host_positive_circuit_negative",
                                      "host_negative_circuit_positive", "both_negative"])
        self.assertEqual(labels.count("at_least_one_tie"), 5)
        self.assertTrue(all(not row["both_mae_positive"] for row in rows))

    def test_joint_hc_loss_must_be_shared_without_averaging_disagreement(self):
        a, _ = self.audit.summarize_group([fold("1", .3, .2)])
        b, _ = self.audit.summarize_group([fold("1", .3, .21, comparison="H_to_HC")])
        with self.assertRaisesRegex(ValueError, "HC"):
            self.audit.conditional_status(a+b)

    def test_dateline_point_and_example_port_use_all_original_ports(self):
        islands, ports, checks = self.audit.reconstruct_geography(*geography_fixture())
        first = next(row for row in islands if row["island_id"] == "1")
        self.assertAlmostEqual(abs(first["representative_longitude"]), 180)
        self.assertEqual(first["source_port_ids"], ["2", "10"])
        self.assertEqual(first["source_port_names"], ["Port Two", "Port Ten"])
        self.assertEqual(first["sampled_country_codes"], ["AA", "BB"])
        self.assertEqual(first["example_port_id"], "2")
        self.assertEqual(first["display_label"], "ID 1 | example port: Port Two")
        self.assertEqual((first["original_ports"], first["marine_ports"]), (2, 1))
        excluded = next(row for row in ports if row["port_id"] == "2")
        self.assertIsNone(excluded["marine_longitude"])

    def test_absent_marine_location_cannot_be_replaced_by_observed_zero(self):
        values = list(geography_fixture())
        next(row for row in values[2] if row["port_id"] == "2")["marine_longitude"] = 0.
        with self.assertRaisesRegex(ValueError, "excluded.*blank"):
            self.audit.reconstruct_geography(*values)

    def test_port_names_are_joined_by_objectid_and_must_match_the_ledger(self):
        values = list(geography_fixture())
        values[2][0]["source_port_name"] = "Wrong place"
        with self.assertRaisesRegex(ValueError, "port name"):
            self.audit.reconstruct_geography(*values)

    def test_missing_adverse_table_record_is_rejected(self):
        expected = [dict(island_id="1", mse_gain=.1), dict(island_id="2", mse_gain=-.1)]
        with self.assertRaisesRegex(ValueError, "coverage"):
            self.audit.compare_table(expected, expected[:1], ("island_id",), "gain")

    def test_undefined_saved_percentage_is_not_accepted_as_zero(self):
        expected = [dict(island_id="1", mse_relative_gain_percent=None)]
        with self.assertRaisesRegex(ValueError, "undefined"):
            self.audit.compare_table(expected, [dict(island_id="1", mse_relative_gain_percent="0")], ("island_id",), "gain")

    def test_changed_frozen_file_fails_hash_gate(self):
        self.assertTrue(hasattr(self.audit, "verify_manifest"), "Frozen-file audit is missing")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = b"locked bytes"
            (root/"source.csv").write_bytes(b"altered data")
            entries = [dict(path="source.csv", bytes=len(original), sha256=hashlib.sha256(original).hexdigest())]
            with self.assertRaisesRegex(ValueError, "hash"):
                self.audit.verify_manifest(root, entries)

    def test_verified_manifest_records_actual_input_identity(self):
        self.assertTrue(hasattr(self.audit, "verify_manifest"), "Frozen-file audit is missing")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = b"locked bytes"
            (root/"source.csv").write_bytes(original)
            sha = hashlib.sha256(original).hexdigest()
            actual = self.audit.verify_manifest(root, [dict(path="source.csv", bytes=len(original), sha256=sha)])
            self.assertEqual(actual, {"source.csv": sha})

    def test_cli_refuses_to_write_even_failure_receipt_into_parent_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            current = root/"results/island51_host_circuit_atlas_20261106"
            parent = root/"results/island51_host_circuit_joint_20261106"
            current.mkdir(parents=True)
            parent.mkdir(parents=True)
            script = current/"verify_atlas.py"
            script.write_text(SCRIPT.read_text())
            result = subprocess.run([sys.executable, str(script), "--results-dir", str(parent)],
                                    capture_output=True, text=True, check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((parent/"independent_atlas_audit.json").exists())


if __name__ == "__main__":
    unittest.main()
