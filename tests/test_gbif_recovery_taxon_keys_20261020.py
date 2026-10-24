"""Source-only preflight: synthetic receipts test contracts, not ecology."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/validate_gbif_recovery_taxon_keys_20261020.py"
KEYS = ROOT / "data/exploratory/gbif_239_partial96_exact_taxon_keys_20261019.csv"
FIELDS = ("species", "accepted_gbif_species_key")
SCHEMA = "chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01"


def original_source_fixture():
    """Retain the source's 135/96/8 design without committing its manifest."""
    receipts = []
    rows = []
    for i in range(239):
        species = "Chlosyne lacinia" if i == 135 else f"Butterfly species{i:03d}"
        key = 1914465 if i == 135 else 10000 + i
        status = (
            "COMPLETE" if i < 135 else
            "PARTIAL_FETCH_NOT_ELIGIBLE" if i < 231 else
            "TAXON_UNRESOLVED"
        )
        receipts.append({
            "species": species,
            "status": status,
            "taxon": None if status == "TAXON_UNRESOLVED" else {"key": key},
        })
        if status == "PARTIAL_FETCH_NOT_ELIGIBLE":
            rows.append([species, key])
    return {
        "schema": SCHEMA,
        "species_requested": 239,
        "species_statuses": {
            "COMPLETE": 135, "PARTIAL_FETCH_NOT_ELIGIBLE": 96,
            "TAXON_UNRESOLVED": 8,
        },
        "species_receipts": receipts,
    }, rows


class TestRecoveryTaxonKeyPreflight(unittest.TestCase):
    def run_preflight(self, source, rows, fields=FIELDS):
        with TemporaryDirectory() as folder:
            root = Path(folder)
            manifest = root / "original.json"
            keys = root / "recovery.csv"
            manifest.write_text(json.dumps(source), encoding="utf-8")
            with keys.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(fields)
                writer.writerows(rows)
            before = (manifest.read_bytes(), keys.read_bytes())
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--original-manifest", str(manifest),
                 "--keys", str(keys)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual((manifest.read_bytes(), keys.read_bytes()), before)
            return result, tuple(hashlib.sha256(b).hexdigest() for b in before)

    def assert_rejected(self, source, rows, *messages):
        result, _ = self.run_preflight(source, rows)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        for message in messages:
            self.assertIn(message, result.stderr)

    def test_exact_96_source_pairs_pass_and_inputs_are_unchanged(self):
        source, rows = original_source_fixture()
        result, hashes = self.run_preflight(source, rows)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt = json.loads(result.stdout)
        self.assertEqual(receipt["status"], "PASS_ORIGINAL_96_TAXON_KEYS")
        self.assertEqual(receipt["recovery_species_verified"], 96)
        self.assertEqual(receipt["original_species_verified"], 239)
        self.assertEqual(receipt["original_manifest_sha256"], hashes[0])
        self.assertEqual(receipt["recovery_keys_sha256"], hashes[1])
        self.assertEqual(result.stderr, "")

    def test_prior_chlosyne_key_is_rejected_before_acquisition(self):
        source, rows = original_source_fixture()
        rows[0][1] = 1914469
        self.assert_rejected(source, rows, "Chlosyne lacinia", "1914465", "1914469")

    def test_duplicate_species_does_not_hide_a_missing_species(self):
        source, rows = original_source_fixture()
        rows[1] = rows[0][:]
        self.assert_rejected(source, rows, "Duplicate recovery species", "Chlosyne lacinia")

    def test_duplicate_taxon_key_is_rejected(self):
        source, rows = original_source_fixture()
        rows[1][1] = rows[0][1]
        self.assert_rejected(source, rows, "Duplicate recovery taxon key", "1914465")

    def test_missing_species_is_named(self):
        source, rows = original_source_fixture()
        self.assert_rejected(source, rows[1:], "Missing recovery species", "Chlosyne lacinia")

    def test_original_complete_species_cannot_replace_a_partial_species(self):
        source, rows = original_source_fixture()
        rows[0] = ["Butterfly species000", 10000]
        self.assert_rejected(source, rows, "Unexpected recovery species", "Butterfly species000")

    def test_malformed_or_nonpositive_taxon_keys_fail_closed(self):
        for bad_key in ("", "not-a-key", "1914465.0", "0", "-1"):
            with self.subTest(key=bad_key):
                source, rows = original_source_fixture()
                rows[0][1] = bad_key
                self.assert_rejected(source, rows, "Invalid recovery taxon key", "Chlosyne lacinia")

    def test_duplicate_original_species_is_rejected(self):
        source, rows = original_source_fixture()
        source["species_receipts"][136] = dict(source["species_receipts"][135])
        self.assert_rejected(source, rows, "Duplicate original species", "Chlosyne lacinia")

    def test_missing_original_taxon_key_fails_closed(self):
        source, rows = original_source_fixture()
        source["species_receipts"][135]["taxon"] = None
        self.assert_rejected(source, rows, "Invalid original taxon key", "Chlosyne lacinia")

    def test_original_status_counts_must_match_receipts(self):
        source, rows = original_source_fixture()
        source["species_receipts"][135]["status"] = "COMPLETE"
        self.assert_rejected(source, rows, "Original receipt status counts")

    def test_wrong_original_source_schema_is_rejected(self):
        source, rows = original_source_fixture()
        source["schema"] = "other_source"
        self.assert_rejected(source, rows, "Original source schema")

    def test_csv_schema_mismatch_is_rejected(self):
        source, rows = original_source_fixture()
        result, _ = self.run_preflight(source, rows, ("species", "taxonKey"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Recovery CSV schema", result.stderr)

    def test_repository_chlosyne_key_matches_original_source_receipt(self):
        with KEYS.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        matching = [r for r in rows if r["species"] == "Chlosyne lacinia"]
        self.assertEqual(len(matching), 1)
        self.assertEqual(int(matching[0]["accepted_gbif_species_key"]), 1914465)


if __name__ == "__main__":
    unittest.main()
