#!/usr/bin/env python3
"""Regression tests for immutable NC temporal source screen; stdlib only."""
import copy
import importlib.util
import math
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("audit_nc",ROOT/"scripts/audit_ncomms_temporal_gate.py")
nc=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(nc)


class SourceAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data,cls.pins=nc.read_pinned(ROOT)

    def test_source_identity_and_coverage(self):
        self.assertEqual(len(self.data),10)
        self.assertEqual(len(self.pins),10)
        for key,pin in self.pins.items():
            self.assertEqual(pin["source_git_blob_sha1"],nc.SOURCE[key][1])

    def test_real_archived_225_species(self):
        values=nc.summarize_species(self.data["species"])
        self.assertEqual(values["original_candidate_butterflies"],225)
        self.assertEqual(values["butterflies_with_later_source_pairs"],82)
        self.assertEqual(values["later_indexed_butterfly_species_island_pairs"],1163)
        self.assertEqual(values["later_pairs_with_recorded_earlier_matching_host"],5)
        self.assertEqual(len(values["exact_butterflies_with_earlier_recorded_host_overlap"]),5)

    def test_independent_archived_counts_and_readiness(self):
        result=nc.audit(self.data,self.pins)
        self.assertEqual(result["status"],"ARCHIVED_SOURCE_COUNTS_RECOMPUTED_NO_ECOLOGICAL_MODEL")
        self.assertFalse(result["screening_passed"])
        self.assertIn("BLOCKED",result["primary_test_readiness"])
        self.assertEqual(result["source_coverage"]["original_source_ports"],194)
        self.assertEqual(result["source_coverage"]["historical_general_plant_islands_with_indexed_speciesKeys"],45)
        self.assertEqual(result["retrospective_vs_prior_comparison"]["later_record_pairs_with_matching_host_2001_2024"],208)
        self.assertEqual(result["retrospective_vs_prior_comparison"]["later_record_pairs_with_matching_host_2001_2010"],5)
        self.assertEqual(result["dated_host_family_preserving_null"]["draws_matching_or_exceeding_observed_pair_fraction"],44)

    def test_existing_155_species_temporal_spatial_host_positive_is_separate(self):
        result=nc.audit(self.data,self.pins)
        models=result["separate_existing_155_species_host_forecast"]
        self.assertEqual(set(models),{"temporal","spatial"})
        self.assertAlmostEqual(models["temporal"]["percentage_relative_species_equal_MSE_improvement"],
                               4.120175607931823,places=9)
        self.assertAlmostEqual(models["spatial"]["percentage_relative_species_equal_MSE_improvement"],
                               2.9952059516988134,places=9)
        for item in models.values():
            self.assertEqual(item["species"],155)
            self.assertEqual(item["test_rows"],208482)
            self.assertFalse(item["time_prior_host_arrivals_verified"])
            self.assertFalse(item["same_source_as_225_species_island_pairs"])

    def test_duplicate_species_rejected(self):
        rows=copy.deepcopy(self.data["species"])
        rows[1]["species"]=rows[0]["species"]
        with self.assertRaises(ValueError):
            nc.summarize_species(rows)

    def test_fraction_poison_rejected(self):
        rows=copy.deepcopy(self.data["species"])
        index=next(i for i,z in enumerate(rows) if z["comparison_evaluable"]=="True")
        rows[index]["fraction_matching_early"]="0.987654"
        with self.assertRaises(ValueError):
            nc.summarize_species(rows)

    def test_unsupported_early_match_rejected(self):
        rows=copy.deepcopy(self.data["species"])
        index=next(i for i,z in enumerate(rows) if z["comparison_evaluable"]=="True")
        rows[index]["matching_early_pairs"]="100000"
        with self.assertRaises(ValueError):
            nc.summarize_species(rows)

    def test_corrupt_source_validation_rejected(self):
        data=copy.deepcopy(self.data)
        data["dated_validation"]["matched_pairs"]=12
        with self.assertRaises(ValueError):
            nc.audit(data,self.pins)

    def test_null_is_a_descriptive_structural_comparator(self):
        result=nc.audit(self.data,self.pins)
        self.assertFalse(result["dated_host_family_preserving_null"]["causal_test"])
        self.assertTrue(result["archived_summary_cannot_identify_the_five_matching_islands"])
        text=nc.render_markdown(result)
        self.assertIn("NOT",text.upper())
        self.assertIn("5/1163",text.replace(" ",""))
        self.assertIn("not a statistical power calculation",text.lower())

    def test_finite_nonnegative_integer_guard(self):
        for bad in ["-3","nan","inf","4.5",""]:
            with self.assertRaises(ValueError):
                nc.fixed_integer(bad,label="test")
        self.assertEqual(nc.fixed_integer("0.0",label="test"),0)
        self.assertTrue(math.isfinite(nc.fixed_integer("200",label="test")))


if __name__=="__main__":
    unittest.main()
