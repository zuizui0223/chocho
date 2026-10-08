"""Unit regression for the same-continent fixed-margin swap kernel."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path

SCRIPTS=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(SCRIPTS))
from analyze_bce_clarke_homogenization_level1_null import within_level1_null


class TestBceLevelOneNull(unittest.TestCase):
    def setUp(self):
        # Two continents with three regions each and 4 synthetic butterflies.
        # Every species has one native and one introduced-only region per Level1.
        self.native=[
           {0,3},{1,4},{2,5},{0,4}
        ]
        self.added=[
           {1,4},{2,5},{0,3},{2,5}
        ]
        self.labels=["A","A","A","B","B","B"]
        self.native_masks=[0]*6
        for i,row in enumerate(self.native):
            for r in row:
                self.native_masks[r]|=1<<i

    def test_preservation_and_repeatability(self):
        x=within_level1_null(
             self.native,self.added,self.native_masks,self.labels,19,20261007)
        y=within_level1_null(
             self.native,self.added,self.native_masks,self.labels,19,20261007)
        self.assertEqual(x,y)
        self.assertTrue(x["sampled_margins_verified"])
        self.assertEqual(x["added_edges_in_original_region_domain"],8)
        self.assertEqual(x["eligible_Level1_groups"],2)
        self.assertGreater(x["accepted_swaps_total"],0)
        self.assertEqual(len(x["regional_null_samples"]),19)
        self.assertTrue(all(0<=t<=1 for t in x["regional_null_samples"]))
        # Original input structures must be immutable.
        self.assertEqual(self.added,[{1,4},{2,5},{0,3},{2,5}])

    def test_refuses_incorrect_native_overlap(self):
        bad=[set(x) for x in self.added]
        bad[0].add(0)
        with self.assertRaisesRegex(ValueError,"disjoint"):
            within_level1_null(
                self.native,bad,self.native_masks,self.labels,4,0)

    def test_refuses_cross_domain_index(self):
        bad=[set(x) for x in self.added]
        bad[0].add(7)
        with self.assertRaisesRegex(ValueError,"outside"):
            within_level1_null(
                self.native,bad,self.native_masks,self.labels,4,0)


if __name__=="__main__":
    unittest.main()
