"""Tiny synthetic guardrails for E. editha host-gap sensitivity logic.

Independent of the actual post-hoc result; tests structural set operations,
not host-choice ecology, demographic substitution, or global inference.
"""
from __future__ import annotations

import runpy
import unittest
from pathlib import Path
from collections import defaultdict

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_euphydryas_missing_hosts_augmentation.py"
NS = runpy.run_path(str(SCRIPT))
footprint = NS["footprint"]
describe = NS["describe"]
plantago_loss = NS["plantago_loss"]


class HostCompletenessSetTests(unittest.TestCase):
    def test_native_missing_host_reclassifies_added_region(self):
        # Plantago contributes r3 alone and overlaps B in r4.
        n = defaultdict(set, {"P": {"r1"}, "B": {"r2"}, "Castilleja": {"r3"}})
        c = defaultdict(set, {
            "P": {"r1", "r3", "r4"},
            "B": {"r2", "r4"},
            "Castilleja": {"r3"},
        })
        base = footprint({"P", "B"}, n, c)
        new = footprint({"P", "B", "Castilleja"}, n, c)
        self.assertEqual(describe(base), {"native": 2, "contemporary": 4, "introduced_added": 2})
        self.assertEqual(describe(new), {"native": 3, "contemporary": 4, "introduced_added": 1})
        self.assertEqual(plantago_loss({"P", "B"}, "P", n, c), {"r3"})
        self.assertEqual(plantago_loss({"P", "B", "Castilleja"}, "P", n, c), set())

    def test_missing_hosts_can_change_link_count_without_region_change(self):
        n = defaultdict(set, {"P": {"r1"}, "B": {"r2"}, "Castilleja": {"r2"}})
        c = defaultdict(set, {
            "P": {"r1", "r3"},
            "B": {"r2"},
            "Castilleja": {"r2"},
        })
        base = footprint({"P", "B"}, n, c)
        aug = footprint({"P", "B", "Castilleja"}, n, c)
        self.assertEqual(describe(base), describe(aug))
        self.assertEqual(describe(base)["introduced_added"], 1)
        self.assertEqual(plantago_loss({"P", "B"}, "P", n, c), {"r3"})
        self.assertEqual(plantago_loss({"P", "B", "Castilleja"}, "P", n, c), {"r3"})

    def test_new_host_never_increases_plantago_sole_dependence(self):
        n = defaultdict(set, {"P": {"r1"}, "B": {"r2"}, "C": {"r3"}})
        c = defaultdict(set, {"P": {"r1","r3"}, "B": {"r2"}, "C": {"r3","r4"}})
        base = plantago_loss({"P", "B"}, "P", n, c)
        augmented = plantago_loss({"P", "B", "C"}, "P", n, c)
        self.assertTrue(augmented <= base)

    def test_within_region_composition_can_change_without_area_gain(self):
        n = defaultdict(set, {
            "P": {"r1"}, "B": {"r2"},
            "Castilleja": {"r2"},
        })
        before = footprint({"P", "B"}, n, n)
        after = footprint({"P", "B", "Castilleja"}, n, n)
        self.assertEqual(describe(before), describe(after))
        x = NS["composition_difference"]({"P", "B"}, {"P", "B", "Castilleja"}, n)
        self.assertEqual(x["regional_host_incidences_added"], 1)
        self.assertEqual(x["regions_with_changed_recorded_host_count"], 1)
        self.assertEqual(x["regions_one_to_two_or_more_recorded_hosts"], 1)
        self.assertEqual(x["baseline_single_recorded_host_regions"], 2)
        self.assertEqual(x["augmented_single_recorded_host_regions"], 1)
        self.assertEqual(x["affected_region_codes_and_host_counts"], [{
            "region": "r2", "baseline_host_species": 1,
            "augmented_host_species": 2, "additional_species": 1
        }])

    def test_both_new_hosts_same_region_change_richness_by_two(self):
        n = defaultdict(set, {
            "old": {"rA","rB"}, "c1": {"rB"}, "c2": {"rB"}
        })
        x = NS["composition_difference"]({"old"}, {"old","c1","c2"}, n)
        self.assertEqual(x["regional_host_incidences_added"], 2)
        self.assertEqual(x["regions_with_changed_recorded_host_count"], 1)
        self.assertEqual(x["regions_by_count_increment"], {"plus_one": 0, "plus_two": 1})
        self.assertEqual(x["regions_one_to_two_or_more_recorded_hosts"], 1)


if __name__ == "__main__":
    unittest.main()
