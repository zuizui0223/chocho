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


if __name__ == "__main__":
    unittest.main()
