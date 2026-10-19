"""Source-checked tests for an exploratory Melitaea reaction-norm audit."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_melitaea_ambiguity_and_crossovers_20261019 import (
    source_rows,
    effect_by_family,
    all_binary_ambiguity_scenarios,
    bootstrap_family_mean,
    symmetric_signflip_diagnostic,
    fit_illustrative_constant_family_host_effect,
    simulate_reversal_count,
)

class TestMelitaeaStageAndSensitivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.families, cls.ambiguous = source_rows()

    def test_original_counts_and_ambiguous_records(self):
        self.assertEqual(len(self.rows), 1680)
        self.assertEqual(len(self.families), 14)
        self.assertEqual(len(self.ambiguous), 4)
        self.assertEqual(
            {(self.rows[i]["Family"], self.rows[i]["Temperature"], self.rows[i]["Plant"])
             for i in self.ambiguous},
            {("6", "28", "Vs"), ("44", "34", "Vs"),
             ("13", "32", "Pl"), ("13", "28", "Vs")}
        )

    def test_all_sixteen_recodings_maintain_fixed_initial_denominators(self):
        outcomes = all_binary_ambiguity_scenarios(self.rows, self.families, self.ambiguous)
        self.assertEqual(len(outcomes), 16)
        self.assertEqual(len({r["ambiguous_success_pattern"] for r in outcomes}), 16)
        lo = min(x["estimate_Vs_minus_Pl"] for x in outcomes)
        hi = max(x["estimate_Vs_minus_Pl"] for x in outcomes)
        self.assertAlmostEqual(lo, 0.05)
        self.assertAlmostEqual(hi, 46/840)
        self.assertEqual({x["family_crossovers"] for x in outcomes}, {7, 8})
        self.assertLess(min(x["symmetric_signflip_two_sided"] for x in outcomes), 0.05)
        self.assertGreater(max(x["symmetric_signflip_two_sided"] for x in outcomes), 0.05)

    def test_primary_and_full_ambiguity_contrasts(self):
        primary, by_t, nflip = effect_by_family(
            self.rows, self.families, frozenset(self.ambiguous))
        all_unknown_fail, _, nflip_strict = effect_by_family(
            self.rows, self.families, frozenset())
        self.assertAlmostEqual(sum(primary)/len(primary), 45/840)
        self.assertAlmostEqual(sum(all_unknown_fail)/len(all_unknown_fail), 43/840)
        self.assertEqual(nflip, 7)
        self.assertEqual(nflip_strict, 8)
        self.assertEqual(set(by_t), {"28","30","32","34"})
        self.assertAlmostEqual(symmetric_signflip_diagnostic(primary), 0.0439453125)
        a, b = bootstrap_family_mean(primary, n=1000)
        self.assertLess(a, sum(primary)/len(primary))
        self.assertGreater(b, sum(primary)/len(primary))

    def test_null_simulation_is_unambiguously_hypothetical(self):
        p = fit_illustrative_constant_family_host_effect(self.rows, self.families)
        self.assertEqual(len(p), 14)
        self.assertTrue(all(len(f)==4 for f in p))
        self.assertTrue(all(0<=x<=1 for f in p for pair in f for x in pair))
        trials = simulate_reversal_count(p, correlations=(0, 0.05, 0.15), n=149)
        self.assertEqual(len(trials), 3)
        self.assertEqual({z["simulations"] for z in trials}, {149})
        self.assertTrue(all(0 <= z["simulated_fraction_ge_7"] <= 1 for z in trials))
        self.assertEqual(
            [z["assumed_larvae_within_group_correlation"] for z in trials],
            [0, 0.05, 0.15]
        )

if __name__ == "__main__":
    unittest.main()
