"""Independent reproducibility checks for 239-species macro habitat filter."""
import sys
import math
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_global_butterfly_habitat_filter_20261019 import (
    read_inputs, annotate, mean, model_matrices, _coefficient,
    wild_genus_signflip, grouped_out_of_genus_cv,
    STRONG_ASSOC, STRONG_AVOID
)

class TestMacroHabitatFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.habitats, cls.digests = read_inputs()
        cls.strong = annotate(cls.rows)
        cls.fams = sorted({r["family"] for r in cls.rows})

    def test_original_data_and_species_numbers(self):
        r=self.rows
        self.assertEqual(len(r),239)
        self.assertEqual(len(self.habitats),339)
        self.assertEqual(sum(int(x["native"]) for x in r),26530)
        self.assertEqual(sum(int(x["contemporary"]) for x in r),41083)
        self.assertEqual(sum(int(x["added"]) for x in r),14553)
        self.assertTrue(all(math.isclose(x["log_expansion"],
              math.log(x["contemporary"]/x["native"]),abs_tol=1e-10) for x in r))
        self.assertEqual(len(self.fams),6)
        self.assertTrue(all(len(s)==64 for s in self.digests.values()))

    def test_exposure_definition_not_reselected_after_fit(self):
        rows=self.strong
        self.assertEqual(len(rows),183)
        self.assertEqual(sum(x["exposure"] for x in rows),150)
        self.assertEqual(sum(not x["exposure"] for x in rows),33)
        self.assertEqual({x["disturbance"] for x in rows},
                         {STRONG_ASSOC,STRONG_AVOID})
        full_trait_counts = {
          t:sum(r["disturbance_affinity"]==t for r in self.habitats)
          for t in (STRONG_ASSOC,STRONG_AVOID)
        }
        self.assertEqual(full_trait_counts[STRONG_ASSOC],193)
        self.assertEqual(full_trait_counts[STRONG_AVOID],68)

    def test_linear_synthetic_effect_identified_with_controls(self):
        # Synthetic algebra sanity-check ONLY: never report as ecological
        # observed response. The original real-data estimator is separate.
        Z, X, real_y = model_matrices(self.strong,self.fams)
        artificial_y = 0.47*Z[:,0] + .03*Z[:,1] + 1.37*X[:,-1]
        synthetic=[{**row,"artificial":float(y)} for row,y in zip(self.strong,artificial_y)]
        b,_,_,_,_,_ = _coefficient(synthetic,self.fams,y_key="artificial")
        self.assertAlmostEqual(b,1.37,places=8)

    def test_wild_grouped_signflip_cannot_claim_exact_randomization(self):
        Z,X,y=model_matrices(self.strong,self.fams)
        val=wild_genus_signflip(Z,y,X,[r["genus"] for r in self.strong],n=99)
        self.assertEqual(val["simulation_draws"],99)
        self.assertEqual(val["status"],"EXPLORATORY_GENUS_WILD_SIGNFLIP_REDUCED_RESIDUAL")
        self.assertGreaterEqual(val["two_sided_tail_fraction_not_randomization_p"],0.01)
        self.assertLessEqual(val["two_sided_tail_fraction_not_randomization_p"],1)

    def test_cv_groups_are_complete_and_score_is_not_resubstitution(self):
        cv=grouped_out_of_genus_cv(self.strong,self.fams,repeat=2)
        self.assertEqual(cv["repeats"],2)
        self.assertGreater(cv["mean_MSE_baseline"],0)
        self.assertGreater(cv["mean_MSE_plus_disturbance"],0)
        self.assertEqual(cv["heldout_grouping"],"butterfly_genus_5fold")
        self.assertTrue(all(0<=x["repeat"]<2 for x in cv["individual_repeats"]))

if __name__=="__main__":
    unittest.main()
