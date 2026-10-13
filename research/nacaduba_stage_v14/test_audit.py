"""Regression and overclaim guards for v14 historical and stage synthesis."""
import unittest
from analyze import audit

class StageAudit(unittest.TestCase):
    def setUp(self): self.d=audit()
    def test_time_documents_not_arrival(self):
        self.assertEqual(self.d["documentation_gaps_years"],{"plant_to_retrospective_photo":114,"plant_to_checklist":118,"plant_to_larval_record":120})
        self.assertIsNone(self.d["estimated_biological_arrival_lag_years"])
    def test_stage_duration_sources(self):
        self.assertEqual(self.d["larval_floral_feeding_days"],[9.5,11])
        self.assertEqual(self.d["egg_to_larval_feeding_end_days"],[11.5,13.5])
        self.assertEqual(self.d["egg_to_adult_reared_days"],[17.5,19.5])
    def test_field_host_and_captive_food_distinct(self):
        self.assertEqual(self.d["field_larval_host_plants"],["Acacia auriculiformis","Acacia plectocarpa","Pithecellobium dulce"])
        self.assertIn("Acacia mangium",self.d["not_field_larval_hosts"])
        self.assertIn("Vachellia farnesiana",self.d["not_field_larval_hosts"])
    def test_no_unobserved_biological_claims(self):
        for k in ("field_larval_survival_rate","flower_phenology_effect","colonization_causal_effect"): self.assertIsNone(self.d[k])
    def test_source_sensitivity_is_arithmetic_only(self):
        self.assertEqual(self.d["source_hypothetical_resource_delta"],{"native":4,"contemporary":26,"introduced_only_added":22})

if __name__=="__main__":unittest.main()
