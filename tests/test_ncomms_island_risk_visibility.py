#!/usr/bin/env python3
"""Regression tests for original 225x51 temporal risk visibility, stdlib only."""
import copy
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("nc_visibility",ROOT/"scripts/build_ncomms_island_risk_visibility.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class VisibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs,cls.pins=m.load_original()
        cls.grid,cls.islands,cls.species,cls.audit=m.build_visibility(cls.inputs,cls.pins)

    def test_input_sha_and_size(self):
        self.assertEqual(len(self.inputs),6)
        self.assertEqual(set(self.inputs),set(m.FILES))
        for key,pin in self.pins.items():
            self.assertEqual(pin["git_blob_sha1"],m.FILES[key][1])

    def test_full_rectangular_visibility_not_risk_cohort(self):
        self.assertEqual(len(self.grid),11475)
        self.assertEqual(len(self.islands),51)
        self.assertEqual(len(self.species),225)
        self.assertEqual(len({(z["species"],z["landmass_id"]) for z in self.grid}),11475)
        self.assertTrue(all(z["source_scope"].endswith("NOT_A_RISK_PREDICTION_TABLE") for z in self.grid))

    def test_exact_early_late_aggregates(self):
        r=self.audit
        self.assertEqual(r["published_pre2011_butterfly_indexed_species_island_cells"],524)
        self.assertEqual(r["published_first2011_2022_indexed_species_island_cells"],430)
        self.assertEqual(r["archived_row_identified_first2011_2017_species_island_cells"],302)
        self.assertEqual(r["unassigned_2018_2022_positive_identities_but_known_marginal_total"],128)
        self.assertEqual(r["all_cells_not_indexed_in_original_butterfly_2001_2010_by_count"],10951)
        self.assertEqual(r["all_cells_without_butterfly_index_in_original_2001_2022_by_count"],10521)
        self.assertEqual(sum(z["unassigned_new_species_2018_2022"] for z in self.islands),128)
        self.assertEqual(sum(z["unassigned_new_islands_2018_2022"] for z in self.species),128)

    def test_positive_row_known_nonpositives_are_undetermined(self):
        known=[z for z in self.grid if z["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX"]
        unknown=[z for z in self.grid if z["butterfly_first_index_2011_2017_at_THIS_species_island"]!="CONFIRMED_FIRST_INDEX"]
        self.assertEqual((len(known),len(unknown)),(302,11173))
        self.assertTrue(all(z["butterfly_pre2011_index_at_THIS_species_island"]=="KNOWN_NOT_PRE2011" for z in known))
        self.assertTrue(all(z["butterfly_pre2011_index_at_THIS_species_island"].startswith("NOT_RECOVERED") for z in unknown))
        self.assertTrue(all(z["full_index_2011_2022_at_THIS_species_island"].startswith("NOT_RECOVERED") for z in unknown))
        self.assertTrue(all(z["historical_2001_2010_HOSTS_record_at_THIS_species_island"].startswith("NOT_RECOVERED") for z in self.grid))
        self.assertEqual(self.audit["first2011_2017_species"],106)
        self.assertEqual(self.audit["first2011_2017_islands"],48)

    def test_source_17_island_case_null_is_preexisting_and_not_causal(self):
        s=self.audit["original_17island_prehost_risk_subset_aggregate"]
        self.assertEqual((s["population_islands"],s["population_butterfly_species"]),(17,60))
        early=s["periods"]["2011_2017_first_index_among_no_pre2011_butterfly_index"]
        late=s["periods"]["2018_2022_first_index_among_still_unindexed_after_2017"]
        self.assertEqual((early["risk_cells"],early["indexed_first_events"],early["historical_host_indexed_first_events"]),(905,113,34))
        self.assertEqual((late["risk_cells"],late["indexed_first_events"],late["historical_host_indexed_first_events"]),(792,33,20))
        self.assertAlmostEqual(early["conditional_margin_null_mean_matched_events"],33.91983967935872)
        self.assertAlmostEqual(late["conditional_margin_null_mean_matched_events"],19.50367401469606)
        self.assertTrue(early["NOT_CAUSAL"])
        self.assertTrue(late["NONINDEXED_PLANT_IS_NOT_ABSENT"])

    def test_no_event_leakage_in_nonrecord_cells(self):
        # The observed 2011–2017 event cohort is data-derived but non-event
        # cells are deliberately not encoded as negative outcomes.
        for z in self.grid:
            if z["first_index_year_IF_2011_2017_case"]=="":
                self.assertEqual(z["modern_host_group_IF_2011_2017_case"],"NOT_AN_INDEXED_2011_2017_EVENT_ROW")

    def test_malformed_source_key_rejected(self):
        data=copy.deepcopy(self.inputs)
        data["event_2011_2017"][0]["first_record_year"]="2022"
        with self.assertRaises(ValueError):
            m.build_visibility(data,self.pins)

    def test_manipulated_species_margin_rejected(self):
        data=copy.deepcopy(self.inputs)
        row=next(z for z in data["species_roster"] if int(z["new_records_on_eligible_islands"])>0)
        row["new_records_on_eligible_islands"]="0"
        with self.assertRaises(ValueError):
            m.build_visibility(data,self.pins)

    def test_modified_17island_null_rejected(self):
        data=copy.deepcopy(self.inputs)
        item=data["original_17island_risk"]["periods"]["2011_2017_first_index_among_no_pre2011_butterfly_index"]["geographic_margin_null"]
        item["null_mean_matched_events"]=40.0
        with self.assertRaises(ValueError):
            m.build_visibility(data,self.pins)

    def test_technical_integer_guard(self):
        for bad in ["-1","notnumber","2.5","",None]:
            with self.assertRaises(ValueError):
                m.get_positive_int(bad,"test")
        self.assertEqual(m.get_positive_int("0","zero"),0)


if __name__=="__main__":
    unittest.main()
