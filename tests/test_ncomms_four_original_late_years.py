"""Tests for four independently archived original-year late-indexed butterfly cases."""
import copy
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("nc_four_years",ROOT/"scripts/recover_ncomms_four_original_late_year_cases.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class OriginalFiveLinkYearTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.builder=m.load_builder()
        src,pins=cls.builder.load_original()
        cls.grid,cls.islands,cls.species,cls.orig_audit=cls.builder.build_visibility(src,pins)
        cls.chain,cls.receipt=m.load_chain()
        cls.rows,cls.late,cls.audit=m.recover_visibility(
            copy.deepcopy(cls.grid),cls.species,cls.islands,cls.orig_audit,
            cls.chain,cls.receipt
        )

    def test_exact_historical_sha_and_scope(self):
        self.assertEqual(self.receipt["git_blob_sha1"],m.CHAIN_GIT_SHA1)
        self.assertEqual(self.chain["local_original_year_parquet_sha256"],m.SOURCE_OLD_PRIVATE_SHA256)
        self.assertEqual(len(self.chain["rows"]),5)

    def test_exact_original_cases_recovered(self):
        self.assertEqual(len(self.rows),11475)
        self.assertEqual(self.audit["input_visibility_grid_original_positive_2011_2017"],302)
        self.assertEqual(self.audit["all_year_resolved_original_source_endpoints"],10)
        self.assertEqual(self.audit["new_original_source_year_2018_2022_case_IDs"],4)
        self.assertEqual(self.audit["remaining_2018_2022_positive_identity_deficit"],124)
        expected={("Colias croceus","169",2019),("Colias hyale","169",2019),
                  ("Pieris rapae","951",2019),("Nymphalis antiopa","715",2021)}
        self.assertEqual({(z["butterfly"],z["landmass_id"],z["first_index_year"]) for z in self.late},expected)

    def test_early_year_source_agrees_with_already_known_positive_events(self):
        origin={(z["species"],z["landmass_id"]):z for z in self.grid}
        recovered={(z["species"],z["landmass_id"]):z for z in self.rows}
        for key,z in recovered.items():
            if origin[key]["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX":
                self.assertEqual(z["butterfly_first_index_2011_2017_at_THIS_species_island"],"CONFIRMED_FIRST_INDEX")
                self.assertEqual(z["first_index_year_IF_2011_2017_case"],origin[key]["first_index_year_IF_2011_2017_case"])
        self.assertEqual(sum(z["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX" for z in self.rows),302)

    def test_other_unidentified_rows_not_used_as_negatives(self):
        not_recovered=[z for z in self.rows if not z["full_index_2011_2022_at_THIS_species_island"].startswith("CONFIRMED_")]
        self.assertEqual(len(not_recovered),11475-302-4)
        self.assertTrue(all(z["full_index_2011_2022_at_THIS_species_island"].startswith("NOT_RECOVERED_") for z in not_recovered))
        self.assertTrue(all(z["historical_2001_2010_HOSTS_record_at_THIS_species_island"].startswith("NOT_RECOVERED_") for z in self.rows))
        self.assertTrue(self.audit["source_named_prior_host_is_NOT_a_valid_complete_population_exposure"])

    def test_year_conflict_fails_closed(self):
        chain=copy.deepcopy(self.chain)
        target=chain["rows"][0]["endpoints"][0]
        target["first_indexed_year_from_original_year_table"]=2017
        with self.assertRaises(ValueError):
            m.recover_visibility(copy.deepcopy(self.grid),self.species,self.islands,
                                 self.orig_audit,chain,self.receipt)

    def test_duplicate_endpoint_fails_closed(self):
        chain=copy.deepcopy(self.chain)
        chain["rows"][1]["endpoints"][0]["island"]=chain["rows"][0]["endpoints"][0]["island"]
        chain["rows"][1]["butterfly"]=chain["rows"][0]["butterfly"]
        with self.assertRaises(ValueError):
            m.recover_visibility(copy.deepcopy(self.grid),self.species,self.islands,
                                 self.orig_audit,chain,self.receipt)

    def test_missing_real_first_year_fails_closed(self):
        chain=copy.deepcopy(self.chain)
        chain["rows"][0]["endpoints"][0]["first_indexed_year_from_original_year_table"]=None
        with self.assertRaises((ValueError,TypeError)):
            m.recover_visibility(copy.deepcopy(self.grid),self.species,self.islands,
                                 self.orig_audit,chain,self.receipt)


if __name__=="__main__":
    unittest.main()
