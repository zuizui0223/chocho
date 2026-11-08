"""Synthetic-unit verification of SHA-required original year Parquet importer logic.

These are small, explicitly synthetic test rows, NOT original data or an
empirical recovery. A true 83,573-row source read is intentionally absent from CI.
"""
import copy
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("ncomms_exact_source_importer",ROOT/"scripts/import_ncomms_exact_original_year_parquet.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class SyntheticOriginalYearTests(unittest.TestCase):
    def setUp(self):
        self.rows=[
            {"speciesKey":1,"landmass_id":10,"year":2005,"records":2,"precise1km_records":2,"landmass_area_km2":123.5},
            {"speciesKey":1,"landmass_id":10,"year":2007,"records":3,"precise1km_records":3,"landmass_area_km2":123.5},
            {"speciesKey":1,"landmass_id":20,"year":2019,"records":1,"precise1km_records":1,"landmass_area_km2":144.2},
            {"speciesKey":2,"landmass_id":10,"year":2013,"records":5,"precise1km_records":5,"landmass_area_km2":123.5},
            {"speciesKey":2,"landmass_id":10,"year":2023,"records":3,"precise1km_records":3,"landmass_area_km2":123.5},
        ]
        self.species=[
            {"species":"Species alpha","taxon_key":"1",
             "early_indexed_islands_2001_2010":1,"new_indexed_islands_2011_2022":1},
            {"species":"Species beta","taxon_key":"2",
             "early_indexed_islands_2001_2010":0,"new_indexed_islands_2011_2022":1},
        ]
        self.islands=[
            {"landmass_id":"10","early_butterfly_species_recorded_2001_2010":1,
             "new_butterfly_species_indexed_2011_2022":1},
            {"landmass_id":"20","early_butterfly_species_recorded_2001_2010":0,
             "new_butterfly_species_indexed_2011_2022":1},
        ]
        self.early={("Species beta","10"):2013}
        self.late={("Species alpha","20"):2019}

    def test_source_synthetic_basic(self):
        q=m.source_spec_check(self.rows,full_original=False)
        self.assertEqual(q["source_rows"],5)
        self.assertEqual(q["occurrence_records_sum"],14)
        self.assertEqual(q["species"],2)
        self.assertEqual(q["island_groups"],2)

    def test_first_year_rows_and_missing_are_distinct(self):
        r,s=m.build_original_visibility(self.rows,self.species,self.islands,self.early,self.late,require_original_population=False)
        self.assertEqual(len(r),4)
        t={(x["species"],x["landmass_id"]):x for x in r}
        self.assertEqual(t["Species alpha","10"]["original_first_indexed_year"],2005)
        self.assertEqual(t["Species alpha","20"]["original_first_indexed_year"],2019)
        self.assertEqual(t["Species beta","10"]["original_first_indexed_year"],2013)
        self.assertEqual(t["Species beta","20"]["original_first_indexed_year"],"")
        self.assertIn("NOT_ECOLOGICAL_ABSENCE",t["Species beta","20"]["not_indexed_2001_2024_state"])
        self.assertEqual(s["source_2001_2010_indexed"],1)
        self.assertEqual(s["source_2011_2017_first_indexed"],1)
        self.assertEqual(s["source_2018_2022_first_indexed"],1)

    def test_conflict_2011_2017_source_rejected(self):
        early={("Species beta","10"):2012}
        with self.assertRaises(ValueError):
            m.build_original_visibility(self.rows,self.species,self.islands,early,self.late,require_original_population=False)

    def test_conflict_four_verified_late_sources_rejected(self):
        late={("Species alpha","20"):2020}
        with self.assertRaises(ValueError):
            m.build_original_visibility(self.rows,self.species,self.islands,self.early,late,require_original_population=False)

    def test_duplicate_input_rows_rejected(self):
        with self.assertRaises(ValueError):
            m.source_spec_check(self.rows+[self.rows[0]],full_original=False)

    def test_conflicting_precise_records_rejected(self):
        rows=copy.deepcopy(self.rows)
        rows[1]["precise1km_records"]=50
        with self.assertRaises(ValueError):
            m.source_spec_check(rows,full_original=False)

    def test_fake_sha_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"original.parquet"
            p.write_bytes(b"not a real Parquet")
            self.assertNotEqual(m.sha256_file(p),m.ORIGINAL_SOURCE_SHA256)

    def test_full_original_cannot_be_spoofed_by_small_example(self):
        with self.assertRaises(ValueError):
            m.source_spec_check(self.rows,full_original=True)


if __name__=="__main__":
    unittest.main()
