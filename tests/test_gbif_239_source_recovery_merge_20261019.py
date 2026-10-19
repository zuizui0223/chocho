"""Unit-tests for the lossless initial-135 + recovery-96 GBIF source merge.

No simulated butterfly observations are interpreted as biological evidence.
Artificial tiny count rows test logic, not ecological associations.
"""
from __future__ import annotations
import sys
import csv
import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from merge_gbif_239_source_with_transport_recovery_20261019 import assemble


def sample_source(recovered_subset=5):
    old=[]
    # Original source exactly 135 complete, 96 transport-partial, 8 unresolved.
    for i in range(239):
        name=f"Butterfly species{i:03}"
        if i<135:
            status="COMPLETE"
            taxon={"key":1000+i}
        elif i<231:
            status="PARTIAL_FETCH_NOT_ELIGIBLE"
            taxon={"key":1000+i}
        else:
            status="TAXON_UNRESOLVED"
            taxon=None
        old.append({
            "species":name,"status":status,"taxon":taxon,
            "windows":{
                "2010_2017":{"complete":status=="COMPLETE"},
                "2018_2025":{"complete":status=="COMPLETE"}
            },
            "errors":[{"msg":"transport failure"}] if status!="COMPLETE" else []
        })
    original={
        "schema":"chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01",
        "species_requested":239,
        "species_taxon_resolved_and_complete":135,
        "species_statuses":{
            "COMPLETE":135,"PARTIAL_FETCH_NOT_ELIGIBLE":96,"TAXON_UNRESOLVED":8
        },
        "species_receipts":old,
        "outputs":{"aggregated_region_basis_rows":1,"source_records_mapped":3},
        "queries":{"windows":[["2010_2017","2010,2017"],["2018_2025","2018,2025"]]}
    }
    recovered=[]
    for i in range(135,231):
        complete=(i-135)<recovered_subset
        recovered.append({
            "species":f"Butterfly species{i:03}",
            "accepted_species_key":1000+i,
            "status":"COMPLETE" if complete else "PARTIAL_TRANSPORT_NOT_ABSENCE",
            "windows":{
                "2010_2017":{"complete":complete},
                "2018_2025":{"complete":complete}
            },
            "errors":[] if complete else [{"msg":"API failed"}],
            "source_validity_counters":{}
        })
    rescue={
        "schema":"chocho_GBIF_96_partial_species_transportreread_20261019_v01",
        "species_attempted":96,
        "original_ordinal_size":300,
        "transport_only_subrequest_size":75,
        "species_receipts":recovered
    }
    basecounts={("Butterfly species000","2010_2017","FRA","HUMAN_OBSERVATION"):3}
    recovercounts={
        (f"Butterfly species{135+i:03}","2018_2025","CAN","PRESERVED_SPECIMEN"):i+1
        for i in range(recovered_subset)
    }
    return original,basecounts,rescue,recovercounts


def run_assemble(raw,old,rec,added):
    # Production merge intentionally requires exact source *file paths*,
    # not Python dicts. Build tiny source-shaped CSV fixtures in a sandbox,
    # then test that the actual file-reading merge enforces its contracts.
    with TemporaryDirectory() as folder:
        root=Path(folder)
        original=root/"initial_135.csv"
        recovered=root/"recovered_96.csv"
        headers=("species","window","wgsrpd3","basis","sampled_record_count")
        for path,values in ((original,old),(recovered,added)):
            with path.open("w",newline="",encoding="utf-8") as fh:
                writer=csv.writer(fh)
                writer.writerow(headers)
                for key,n in sorted(values.items()):
                    writer.writerow([*key,n])
        return assemble(raw,original,rec,recovered)


class TestMerge(unittest.TestCase):
    def test_exact_source_promotion_and_counts(self):
        raw,old,rec,added=sample_source(5)
        result,metadata=run_assemble(raw,old,rec,added)
        self.assertEqual(metadata["species_taxon_resolved_and_complete"],140)
        self.assertEqual(metadata["species_statuses"],{
            "COMPLETE":140,"PARTIAL_FETCH_NOT_ELIGIBLE":91,"TAXON_UNRESOLVED":8
        })
        self.assertEqual(metadata["merge_metadata"]["recovered_new_complete"],5)
        self.assertEqual(metadata["merge_metadata"]["still_transport_incomplete"],91)
        self.assertEqual(metadata["outputs"]["source_records_mapped"],3+sum(range(1,6)))
        self.assertEqual(set(result),set(old)|set(added))
        species={r["species"]:r for r in metadata["species_receipts"]}
        self.assertEqual(species["Butterfly species135"]["status"],"COMPLETE")
        self.assertTrue(species["Butterfly species135"]["originally_incomplete_source_recovered"])
        self.assertEqual(species["Butterfly species140"]["status"],"PARTIAL_FETCH_NOT_ELIGIBLE")
        self.assertEqual(species["Butterfly species231"]["status"],"TAXON_UNRESOLVED")

    def test_reject_wrong_taxon_and_duplicate_record(self):
        raw,old,rec,added=sample_source(1)
        rec["species_receipts"][0]["accepted_species_key"]=123456
        with self.assertRaises(ValueError):
            run_assemble(raw,old,rec,added)
        raw,old,rec,added=sample_source(1)
        added.update(old)
        with self.assertRaises(ValueError):
            run_assemble(raw,old,rec,added)

    def test_reject_partial_recovery_leaking_into_observed(self):
        raw,old,rec,added=sample_source(0)
        added[("Butterfly species135","2018_2025","CAN","HUMAN_OBSERVATION")]=1
        with self.assertRaises(ValueError):
            run_assemble(raw,old,rec,added)

    def test_reject_wrong_recovery_species_set_and_unready_windows(self):
        raw,old,rec,added=sample_source(1)
        rec["species_receipts"][0]["species"]="Wrong butterfly"
        with self.assertRaises(ValueError):
            run_assemble(raw,old,rec,added)
        raw,old,rec,added=sample_source(1)
        rec["species_receipts"][0]["windows"]["2018_2025"]["complete"]=False
        with self.assertRaises(ValueError):
            run_assemble(raw,old,rec,added)

    def test_zero_recovery_is_explicit_missingness(self):
        raw,old,rec,added=sample_source(0)
        merged,metadata=run_assemble(raw,old,rec,added)
        self.assertEqual(metadata["species_statuses"],{
            "COMPLETE":135,"PARTIAL_FETCH_NOT_ELIGIBLE":96,"TAXON_UNRESOLVED":8})
        self.assertEqual(merged,old)
        self.assertTrue(metadata["sampling_is_systematic_and_not_random_or_exhaustive"])


if __name__=="__main__":
    unittest.main()
