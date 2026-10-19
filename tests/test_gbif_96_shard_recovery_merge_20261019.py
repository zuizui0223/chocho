"""Prevent pseudo-absence, species duplication and GBIF record loss in 4x24 recovery."""
from __future__ import annotations
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from merge_gbif_partial96_shard_receipts_20261019 import merge, SHARDS, SCHEMA, FIELDS


class TestGBIF96ShardMerging(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.taxonomy=self.root/"original_96.csv"
        self.shards=self.root/"shards"
        self.taxonomy.parent.mkdir(exist_ok=True,parents=True)
        with self.taxonomy.open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f);w.writerow(["species","accepted_gbif_species_key"])
            for i in range(96):w.writerow([f"Species_{i:03d}",10000+i])
        for shard,(start,end) in enumerate(SHARDS):
            folder=self.shards/f"gbif-96-shard-{shard}"
            folder.mkdir(exist_ok=True,parents=True)
            receipt=[]
            present=[]
            for i in range(start,end):
                sp=f"Species_{i:03d}"
                complete=i%3==0
                receipt.append({
                    "species":sp,"accepted_species_key":10000+i,
                    "status":"COMPLETE" if complete else "PARTIAL_TRANSPORT_NOT_ABSENCE",
                    "windows":{
                        "2010_2017":{"complete":True},
                        "2018_2025":{"complete":True if complete else False}
                    }
                })
                if complete:present.append([sp,"2010_2017","ABC","HUMAN_OBSERVATION",1])
            self._write_shard(shard,start,end,receipt,present)

    def _write_shard(self,shard,start,end,receipts,counts):
        folder=self.shards/f"gbif-96-shard-{shard}"
        meta={
            "schema":SCHEMA,"species_attempted":end-start,
            "shard_start_index":start,"shard_end_exclusive":end,
            "original_ordinal_size":300,"transport_only_subrequest_size":75,
            "retained_original_source_sampling_protocol":True,
            "species_receipts":receipts,
            "complete_species":sum(r["status"]=="COMPLETE" for r in receipts),
            "not_completed_species":sum(r["status"]!="COMPLETE" for r in receipts),
            "mapped_complete_records":sum(int(r[-1]) for r in counts),
            "regional_basis_count_rows":len(counts)
        }
        (folder/"gbif_96_recovery_receipts.json").write_text(
            json.dumps(meta),encoding="utf-8"
        )
        with (folder/"gbif_96_recovered_region_basis_counts.csv").open(
                "w",encoding="utf-8",newline="") as f:
            writer=csv.writer(f);writer.writerow(FIELDS);writer.writerows(counts)

    def test_all_96_present_no_pseudo_absence(self):
        counts,metadata=merge(self.shards,self.taxonomy)
        self.assertEqual(metadata["species_attempted"],96)
        self.assertEqual(metadata["complete_species"],32)
        self.assertEqual(metadata["not_completed_species"],64)
        self.assertEqual(metadata["mapped_complete_records"],32)
        self.assertEqual(len(counts),32)
        self.assertEqual([r["species"] for r in metadata["species_receipts"]],
                         [f"Species_{i:03d}" for i in range(96)])
        self.assertTrue(metadata["partial_taxa_are_not_absent"])
        self.assertTrue(metadata["no_original_135_records_refetched"])

    def test_unexpected_counts_for_incomplete_species_fail(self):
        p=self.shards/"gbif-96-shard-0"
        with (p/"gbif_96_recovered_region_basis_counts.csv").open("a") as f:
            f.write("Species_001,2010_2017,ABC,HUMAN_OBSERVATION,1\n")
        with self.assertRaisesRegex(ValueError,"partial-species"):
            merge(self.shards,self.taxonomy)

    def test_wrong_shard_cannot_steal_a_species(self):
        p=self.shards/"gbif-96-shard-1"/"gbif_96_recovery_receipts.json"
        data=json.loads(p.read_text())
        data["species_receipts"][0]["species"]="Species_000"
        p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,"non-frozen species"):
            merge(self.shards,self.taxonomy)

    def test_missing_shard_fails_closed(self):
        p=self.shards/"gbif-96-shard-3"/"gbif_96_recovery_receipts.json"
        p.unlink()
        with self.assertRaises(FileNotFoundError):
            merge(self.shards,self.taxonomy)

    def test_source_key_drift_fails_closed(self):
        p=self.shards/"gbif-96-shard-2"/"gbif_96_recovery_receipts.json"
        data=json.loads(p.read_text())
        data["species_receipts"][0]["accepted_species_key"]=1
        p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,"taxon key drift"):
            merge(self.shards,self.taxonomy)


if __name__=="__main__":
    unittest.main()
