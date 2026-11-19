#!/usr/bin/env python3
"""Transcription-only check of Fall & Drezner (2011) Table 1, Tonga.

The source's 'water' category must NOT be renamed seawater thalassochory.
This is not a new 2020 Tonga species-level analysis or independent replication
of EuDiS v5 native-minus-introduced 22-island thalassochory result.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"results/tropical_dispersal_source_bridge_20261118/TONGA2011_PUBLISHED_TABLE1_COUNTS.json"

def main():
    d=json.loads(DATA.read_text(encoding="utf-8"))
    assert d["source"]["doi"]=="10.2984/65.2.143"
    rows={r["key"]:r for r in d["origin_classes"]}
    expected={"native":(364,137),"polynesian_introduction":(64,35),
              "post_european_introduction":(309,137)}
    for key,(n,w) in expected.items():
        assert (rows[key]["total_species"],rows[key]["water_guild_species"])==(n,w)
    props={key:100*w/n for key,(n,w) in expected.items()}
    diff=props["native"]-props["post_european_introduction"]
    assert -6.71<diff<-6.69
    assert "not specific EuDiS" in d["traits"],"Protect distinct trait definitions"
    print(json.dumps({
        "status":"PASS",
        "published_aggregate_source":"Fall and Drezner 2011, Table 1",
        "trait":"multiple-guild broad water dispersal, NOT thalassochory",
        "percent_by_origin":props,
        "native_minus_post_European_pp":diff,
        "scope":"one archipelago, 2011 published aggregates, not global matched WCVP",
        "2020_1020_species_source_downloaded":False},
        ensure_ascii=False,sort_keys=True,indent=2))

if __name__=="__main__":main()
