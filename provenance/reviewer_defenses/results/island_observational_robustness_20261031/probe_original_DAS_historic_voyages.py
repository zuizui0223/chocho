#!/usr/bin/env python3
"""Audit documented original Dutch Asiatic Shipping VOC voyage file identity.

Source: Huygens Dutch Asiatic Shipping, downloaded/normalized by
IBM/chuk-mcp-maritime-archives scripts/download_das.py.
Frozen IBM repository at 33ae4ef8dfa6acd5f35e619e079ceec518306c9f,
data/voyages.json blob f7a55fcc91fec1e845c955cebd3355eb5bfa8a58.
This pilot DOES NOT pretend it is the IISG 57,669 Asia port-to-port journeys:
DAS 8,194 voyage records are a distinct source with different temporal and
geographical coverage, and may contain coarse Netherlands-to-Asia voyages.

The original individual vessel records and ship names are not republished.
"""
from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path

SOURCE=Path("_historical_DAS/data/voyages.json")
OUT=Path("provenance/reviewer_defenses/results/island_observational_robustness_20261031/"
         "das_1595_1795_source_schema_audit.json")
SOURCE_COMMIT="33ae4ef8dfa6acd5f35e619e079ceec518306c9f"
SOURCE_BLOB="f7a55fcc91fec1e845c955cebd3355eb5bfa8a58"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(path=SOURCE):
    raw=json.loads(path.read_text(encoding="utf8"))
    if isinstance(raw,dict) and isinstance(raw.get("voyages"),list):
        voyages=raw["voyages"]
    elif isinstance(raw,list):
        voyages=raw
    else:
        raise RuntimeError("Unknown DAS normalized voyage source shape")
    if not (7000<=len(voyages)<=10000):
        raise RuntimeError("Wrong DAS voyage source count; no trade network inferred")
    if not all(isinstance(x,dict) for x in voyages):
        raise RuntimeError("Expected structured voyage dictionaries")
    columns=sorted({k for x in voyages for k in x})
    profiles={}
    for k in columns:
        present=[x[k] for x in voyages if x.get(k) not in (None,"",[],{})]
        values=[str(v) for v in present if not isinstance(v,(list,dict))]
        counts=Counter(values)
        profile={
            "row_coverage":len(present),
            "scalar_distinct_values":len(counts),
            "non_scalar_values":len(present)-len(values),
        }
        low=k.lower()
        if any(token in low for token in ("port","depart","arrival","destination","route","call","place","origin")) and not any(t in low for t in ("date","ship","captain","vessel","person")):
            profile["frequent_source_place_names"]=[
                {"name":name,"n_rows":count}
                for name,count in counts.most_common(20)
                if count>=2
            ]
        profiles[k]=profile
    # Port destination/departure PLACE name vocabulary is public, low
    # cardinality and needed for auditable historical gazetteer matching.
    # Do not release names of ships or individuals.
    port_vocabulary={
        k:sorted({str(r[k]).strip() for r in voyages if r.get(k)})
        for k in ("departure_port","destination_port")
    }
    return {
        "historical_port_name_vocabulary":port_vocabulary,
        "status":"HISTORICAL_DAS_SOURCE_SHAPE_ONLY_NOT_PORT_PAIR_MATCH",
        "source_citation":"Dutch Asiatic Shipping (DAS), Huygens Institute, normalized in IBM/chuk-mcp-maritime-archives",
        "source_URL":"https://resources.huygens.knaw.nl/das",
        "derived_repository":"IBM/chuk-mcp-maritime-archives",
        "normalized_repo_commit_SHA":SOURCE_COMMIT,
        "normalized_original_Git_blob_SHA":SOURCE_BLOB,
        "source_file_SHA256":sha256(path),
        "record_count":len(voyages),
        "source_keys":columns,
        "column_profiles":profiles,
        "limits":[
            "The Dutch Asiatic Shipping database is Dutch-only, not a global historic maritime network.",
            "The IBM normalized source may not preserve every intra-Asian stop; check fields before converting to port-to-port edges.",
            "Historic place names and coordinates must be matched to original WPI island locations without ambiguity.",
            "Current port name and 2019 coordinate do not guarantee historic port continuous operation.",
            "The published voyages were not sampled based on butterfly and host outcomes.",
            "No original cargo entries or transported biological taxa are identified by this schema.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,default=SOURCE)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    result=audit(args.input)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="column_profiles"},indent=2))


if __name__=="__main__":
    main()
