#!/usr/bin/env python3
"""Source-checked Zenodo colonial-flora archive structure only.

Lenzner et al. 2022 DOI 10.5281/zenodo.6640264. Full zip is downloaded to
a temporary runner directory. Only relative member paths, file sizes and
qualified source-file roles are persisted. No data matrices rehosted.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json
import hashlib
import zipfile

ZIP=Path("_temp_historical_sources/lenzner_6640264.zip")
EXPECTED_MD5="16f5585d2500676951a4b664bbf0ffde"
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/lenzner_archive_structure.json")


def audit(path=ZIP):
    if not path.exists():
        raise RuntimeError("Missing fixed source file")
    h=hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(2**20),b""):
            h.update(chunk)
    if h.hexdigest()!=EXPECTED_MD5:
        raise RuntimeError("Source MD5 mismatch; halt rather than use unrelated archive")
    with zipfile.ZipFile(path) as z:
        members=[]
        for m in z.infolist():
            lower=m.filename.lower()
            matches=[s for s in ("empire","colon","region","spatial","glonaf",
                                 "distance","duration","shapefile","geo","centroid",
                                 "coordinate","occup","island","table","license",
                                 "readme","data","script")
                     if s in lower]
            members.append({"path":m.filename,"bytes":m.file_size,
                            "compressed_bytes":m.compress_size,
                            "is_directory":m.is_dir(),
                            "matching_roles":matches})
    kinds=Counter(Path(x["path"]).suffix.lower() for x in members if not x["is_directory"])
    return {
        "status":"SOURCE_ARCHIVE_STRUCTURE_ONLY_NOT_EMPIRICAL_COLONIAL_JOIN",
        "doi":"10.5281/zenodo.6640264",
        "md5":EXPECTED_MD5,
        "size_bytes":path.stat().st_size,
        "entries":len(members),
        "extensions":dict(sorted(kinds.items())),
        "files":members,
        "interpretation_limits":[
            "File names do not prove colonial history regional coverage.",
            "Any source data reuse requires inspecting actual terms; Zenodo metadata reports other-open.",
            "GloNAF region geometries must be crosswalked to WGSRPD3 rather than matched by modern country name.",
            "Results must be distinguished from the already-published legacy of colonial history on plant similarity.",
        ],
    }


def main():
    d=audit()
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in d.items() if k!="files"},indent=2))


if __name__=="__main__":
    main()
