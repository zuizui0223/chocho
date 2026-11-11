#!/usr/bin/env python3
"""Build the ONE island-and-port butterfly Research Article as coauthor review DOCX.

This builds the ONE integrated manuscript v0.3, not the frozen global-only GEB
source. The two main figures appear in Results with their
source-audited, noncausal captions. Document is DRAFT; not submitted, journal
typeset, bibliographically verified or coauthor approved.
"""
from __future__ import annotations
import argparse
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/"manuscript/integrated_butterfly_resource_access_20261109"
SOURCE=M/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
CAPTIONS=M/"FIGURE_CAPTIONS_V03.md"
STATUS="DRAFT_FOR_COAUTHOR_REVIEW_NOT_SUBMITTED"

def build(output:Path):
    if not shutil.which("pandoc"):
        raise RuntimeError("pandoc is needed to build the formatted review manuscript")
    original=SOURCE.read_text(encoding="utf8")
    if not original.startswith("# Human-redistributed host plants"):raise ValueError("Wrong canonical manuscript")
    # The DOCX title comes from Pandoc metadata; avoid duplicating the same
    # heading immediately below it in the manuscript body.
    original=original.partition("\n")[2].lstrip("\n")
    if original.count("## Structured abstract")!=1 or original.count("## 3. Results")!=1:
        raise ValueError("Not one canonical manuscript")
    draft=["---",
           "title: 'Human-redistributed host plants disproportionately expand potential larval resources in island regions'",
           "subtitle: 'Single integrated Research Article v0.3 - coauthor review draft, not for submission'",
           "lang: en",
           "---",
           "",
           "> **REVIEW COPY (2026-11-09):** One combined study. No coauthor approval, publisher style or transport/colonisation mechanism is claimed.",
           "",
           original.replace("../../figures/", "figures/").strip(), ""]
    supplement = (M / "TABLE_S1_six_island_cases.md").read_text(encoding="utf8")
    draft.extend(["", supplement])
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        md=Path(tmp)/"review_package.md"
        md.write_text("\n".join(draft),encoding="utf8")
        cmd=["pandoc",str(md),"--from=markdown+tex_math_dollars+link_attributes",
             "--to=docx","--standalone","--resource-path="+str(ROOT),
             "--output="+str(output)]
        p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
        if p.returncode:
            raise RuntimeError(f"Pandoc failed: {p.stderr[-2500:]}")
    with zipfile.ZipFile(output) as z:
        bad=z.testzip()
        if bad:raise ValueError("DOCX media CRC failed: "+bad)
        media=[n for n in z.namelist() if n.startswith("word/media/")]
        if len(media)<15:raise ValueError(f"Only {len(media)} embedded butterfly figures")
        xml=z.read("word/document.xml").decode("utf-8")
        for term in ("Human-redistributed host plants","Structured abstract","Figure 1","Figure 2",
                     "butterfly"):
            if term not in xml:raise ValueError("Missing canonical integrated review manuscript/figure text: "+term)
    receipt={
      "schema":"chocho_single_island_port_butterfly_coauthor_review_DOCX_v03",
      "status":STATUS,
      "source_manuscript":SOURCE.relative_to(ROOT).as_posix(),
      "source_figure_captions":CAPTIONS.relative_to(ROOT).as_posix(),
      "figure_count":2,
      "in_text_case_maps":1,
      "supplementary_figures":12,
      "supplementary_tables":1,
      "embedded_media_files":len(media),
      "not_yet_checked":["Word/PDF full-page visual rendering by LibreOffice",
         "journal formatting compliance","bibliography and data licenses","coauthor approval",
         "historical butterfly host use, established colonisation, ship-mediated transfer"],
    }
    r=output.with_suffix(".source_receipt.json")
    r.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"file":str(output),"bytes":output.stat().st_size,
                      "embedded_figures":len(media),"status":STATUS},indent=2))
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    build(a.output)
