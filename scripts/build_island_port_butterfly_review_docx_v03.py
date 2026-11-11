#!/usr/bin/env python3
"""Build the ONE island-and-port butterfly Research Article as coauthor review DOCX.

This builds the ONE integrated manuscript v0.3, not the frozen global-only GEB
source. Full five figure files are embedded as a review appendix with their
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
from docx import Document

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/"manuscript/integrated_butterfly_resource_access_20261109"
SOURCE=M/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
CAPTIONS=M/"FIGURE_CAPTIONS_V03.md"
F=ROOT/"figures/integrated_island_port_butterfly_v02"
F5=ROOT/"figures/integrated_island_port_butterfly_v03"
FIG_STEMS={
  1:"Fig1_resource_geography_scales_v03",
  2:"Fig2_islands_resource_butterflies_v02",
  3:"Fig3_ports_aliens_known_butterfly_hosts_v02",
  4:"Fig4_butterflies_and_port_host_geography_v02",
  5:"Fig5_hosts_marine_conditional_vs_transfer_v03",
}
STATUS="DRAFT_FOR_COAUTHOR_REVIEW_NOT_SUBMITTED"

def caption_sections():
    text=CAPTIONS.read_text(encoding="utf8")
    chunks=re.split(r"(?m)^## Figure (\d+) \| ",text)
    if len(chunks)!=11:raise ValueError("Expected five full scientific figure captions")
    # Each figure consumes id, title and paragraph until next heading.
    out={}
    for i in range(1,len(chunks),2):
        n=int(chunks[i])
        title_body=chunks[i+1]
        title,_,body=title_body.partition("\n")
        # The source and CI links below the five scientific captions are
        # internal author QA, not a Figure 5 legend. Keep them in the GitHub
        # scientific record, not in the coauthor manuscript appendix.
        body=body.split("\n## Shared scientific source and QA",1)[0].strip()
        out[n]={"title":title.strip(),"body":body}
    if sorted(out)!=list(range(1,6)):raise ValueError("Figure number/order changed")
    return out

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
    caps=caption_sections()
    for i,stem in FIG_STEMS.items():
        p=(F5 if i in (1,5) else F)/(stem+".png")
        if not p.exists():
            raise ValueError(f"Missing source-reviewed Figure {i} image")

    draft=["---",
           "title: 'Human-redistributed host plants amplify butterfly resources on islands: ports and marine geography'",
           "subtitle: 'Single integrated Research Article v0.3 - coauthor review draft, not for submission'",
           "lang: en",
           "---",
           "",
           "> **REVIEW COPY (2026-11-09):** One combined study. No coauthor approval, publisher style or transport/colonisation mechanism is claimed.",
           "",
           original.replace("../../figures/", "figures/").strip(), "",
           "# Five verified scientific figures - review appendix", "",
           "Figure 1 summarizes three scales of host availability and recording. Figures 2–4 retain their source plots. Figure 5 compares marine and host predictors. Captions distinguish their butterfly, botanical and geographic populations. These are *not* biological observations of shipping, host feeding or butterfly arrival.", ""]
    for i,stem in FIG_STEMS.items():
        name=f"Figure {i} | {caps[i]['title']}"
        p=(F5 if i in (1,5) else F)/(stem+".png")
        # The embedded figure is review-only; the five-page vector PDF separately
        # preserves the original vector plot and exact source receipt hashes.
        draft.extend([f"## {name}","",f"![{name}]({p.relative_to(ROOT).as_posix()}){{width=15cm}}",
                      "",caps[i]["body"],""])
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
    # Give the concluding multi-panel Figure 5 a clear page start.
    # Its expanded v0.3 legend may continue to a second page; that is
    # acceptable for a review appendix but MUST be checked visually.
    # Original research figures and their source SHA receipts are unmodified.
    review=Document(output)
    final_headings=[
        p for p in review.paragraphs
        if p.style.name.startswith("Heading")
        and p.text.startswith("Figure 5 |")
    ]
    if len(final_headings)!=1:
        raise ValueError(f"Cannot locate unique Figure 5 heading: {len(final_headings)}")
    final_headings[0].paragraph_format.page_break_before=True
    review.save(output)
    with zipfile.ZipFile(output) as z:
        bad=z.testzip()
        if bad:raise ValueError("DOCX media CRC failed: "+bad)
        media=[n for n in z.namelist() if n.startswith("word/media/")]
        if len(media)<11:raise ValueError(f"Only {len(media)} embedded butterfly figures")
        xml=z.read("word/document.xml").decode("utf-8")
        for term in ("Human-redistributed host plants","Structured abstract","Figure 1","Figure 5",
                     "butterfly"):
            if term not in xml:raise ValueError("Missing canonical integrated review manuscript/figure text: "+term)
    receipt={
      "schema":"chocho_single_island_port_butterfly_coauthor_review_DOCX_v03",
      "status":STATUS,
      "source_manuscript":SOURCE.relative_to(ROOT).as_posix(),
      "source_figure_captions":CAPTIONS.relative_to(ROOT).as_posix(),
      "figure_count":5,
      "in_text_case_maps":1,
      "supplementary_figures":5,
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
