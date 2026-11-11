#!/usr/bin/env python3
"""Make the source-verified ONE-PAPER five-figure vector review PDF.

This is not a journal submission package. Figures 1/5 are sourced from
their original source-locked PDFs; Figures 2/3/4 are converted from their
source-locked original SVGs to vector PDF, then joined in figure order.
The SHA256 of each source raster is crosschecked against archived original
render receipts to prevent a misleading set of mixed revisions.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import cairosvg
import fitz

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"figures/integrated_island_port_butterfly_v02"
STEMS=[
 "Fig1_island_port_butterfly_framework_v02",
 "Fig2_islands_resource_butterflies_v02",
 "Fig3_ports_aliens_known_butterfly_hosts_v02",
 "Fig4_butterflies_and_port_host_geography_v02",
 "Fig5_historical_hosts_butterfly_access_falsification_v02",
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def build(outdir):
    a=json.loads((DIR/"FIGURES_234_SOURCE_RECEIPT_V02.json").read_text())
    b=json.loads((DIR/"FIGURES_15_SOURCE_RECEIPT_V02.json").read_text())
    expected="manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md"
    assert a["main_manuscript"]==b["canonical_manuscript"]==expected
    outdir.mkdir(parents=True,exist_ok=True)
    merged=fitz.open();docs=[]
    for i,stem in enumerate(STEMS,1):
        png=DIR/(stem+".png")
        svg=DIR/(stem+".svg")
        receipt=(b["files"] if i in (1,5) else a["visual_files"])
        for path in (png,svg):
            assert path.exists() and sha(path)==receipt[path.name]["sha256"],f"Wrong original rendered figure: {path.name}"
        if i in (1,5):
            pdf=DIR/(stem+".pdf")
            assert pdf.exists() and sha(pdf)==receipt[pdf.name]["sha256"]
        else:
            pdf=outdir/(stem+".pdf")
            cairosvg.svg2pdf(bytestring=svg.read_bytes(),write_to=str(pdf))
        d=fitz.open(pdf)
        if d.page_count!=1 or d[0].rect.width<300:raise ValueError("Invalid 1-page figure "+stem)
        merged.insert_pdf(d)
        d.close()
        docs.append({"figure":i,"source_png_sha256":sha(png),"source_svg_sha256":sha(svg),
                     "vector_pdf_sha256":sha(pdf),"vector_pdf_bytes":pdf.stat().st_size})
    assert merged.page_count==5
    dest=outdir/"ISLAND_PORT_BUTTERFLIES_MAIN_FIGURES_1_TO_5_REVIEW_V02.pdf"
    merged.set_metadata({"title":"Islands and ports reveal introduced host plants and butterfly records - Figures 1–5",
                         "subject":"Scientific review of five noncausal figures; not journal submission"})
    merged.save(dest,garbage=4,deflate=True)
    merged.close()
    receipt={
       "schema":"chocho_island_port_one_paper_five_figure_review_pdf_v02",
       "status":"ALL_FIVE_SOURCE_LOCKED_FIGURES_VISUALLY_READY_FOR_COAUTHOR_REVIEW_NOT_JOURNAL_SUBMISSION",
       "canonical_manuscript":expected,
       "figures":docs,"review_pdf_sha256":sha(dest),
       "page_count":5,"review_pdf_bytes":dest.stat().st_size,
       "source_15_receipt":b["schema"],"source_234_receipt":a["schema"],
       "boundary":["Not an identified plant-introduction or shipborne butterfly colonisation pathway",
          "Actual journal display compliance, bibliography, specimen licensing and coauthor approval are separate gates"],
    }
    path=outdir/"ISLAND_PORT_BUTTERFLIES_FIVE_FIGURES_PDF_RECEIPT_V02.json"
    path.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"pdf":str(dest),"sha256":sha(dest),"pages":5},ensure_ascii=False))
    return dest

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=DIR)
    args=p.parse_args()
    build(args.out_dir)
