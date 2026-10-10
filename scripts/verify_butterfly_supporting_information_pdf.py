#!/usr/bin/env python3
"""Post-layout QC for the editable GEB Supporting Information Word artifact.

Convert its DOCX with LibreOffice before calling this checker. This checker
compares PDF extracted text to the source's required headings/key results and
reports the pages present. It does NOT substitute for visual human inspection.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


def call(*args: str) -> str:
    return subprocess.check_output(args, text=True, encoding="utf-8", errors="replace", timeout=45)


def inspect(pdf: Path) -> dict:
    if not pdf.is_file() or pdf.stat().st_size < 10_000:
        raise RuntimeError("Supporting PDF is missing or implausibly small")
    info=call("pdfinfo", str(pdf))
    pages_match=re.search(r"(?m)^Pages:\s+(\d+)$", info)
    if not pages_match:
        raise RuntimeError("Cannot parse Supporting Information PDF page count")
    pages=int(pages_match.group(1))
    if not 3 <= pages <= 80:
        raise RuntimeError(f"Unlikely supplementary page count: {pages}")
    content=call("pdftotext", "-layout", str(pdf), "-")
    # Ligatures, non-breaking spaces and ordinary line breaks may vary between
    # document renderers; normalize whitespace but not scientific values.
    text=re.sub(r"\s+", " ", content).lower()
    required=[
        "supplementary information",
        *[f"supplementary table s{k}." for k in range(1,10)],
        "host-interaction knowledge sensitivity",
        "0.003616",
        "0.003436",
        "50.1%",
        "54.9%",
        "1,027",
        "competition",
    ]
    missing=[phrase for phrase in required if phrase.lower() not in text]
    if missing:
        raise RuntimeError("Page-rendered SI is missing expected scientific text: "+repr(missing))
    # A separate binary check catches corrupt output with duplicate/blank pages.
    page_texts=content.split("\f")
    meaningful=[z for z in page_texts if len(re.sub(r"\s+","",z))>15]
    if len(meaningful)!=pages:
        raise RuntimeError(f"PDF contains apparently blank pages: {len(meaningful)}/{pages}")
    return {
        "schema":"chocho_geb_supporting_info_pdf_layout_text_check_v0.1",
        "status":"RENDERED_TEXT_COMPLETE_NOT_MANUAL_VISUAL_APPROVAL",
        "pages":pages,
        "pdf_bytes":pdf.stat().st_size,
        "expected_table_headings":9,
        "critical_source_numbers_present":True,
        "no_blank_pages_by_text_length":True,
        "limitation":"This test cannot reliably detect clipped table columns, overlapping glyphs, or visual page quality; visual inspection of the page-rendered PDF remains necessary.",
    }


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    result=inspect(args.pdf)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
