#!/usr/bin/env python3
"""Build an actual v0.4 coauthor-review DOCX/PDF-source image bundle.

This differs from archived GEB blinded DOCX CI, which renders the unrelated
v0.2 butterfly specialization manuscript. This consumes the live one-paper
v0.4 main text and its four verified SVGs, not source-archive supplement text.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import cairosvg
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

REL="manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
SPLIT="## Supplementary results: island record histories and maritime comparisons"
FIGRE=r"!\[Figure ([1-4])\. ([^\]]+)\]\(([^)]+\.svg)\)"
TOKEN=re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|" + r"\x60[^\x60]+\x60|" + r"\[[^\]]+\]\([^)]+\))")

def add_inline(p, raw):
    cursor=0
    for match in TOKEN.finditer(raw):
        if match.start()>cursor:
            p.add_run(raw[cursor:match.start()])
        tok=match.group()
        if tok.startswith("**"):
            run=p.add_run(tok[2:-2]);run.bold=True
        elif tok.startswith("*"):
            run=p.add_run(tok[1:-1]);run.italic=True
        elif tok.startswith("\x60"):
            run=p.add_run(tok[1:-1])
            run.font.name="Courier New"
        elif tok.startswith("["):
            m=re.match(r"\[([^\]]+)\]\(([^)]+)\)",tok)
            if m:
                r=p.add_run(m.group(1))
                if m.group(2).startswith(("https://","http://")):
                    r.font.color.rgb=RGBColor(29,89,124)
        cursor=match.end()
    if cursor<len(raw):
        p.add_run(raw[cursor:])

def set_line_numbering(doc):
    sect=doc.sections[0]
    secpr=sect._sectPr
    ln=OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"),"1")
    ln.set(qn("w:start"),"0")
    ln.set(qn("w:restart"),"continuous")
    secpr.append(ln)

def add_doc_paragraph(doc,line):
    match=re.match(r"^(#{1,4})\s+(.*)$",line)
    if match:
        level=min(len(match.group(1)),3)
        doc.add_heading(match.group(2),level=level)
        return
    if line.startswith("- "):
        p=doc.add_paragraph(style="List Bullet")
        add_inline(p,line[2:])
        return
    if line.startswith("> "):
        p=doc.add_paragraph()
        p.style="Quote"
        add_inline(p,line[2:])
        return
    p=doc.add_paragraph()
    add_inline(p,line)

def build(args):
    source=args.repo_root.resolve()
    text=(source/REL).read_text(encoding="utf8")
    if SPLIT not in text:
        raise RuntimeError("Cannot locate v0.4 original supplementary boundary")
    main=text.split(SPLIT,1)[0]
    if main.count("## Structured abstract")!=1 or main.count("## 1. Introduction")!=1:
        raise RuntimeError("Source article structure invalid")
    if not main.startswith("# Including introduced plants weakens the island isolation gradient"):
        raise RuntimeError("Not the integrated island-port-sea v0.4 manuscript")
    figs=list(re.finditer(FIGRE,main))
    if len(figs)!=4 or [m.group(1) for m in figs]!=list("1234"):
        raise RuntimeError("Main article must contain exactly four correct figure links")
    args.output_docx.parent.mkdir(parents=True,exist_ok=True)
    assetdir=args.output_docx.parent/"assets"
    assetdir.mkdir(exist_ok=True)
    doc=Document()
    sec=doc.sections[0]
    sec.top_margin=Inches(.85)
    sec.bottom_margin=Inches(.8)
    sec.left_margin=Inches(.85)
    sec.right_margin=Inches(.8)
    doc.styles["Normal"].font.name="Aptos"
    doc.styles["Normal"].font.size=Pt(10)
    doc.styles["Normal"].paragraph_format.space_after=Pt(5)
    doc.styles["Normal"].paragraph_format.line_spacing=1.13
    for style in ["Heading 1","Heading 2","Heading 3"]:
        doc.styles[style].font.name="Aptos Display"
        doc.styles[style].font.color.rgb=RGBColor(26,58,80)
    doc.styles["Heading 1"].font.size=Pt(15)
    doc.styles["Heading 2"].font.size=Pt(12)
    doc.styles["Heading 3"].font.size=Pt(11)
    header=sec.header.paragraphs[0]
    header.text="GEB · ISLAND BUTTERFLY HOST RESOURCES · COAUTHOR REVIEW DRAFT · NOT SUBMITTED"
    header.style="Caption"
    footer=sec.footer.paragraphs[0]
    footer.alignment=WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run("Source v0.4 · 19 Nov 2026 · page ")
    fld=OxmlElement("w:fldSimple");fld.set(qn("w:instr"),"PAGE");footer._p.append(fld)
    set_line_numbering(doc)
    doc.core_properties.title=main.splitlines()[0].removeprefix("# ")
    doc.core_properties.subject="Unsubmitted GEB coauthor review draft (v0.4)"
    doc.core_properties.author=""
    doc.core_properties.last_modified_by=""
    images=[]
    rendered=[]
    for raw in main.splitlines():
        line=raw.strip()
        if not line:
            continue
        if line.startswith("**ONE integrated Research Article") or line.startswith("**Status:**"):
            add_doc_paragraph(doc,line);continue
        fm=re.fullmatch(FIGRE,line)
        if fm:
            n,title,link=fm.groups()
            orig=(source/"manuscript/integrated_butterfly_resource_access_20261109"/link).resolve()
            if not orig.is_relative_to(source) or not orig.is_file():
                raise RuntimeError("Invalid main figure path: "+str(orig))
            png=assetdir/("Figure_"+n+".png")
            cairosvg.svg2png(url=str(orig),write_to=str(png),output_width=1700,background_color="white")
            rendered.append(png)
            p=doc.add_paragraph()
            p.alignment=WD_ALIGN_PARAGRAPH.CENTER
            available=sec.page_width-sec.left_margin-sec.right_margin
            p.add_run().add_picture(str(png),width=available)
            cap=doc.add_paragraph("Figure "+n+". "+title,style="Caption")
            cap.alignment=WD_ALIGN_PARAGRAPH.CENTER
            images.append({"figure":int(n),"svg":str(orig.relative_to(source)),"png":str(png.name),"size_bytes":png.stat().st_size})
            continue
        # Markdown horizontal rules and table pipe rows should not become prose.
        if line=="---":
            continue
        add_doc_paragraph(doc,line)
    doc.save(args.output_docx)
    if args.output_docx.stat().st_size<30000:
        raise RuntimeError("Generated review DOCX implausibly small")
    if len(images)!=4:
        raise RuntimeError("Actual integrated manuscript did not embed four main figures")
    from zipfile import ZipFile
    with ZipFile(args.output_docx) as z:
        text_xml=z.read("word/document.xml").decode("utf8")
        names=z.namelist()
        embedded=[name for name in names if name.startswith("word/media/")]
        if len(embedded)!=4:
            raise RuntimeError("DOCX does not physically embed four figure images")
        for phrase in ("82 of 97","54.3%","40%","11%","0.29","Pieris cheiranthi","0.36"):
            if phrase not in text_xml:
                raise RuntimeError("Review DOCX missing source-verified positive: "+phrase)
    receipt={
        "status":"COAUTHOR_REVIEW_DOCX_BUILT_NOT_SUBMITTED",
        "actual_source_manuscript":REL,
        "source_main_excludes_archived_S_A_to_F_by_explicit_section_boundary":True,
        "main_markdown_words_approx":len(main.split()),
        "figures":images,
        "docx_bytes":args.output_docx.stat().st_size,
        "source_note":"Old v0.2 blinded-review workflow is not a v0.4 article visual audit"
    }
    (args.output_docx.parent/"review_build_receipt.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GEB_V04_CURRENT_ARTICLE_DOCX_BUILT",json.dumps(receipt,ensure_ascii=False),flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--output-docx",type=Path,required=True)
    build(ap.parse_args())

if __name__=="__main__":
    main()
