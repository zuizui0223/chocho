#!/usr/bin/env python3
"""GEB current-v0.4 coauthor REVIEW build, NOT proof of publication data access."""
import argparse
import json
import re
from pathlib import Path
from zipfile import ZipFile

import cairosvg
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from build_integrated_geb_v04_coauthor_review import add_inline, set_line_numbering

B=Path("manuscript/integrated_butterfly_resource_access_20261109")
M=B/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
S=B/"GEB_SUPPORTING_INFORMATION_V05.md"
F=B/"FIGURE_CAPTIONS_V03.md"
START="## References (DOI identities checked; final source and style audit outstanding)"
FIG=re.compile(r"!\[([^\]]+)\]\(([^)]+\.svg)\)")
KEYWORDS=("alien plants","butterfly–plant interactions","island biogeography",
          "larval host plants","ports","resource geography","spatial prediction")
RUNNING="Island butterfly host resources"

def neutral(text):
    # Remove inaccessible author-owned and repository-relative Markdown links
    # from the anonymous editorial copy. Record the original URLs in the
    # private source audit; NEVER substitute a fictional reviewer DOI or URL.
    text=re.sub(r"\[([^\]]+)\]\(https?://github\.com/zuizui0223/chocho[^)]*\)",
                r"\1",text)
    # Keep ![figure](../../figures/...) intact; only strip ordinary links.
    text=re.sub(r"(?<!!)\[([^\]]+)\]\((?:\.\./)[^)]*\)",r"\1",text)
    return re.sub(r"https?://github\.com/zuizui0223/chocho[^\s)]*",
                  lambda m: "archived source material" +
                  re.search(r"[.,;:!?]*$", m.group(0)).group(0), text)

def journal_si(text):
    # Retain the scientific SI numbers/sections while excluding internal
    # editorial directions that belong only in the Git source-audit record.
    substitutions={
      "**20 November 2026; journal-facing coauthor-review working version, unsubmitted.** Eight consecutively numbered appendices; figures S3.1 and S5.1. The underlying complete old v0.3 source archives remain versioned and separate; this is a source-qualified summary, not a replacement for downloadable data and code.":
      "Eight appendices provide data definitions and sensitivity analyses, with supporting Figures S3.1 and S5.1.",
      "**Submission source gate:** HOSTS and WCVP v13 source DOI/version, GSHHG and original 224-butterfly GBIF DOI are fixed. The *other* 2026 butterfly first-index/plant port-acquisition sources still lack an exact, reviewer-accessible original derived GBIF citation/DOI and some original files have distribution limitations. An anonymous read-only repository/archive containing analysis code and redistributable derived data must be provided at **peer review** and deposited in a stable repository for **publication**. This SI does not magically create those original records or a DOI.":
      "**Data-source limits.** HOSTS, WCVP v13, GSHHG and the original 224-butterfly GBIF download have fixed source identifiers. The separate historical first-index and port-buffer GBIF queries did not preserve original occurrence-level contributing publisher identifiers. Reproduction therefore depends on archived query receipts, permitted derived tables and source-specific access terms; these source limits are distinct from ecological evidence.",
      "These two figures must be embedded as actual media in any uploaded SI DOCX/PDF, not only linked. ":
      "",
      "### S5.3. Crucial ORIGINAL eligible-region sensitivity (PR #183)":
      "### S5.3. Original eligible-region sensitivity",
      "The [complete R01–R37 and PR #181–#187 evidence coverage table](GEB_FULL_EVIDENCE_AUDIT_20261120.md) grades the 37 historically sourced findings by main-text, source-coupled SI and exploratory/archive tiers.":
      "The source-evidence coverage table classifies 37 historically sourced findings into main-text, source-coupled Supporting Information and exploratory/archive tiers.",
    }
    for old,new in substitutions.items():
        if text.count(old)!=1:
            raise RuntimeError("SI editorial export source text drift")
        text=text.replace(old,new)
    return neutral(text)

def docbase(title,main=False):
    d=Document()
    sec=d.sections[0]
    sec.top_margin=Inches(.8)
    sec.bottom_margin=Inches(.8)
    sec.left_margin=Inches(.85)
    sec.right_margin=Inches(.85)
    normal=d.styles["Normal"]
    normal.font.name="Aptos"
    normal.font.size=Pt(10)
    normal.paragraph_format.space_after=Pt(5)
    d.core_properties.title=title
    d.core_properties.author=""
    d.core_properties.last_modified_by=""
    sec.header.paragraphs[0].text="GEB · "+("BLINDED RESEARCH ARTICLE" if main else "SUPPORTING INFORMATION")
    if main:set_line_numbering(d)
    return d

def markdown(doc,txt,root,out,with_images=True,fig_prefix="main"):
    lines=neutral(txt).splitlines()
    imgs=[]
    i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line or line=="---":i+=1;continue
        fm=FIG.fullmatch(line)
        if fm:
            if not with_images:i+=1;continue
            cap,href=fm.groups()
            p=(root/B/href).resolve()
            if not p.is_relative_to(root) or not p.is_file():
                raise RuntimeError("Missing/unsafe figure "+href)
            dest=out/("figure_"+fig_prefix+"_"+str(len(imgs)+1)+".png")
            cairosvg.svg2png(url=str(p),write_to=str(dest),output_width=1500,background_color="white")
            para=doc.add_paragraph()
            para.add_run().add_picture(str(dest),width=Inches(6.6))
            para=doc.add_paragraph(style="Caption")
            add_inline(para,cap)
            imgs.append({"file":str(dest),"source":href,"bytes":dest.stat().st_size})
            i+=1;continue
        head=re.match(r"^(#{1,4}) (.*)",line)
        if head:
            doc.add_heading(head.group(2),level=min(len(head.group(1)),3));i+=1;continue
        if line.startswith("|") and i+1<len(lines) and re.fullmatch(r"[\s:|\-]+",lines[i+1].strip()):
            vals=[x.strip() for x in line.strip("|").split("|")]
            rows=[]
            i+=2
            while i<len(lines) and lines[i].lstrip().startswith("|"):
                z=[x.strip() for x in lines[i].strip().strip("|").split("|")]
                if len(z)==len(vals):rows.append(z)
                i+=1
            table=doc.add_table(rows=1,cols=len(vals))
            table.style="Table Grid"
            # Repeat table headers on every rendered page (w:tblHeader).
            trpr=table.rows[0]._tr.get_or_add_trPr()
            repeat=OxmlElement("w:tblHeader")
            repeat.set(qn("w:val"),"true")
            trpr.append(repeat)
            for j,v in enumerate(vals):
                table.rows[0].cells[j].text=v.replace("**","").replace("`","")
            for row in rows:
                cells=table.add_row().cells
                for j,v in enumerate(row):
                    cells[j].text=v.replace("**","").replace("`","")
                    for p in cells[j].paragraphs:
                        for run in p.runs:run.font.size=Pt(8)
            # Keep each short scientific table row on one page. In particular,
            # an S2 denominator label must stay beside its numerical values.
            for row in table.rows:
                row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
            continue
        if line.startswith("- "):
            pa=doc.add_paragraph(style="List Bullet")
            add_inline(pa,line[2:])
        else:
            pa=doc.add_paragraph()
            add_inline(pa,line)
        i+=1
    return imgs

def build(root,out):
    main=(root/M).read_text(encoding="utf8")
    si=(root/S).read_text(encoding="utf8")
    figcaps=(root/F).read_text(encoding="utf8")
    if main.count(START)!=1 or len(re.findall(r"^## Appendix S[1-8]\.",si,re.M))!=8:
        raise RuntimeError("Current source article or SI structure missing")
    abstract=main.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    a_words=len(re.sub(r"[*#]","",abstract).split())
    mainbody="## 1. Introduction"+main.split("## 1. Introduction",1)[1].split("## Data and code availability",1)[0]
    # Per GEB: explicitly cite each separate Supporting Information appendix.
    appendices="Detailed supporting methods and evidence appear in Appendix S1, Appendix S2, Appendix S3, Appendix S4, Appendix S5, Appendix S6, Appendix S7 and Appendix S8 in Supporting Information.\n"
    if mainbody.count("## 3. Results")!=1:raise RuntimeError("Results section not unique")
    mainbody=mainbody.replace("## 3. Results",appendices+"\n## 3. Results")
    # The legacy working paper's old S.F/S.G pointers are not valid appendix
    # identifiers in the new, separate journal-facing S1–S8 SI.
    converted=[]
    for line in mainbody.splitlines():
        if "Supplement S.F" in line:
            if "species-equal" in line or "contributors were collapsed" in line:
                appendix="Appendix S2"
            elif "source-built" in line or "potential-resource composition" in line:
                appendix="Appendix S4"
            else:
                appendix="Appendix S3"
            line=line.replace("Supplement S.F",appendix)
        if "Supplement S.B–S.D" in line:
            line=line.replace("Supplement S.B–S.D","Appendix S7")
        if "Supplement S.G.2–S.G.5" in line:
            line=line.replace("Supplement S.G.2–S.G.5","Appendix S5")
        converted.append(line)
    mainbody="\n".join(converted)+"\n"
    if re.search(r"Supplement S\.[A-G]",mainbody):
        raise RuntimeError("Legacy source appendix labels remain in GEB independent SI export")
    bodywords=len(re.sub(r"[*#]","",mainbody).split())
    if a_words>300 or bodywords>5000:
        raise RuntimeError("Current source exceeds GEB abstract/main word goals")
    figures=[m.group() for m in FIG.finditer(mainbody)]
    if len(figures)!=4:
        raise RuntimeError("Four actual source main figures are required")
    reflist=[s for s in main.split(START,1)[1].split("**Reference scope:",1)[0].splitlines() if s.startswith("- ")]
    if len(reflist)!=39:
        raise RuntimeError("Unexpected source bibliography size")
    data=main.split("## Data and code availability",1)[1].split("## Figure legends",1)[0]
    data=re.sub(r"^\s*\(draft; final rights and archive check required\)\s*", "", data, flags=re.I)
    data=data.replace("not publicly archived in this Draft", "not included in the present source archive")
    data=data.replace("Original source/data sharing permissions, final citations, rights for individual larval photographs, independent GBIF index DOIs, figure accessibility and coauthor approval must be resolved before any submission. Draft PR #179 remains unmerged.",
                      "Some independent historical GBIF acquisition sources do not retain original occurrence-level publisher identities; original source access and material permissions are described in the accompanying source provenance and Supporting Information.")
    figs_legend=figcaps.split("## Figure source and rights record",1)[0]
    figs_legend=figs_legend[figs_legend.index("## Figure 1 |"):]
    # The archived Figure 4 inset retains its historical S.G footer; give
    # readers the actual location in the separate S1–S8 Supporting Information.
    if figs_legend.count("Supplementary Figure S.G-D") != 1:
        raise RuntimeError("Original full South America figure locator drift")
    figs_legend=figs_legend.replace("Supplementary Figure S.G-D",
        "Figure S5.1 (Appendix S5; called Supplement S.G in the inset)")
    mainclean=(main.splitlines()[0]+"\n\n**Running title:** "+RUNNING+
             "\n\n## Structured abstract"+abstract+"\n**Keywords:** "+
             ", ".join(KEYWORDS)+"\n\n"+mainbody+"\n## References\n"+
             "\n".join(reflist)+"\n\n## Data and Code Availability Statement\n"+
             data+"\n## Figure legends\n"+figs_legend)
    mainclean=neutral(mainclean)
    for forbidden in ("COAUTHOR PREFLIGHT", "Journal submission SOURCE GATE OPEN",
                      "reviewer-access link required", "Draft PR #179"):
        if forbidden in mainclean:
            raise RuntimeError("Internal editorial status leaked into main export: "+forbidden)
    if "zuizui0223" in mainclean.lower():
        raise RuntimeError("Non-blinded author repository URL in clean manuscript")
    out.mkdir(parents=True,exist_ok=True)
    (out/"GEB_v05_blinded_main.md").write_text(mainclean,encoding="utf8")
    siclean=journal_si(si)
    for forbidden in ("COAUTHOR PREFLIGHT","coauthor-review working version",
                      "This SI does not magically", "must be embedded as actual media",
                      "reviewer-access link required", "PR #"):
        if forbidden in siclean:
            raise RuntimeError("Internal editorial status leaked into SI export: "+forbidden)
    (out/"GEB_v05_separate_SI.md").write_text(siclean,encoding="utf8")
    mdoc=docbase("GEB island-port-sea butterfly potential resource opportunities",main=True)
    imagesmain=markdown(mdoc,mainclean,root,out,with_images=False)
    for i,raw in enumerate(figures,1):
        imagesmain.extend(markdown(mdoc,raw,root,out,with_images=True,fig_prefix="main_"+str(i)))
    if len(imagesmain)!=4:raise RuntimeError("Main four SVGs not embedded")
    sf=out/"GEB_v05_BLINDED_MAIN_4figures.docx"
    mdoc.save(sf)
    sdoc=docbase("GEB S1–S8 separate Supporting Information",main=False)
    imagessi=markdown(sdoc,siclean,root,out,with_images=True,fig_prefix="SI")
    if len(imagessi)!=2:raise RuntimeError("GEB SI S3.1/S5.1 must be embedded")
    si_file=out/"GEB_v05_SUPPORTING_INFORMATION_S1-S8.docx"
    sdoc.save(si_file)
    for p,n in [(sf,4),(si_file,2)]:
        with ZipFile(p) as z:
            actual=sum(v.startswith("word/media/") for v in z.namelist())
            if actual!=n:raise RuntimeError("Missing embedded figure in "+p.name)
            if "zuizui0223" in z.read("word/document.xml").decode("utf8").lower():
                raise RuntimeError("Potential author-identity leak in "+p.name)
    status={
       "state":"DOCUMENT_STRUCTURE_VALID; JOURNAL_SUBMISSION_BLOCKED",
       "abstract_words":a_words,"main_body_words":bodywords,
       "references":len(reflist),"running_title":RUNNING,
       "alphabetical_keywords":list(KEYWORDS),
       "main_embedded_figures":len(imagesmain),"SI_appendices":8,
       "SI_embedded_figures":len(imagessi),
       "blocking_gates":[
         "Independent original GBIF/historical and 1596-port acquisition contributor DOI",
         "Reviewer-accessible anonymized archived source data/code (GitHub alone not stable)",
         "Source rights/third-party photo approvals and WPI original conversion provenance",
         "Coauthor-approved separate identifying title page and declarations",
         "APA6 first in-text citation and complete author lists"
       ]
    }
    (out/"GEB_v05_REVIEW_PREFLIGHT.json").write_text(
        json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GEB_V05_BLINDED_MAIN_AND_SEPARATE_SI_PASS",
          json.dumps(status,ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    z=ap.parse_args()
    build(z.repo_root.resolve(),z.outdir)
