#!/usr/bin/env python3
"""Export the current editable manuscript Markdown and regenerated PNGs to Word.

This is an editable, content-checked preview; original author-local pagination is
not guaranteed. No access to private source records is required.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from docx.table import Table
from docx.text.paragraph import Paragraph

ROOT=Path(__file__).resolve().parent
WIDTHS=[[2,1.9,2.9],[2.9,1.95,1.95],[1.25,.73,.87,1.12,.85,.96,1.02],
        [2.85,1.97,1.98],[2.4,2.2,2.2],[2,.79,.77,.79,.79,.81,.85],
        [2.35,.64,.89,1.37,1.55],[3.47,.65,.9,1.78],[2.27,2.27,2.26]]
IMAGE=re.compile(r'^!\[Figure ([\w.]+)\]\(([\w]+\.png)\)$')
EMPHASIS=re.compile(r'(\*\*|\*)([^*]+)\1')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def plain(text):return EMPHASIS.sub(lambda m:m[2],text)

def add_runs(p,text):
    pos=0
    for m in EMPHASIS.finditer(text):
        if m.start()>pos:p.add_run(text[pos:m.start()])
        r=p.add_run(m[2]);r.bold=m[1]=='**';r.italic=m[1]=='*';pos=m.end()
    if pos<len(text):p.add_run(text[pos:])
    assert p.text==plain(text)

def parse(path):
    out=[]
    for b in path.read_text(encoding='utf-8').strip().split('\n\n'):
        if (m:=IMAGE.fullmatch(b)):
            out.append(('figure',m[1],m[2]))
        elif b.startswith('| '):
            rows=[line[2:-2].split(' | ') for line in b.splitlines()]
            assert rows[1]==['---']*len(rows[0]) and all(len(r)==len(rows[0]) for r in rows)
            out.append(('table',[rows[0],*rows[2:]]))
        elif (m:=re.match(r'^(#{1,3}) (.+)$',b)):
            out.append(('paragraph',f'Heading {len(m[1])}',m[2]))
        else:
            assert '\n' not in b, f'Unexpected multiline block: {b[:60]}'
            style='Caption' if b.startswith('Figure ') else 'Supplementary Table Caption' if b.startswith('**Table S') else 'Normal'
            out.append(('paragraph',style,b))
    return out

def configure(doc,main):
    sec=doc.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=Inches(.8);sec.left_margin=sec.right_margin=Inches(.85)
    sec.header_distance=Inches(.36);sec.footer_distance=Inches(.37)
    for name,size,spacing,before,after in [('Normal',12,1.5,0,6),('Heading 1',16,1.15,0,10),('Heading 2',14,1.15,12,5),('Heading 3',12.5,1.15,9,4),('Caption',10.5,1.05,3,0),('Supplementary Table Caption',10.5,1.15,6,4)]:
        if name not in doc.styles:doc.styles.add_style(name,WD_STYLE_TYPE.PARAGRAPH)
        st=doc.styles[name];st.font.name='Times New Roman';st.font.size=Pt(size);st.font.color.rgb=RGBColor(0,0,0);st.font.italic=False;st.font.bold=name.startswith('Heading')
        pf=st.paragraph_format;pf.line_spacing=spacing;pf.space_before=Pt(before);pf.space_after=Pt(after);pf.keep_with_next=name.startswith('Heading') or name=='Supplementary Table Caption'
    p=sec.footer.paragraphs[0];p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    f=OxmlElement('w:fldSimple');f.set(qn('w:instr'),'PAGE');p._p.append(f)
    if main:
        line=OxmlElement('w:lnNumType');line.set(qn('w:countBy'),'1');line.set(qn('w:start'),'0');line.set(qn('w:restart'),'continuous');sec._sectPr.insert_element_before(line,'w:pgNumType','w:cols','w:docGrid')
    doc.core_properties.author='';doc.core_properties.last_modified_by='';doc.core_properties.title='';doc.core_properties.subject='';doc.core_properties.comments=''

def expected_blocks(blocks):
    result=[]
    for b in blocks:
        if b[0]=='figure':result.append(('figure',))
        elif b[0]=='table':result.append(('table',[[plain(c) for c in r] for r in b[1]]))
        else:result.append(('paragraph',plain(b[2])))
    return result

def actual_blocks(doc):
    out=[]
    for e in doc.element.body:
        if e.tag==qn('w:p'):
            p=Paragraph(e,doc)
            if p.text:out.append(('paragraph',p.text))
            elif p._element.xpath('.//w:drawing'):out.append(('figure',))
        elif e.tag==qn('w:tbl'):
            t=Table(e,doc);out.append(('table',[[c.text for c in r.cells] for r in t.rows]))
    return out

def build(source,out,main):
    doc=Document();configure(doc,main);blocks=parse(source);nfig=ntab=npara=0;ref=False
    figures=[]
    for b in blocks:
        if b[0]=='figure':
            _,number,image=b;p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.keep_with_next=True;p.paragraph_format.page_break_before=main;p.paragraph_format.space_after=Pt(5)
            imagepath=source.parent/image;p.add_run().add_picture(str(imagepath),width=Inches(6.55 if main else 6.45));nfig+=1;figures.append({'number':number,'sha256':sha(imagepath)})
        elif b[0]=='table':
            rows=b[1];widths=WIDTHS[ntab];assert len(widths)==len(rows[0]);ntab+=1
            table=doc.add_table(rows=len(rows),cols=len(rows[0]));table.style='Table Grid';table.autofit=False
            for c,w in zip(table.columns,widths):c.width=Inches(w)
            for i,row in enumerate(table.rows):
                pr=row._tr.get_or_add_trPr();pr.append(OxmlElement('w:cantSplit'))
                if i==0:pr.append(OxmlElement('w:tblHeader'))
                for j,cell in enumerate(row.cells):
                    cell.width=Inches(widths[j]);p=cell.paragraphs[0];add_runs(p,rows[i][j]);p.paragraph_format.line_spacing=1.1;p.paragraph_format.space_after=Pt(3);p.paragraph_format.space_before=Pt(2)
                    for r in p.runs:
                        r.font.name='Times New Roman';r.font.size=Pt(10 if len(widths)>=7 else 10.5)
                        if i==0:r.bold=True
            # No additional empty paragraph: source block order is preserved.
        else:
            _,style,text=b;p=doc.add_paragraph(style=style);add_runs(p,text);npara+=1
            if text in ('References','Supporting Information references and source citations'):ref=True
            elif style.startswith('Heading'):ref=False
            elif ref and style=='Normal' and not text.startswith('The main text contains'):
                pf=p.paragraph_format;pf.left_indent=Inches(.25);pf.first_line_indent=Inches(-.25);pf.line_spacing=1.15;pf.space_after=Pt(5)
    expected=(116,5,0) if main else (155,3,9)
    # Counts describe the frozen initial text; future intentional changes require a corresponding version update.
    assert (npara,nfig,ntab)==expected, (source.name,npara,nfig,ntab,expected)
    out.parent.mkdir(parents=True,exist_ok=True);doc.save(out);reopened=Document(out)
    assert actual_blocks(reopened)==expected_blocks(blocks), f'Word round-trip mismatch: {source.name}'
    embedded=[hashlib.sha256(r.target_part.blob).hexdigest() for r in reopened.part.rels.values() if r.reltype.endswith('/image')]
    assert sorted(embedded)==sorted(f['sha256'] for f in figures)
    cells=sum(len(t.rows)*len(t.columns) for t in reopened.tables)
    assert cells==(0 if main else 314)
    return {'source':source.name,'source_sha256':sha(source),'docx':out.name,'docx_sha256':sha(out),'paragraphs':npara,'figures':figures,'tables':ntab,'table_cells':cells,'text_and_cells_roundtrip':'exact','image_bytes_roundtrip':'exact'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path,default=ROOT);p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
    rows=[build(a.source_dir/(n+'.md'),a.outdir/(n+'.docx'),n=='Manuscript') for n in ['Manuscript','Supporting_Information']]
    report={'status':'EDITABLE_SOURCE_WORD_EXPORT_PASS','scope':'Exact text, table cells and embedded image bytes; no claim of identical pagination to the author-local Word candidate or journal submission readiness.','documents':rows}
    (a.outdir/'word_export_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
