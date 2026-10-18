#!/usr/bin/env python3
"""Pre-response eligibility audit of public IMMP butterfly life-stage records.

Input: a SHA-pinned public Dryad XLSX. Read column names and sheet dimensions
only, WITHOUT analysing outcome values, testing ecological sign hypotheses,
reidentifying monitoring sites or claiming larval survival. This guard exists
to prevent opportunistic use of an incomplete dataset as butterfly demography.

Only standard-library OOXML inspection, not spreadsheet editing/export.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

EXPECTED_SHA256="d9ed67778e96643ea4414da8971340b9c769a0a0aaaf70871700ecce36016a3d"
SRC_DOI="10.5061/dryad.mkkwh718x"
SRC_VERSION=6
SRC_DATE="2026-09-02"
MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P_REL="http://schemas.openxmlformats.org/package/2006/relationships"

VARIABLE_CLASSES={
    "egg":re.compile(r"egg|ovipos",re.I),
    "juvenile":re.compile(r"larva|larvae|instar|caterpillar|pupa|chrysal",re.I),
    "adult":re.compile(r"adult|adultcount",re.I),
    "milkweed":re.compile(r"milkweed|asclepias|host.?plant|plant.?species",re.I),
    "site":re.compile(r"site|plot|transect|location",re.I),
    "survey_date":re.compile(r"date|month|year|visit",re.I),
    "effort":re.compile(r"effort|observed|searched|surveyed|plants.?checked",re.I),
}


def col_to_idx(c:str)->int:
    num=0
    for x in c:
        if not "A"<=x<="Z":
            break
        num=num*26+ord(x)-ord("A")+1
    return num


def sheet_summaries(path:Path)->list[dict]:
    with zipfile.ZipFile(path) as archive:
        book=ET.fromstring(archive.read("xl/workbook.xml"))
        rels=ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rid_map={}
        for elem in rels.findall(f"{{{P_REL}}}Relationship"):
            target=elem.attrib.get("Target","")
            if target.startswith("/"):
                resolved=target.lstrip("/")
            elif target.startswith("xl/"):
                resolved=target
            else:
                resolved="xl/"+target
            rid_map[elem.attrib.get("Id","")]=resolved
        strings=[]
        if "xl/sharedStrings.xml" in archive.namelist():
            x=ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for e in x.findall(f"{{{MAIN}}}si"):
                strings.append("".join(n.text or "" for n in e.iter(f"{{{MAIN}}}t")))
        result=[]
        for e in book.findall(f".//{{{MAIN}}}sheet"):
            name=e.attrib["name"]
            rid=e.attrib.get(f"{{{REL}}}id","")
            resolved=rid_map.get(rid)
            if not resolved or resolved not in archive.namelist():
                raise ValueError("workbook sheet relationships cannot be resolved")
            x=ET.fromstring(archive.read(resolved))
            body=x.find(f".//{{{MAIN}}}sheetData")
            if body is None:
                continue
            count=0
            preview=[]
            # row-level metadata only; NO observations or data values leave
            # the audit (other than checking the presence of candidate headers).
            for row in body.findall(f"{{{MAIN}}}row"):
                count+=1
                if len(preview)<6:
                    items=[]
                    for cell in row.findall(f"{{{MAIN}}}c"):
                        cell_type=cell.attrib.get("t","")
                        v=cell.find(f"{{{MAIN}}}v")
                        inline=cell.find(f"{{{MAIN}}}is")
                        if cell_type=="s" and v is not None and v.text:
                            text=strings[int(v.text)]
                        elif cell_type=="inlineStr" and inline is not None:
                            text="".join(k.text or "" for k in inline.iter(f"{{{MAIN}}}t"))
                        else:
                            text=(v.text or "") if v is not None else ""
                        if text.strip():
                            items.append((col_to_idx(cell.attrib.get("r","")),text.strip()))
                    if items:
                        preview.append(items)
            # Identify candidate headers from early lines by column text,
            # but do not show raw values from outcome rows in the report.
            headers=[]
            for line in preview[:4]:
                found=sum(
                    bool(pattern.search(txt))
                    for _,txt in line
                    for pattern in VARIABLE_CLASSES.values()
                )
                if found>=2 or len(headers)==0 and count>0:
                    if found >= len([1 for _,c in headers]):
                        headers=line
            clean_headers=[x[:90] for _,x in sorted(headers)]
            candidates={
                key:[c for c in clean_headers if regex.search(c)]
                for key,regex in VARIABLE_CLASSES.items()
            }
            result.append({
                "worksheet":name,
                "nonempty_row_count_upper_bound":count,
                "header_candidates":clean_headers[:120],
                "candidate_variable_classes":candidates,
            })
        return result


def eligibility(sheets:list[dict])->dict:
    adult=[s for s in sheets if s["candidate_variable_classes"]["adult"]]
    egg=[s for s in sheets if s["candidate_variable_classes"]["egg"]]
    juvenile=[s for s in sheets if s["candidate_variable_classes"]["juvenile"]]
    host=[s for s in sheets if s["candidate_variable_classes"]["milkweed"]]
    both=[s for s in sheets if s in egg and s in juvenile]
    return {
        "stage_count_columns_exist_in_same_worksheet":bool(both),
        "egg_worksheets":[x["worksheet"] for x in egg],
        "juvenile_worksheets":[x["worksheet"] for x in juvenile],
        "host_or_milkweed_worksheets":[x["worksheet"] for x in host],
        "adult_worksheets":[x["worksheet"] for x in adult],
        "paired_native_versus_introduced_species_at_same_site":
            "NOT_IDENTIFIABLE_FROM_COLUMN_NAMES",
        "cohort_linked_egg_to_recruit_survival":
            "NOT_IDENTIFIABLE_FROM_COLUMN_NAMES",
        "test_of_negative_butterfly_recruitment":
            "NOT_PERFORMED",
        "action":
            "Requires explicit plant-species-native-status, same-site/season sampling comparability and genuine stage-followup before any positive or negative ecological fitness statement.",
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--public-xlsx",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    a=p.parse_args()
    sha=hashlib.sha256(a.public_xlsx.read_bytes()).hexdigest()
    if sha!=EXPECTED_SHA256:
        raise ValueError("2026 Dryad data version/byte digest changed")
    sheets=sheet_summaries(a.public_xlsx)
    if not any(s["worksheet"]=="Act2" for s in sheets):
        raise ValueError("published workbook lacks stated Activity 2 egg/larval sheet")
    obj={
        "schema":"chocho_external_immp_monarch_lifestage_schema_gate_v0.1",
        "status":"OBSERVATIONAL_SCHEMA_AUDIT_ONLY_NO_BUTTERFLY_FITNESS_TEST",
        "source":{
            "doi":SRC_DOI,
            "version":SRC_VERSION,
            "published":SRC_DATE,
            "workbook_sha256":sha,
            "source_license":"CC0-1.0",
            "record_period":"2021-2022 California public IMMP subset; not global Lepidoptera",
        },
        "worksheets":sheets,
        "eligibility":eligibility(sheets),
        "no_inference":[
            "Number of juveniles divided by number of eggs observed in the same visit is NOT survival of an identified cohort.",
            "Naturalized versus native milkweed comparison requires actually identified host species and site-season support, not mere presence of an adult count column.",
            "Stage-specific count presence cannot identify oviposition diversion, adaptive preference, per-egg fitness or causal invasive-plant impacts.",
            "This one-species dataset cannot show global nonmonotonicity of butterfly recruitment without an independent cross-species sample.",
        ],
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in obj.items() if k!="worksheets"},indent=2))


if __name__=="__main__":
    main()
