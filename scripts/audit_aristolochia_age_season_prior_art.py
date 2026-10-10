#!/usr/bin/env python3
"""Source-check 2023 Aristolochia host ontogeny × season butterfly performance.

Prior published butterfly data, NOT a new biological result. Fetches article
XML/HTML only, never invents per-plant raw records. Prevents overclaiming
pure food-mass mechanisms on supplementing fresh leaves.
"""
from __future__ import annotations
import argparse,hashlib,json,re
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.request import Request,urlopen
from html.parser import HTMLParser

URLS=[
 "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10265686/fullTextXML",
 "https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2023.1145363/full"
]
DOI="10.3389/fpls.2023.1145363"

class Plain(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts=[]
        self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag in ("script","style","nav"):
            self.skip+=1
    def handle_endtag(self,tag):
        if tag in ("script","style","nav") and self.skip:
            self.skip-=1
    def handle_data(self,data):
        if not self.skip:self.parts.append(data)

def plain_text(raw:bytes):
    stripped=raw.lstrip()
    if stripped.startswith(b"<?xml") or stripped.startswith(b"<article"):
        root=ET.fromstring(raw)
        return " ".join(root.itertext())
    parser=Plain()
    parser.feed(raw.decode("utf-8",errors="replace"))
    return " ".join(parser.parts)

def normalized(text:str):
    return re.sub(r"\s+"," ",text.replace("\u00a0"," ")).lower()

def audit_text(raw:bytes,url:str):
    t=normalized(plain_text(raw))
    gate={
       "correct_DOI":DOI in t,
       "original_species_Sericinus": "sericinus montela" in t,
       "second_species_Spodoptera": "spodoptera exigua" in t,
       "host_plants_Aristolochia_contorta": "aristolochia contorta" in t,
       "published_age_factor": "1st-year" in t or "1 st -year" in t or "1 st - and 3 rd" in t,
       "published_two_seasons":"july" in t and "september" in t,
       "published_leaf_C_N":"c/n" in t,
       "published_two_acids":"aristolochic acid 1" in t and "aristolochic acid 2" in t,
       "reported_short_test_duration":"six days" in t or "6 days" in t,
       "reported_temperature_24_C":"24°c" in t or "24 °c" in t or "24° c" in t or "24 c" in t,
       "reported_leaf_replenishment":"50%" in t and "replaced" in t,
    }
    return {"status":"AUTHOR_SOURCE_ARTICLE_VERIFIED" if all(gate.values()) else "RETRIEVED_SOURCE_NOT_FULLY_MATCHED",
            "matched_checks":gate,"source_url":url,
            "retrieved_bytes_sha256":hashlib.sha256(raw).hexdigest(),
            "source_kind":"published_article_fulltext_NOT_original_individual_records",
            "source_text_chars":len(t),
            "has_original_rowlevel_per_plant_data":False}

def download(url):
    req=Request(url,headers={"User-Agent":"chocho-research-public-article-audit/1.0","Accept":"application/xml, text/html, */*"})
    with urlopen(req,timeout=35) as resp:
        content=resp.read(3_000_001)
    if not content or len(content)>3_000_000:
        raise ValueError("article source empty or over limit")
    return content

def run():
    logs=[]
    for url in URLS:
        try:
            result=audit_text(download(url),url)
            if result["status"]=="AUTHOR_SOURCE_ARTICLE_VERIFIED":
                break
            logs.append({"source":url,"status":result["status"],
                         "failed_keys":[k for k,v in result["matched_checks"].items() if not v]})
        except (OSError,ValueError,ET.ParseError) as exc:
            logs.append({"source":url,"status":"SOURCE_ACCESS_OR_PARSE_BLOCKED",
                         "detail":f"{type(exc).__name__}: {str(exc)[:180]}"})
    else:
        result={"status":"ORIGINAL_ARTICLE_ACCESS_NOT_VERIFIED_IN_AUTOMATED_RUN",
                "source_kind":"public literature without original raw individual observations"}
    return {"schema":"chocho_aristolochia_ontogeny_season_2023_source_gate_v01",
            "original_source":DOI,"paper_title":"Age-dependent resistance of a perennial herb, Aristolochia contorta against specialist and generalist leaf-chewing herbivores",
            "result":result,"fallback_log":logs,
            "metadata_facts_as_published":{
               "two_seasons":["July","September"],
               "plant_year_classes":[1,2,3],
               "butterfly_feeding_start_instar":2,
               "specialist":"Sericinus montela",
               "nonbutterfly_comparator":"Spodoptera exigua (moth)",
               "study_duration_days":6,
               "reported_fixed_24C_in_July_and_September":True,
               "leaves_replaced_when_over_half_consumed":True,
               "AAI_July_age_year1_leaf_ng_per_mg":83.80,
               "AAI_July_age_year2_leaf_ng_per_mg":7.61,
               "AAI_July_age_year3_leaf_ng_per_mg":2.26,
               "statistical_age_by_season_interaction_for_Sericinus_RGR_F":8.952
            },
            "metadata_facts_not_new_discovery":True,
            "available_2026_Korean_same_ramet_original_data":False,
            "can_identify_2023_Kyoto_interspecific_competition_mechanism":False,
            "GEB_PR38_untouched":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receipt",type=Path,required=True)
    args=ap.parse_args()
    r=run()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False,indent=2))

if __name__=="__main__":main()
