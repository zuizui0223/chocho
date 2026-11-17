#!/usr/bin/env python3
"""Source audit of 1987 Cape Breton Highlands National Park butterfly checklist.

Reads local pdftotext output from original government PDF (no OCR). Records
exact modern-name hits only, with physical PDF page; none of this is an
absence, native-range, colonization, or current endemism assignment.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
COHORT="results/butterfly_GBIF_effort_origin_gate_20261117/LATER_FIRST_INDICES_TO_REVIEW_20261117.csv"
PIN="f439350de7e2685a1af7197e71df8e6f2b739601"
SOURCE="https://publications.gc.ca/collections/collection_2016/aac-aafc/agrhist/A42-121-1987-eng.pdf"

def blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--text",type=Path,required=True)
    parser.add_argument("--pdf-sha256",type=str,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    p=ROOT/COHORT
    if blob(p)!=PIN:raise RuntimeError("frozen 317-source drift")
    with p.open(newline="",encoding="utf-8") as f:
        orig=list(csv.DictReader(f))
    target=[r for r in orig if r["landmass_id"]=="79"]
    if len(target)!=17:raise RuntimeError("Cape Breton exact study cohort != 17")
    document=args.text.read_text(encoding="utf-8",errors="replace")
    pages=document.split("\f")
    results=[]
    for item in target:
        taxon=item["species"]
        genus,epithet=taxon.split()
        pattern=re.compile(r"\b"+re.escape(genus)+r"\s+"+re.escape(epithet)+r"\b",flags=re.IGNORECASE)
        hits=[]
        for i,content in enumerate(pages):
            if i<180:continue
            for m in list(pattern.finditer(content))[:4]:
                fragment=" ".join(content[max(0,m.start()-150):min(len(content),m.end()+230)].split())
                hits.append({"physical_PDF_page_1based":i+1,"context":fragment})
        results.append({"species":taxon,"speciesKey":item["speciesKey"],
          "original_GBIF_first_index":item["first_GBIF_index_window"],
          "exact_binomial_hits_in_original_government_1987_pdf":len(hits),
          "first_four_contexts":hits[:4]})
    summaries=[]
    for i,page in enumerate(pages):
        if i>=200 and re.search(r"\bButterflies\b|\bHesperiidae\b|\bPapilionidae\b",page,re.IGNORECASE):
            summaries.append({"physical_PDF_page_1based":i+1,"first_420_chars":" ".join(page.split())[:420]})
    result={
      "schema":"chocho_government_of_Canada_1987_cape_breton_17_candidate_taxon_search_v1",
      "status":"RAW_PUBLICATION_PDF_TEXT_SEARCH_NOT_TAXONOMIC_RECONCILIATION",
      "official_source_URL":SOURCE,
      "original_PDF_SHA256":args.pdf_sha256,
      "pages_extracted":len(pages),
      "total_text_characters":len(document),
      "indexed_candidate_species":len(target),
      "exact_binomial_present":sum(bool(x["exact_binomial_hits_in_original_government_1987_pdf"]) for x in results),
      "taxa":results,
      "butterfly_page_locator":summaries[:20],
      "interpretation":"Exact names in original government text indicate possible historically documented taxa, but a line may be hypothetical, expected, or cited outside the island: verify individual species account before promoting to independent exact-island historical presence."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"source":SOURCE,"sha256":args.pdf_sha256,"pages":len(pages),
      "exact_binomial_present":result["exact_binomial_present"],
      "taxa":[{"species":z["species"],"hits":z["exact_binomial_hits_in_original_government_1987_pdf"],"context":z["first_four_contexts"][:1]} for z in results],
      "butterfly_locator":summaries[-12:]},ensure_ascii=False,indent=2))

if __name__=="__main__":main()
