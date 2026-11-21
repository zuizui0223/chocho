#!/usr/bin/env python3
"""Fetch SOURCE-IDENTIFIED Crossref author-list candidates for GEB reference copyedit.

Crossref metadata alone is NOT author-list acceptance. Keep the original
reference unchanged unless the DOI, lead author, year and publisher page agree.
This script neither edits the manuscript nor submits/publishes a dataset.
"""
from __future__ import annotations

import argparse
import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PAPER = "manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
DOI_RE = re.compile(r"https://doi.org/(\S+)", re.I)
YEAR_RE = re.compile(r"\((20\d{2}|19\d{2})\)")
FETCH_URL = "https://api.crossref.org/works/"
USER_AGENT = "chocho-GEB-v05-reference-metadata-preflight/1.0 (public DOI metadata check)"


def candidates(text: str) -> list[dict]:
    parts = text.split("## References (", 1)
    if len(parts) != 2:
        raise ValueError("Reference marker missing")
    refs = parts[1].split("## Supplementary results:", 1)[0]
    output = []
    for line in refs.splitlines():
        if not line.startswith("- ") or "et al." not in line:
            continue
        doi = DOI_RE.search(line)
        year = YEAR_RE.search(line)
        if not doi or not year:
            raise ValueError("Abbreviated reference lacks year or resolvable DOI: " + line)
        output.append({"doi": doi.group(1).rstrip(".,;"), "paper_year": int(year.group(1)),
                       "original_reference": line})
    if not output:
        raise ValueError("Expected abbreviated author references not found")
    if len({item["doi"].casefold() for item in output}) != len(output):
        raise ValueError("DOI duplicate in abbreviated references")
    return output


def fetch(doi: str) -> dict:
    url = FETCH_URL + urllib.parse.quote(doi, safe="/")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=25) as reply:
        data = json.load(reply)
    entry = data.get("message", {})
    if entry.get("DOI", "").casefold() != doi.casefold():
        raise ValueError("Crossref returned a different DOI for " + doi)
    return entry


def meta_years(msg: dict) -> list[int]:
    years = set()
    for k in ("published", "published-print", "published-online", "issued", "created"):
        x = msg.get(k, {}).get("date-parts") or []
        if x and x[0] and isinstance(x[0][0], int):
            years.add(x[0][0])
    return sorted(years)


def normalize(msg: dict, item: dict) -> dict:
    authors = [{"given": a.get("given", ""), "family": a.get("family", ""),
                "sequence": a.get("sequence", "")}
               for a in msg.get("author", [])]
    lead = authors[0]["family"].lower() if authors else ""
    paper = item["original_reference"].lower()
    lead_ok = bool(lead and lead in paper[:max(85, len(lead) + 25)])
    years = meta_years(msg)
    return {**item, "crossref_title": html.unescape(re.sub(r"<[^>]+>", "", (msg.get("title") or [""])[0])),
            "crossref_container": (msg.get("container-title") or [""])[0],
            "crossref_years": years, "paper_year_present_in_metadata": item["paper_year"] in years,
            "crossref_authors": authors, "crossref_author_count": len(authors),
            "lead_family_text_check": lead_ok,
            "metadata_only_NOT_APPROVED": True}


def run(text: str) -> dict:
    data = []
    for item in candidates(text):
        try:
            err = None
            for attempt in range(3):
                try:
                    msg = fetch(item["doi"])
                    data.append(normalize(msg, item))
                    break
                except (TimeoutError, urllib.error.URLError, ValueError) as exc:
                    err = str(exc)
                    if attempt < 2:
                        time.sleep(1.5 * (attempt + 1))
            else:
                data.append({**item, "error": err, "metadata_only_NOT_APPROVED": True})
        except Exception as exc:
            data.append({**item, "error": str(exc), "metadata_only_NOT_APPROVED": True})
    return {"schema": "chocho_GEB_v05_reference_metadata_candidates_v1",
            "scope": "original reference lines containing et al.; no manuscript changes",
            "source": "Crossref read-only DOI public metadata; requires human/source verification",
            "candidate_count": len(data),
            "fetched_count": sum("crossref_authors" in x for x in data),
            "review_required_before_substitution": True,
            "items": data}


def self_test() -> None:
    s = ("## References (DOI identities checked; final source and style audit outstanding)\n"
         "- Bach et al. (2022). Placeholder. https://doi.org/10.1111/ecog.06227\n"
         "## Supplementary results: anything\n")
    assert candidates(s)[0]["doi"] == "10.1111/ecog.06227"
    d = normalize({"DOI": "10.1111/ecog.06227", "author": [{"given":"A","family":"Bach"}],
                   "title":["Hello"], "issued":{"date-parts":[[2022]]}},
                  candidates(s)[0])
    assert d["paper_year_present_in_metadata"] and d["lead_family_text_check"]
    print("CROSSREF_AUTHOR_CANDIDATE_OFFLINE_SELFTEST_PASS")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-root", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    if a.paper_root is None or a.output is None:
        ap.error("--paper-root and --output are required")
    r = run((a.paper_root / PAPER).read_text(encoding="utf-8"))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"candidate_count":r["candidate_count"],"fetched_count":r["fetched_count"]}))
    if r["candidate_count"] != 14 or r["fetched_count"] < 12:
        raise SystemExit("Insufficient Crossref metadata; do not cite as complete")


if __name__ == "__main__":
    main()
