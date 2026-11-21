#!/usr/bin/env python3
"""Source-guarded citation copyedit for the separate journal-facing GEB files.

The historical article and complete source bibliography remain unchanged.
The curated manifest records the primary metadata and all allowed replacements.
"""
import hashlib
import json
import re
from pathlib import Path

MANIFEST = Path("manuscript/integrated_butterfly_resource_access_20261109/GEB_REFERENCE_COPYEDIT_V01.json")


def load_manifest(root, historical_references):
    manifest = json.loads((root / MANIFEST).read_text(encoding="utf8"))
    if manifest["schema_version"] != 1:
        raise RuntimeError("Unsupported GEB citation manifest")
    for path, expected in manifest["source_sha256"].items():
        if hashlib.sha256((root / path).read_bytes()).hexdigest() != expected:
            raise RuntimeError("Citation source changed; review manifest: " + path)
    if historical_references != manifest["historical_references"]:
        raise RuntimeError("Historical bibliography differs from audited source")
    records = manifest["references"]
    ids = [record["id"] for record in records]
    if len(ids) != len(set(ids)):
        raise RuntimeError("Duplicate reference identity")
    for record in records:
        authors = record["authors"]
        if not authors or record["reference"].count(" ... ") != (len(authors) >= 8):
            raise RuntimeError("APA6 reference author-count mismatch: " + record["id"])
        if "et al." in record["reference"]:
            raise RuntimeError("Abbreviated reference author list: " + record["id"])
        if str(record["year"]) not in record["reference"]:
            raise RuntimeError("Reference date mismatch: " + record["id"])
        if record.get("doi") and "https://doi.org/" + record["doi"] not in record["reference"]:
            raise RuntimeError("Reference DOI mismatch: " + record["id"])
    return manifest


def apply_copyedit(text, scope, manifest):
    for edit in manifest["replacements"][scope]:
        actual = text.count(edit["before"])
        if actual != edit["expected_count"]:
            raise RuntimeError(f"Citation source drift ({scope}/{edit['id']}): {actual}")
        text = text.replace(edit["before"], edit["after"])
    return text


def reference_lines(manifest, scope):
    records = [r for r in manifest["references"] if scope in r["documents"]]
    records.sort(key=lambda r: r["sort_key"].casefold())
    return ["- " + r["reference"] for r in records]


def si_bibliography(text, manifest):
    old = manifest["si_reference_block"]
    if text.count(old) != 1:
        raise RuntimeError("Supporting Information bibliography drift")
    return text.replace(old, "\n".join(reference_lines(manifest, "si")))


def verify_scientific_numbers(before, after, manifest):
    # Only remove registered author-date strings. Every remaining numerical
    # token (including signs, dates, cohorts, intervals and figure paths) must
    # survive the copyedit in the same order.
    def tokens(text):
        citations = {c for r in manifest["references"]
                     for c in (r["first_citation"], r["later_citation"])}
        for citation in sorted(citations, key=len, reverse=True):
            text = text.replace(citation, "")
        return re.findall(r"[+−-]?\d+(?:[.,]\d+)*%?", text)
    if tokens(before) != tokens(after):
        raise RuntimeError("Scientific numerical tokens changed during citation copyedit")
    return len(tokens(after))


def audit_citations(text, manifest, scope):
    records = [r for r in manifest["references"] if scope in r["documents"]]
    all_forms = {}
    counts = {}
    for record in records:
        forms = {record["first_citation"], record["later_citation"]}
        for form in forms:
            all_forms[form] = record["id"]
        hits = sorted((m.start(), form) for form in forms
                      for m in re.finditer(re.escape(form), text))
        if not hits:
            raise RuntimeError(f"Uncited {scope} reference: {record['id']}")
        if hits[0][1] != record["first_citation"]:
            raise RuntimeError(f"APA6 first citation incomplete ({scope}): {record['id']}")
        counts[record["id"]] = len(hits)
    # Source observation dates are not bibliographic citations. Every other
    # parenthetical author-year segment must resolve to an audited reference.
    allowed_dates = {"wild, 2015", "garden, 2011"}
    for group in re.findall(r"\(([^()]*)\)", text):
        cited = []
        for segment in group.split(";"):
            segment = segment.strip()
            if not re.fullmatch(r"[^\d]+, (?:19|20)\d{2}[a-z]?", segment):
                continue
            if segment in allowed_dates:
                continue
            if segment not in all_forms:
                raise RuntimeError(f"Unresolved {scope} author-date citation: {segment}")
            cited.append(all_forms[segment])
        if len(cited) > 1:
            keys = {r["id"]: r["sort_key"].casefold() for r in records}
            if cited != sorted(cited, key=keys.get):
                raise RuntimeError(f"Unalphabetized {scope} citation group: {group}")
    return counts


def citation_receipt(manifest, main_text, si_text, numerical_tokens):
    main_counts = audit_citations(main_text, manifest, "main")
    si_counts = audit_citations(si_text, manifest, "si")
    return {
        "status": "PASS: SOURCE-GUARDED APA6 CITATION AND REFERENCE ALIGNMENT",
        "manifest": str(MANIFEST),
        "metadata_checked_date": manifest["checked_date"],
        "historical_references_preserved": len(manifest["historical_references"]),
        "historical_references_omitted_from_journal_export": len(manifest["omitted_historical_references"]),
        "main_reference_count": len(main_counts),
        "SI_reference_count": len(si_counts),
        "main_citation_counts": main_counts,
        "SI_citation_counts": si_counts,
        "main_scientific_numeric_tokens_preserved": numerical_tokens,
        "remaining_source_limits": manifest["remaining_source_limits"],
    }
