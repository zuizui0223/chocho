#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
from html.parser import HTMLParser
from pathlib import Path


class TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_cell = False
        self.cell_parts: list[str] = []
        self.row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag.lower() == "tr":
            self.row = []
        elif tag.lower() in {"td", "th"}:
            self.in_cell = True
            self.cell_parts = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in {"td", "th"} and self.in_cell:
            text = html.unescape(" ".join(self.cell_parts))
            text = re.sub(r"\s+", " ", text).strip()
            self.row.append(text)
            self.in_cell = False
            self.cell_parts = []
        elif tag == "tr":
            if self.row:
                self.rows.append(self.row)
            self.row = []


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_botanical(value: str) -> tuple[set[str], set[str]]:
    value = value.replace("×", "x")
    exact: set[str] = set()
    genera: set[str] = set()

    # FAO sometimes writes "Phaseolus and Vigna spp.".
    for a, b in re.findall(
        r"\b([A-Z][A-Za-z-]+)\s+and\s+([A-Z][A-Za-z-]+)\s+spp\.?",
        value,
    ):
        genera.update((a, b))

    for genus in re.findall(r"\b([A-Z][A-Za-z-]+)\s+spp\.?", value):
        genera.add(genus)

    # Extract binomials even when followed by var./ssp. or parenthetical synonyms.
    for genus, epithet in re.findall(
        r"\b([A-Z][A-Za-z-]+)\s+([a-z][a-z-]+)\b",
        value,
    ):
        if epithet in {"spp", "sp"}:
            continue
        if genus in {"See", "Hybrid", "Mixture", "Various"}:
            continue
        exact.add(f"{genus} {epithet}")

    return exact, genera


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    raw = args.html.read_text(encoding="utf-8", errors="replace")
    parser = TableParser()
    parser.feed(raw)

    exact: set[str] = set()
    genera: set[str] = set()
    source_rows = 0
    for row in parser.rows:
        if len(row) < 3:
            continue
        if row[0].casefold() == "crop name":
            continue
        botanical = row[1].strip()
        if not botanical:
            continue
        # Crop-code cells are numeric-ish. This avoids classification tables
        # elsewhere on the page being mistaken for Appendix 4.
        code = row[2].strip()
        if not re.match(r"^\d", code):
            continue
        e, g = parse_botanical(botanical)
        if e or g:
            source_rows += 1
            exact.update(e)
            genera.update(g)

    if len(exact) < 80:
        raise RuntimeError(f"too few FAO crop binomials parsed: {len(exact)}")
    if len(genera) < 10:
        raise RuntimeError(f"too few FAO crop genera parsed: {len(genera)}")
    required = {"Medicago sativa", "Zea mays", "Oryza sativa", "Dactylis glomerata"}
    if not required <= exact:
        raise RuntimeError(f"required crop controls missing: {sorted(required - exact)}")
    # Avena and Trifolium are explicitly genus-level in the FAO source.
    for required_genus in ("Avena", "Trifolium", "Lotus"):
        if required_genus not in genera:
            raise RuntimeError(f"required FAO crop genus missing: {required_genus}")

    payload = {
        "schema": "chocho_fao_crop_taxa_v0.1",
        "source": "FAO ICC Appendix 4",
        "source_url": "https://www.fao.org/4/a0135e/A0135E10.htm",
        "source_html_sha256": sha256_path(args.html),
        "parsed_crop_rows": source_rows,
        "exact_binomials": sorted(exact),
        "genus_wildcards": sorted(genera),
        "controls": {
            "Medicago_sativa": "Medicago sativa" in exact,
            "Zea_mays": "Zea mays" in exact,
            "Oryza_sativa": "Oryza sativa" in exact,
            "Avena_genus": "Avena" in genera,
            "Trifolium_genus": "Trifolium" in genera,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "exact_binomials": len(exact),
        "genus_wildcards": len(genera),
        "source_rows": source_rows,
        "source_sha256": payload["source_html_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
