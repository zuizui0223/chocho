#!/usr/bin/env python3
"""Audit 239-species selection from independently listed Macaronesian butterflies.

Run after audit_independent_checklist_hosts.py, or against the frozen CSVs
stored with this result. Only exact butterfly binomials are counted. Does
not infer the reasons for sample exclusion.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
SOURCE = HERE / "species_island_coverage.csv"
FROZEN = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def audit():
    pairs = read_csv(SOURCE)
    frozen = read_csv(FROZEN)
    original_names = {r["species"] for r in frozen}
    taxa = {}
    for row in pairs:
        name = row["scientific_name"]
        if name not in taxa:
            taxa[name] = {
                "scientific_name": name, "checklist_family": row["checklist_family"],
                "island_count": 0, "archipelagos": set(),
                "hosts_named_genus": int(row["hosts_named_host_genus"]),
                "hosts_named_species": int(row["hosts_named_host_species"]),
                "in_frozen239": int(name in original_names), "madeira_origin": "",
            }
        item = taxa[name]
        item["island_count"] += 1
        item["archipelagos"].add(row["archipelago"])
        if row["island"] == "Madeira, Madeira":
            item["madeira_origin"] = row["madeira_origin"]

    rows = []
    for item in sorted(taxa.values(), key=lambda r: r["scientific_name"]):
        archipelagos = sorted(item.pop("archipelagos"))
        rows.append({
            "scientific_name": item["scientific_name"],
            "checklist_family": item["checklist_family"],
            "island_count": item["island_count"],
            "archipelago_count": len(archipelagos),
            "archipelagos": "; ".join(archipelagos),
            "hosts_named_genus": item["hosts_named_genus"],
            "hosts_named_species": item["hosts_named_species"],
            "in_frozen239": item["in_frozen239"],
            "madeira_origin": item["madeira_origin"],
        })

    in_model = [r for r in rows if r["in_frozen239"]]
    host_known = [r for r in rows if r["hosts_named_genus"]]
    outside_but_host_known = [
        r for r in rows if r["hosts_named_genus"] and not r["in_frozen239"]
    ]
    endemic = [r for r in rows if r["madeira_origin"] == "island_endemic_species"]
    native = [r for r in rows if r["madeira_origin"] == "native_nonendemic"]
    assert len(frozen) == 239 and len(rows) == 42
    assert len(host_known) == 19 and len(in_model) == 12
    assert len(outside_but_host_known) == 7
    assert len(endemic) == 3 and sum(r["in_frozen239"] for r in endemic) == 0
    assert len(native) == 7 and sum(r["in_frozen239"] for r in native) == 6

    outfile = HERE / "independent_42_species_vs_frozen239.csv"
    with outfile.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "scientific_name", "checklist_family", "island_count",
            "archipelago_count", "archipelagos", "hosts_named_genus",
            "hosts_named_species", "in_frozen239", "madeira_origin",
        ])
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "source_denominator": len(rows),
        "source_with_named_host_genus": len(host_known),
        "in_frozen239": len(in_model),
        "outside_frozen239": len(rows) - len(in_model),
        "nonselected_named_host": len(outside_but_host_known),
        "nonselected_no_named_host": len(rows)-len(host_known),
        "model239_size": len(frozen),
        "selected_species": [r["scientific_name"] for r in in_model],
        "nonselected_despite_HOSTS": [
            r["scientific_name"] for r in outside_but_host_known
        ],
        "madeira_endemic_in_frozen": sum(r["in_frozen239"] for r in endemic),
        "madeira_native_nonendemic_in_frozen": sum(
            r["in_frozen239"] for r in native
        ),
    }
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    audit()
