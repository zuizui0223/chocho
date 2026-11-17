#!/usr/bin/env python3
"""Independent butterfly-checklist denominator vs frozen NHM HOSTS coverage.

Inputs are archived CSV files, NOT species pre-selected from HOSTS. Uses
standard-library Python only. No missing host match is interpreted as a
biological absence or a true absence of known hosts in other sources.
"""
import argparse
import csv
import hashlib
import json
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

CHECKLIST_TAXON_BLOB = "67f3126314c7facd5b133ca06f18f397b8538549"
CHECKLIST_DISTRIBUTION_BLOB = "1e57f4cbf2ec78e84a7ba2d21b112eda49396130"
HOSTS_RESOURCE_BLOB = "06749c04bd5b7a754d6f7717a31c01c51357cae8"
MADEIRA_STATUS_SOURCE = (
    "Borges et al. (eds.) 2008, A list of the terrestrial fungi, flora "
    "and fauna of Madeira and Selvagens archipelagos, pp. 342–343"
)

# Status applies to Madeira archipelago only; copied from an independent
# published faunal checklist. n=native; END=island endemic; MAC=Macaronesian;
# i=introduced; m=migrant; dagger=historically extinct.
# Two END cases are described at subspecific rank in the older checklist.
MADEIRA_STATUS = {
    "Colias crocea": ("native_nonendemic", "n; checklist Colias croceus"),
    "Danaus plexippus": ("native_nonendemic", "n"),
    "Lampides boeticus": ("native_nonendemic", "n"),
    "Leptotes pirithous": ("native_nonendemic", "n"),
    "Pararge aegeria": ("native_nonendemic", "n"),
    "Vanessa atalanta": ("native_nonendemic", "n"),
    "Vanessa cardui": ("native_nonendemic", "n"),
    "Gonepteryx maderensis": ("island_endemic_species", "END"),
    "Hipparchia maderensis": ("island_endemic_species", "END"),
    "Pararge xiphia": ("island_endemic_species", "END"),
    "Pieris wollastoni": (
        "historic_extinct_endemic_taxon",
        "END dagger; older taxon Pieris brassicae wollastoni",
    ),
    "Lycaena phlaeas": (
        "endemic_subspecies_only",
        "END applies to Lycaena phlaeas phlaeoides, not all L. phlaeas",
    ),
    "Vanessa vulcania": ("macaronesian", "MAC"),
    "Pieris rapae": ("introduced", "i"),
    "Issoria lathonia": ("migrant", "m; INBO both present and irregular on Madeira"),
}
ALIAS_CANDIDATES_UNVERIFIED = {
    "Cyclyrius webbianus": "Leptotes webbianus",
    "Pontia daplidice": "Pieris daplidice",
}
OUTPUT_COLUMNS = [
    "archipelago", "island", "scientific_name", "checklist_family",
    "occurrence_status", "has_irregular_duplicate", "madeira_origin",
    "madeira_source_note", "hosts_exact_insect_binomial",
    "hosts_named_host_genus", "hosts_named_host_species",
    "hosts_exact_rows", "hosts_named_genus_rows",
    "hosts_named_species_rows", "hosts_family_match",
    "unverified_alias_candidate",
]


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def normalize_name(name):
    text = unicodedata.normalize("NFKD", name or "").lower()
    return " ".join(re.findall(r"[a-z]+", text))


def valid_plant_genus(name):
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z-]*", (name or "").strip())) and (
        (name or "").lower() not in {"unknown", "unidentified", "indet", "fungi"}
    )


def valid_plant_species(name):
    return bool(re.fullmatch(r"[a-z][a-z-]*", (name or "").strip())) and (
        (name or "").lower() not in {"sp", "spp", "indet", "unknown", "cf", "aff"}
    )


def audit(taxon_csv, distribution_csv, hosts_csv, out_dir, strict=True):
    taxon = read_csv(taxon_csv)
    distribution = read_csv(distribution_csv)
    hosts = read_csv(hosts_csv)

    taxonomy = {r["taxonID"]: r for r in taxon}
    hosts_index = defaultdict(lambda: {
        "rows": 0, "named_genus": 0, "named_species": 0, "families": set()
    })
    for r in hosts:
        genus = normalize_name(r.get("Insect Genus", ""))
        epithet = normalize_name(r.get("Insect Species", ""))
        if len(genus.split()) != 1 or len(epithet.split()) != 1:
            continue
        entry = hosts_index[genus + " " + epithet]
        entry["rows"] += 1
        if valid_plant_genus(r.get("Hostplant Genus", "")):
            entry["named_genus"] += 1
            if valid_plant_species(r.get("Hostplant Species", "")):
                entry["named_species"] += 1
        entry["families"].add((r.get("Insect Family") or "").strip())

    # Build denominator exclusively from independent island presence status.
    island_status = defaultdict(set)
    for r in distribution:
        locality = r.get("locality") or ""
        if not locality.startswith(("Azores, ", "Canary Islands, ", "Madeira, ")):
            continue  # Do not count archipelago aggregate or mainland countries.
        taxon_id = r.get("taxonID")
        if taxon_id not in taxonomy:
            raise ValueError("Unknown independent taxonID: " + str(taxon_id))
        island_status[(locality, taxon_id)].add(r.get("occurrenceStatus", ""))

    audited = []
    for (island, taxon_id), states in sorted(island_status.items()):
        if "present" not in states:
            continue
        t = taxonomy[taxon_id]
        name = t["scientificName"]
        matched = hosts_index.get(normalize_name(name), {
            "rows": 0, "named_genus": 0, "named_species": 0, "families": set()
        })
        status, source_note = ("not_independently_classified", "")
        if island.startswith("Madeira, "):
            status, source_note = MADEIRA_STATUS.get(
                name, ("not_in_2008_status_crosswalk", "")
            )
        audited.append({
            "archipelago": island.split(", ")[0],
            "island": island,
            "scientific_name": name,
            "checklist_family": t["family"],
            "occurrence_status": "present",
            "has_irregular_duplicate": int("irregular" in states),
            "madeira_origin": status,
            "madeira_source_note": source_note,
            "hosts_exact_insect_binomial": int(matched["rows"] > 0),
            "hosts_named_host_genus": int(matched["named_genus"] > 0),
            "hosts_named_host_species": int(matched["named_species"] > 0),
            "hosts_exact_rows": matched["rows"],
            "hosts_named_genus_rows": matched["named_genus"],
            "hosts_named_species_rows": matched["named_species"],
            "hosts_family_match": int(t["family"] in matched["families"]),
            "unverified_alias_candidate": ALIAS_CANDIDATES_UNVERIFIED.get(name, ""),
        })

    if len({(r["island"], r["scientific_name"]) for r in audited}) != len(audited):
        raise AssertionError("Multiple independent names share one island pair")
    by_species = {}
    for r in audited:
        by_species.setdefault(r["scientific_name"], r)

    def cover(rows):
        return {
            "n": len(rows),
            "hosts_exact_binomial": sum(r["hosts_exact_insect_binomial"] for r in rows),
            "hosts_named_genus": sum(r["hosts_named_host_genus"] for r in rows),
            "hosts_named_species": sum(r["hosts_named_host_species"] for r in rows),
        }

    per_archipelago = {}
    for archipelago in ("Azores", "Canary Islands", "Madeira"):
        per_archipelago[archipelago] = cover(
            [r for r in audited if r["archipelago"] == archipelago]
        )
    madeira_strata = {
        origin: cover([r for r in audited if r["island"] == "Madeira, Madeira"
                       and r["madeira_origin"] == origin])
        for origin in sorted(set(r["madeira_origin"] for r in audited
                                 if r["island"] == "Madeira, Madeira"))
    }
    source_path = {"taxon": taxon_csv, "distribution": distribution_csv, "hosts": hosts_csv}
    result = {
        "schema": "chocho_independent_butterfly_checklist_host_coverage_v1",
        "scope": "INBO 18 individual Macaronesian islands; occurrenceStatus=present",
        "unit": "checklist butterfly species x actual named island",
        "sources": {
            "inbo_gbif_doi": "10.15468/ye7whj",
            "taxon_git_blob": CHECKLIST_TAXON_BLOB,
            "distribution_git_blob": CHECKLIST_DISTRIBUTION_BLOB,
            "hosts_mirror_git_commit": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "hosts_resource_git_blob": HOSTS_RESOURCE_BLOB,
            "madeira_status_reference": MADEIRA_STATUS_SOURCE,
            "inputs_sha256": {
                k: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                for k, path in source_path.items()
            },
        },
        "source_row_counts": {
            "taxon": len(taxon), "distribution": len(distribution), "hosts": len(hosts)
        },
        "present_islands": len(set(r["island"] for r in audited)),
        "species_island": cover(audited),
        "independent_unique_species": cover(list(by_species.values())),
        "archipelago_species_island": per_archipelago,
        "madeira_island_origin_strata": madeira_strata,
        "duplicates_present_and_irregular": [
            {"island": r["island"], "taxon": r["scientific_name"]}
            for r in audited if r["has_irregular_duplicate"]
        ],
        "unverified_taxonomic_aliases": ALIAS_CANDIDATES_UNVERIFIED,
        "limitations": [
            "Exact binomial coverage is not taxonomically harmonized HOSTS coverage.",
            "No HOSTS match does not imply the larval host is truly unknown.",
            "Independent checklist occurrence is not contemporary breeding verification.",
            "Madeira 2008 origin statuses are historical and are not projected to other islands.",
            "Lycaena phlaeas has an endemic subspecies: species-level HOSTS cannot validate subspecific use.",
            "Pieris wollastoni is extinct historically; older source classified it as P. brassicae wollastoni.",
        ],
    }

    if strict:
        assert (len(taxon), len(distribution), len(hosts)) == (496, 8020, 140485)
        assert result["present_islands"] == 18
        assert result["species_island"]["n"] == 227
        assert result["independent_unique_species"]["n"] == 42
        assert result["independent_unique_species"]["hosts_named_genus"] == 19
        assert result["independent_unique_species"]["hosts_named_species"] == 18
        assert result["species_island"]["hosts_named_genus"] == 168
        assert result["archipelago_species_island"]["Azores"]["n"] == 61
        assert result["archipelago_species_island"]["Canary Islands"]["n"] == 141
        assert result["archipelago_species_island"]["Madeira"]["n"] == 25
        assert madeira_strata["native_nonendemic"]["n"] == 7
        assert madeira_strata["native_nonendemic"]["hosts_named_genus"] == 7
        assert madeira_strata["island_endemic_species"]["n"] == 3
        assert madeira_strata["island_endemic_species"]["hosts_named_genus"] == 0

    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "species_island_coverage.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(audited)
    with (out_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxon", required=True, type=Path)
    parser.add_argument("--distribution", required=True, type=Path)
    parser.add_argument("--hosts", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--no-expected-assertions", action="store_true")
    args = parser.parse_args()
    audit(args.taxon, args.distribution, args.hosts, args.out_dir,
          strict=not args.no_expected_assertions)


if __name__ == "__main__":
    main()
