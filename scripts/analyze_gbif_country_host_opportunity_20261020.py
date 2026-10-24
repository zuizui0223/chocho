#!/usr/bin/env python3
"""Independent country-grain full-index butterfly/resource association.

This does not replace the bounded WGSRPD3 sample. Country co-coverage is not
local co-occurrence, host use, absence, breeding, or a causal introduction test.
The source-blind protocol is GBIF_COUNTRY_FULL_INDEX_MACRO_20261020.md.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np

from acquire_gbif_country_facets_20261020 import (
    SCHEMA as FACET_SCHEMA, load_original_manifest, build_queries, QUERY_IDENTITY,
)
from analyze_butterfly_239_GBIF_occurrence_host_opportunity_20261019 import (
    source_augmented_model, climate_and_plant_candidates, stratum_information_audit,
)

SEED = 20261020
ORIGINAL_SHA = "18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2"
SOURCES = ("all_accepted_bases", "core")
WINDOWS = ("pooled", "2010_2017", "2018_2025")
PRIMARY = ("all_accepted_bases", "pooled", 3, .05, 1)


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collision_species(records):
    bykey = defaultdict(list)
    for row in records:
        if row.get("taxon_key") is not None:
            bykey[row["taxon_key"]].append(row["species"])
    return {sp for names in bykey.values() if len(names) > 1 for sp in names}


def load_observations(count_path, receipt_path, original_path):
    original = load_original_manifest(original_path)
    if original["original_manifest_sha256"] != ORIGINAL_SHA:
        raise ValueError("Original taxonomic source bytes changed")
    receipt = json.loads(Path(receipt_path).read_text())
    if receipt.get("schema") != FACET_SCHEMA or receipt.get("original_manifest_sha256") != ORIGINAL_SHA:
        raise ValueError("Country source schema or original source hash changed")
    if file_sha(count_path) != receipt["outputs"]["country_counts_sha256"]:
        raise ValueError("Country aggregate CSV hash differs from acquisition receipt")
    queries = {q["query_id"]: q for q in build_queries(original)}
    received = receipt["query_receipts"]
    if len(received) != len(queries) or len({q["query_id"] for q in received}) != len(queries):
        raise ValueError("Missing or duplicate query receipt")
    for row in received:
        wanted = queries.get(row["query_id"])
        if wanted is None or any(row.get(k) != wanted[k] for k in QUERY_IDENTITY):
            raise ValueError("Query changed from frozen species/key/filter identity")
    expected_species = {r["species"]: r for r in original["species_receipts"]}
    records = receipt["species_receipts"]
    if len(records) != 239 or {r["species"] for r in records} != set(expected_species):
        raise ValueError("Country source lost original 239-species accounting")
    for row in records:
        old = expected_species[row["species"]]
        key = (old.get("taxon") or {}).get("key") if old["status"] != "TAXON_UNRESOLVED" else None
        if row.get("taxon_key") != key:
            raise ValueError("Country source accepted key differs from original")
    aliases = collision_species(records)
    if aliases != {"Oeneis bore", "Oeneis polixenes"}:
        raise ValueError("Original accepted-key collision set changed")
    complete = {s: set() for s in SOURCES}
    expected = {}
    query_lookup = {(q["species"], q["basis_scope"], q["window"]): q for q in received}
    for row in records:
        for scope in SOURCES:
            detail = row["basis_scopes"][scope]
            windows = [query_lookup.get((row["species"], scope, w)) for w in WINDOWS[1:]]
            actual_complete = all(q and q.get("status") == "COMPLETE" and q.get("complete") is True for q in windows)
            if (detail["status"] == "COMPLETE") != actual_complete:
                raise ValueError("Species completeness disagrees with both window receipts")
            if not actual_complete:
                continue
            if row["species"] not in aliases:
                complete[scope].add(row["species"])
            for q in windows:
                for country, n in q["country_counts"].items():
                    expected[(row["species"], q["window"], country, scope)] = n
    actual = {}
    with Path(count_path).open(newline="") as f:
        for row in csv.DictReader(f):
            key = (row["species"], row["window"], row["country"], row["basis_scope"])
            if key in actual:
                raise ValueError("Duplicate country aggregate row")
            n = int(row["record_count"])
            if n < 0:
                raise ValueError("Negative indexed count")
            actual[key] = n
    if actual != expected:
        raise ValueError("Country CSV differs from complete-window query country facets")
    counts = {(s, w): Counter() for s in SOURCES for w in WINDOWS}
    for (sp, window, country, scope), n in actual.items():
        if sp in complete[scope]:
            counts[scope, window][sp, country] += n
            counts[scope, "pooled"][sp, country] += n
    original_complete = {r["species"] for r in original["species_receipts"] if r["status"] == "COMPLETE"} - aliases
    return counts, complete, original_complete, receipt, sorted(aliases)


def load_geography(summary_path, links_path, regions, atoms_path, atoms_summary_path):
    summary = json.loads(Path(summary_path).read_text())
    if summary.get("schema") != "country_wgsrpd_area_crosswalk_v1":
        raise ValueError("Unexpected country crosswalk schema")
    if summary["output_links_sha256"] != file_sha(links_path):
        raise ValueError("Country crosswalk hash mismatch")
    if summary["thresholds"]["material_min_intersection_over_smaller_polygon"] != .05 or summary["thresholds"]["guard_min_intersection_over_smaller_polygon"] != .005:
        raise ValueError("Crosswalk thresholds differ from the source-blind protocol")
    links = defaultdict(dict)
    with Path(links_path).open(newline="") as f:
        for row in csv.DictReader(f):
            if row["wgsrpd3"] in links[row["country"]]:
                raise ValueError("Duplicate country/WGSRPD3 link")
            links[row["country"]][row["wgsrpd3"]] = row
    index = {r: i for i, r in enumerate(regions)}
    atoms_summary = json.loads(Path(atoms_summary_path).read_text())
    for key, path in (("atoms_csv_sha256", atoms_path),
                      ("crosswalk_summary_sha256", summary_path),
                      ("crosswalk_links_sha256", links_path)):
        if atoms_summary.get(key) != file_sha(path):
            raise ValueError("Disjoint area-atom provenance mismatch: " + key)
    atoms = defaultdict(list)
    with Path(atoms_path).open(newline="") as f:
        seen_atoms = set()
        for row in csv.DictReader(f):
            members = frozenset(row["wgsrpd3_members"].split("|"))
            atom_key = (row["country"], members)
            if atom_key in seen_atoms or not members or float(row["area_km2"]) <= 0:
                raise ValueError("Duplicate or invalid disjoint country area atom")
            seen_atoms.add(atom_key)
            atoms[row["country"]].append((members, float(row["area_km2"])))
    countries, excluded = {}, {}
    for row in summary["countries"]:
        code = row["country"]
        material = set(row["material_regions"])
        guard = set(row["guard_regions"])
        reasons = []
        if not row["iso_eligible"]:
            reasons.append("country_code_or_boundary_ambiguity")
        if row["material_coverage"] < .95:
            reasons.append("material_geographic_coverage_below_95pct")
        if len(row["material_level1_codes"]) != 1:
            reasons.append("not_exactly_one_material_Level1")
        if not material or not material <= guard or not guard <= set(index):
            reasons.append("unknown_or_unavailable_climate_geography")
        if reasons:
            excluded[code] = {"reasons": reasons, "names": row["names"],
                              "material_coverage": row["material_coverage"],
                              "material_regions": sorted(material), "mapping_flags": row["mapping_flags"]}
            continue
        if not atoms[code] or any(not members <= guard for members, _ in atoms[code]):
            raise ValueError("Country atoms do not match the guard-region crosswalk")
        union_area = sum(area for _, area in atoms[code])
        expected_area = row["country_area_km2"] * row["guard_coverage"]
        if not np.isclose(union_area, expected_area, rtol=1e-7, atol=1e-6):
            raise ValueError("Country atom areas do not reproduce exact guard union")
        countries[code] = {
            "country": code, "names": row["names"], "area": row["country_area_km2"],
            "level1": row["material_level1_codes"][0],
            "material": {index[r]: float(links[code][r]["intersection_km2"]) for r in material},
            "guard": {index[r]: float(links[code][r]["intersection_km2"]) for r in guard},
            "atoms": [{"members": {index[r] for r in members}, "area": area} for members, area in atoms[code]],
            "material_coverage": row["material_coverage"],
            "material_overlap_excess_fraction": row["material_overlap_excess_fraction"],
        }
    return countries, {"eligible_country_count": len(countries), "eligible_countries": sorted(countries),
                       "countries_excluded": excluded,
                       "exclusion_reason_counts": dict(Counter(r for v in excluded.values() for r in v["reasons"])),
                       "eligible_country_Level1_counts": dict(Counter(v["level1"] for v in countries.values())),
                       "max_eligible_material_overlap_excess": max((v["material_overlap_excess_fraction"] for v in countries.values()), default=0),
                       "crosswalk_source": summary["source"],
                       "country_area_denominator": "entire dissolved map-unit country area; unmapped remainder retained",
                       "polygon_overlap_rule": "exact disjoint membership atoms; each geographic portion counted once for any selected region subset",
                       "atoms_summary_sha256": file_sha(atoms_summary_path), "atoms_csv_sha256": file_sha(atoms_path)}


def classify_country(country, *, native, added, climate):
    def fraction(mask, weights):
        chosen = {j for j in weights if mask >> j & 1}
        union_area = sum(atom["area"] for atom in country["atoms"] if atom["members"] & chosen)
        value = union_area / country["area"]
        if not 0 <= value <= 1 + 1e-8:
            raise ValueError("Invalid exact country subset union area")
        return min(1.0, value)
    native_present = any(native >> j & 1 for j in country["guard"])
    return {
        "native_guard_present": native_present,
        "native_fraction": fraction(native, country["guard"]),
        "added_country": not native_present and any(added >> j & 1 for j in country["material"]),
        "added_fraction": fraction(added, country["material"]),
        "climate_fraction": fraction(climate, country["material"]),
        "added_climate_fraction": fraction(added & climate, country["material"]),
        "all_material_climate": all(climate >> j & 1 for j in country["material"]),
    }


def prepare_design(source_root, crosswalk_summary, crosswalk_links, crosswalk_atoms, atoms_summary):
    model, regions, climate, provenance = source_augmented_model(source_root)
    bysp, _, _, _ = climate_and_plant_candidates(model, regions, climate)
    countries, geography = load_geography(crosswalk_summary, crosswalk_links, regions, crosswalk_atoms, atoms_summary)
    exposures, audit_rows = {}, []
    for sp, features in bysp.items():
        for q in (.05, .10):
            mask = features["masks"][q]
            if mask is None:
                continue
            for code, country in countries.items():
                result = classify_country(country, native=features["native"], added=features["added"], climate=mask)
                exposures[sp, code, q] = result
                if q == .05:
                    audit_rows.append({"species": sp, "country": code, "Level1": country["level1"],
                                       "country_area_km2": country["area"], **result})
    design = {"executed_utc": datetime.now(timezone.utc).isoformat(), "geography": geography,
              "host_climate_source_species": len(bysp),
              "host_climate_eligible_species": len({k[0] for k in exposures}),
              "source_provenance": provenance,
              "source_blind_candidate_counts": {str(q): {
                  "outside_native_and_90pct_analogue_cells": sum(not x["native_guard_present"] and x["climate_fraction"] >= .9 - 1e-12 for k, x in exposures.items() if k[2] == q),
                  "added_outside_native_and_90pct_analogue_cells": sum(x["added_country"] and x["climate_fraction"] >= .9 - 1e-12 for k, x in exposures.items() if k[2] == q),
              } for q in (.05, .10)}}
    return countries, exposures, design, audit_rows


def build_strata(countries, exposures, counts, *, targets, effort_species, gate,
                 q, record_minimum, geography_mode="primary", stratification="level1"):
    effort_names = defaultdict(set)
    effort_records = Counter()
    for (sp, code), n in counts.items():
        if n > 0 and sp in effort_species and code in countries:
            effort_names[code].add(sp)
            effort_records[code] += n
    groups = defaultdict(list)
    audit = Counter()
    for sp in sorted(targets):
        for code, country in countries.items():
            x = exposures.get((sp, code, q))
            if x is None or x["native_guard_present"] or x["climate_fraction"] < .9 - 1e-12:
                continue
            if geography_mode == "all_climate" and not x["all_material_climate"]:
                continue
            if geography_mode == "single_region" and len(country["material"]) != 1:
                continue
            if geography_mode == "added90" and x["added_country"] and x["added_fraction"] < .9 - 1e-12:
                audit["excluded_intermediate_added_cells"] += 1
                continue
            richness = len(effort_names[code] - {sp})
            if richness < gate:
                continue
            other_records = effort_records[code] - (counts.get((sp, code), 0) if sp in effort_species else 0)
            group = country["level1"]
            if stratification == "level1_area":
                group += "|area" + str(0 if country["area"] < 10000 else 1 if country["area"] < 100000 else 2)
            elif stratification == "level1_record_intensity":
                group += "|records" + str(0 if other_records < 100 else 1 if other_records < 10000 else 2)
            elif stratification != "level1":
                raise ValueError("Unknown country stratification")
            groups[sp, group].append({"country": code, "added": int(x["added_country"]),
                                     "positive": int(counts.get((sp, code), 0) >= record_minimum),
                                     "area_score": x["added_climate_fraction"],
                                     "other_species_richness": richness, "other_species_records": other_records})
            audit["candidate_country_species_cells"] += 1
    strata = []
    for (sp, group), cells in sorted(groups.items()):
        N = len(cells)
        A = sum(c["added"] for c in cells)
        K = sum(c["positive"] for c in cells)
        O = sum(c["added"] * c["positive"] for c in cells)
        if 0 < A < N:
            strata.append({"species": sp, "genus": sp.split()[0], "stratum": group,
                           "N": N, "A": A, "K": K, "O": O, "expected": A * K / N,
                           "cells": cells})
    return strata, dict(audit)


def summarize_strata(strata, *, label, null_draws, bootstrap_draws, continuous=False):
    information = stratum_information_audit(strata)
    totalK = sum(x["K"] for x in strata)
    species = {x["species"] for x in strata}
    totalO = sum(x["O"] for x in strata)
    expected = sum(x["expected"] for x in strata)
    if continuous:
        totalO = sum(sum(c["area_score"] * c["positive"] for c in x["cells"]) for x in strata)
        expected = sum(x["K"] * np.mean([c["area_score"] for c in x["cells"]]) for x in strata)
    informative = information["nonzero_conditional_variance_strata"]
    if continuous:
        informative = sum(0 < x["K"] < x["N"] and
                          np.ptp([c["area_score"] for c in x["cells"]]) > 1e-12 for x in strata)
    enough = (informative >= 50 and len(species) >= 25 and totalK >= 25)
    out = {
        "status": "EXPLORATORY_MEETS_INFORMATION_THRESHOLD" if enough else "BELOW_PROTOCOL_INFORMATION_THRESHOLD",
        "label": label, "information": information, "species_with_resource_contrast": len(species),
        "species_with_positive_country": len({x["species"] for x in strata if x["K"] > 0}),
        "candidate_country_species_n": sum(x["N"] for x in strata),
        "positive_country_species_n": totalK, "candidate_added_country_species_n": sum(x["A"] for x in strata),
        "observed_overlap": float(totalO), "conditional_expected_overlap": float(expected),
        "observed_fraction": float(totalO / totalK) if totalK else None,
        "expected_fraction": float(expected / totalK) if totalK else None,
        "observed_minus_expected_pp": float(100 * (totalO - expected) / totalK) if totalK else None,
        "statistic": "added_and_analogue_country_area_fraction" if continuous else "binary_added_country_fraction",
        "interpretation": "source-conditional geographical benchmark; not a biological absence, causal, or breeding test",
    }
    if continuous:
        out["continuous_score_informative_strata"] = int(informative)
        out["information_gate_reference"] = ">=50 strata with 0<K<N and nonconstant country area score; information field also retains binary-contrast audit"
    if not totalK or not strata:
        return out
    rng = np.random.default_rng(SEED + int(hashlib.sha256(label.encode()).hexdigest()[:8], 16))
    sim = np.zeros(null_draws, dtype=float)
    bygenus = defaultdict(lambda: np.zeros(3))
    bylevel = defaultdict(lambda: np.zeros(3))
    bycountry = defaultdict(lambda: np.zeros(3))
    for row in strata:
        obs, exp = row["O"], row["expected"]
        if continuous:
            scores = np.array([c["area_score"] for c in row["cells"]])
            obs = sum(c["area_score"] * c["positive"] for c in row["cells"])
            exp = row["K"] * float(scores.mean())
            if row["K"] == row["N"]:
                sim += scores.sum()
            elif row["K"] > 0:
                # Random ranks generate a uniform K-subset without replacement.
                chosen = np.argpartition(rng.random((null_draws, row["N"])), row["K"] - 1, axis=1)[:, :row["K"]]
                sim += scores[chosen].sum(axis=1)
        else:
            sim += rng.hypergeometric(row["K"], row["N"] - row["K"], row["A"], size=null_draws)
        bygenus[row["genus"]] += (obs, exp, row["K"])
        bylevel[row["stratum"].split("|")[0]] += (obs, exp, row["K"])
        for c in row.get("cells", []):
            score = c["area_score"] if continuous else c["added"]
            bycountry[c["country"]] += (score * c["positive"], score * row["K"] / row["N"], c["positive"])
    values = np.array([bygenus[g] for g in sorted(bygenus)])
    picks = rng.integers(0, len(values), (bootstrap_draws, len(values)))
    bootstrap = values[picks].sum(axis=1)
    bootstrap = bootstrap[bootstrap[:, 2] > 0]
    effects = 100 * (bootstrap[:, 0] - bootstrap[:, 1]) / bootstrap[:, 2]
    out.update(
        null_draws=null_draws,
        conditional_null_right_tail=float((np.sum(sim >= totalO - 1e-12) + 1) / (null_draws + 1)),
        null_q025_q500_q975=[float(np.quantile(sim, q)) for q in (.025, .5, .975)],
        genus_bootstrap_clusters=len(values), genus_bootstrap_valid_draws=len(effects),
        genus_bootstrap_excess_pp_95ci=[float(np.quantile(effects, q)) for q in (.025, .975)] if len(effects) else None,
        Level1_components={k: {"observed": float(v[0]), "expected": float(v[1]), "positive": int(v[2])} for k, v in sorted(bylevel.items())},
        country_components={k: {"observed": float(v[0]), "expected": float(v[1]), "positive": int(v[2])} for k, v in sorted(bycountry.items())},
    )
    return out


def save_design(output_dir, design, rows):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "country_design_audit.json").write_text(json.dumps(design, indent=2, sort_keys=True) + "\n")
    if rows:
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        (output_dir / "country_resource_area_audit.csv.gz").write_bytes(gzip.compress(stream.getvalue().encode(), mtime=0))


def run(args):
    countries, exposures, design, rows = prepare_design(args.source_root, args.crosswalk_summary, args.crosswalk_links, args.crosswalk_atoms, args.atoms_summary)
    save_design(args.output_dir, design, rows)
    if args.design_only:
        print(json.dumps({"status": "SOURCE_BLIND_DESIGN_ONLY", "countries": len(countries), "candidates": design["source_blind_candidate_counts"]}))
        return
    if not args.gbif_counts or not args.gbif_manifest:
        raise ValueError("Both GBIF counts and manifest are required outside design-only mode")
    counts, complete, original_complete, receipt, aliases = load_observations(args.gbif_counts, args.gbif_manifest, args.original_manifest)
    results = {}
    primary_strata = None
    def evaluate(scope, window, gate, q, threshold, *, label, targets=None, effort=None, geography_mode="primary", stratification="level1", primary=False):
        strata, audit = build_strata(countries, exposures, counts[scope, window],
            targets=complete[scope] if targets is None else targets,
            effort_species=complete[scope] if effort is None else effort,
            gate=gate, q=q, record_minimum=threshold, geography_mode=geography_mode, stratification=stratification)
        result = summarize_strata(strata, label=label,
            null_draws=args.primary_null_draws if primary else args.sensitivity_null_draws,
            bootstrap_draws=2999 if primary else 499)
        result.update(source_basis=scope, window=window, other_species_gate=gate,
            climate_lower_quantile=q, focal_record_threshold=threshold, geography_mode=geography_mode,
            stratification=stratification, pre_contrast_candidate_audit=audit,
            target_species_count=len(complete[scope] if targets is None else targets),
            effort_species_count=len(complete[scope] if effort is None else effort))
        return result, strata
    for settings in product(SOURCES, WINDOWS, (1, 3, 5), (.05, .10), (1, 5)):
        scope, window, gate, q, threshold = settings
        label = f"{scope}|{window}|E{gate}|q{q:.2f}|records{threshold}"
        result, strata = evaluate(*settings, label=label, primary=settings == PRIMARY)
        results[label] = result
        if settings == PRIMARY:
            primary_strata = strata
            primary_result = result
    diagnostics = {}
    for mode in ("all_climate", "single_region", "added90"):
        diagnostics[mode], _ = evaluate(*PRIMARY, label=mode, geography_mode=mode)
    for stratification in ("level1_area", "level1_record_intensity"):
        diagnostics[stratification], _ = evaluate(*PRIMARY, label=stratification, stratification=stratification)
    diagnostics["continuous_added_analogue_area"] = summarize_strata(primary_strata,
        label="continuous_added_analogue_area", null_draws=args.sensitivity_null_draws,
        bootstrap_draws=499, continuous=True)
    original_available = original_complete & complete["all_accepted_bases"]
    diagnostics["original_targets_and_effort"], _ = evaluate(*PRIMARY,
        label="original_targets_and_effort", targets=original_available, effort=original_available)
    diagnostics["all_targets_original_effort"], _ = evaluate(*PRIMARY,
        label="all_targets_original_effort", effort=original_available)
    record_audit = {}
    for scope in SOURCES:
        c = counts[scope, "pooled"]
        record_audit[scope] = {
            "source_complete_noncolliding_species": len(complete[scope]),
            "indexed_records_with_country": sum(c.values()),
            "indexed_records_in_eligible_countries": sum(n for (sp, country), n in c.items() if country in countries),
            "countries_with_records": len({country for (sp, country), n in c.items() if n > 0}),
            "eligible_countries_with_records": len({country for (sp, country), n in c.items() if n > 0 and country in countries}),
            "record_totals_by_country": dict(sorted(Counter({country: sum(n for (sp, cc), n in c.items() if cc == country) for country in {cc for _, cc in c}}).items())),
        }
    result = {
        "schema": "chocho_full_index_country_host_opportunity_20261020_v01",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "status": primary_result["status"], "primary": primary_result,
        "all_72_prespecified_settings": results, "diagnostics": diagnostics,
        "country_design": design, "source_record_audit": record_audit,
        "source_complete_species_by_basis": {k: sorted(v) for k, v in complete.items()},
        "source_excluded_key_collision_names": aliases,
        "source_unresolved_names": sorted(r["species"] for r in receipt["species_receipts"] if r["status"] == "TAXON_UNRESOLVED"),
        "source_query_statuses": receipt["query_statuses"],
        "provenance": {"original_manifest_sha256": ORIGINAL_SHA, "counts_sha256": file_sha(args.gbif_counts),
            "country_receipt_sha256": file_sha(args.gbif_manifest), "source_response_set_sha256": receipt["source_response_set_sha256"],
            "crosswalk_summary_sha256": file_sha(args.crosswalk_summary), "crosswalk_links_sha256": file_sha(args.crosswalk_links),
            "crosswalk_atoms_sha256": file_sha(args.crosswalk_atoms), "atoms_summary_sha256": file_sha(args.atoms_summary),
            "analysis_script_sha256": file_sha(__file__)},
        "limitations": [
            "Independent country-grain analysis; not recovery of the bounded WGSRPD3 source or a same-grain effect comparison.",
            "Maximum229 non-colliding taxa; original eight unresolved taxa and two shared-key names are not ecological negatives.",
            "Strict country mapping excludes major multi-Level1 countries, territories and ambiguous boundaries; no unrestricted global estimate.",
            "Country co-coverage does not establish local butterfly-host co-occurrence or larval use.",
            "No biological absence, colonisation, population change, reproduction, fitness, or causal introduction effect.",
            "Live GBIF country facets reflect collection and reporting biases, not independent individuals or measured effort.",
            "Host climate is a potential-resource proxy; coarse units spanning countries can misclassify local resources.",
            "Genus bootstrap does not capture all shared host, country, observer, or phylogenetic dependence.",
        ],
    }
    (args.output_dir / "gbif_country_macro_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (args.output_dir / "gbif_country_primary_strata.json").write_text(json.dumps(primary_strata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "primary": primary_result,
        "source": {k: {kk: vv for kk, vv in v.items() if kk != "record_totals_by_country"} for k, v in record_audit.items()},
        "diagnostics": {k: {f: v.get(f) for f in ("status", "observed_minus_expected_pp", "positive_country_species_n", "genus_bootstrap_excess_pp_95ci")} for k, v in diagnostics.items()}}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--original-manifest", type=Path, required=True)
    parser.add_argument("--gbif-counts", type=Path)
    parser.add_argument("--gbif-manifest", type=Path)
    parser.add_argument("--crosswalk-summary", type=Path, required=True)
    parser.add_argument("--crosswalk-links", type=Path, required=True)
    parser.add_argument("--crosswalk-atoms", type=Path, required=True)
    parser.add_argument("--atoms-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--design-only", action="store_true")
    parser.add_argument("--primary-null-draws", type=int, default=9999)
    parser.add_argument("--sensitivity-null-draws", type=int, default=999)
    args = parser.parse_args()
    if min(args.primary_null_draws, args.sensitivity_null_draws) < 1:
        parser.error("Null draw counts must be positive")
    run(args)
