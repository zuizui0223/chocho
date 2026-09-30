#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED = {
    "resource_eligible": 239,
    "host_taxonomy_adequate": 215,
    "expanded": 206,
    "host_taxonomy_adequate_expanded": 191,
    "one_family_host_taxonomy_adequate_expanded": 82,
}
BANDS = ("1_family", "2_families", "3_to_5_families", "6plus_families")


def band(x: float) -> str:
    if x == 1:
        return "1_family"
    if x == 2:
        return "2_families"
    if 3 <= x <= 5:
        return "3_to_5_families"
    if x >= 6:
        return "6plus_families"
    raise ValueError(x)


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: (x[i], i))
    out = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and x[order[j]] == x[order[i]]:
            j += 1
        r = 0.5 * ((i + 1) + j)
        for k in order[i:j]:
            out[k] = r
        i = j
    return out


def corr(x, y):
    if len(x) < 3 or len(x) != len(y):
        return None
    mx, my = sum(x) / len(x), sum(y) / len(y)
    dx, dy = [v - mx for v in x], [v - my for v in y]
    sx, sy = sum(v * v for v in dx), sum(v * v for v in dy)
    if sx <= 0 or sy <= 0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / math.sqrt(sx * sy)


def spearman(x, y):
    return corr(ranks(x), ranks(y)) if len(x) >= 3 else None


def read_descriptors(path):
    out = {}
    with path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        need = {"species", "host_family_count", "host_wgsrpd3_unit_count", "resolved_host_species"}
        if not need <= set(r.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in r:
            s = row["species"].strip()
            out[s] = {
                "host_family_count": float(row["host_family_count"]),
                "native_descriptor_units": int(row["host_wgsrpd3_unit_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
            }
    if len(out) != 339:
        raise RuntimeError(f"expected 339 descriptors, got {len(out)}")
    return out


def read_families(path):
    out = {}
    with path.open(newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        if not {"Species", "Family"} <= set(r.fieldnames or ()):
            raise RuntimeError("LepTraits Species/Family missing")
        for row in r:
            s, fam = row["Species"].strip(), row["Family"].strip()
            if s and fam and s not in out:
                out[s] = fam
    return out


def read_pairs(path):
    out = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if not {"insect_species", "accepted_plant_name_id"} <= set(r.fieldnames or ()):
            raise RuntimeError("interaction schema drift")
        for row in r:
            s, h = row["insect_species"].strip(), row["accepted_plant_name_id"].strip()
            if s and h:
                out[s].add(h)
    return out


def read_units(path):
    out = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if not {"accepted_plant_name_id", "area_code_l3"} <= set(r.fieldnames or ()):
            raise RuntimeError("distribution schema drift")
        for row in r:
            h, u = row["accepted_plant_name_id"].strip(), row["area_code_l3"].strip()
            if h and u:
                out[h].add(u)
    return out


def _comma_decimal_code(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    return value.replace(",", ".")


def level1_crosswalk(level1_path, level2_path, level3_path):
    level1 = {}
    with level1_path.open(encoding="utf-8-sig") as f:
        header = next(f, None)
        for line in f:
            parts = line.rstrip("\r\n").split("*")
            if len(parts) >= 2:
                code = _comma_decimal_code(parts[0])
                name = parts[1].strip()
                if code and name:
                    level1[code] = name

    level2_to_level1 = {}
    with level2_path.open(encoding="utf-8-sig") as f:
        header = next(f, None)
        for line in f:
            parts = line.rstrip("\r\n").split("*")
            if len(parts) >= 3:
                l2 = _comma_decimal_code(parts[0])
                l1_code = _comma_decimal_code(parts[2])
                if l2 and l1_code:
                    level2_to_level1[l2] = l1_code

    out = {}
    with level3_path.open(encoding="utf-8-sig") as f:
        header = next(f, None)
        for line in f:
            parts = line.rstrip("\r\n").split("*")
            if len(parts) >= 3:
                l3 = parts[0].strip()
                l2 = _comma_decimal_code(parts[2])
                l1_code = level2_to_level1.get(l2)
                l1_name = level1.get(l1_code or "")
                if l3 and l1_name:
                    if l3 in out and out[l3] != l1_name:
                        raise RuntimeError(f"WGSRPD crosswalk conflict for {l3}")
                    out[l3] = l1_name
    if len(out) < 300:
        raise RuntimeError(f"WGSRPD hierarchy crosswalk unexpectedly small: {len(out)}")
    return out


def contribution(host_added, target):
    if not target:
        return None, None, 0
    credit = defaultdict(float)
    for u in target:
        hs = [h for h, units in host_added.items() if u in units]
        if not hs:
            raise RuntimeError(f"added unit without host contributor: {u}")
        w = 1 / len(hs)
        for h in hs:
            credit[h] += w
    shares = [v / len(target) for v in credit.values()]
    return 1 / sum(v * v for v in shares), max(shares), len(shares)


def summarize_subset(names, rows):
    rr = [rows[s] for s in sorted(names)]
    fam = Counter(r["butterfly_family"] for r in rr)
    hb = Counter(r["host_breadth_stratum"] for r in rr)
    n = len(rr)
    probs = [v / n for v in fam.values()] if n else []
    return {
        "species": n,
        "butterfly_family_counts": dict(sorted(fam.items())),
        "host_breadth_strata": {k: hb.get(k, 0) for k in BANDS},
        "largest_family_fraction": None if not n else max(fam.values()) / n,
        "effective_butterfly_family_number": None if not probs else 1 / sum(p * p for p in probs),
    }


def grouped_replication(rows, field, rule):
    groups = defaultdict(list)
    for row in rows:
        groups[row[field]].append(row)
    out = []
    for name, rr in sorted(groups.items()):
        use = [r for r in rr if r["host_taxonomy_adequate"] and r["introduced_added_units"] > 0]
        counts = Counter(r["host_breadth_stratum"] for r in use)
        n_good_bands = sum(counts.get(k, 0) >= rule["min_per_stratum"] for k in BANDS)
        x_eff = [r["host_family_count"] for r in use if r["effective_contributor_number"] is not None]
        y_eff = [r["effective_contributor_number"] for r in use if r["effective_contributor_number"] is not None]
        x_max = [r["host_family_count"] for r in use if r["maximum_single_host_fractional_share"] is not None]
        y_max = [r["maximum_single_host_fractional_share"] for r in use if r["maximum_single_host_fractional_share"] is not None]
        out.append({
            field: name,
            "species": len(rr),
            "adequate_expanded_species": len(use),
            "host_breadth_strata": {k: counts.get(k, 0) for k in BANDS},
            "informative_for_directional_replication": (
                len(use) >= rule["min_total"]
                and n_good_bands >= rule["min_informative_strata"]
                and len(set(r["host_family_count"] for r in use)) >= 2
            ),
            "rho_host_family_vs_effective_contributors": spearman(x_eff, y_eff),
            "rho_host_family_vs_maximum_single_host_share": spearman(x_max, y_max),
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--level1-table", type=Path, required=True)
    ap.add_argument("--level2-table", type=Path, required=True)
    ap.add_argument("--level3-table", type=Path, required=True)
    ap.add_argument("--independent-panel-json", type=Path)
    ap.add_argument("--climate-result-json", type=Path)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-species-region-csv", type=Path, required=True)
    a = ap.parse_args()

    protocol = json.loads(a.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "chocho_geb_generality_audit_protocol_v0.1":
        raise RuntimeError("unexpected protocol")
    if protocol.get("status") != "FROZEN_BEFORE_GEB_GENERALITY_AUDIT_RESULT":
        raise RuntimeError("protocol is not frozen")

    desc = read_descriptors(a.descriptors_csv)
    fam = read_families(a.leptraits_csv)
    pairs = read_pairs(a.insect_host_csv)
    native = read_units(a.native_distribution_csv)
    contemp = read_units(a.contemporary_distribution_csv)
    l1 = level1_crosswalk(a.level1_table, a.level2_table, a.level3_table)

    species = {}
    region_rows = []
    missing_codes = set()

    for s, d in desc.items():
        hf = d["host_family_count"]
        if hf <= 0 or d["native_descriptor_units"] <= 0:
            continue
        if s not in fam:
            raise RuntimeError(f"missing LepTraits Family for {s}")
        hosts = sorted(pairs.get(s, ()))
        if not hosts:
            raise RuntimeError(f"missing resolved hosts for {s}")
        nu, cu, host_added = set(), set(), {}
        for h in hosts:
            n, c = set(native.get(h, ())), set(contemp.get(h, ()))
            if not n <= c:
                raise RuntimeError(f"native not subset contemporary: {h}")
            nu |= n
            cu |= c
            host_added[h] = c - n
        if len(nu) != d["native_descriptor_units"]:
            raise RuntimeError(f"native footprint drift for {s}")
        au = cu - nu
        eff, mx, nh = contribution(host_added, au)
        adequate = d["resolved_host_species"] >= hf
        hb = band(hf)
        species[s] = {
            "species": s,
            "butterfly_family": fam[s],
            "host_family_count": hf,
            "host_breadth_stratum": hb,
            "resolved_host_species": d["resolved_host_species"],
            "host_taxonomy_adequate": adequate,
            "native_resource_units": len(nu),
            "contemporary_resource_units": len(cu),
            "introduced_added_units": len(au),
            "effective_contributor_number": eff,
            "maximum_single_host_fractional_share": mx,
            "contributing_host_species": nh,
        }

        for u in nu | cu:
            if u not in l1:
                missing_codes.add(u)
        regions = sorted({l1[u] for u in nu | cu if u in l1})
        for reg in regions:
            nr = {u for u in nu if l1.get(u) == reg}
            cr = {u for u in cu if l1.get(u) == reg}
            ar = cr - nr
            local_added = {h: {u for u in units if l1.get(u) == reg} for h, units in host_added.items()}
            reff, rmx, rnh = contribution(local_added, ar)
            region_rows.append({
                "species": s,
                "butterfly_family": fam[s],
                "level1_region": reg,
                "host_family_count": hf,
                "host_breadth_stratum": hb,
                "resolved_host_species": d["resolved_host_species"],
                "host_taxonomy_adequate": adequate,
                "native_resource_units": len(nr),
                "contemporary_resource_units": len(cr),
                "introduced_added_units": len(ar),
                "regional_log1p_expansion": math.log1p(len(cr)) - math.log1p(len(nr)),
                "effective_contributor_number": reff,
                "maximum_single_host_fractional_share": rmx,
                "contributing_host_species": rnh,
            })

    if missing_codes:
        raise RuntimeError("WGSRPD3 codes missing from official crosswalk: " + ", ".join(sorted(missing_codes)[:20]))

    eligible = set(species)
    adequate = {s for s, r in species.items() if r["host_taxonomy_adequate"]}
    expanded = {s for s, r in species.items() if r["introduced_added_units"] > 0}
    ade_exp = adequate & expanded
    one = {s for s in ade_exp if species[s]["host_family_count"] == 1}
    funnel = {
        "resource_eligible": len(eligible),
        "host_taxonomy_adequate": len(adequate),
        "expanded": len(expanded),
        "host_taxonomy_adequate_expanded": len(ade_exp),
        "one_family_host_taxonomy_adequate_expanded": len(one),
    }
    if funnel != EXPECTED:
        raise RuntimeError(f"known funnel drift: {funnel} != {EXPECTED}")

    subsets = {
        "resource_eligible": eligible,
        "host_taxonomy_adequate": adequate,
        "expanded": expanded,
        "host_taxonomy_adequate_expanded": ade_exp,
        "one_family_host_taxonomy_adequate_expanded": one,
    }
    if a.independent_panel_json:
        p = json.loads(a.independent_panel_json.read_text(encoding="utf-8"))
        ss = set(map(str, p.get("species", [])))
        if len(ss) != 32:
            raise RuntimeError("independent panel != 32")
        subsets["independent_panel"] = ss
    if a.climate_result_json:
        p = json.loads(a.climate_result_json.read_text(encoding="utf-8"))
        ss = {str(r["species"]) for r in p.get("species", [])}
        if len(ss) != 24:
            raise RuntimeError("climate informative != 24")
        subsets["climate_informative"] = ss

    funnel_rep = {k: summarize_subset(v, species) for k, v in subsets.items()}
    reg_rep = grouped_replication(region_rows, "level1_region", protocol["informative_thresholds"]["level1_region"])
    fam_rep = grouped_replication(list(species.values()), "butterfly_family", protocol["informative_thresholds"]["butterfly_family"])
    ir = [r for r in reg_rep if r["informative_for_directional_replication"]]
    iff = [r for r in fam_rep if r["informative_for_directional_replication"]]

    result = {
        "schema": "chocho_geb_generality_audit_result_v0.1",
        "status": "PROSPECTIVE_GENERALITY_AUDIT_COMPLETE",
        "known_funnel_reproduced": funnel,
        "funnel_taxonomic_representation": funnel_rep,
        "level1_region_replication": reg_rep,
        "butterfly_family_replication": fam_rep,
        "replication_summary": {
            "informative_level1_regions": len(ir),
            "regions_expected_positive_effective_contributors": sum((r["rho_host_family_vs_effective_contributors"] or 0) > 0 for r in ir),
            "regions_expected_negative_max_host_share": sum(r["rho_host_family_vs_maximum_single_host_share"] is not None and r["rho_host_family_vs_maximum_single_host_share"] < 0 for r in ir),
            "informative_butterfly_families": len(iff),
            "families_expected_positive_effective_contributors": sum((r["rho_host_family_vs_effective_contributors"] or 0) > 0 for r in iff),
            "families_expected_negative_max_host_share": sum(r["rho_host_family_vs_maximum_single_host_share"] is not None and r["rho_host_family_vs_maximum_single_host_share"] < 0 for r in iff),
        },
        "claim_boundary": protocol["claim_boundary"],
    }

    a.output_json.parent.mkdir(parents=True, exist_ok=True)
    a.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields = [
        "species", "butterfly_family", "level1_region", "host_family_count",
        "host_breadth_stratum", "resolved_host_species", "host_taxonomy_adequate",
        "native_resource_units", "contemporary_resource_units", "introduced_added_units",
        "regional_log1p_expansion", "effective_contributor_number",
        "maximum_single_host_fractional_share", "contributing_host_species",
    ]
    with a.output_species_region_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(region_rows)

    print(json.dumps({
        "known_funnel_reproduced": funnel,
        "replication_summary": result["replication_summary"],
        "resource_eligible_family_counts": funnel_rep["resource_eligible"]["butterfly_family_counts"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
