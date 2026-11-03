#!/usr/bin/env python3
"""Q2-A: trace island gains to introduced larval-hosts' *native source regions*.

Four European metropole references follow Lenzner et al. 2022; these are
plant native-range associations, NOT evidence of historical transport routes,
nor evidence that butterflies themselves originate in those metropoles.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict, Counter
from pathlib import Path

EMPIRE_HOMELAND = {
    "Great_Britain": {"GRB"},
    "Spain": {"SPA"},
    "Portugal": {"POR"},
    "Netherlands": {"NET"},
}
# France is recorded as a comparison, not as one of the four Lenzner et al. empire exposures.
OTHER_HOMELAND = {"France": {"FRA"}}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    for x in ["panel-csv", "host-csv", "native-csv", "contemporary-csv",
              "q1-summary", "level1-json", "output-dir"]:
        p.add_argument("--" + x, required=True, type=Path)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    q1 = json.loads(a.q1_summary.read_text(encoding="utf-8"))
    panel = list(read_csv(a.panel_csv))
    butterflies = {r["species"] for r in panel}
    if len(panel) != 239 or q1["panel_total"]["n_panel_mismatches"]:
        raise RuntimeError("Q1 source-verified panel is not the frozen 239-species matrix")
    if len(q1["classification"]["strict_island_codes"]) != 126 or "CNY" not in q1["classification"]["strict_island_codes"]:
        raise RuntimeError("Q2 v0.2 requires the 126-island corrected Q1, including Canary Islands")
    expected_island_added = int(q1["groups"]["island"]["added_butterfly_region_units"])
    if expected_island_added != 5741:
        raise RuntimeError(f"Q2 v0.2 geographic correction not reconciled: {expected_island_added}")
    # Original validation is unchanged; v0.2 only changes geography labels.
    level1 = json.loads(a.level1_json.read_text(encoding="utf-8"))["units"]
    native_host = defaultdict(set)
    contemporary_host = defaultdict(set)
    for row in read_csv(a.native_csv):
        native_host[row["accepted_plant_name_id"]].add(row["area_code_l3"])
    for row in read_csv(a.contemporary_csv):
        contemporary_host[row["accepted_plant_name_id"]].add(row["area_code_l3"])
    species_hosts = defaultdict(set)
    for row in read_csv(a.host_csv):
        if row["insect_species"] in butterflies:
            species_hosts[row["insect_species"]].add(row["accepted_plant_name_id"])
    missing_native = [(h, code) for h, regs in native_host.items()
                      for code in regs if code not in contemporary_host.get(h, ())]
    if missing_native: raise RuntimeError(f"Native hosts not subset contemporary: {missing_native[:10]}")

    # Reconstruct all 239 x WGSRPD3 resource cells from exact pinned sources.
    # No dependence on expired temporary GitHub Action output artifacts.
    added = []
    for record in panel:
        butterfly = record["species"]
        hosts = species_hosts.get(butterfly, ())
        n = set()
        c = set()
        for host in hosts:
            n.update(native_host.get(host, ()))
            c.update(contemporary_host.get(host, ()))
        if (not n <= c or len(n) != int(record["native_resource_units"])
                or len(c) != int(record["contemporary_resource_units"])):
            raise RuntimeError(f"Source footprint drift for {butterfly}: {len(n)=}, {len(c)=}")
        for code in sorted(c - n):
            added.append({
                "species": butterfly,
                "area_code_l3": code,
                "geographic_group": "island" if code in q1["classification"]["strict_island_codes"] else (
                    "mixed" if code in q1["classification"]["excluded_mixed_codes"] else "mainland"
                )
            })
    if len(added) != 14553:
        raise RuntimeError(f"Added source cells did not reproduce GEB 14,553: {len(added)}")

    # A true new resource cell has no native host in its island region but at
    # least one introduced host present in WCVP contemporary distribution.
    island_cells = [r for r in added if r["geographic_group"] == "island"]
    if len(island_cells) != expected_island_added: raise RuntimeError(f"Unexpected island added count: {len(island_cells)}")
    cell_rows = []
    edge_rows = []
    grouped = defaultdict(list)
    for row in island_cells:
        butterfly, island = row["species"], row["area_code_l3"]
        hosts = species_hosts.get(butterfly, set())
        all_present = [h for h in hosts if island in contemporary_host.get(h, ())]
        introduced = [h for h in all_present if island not in native_host.get(h, ())]
        if not introduced or len(introduced) != len(all_present):
            raise RuntimeError(f"Bad new island resource cell {butterfly} {island}")
        metropole = {emp: sum(bool(native_host.get(h, set()) & regions)
                             for h in introduced)
                     for emp, regions in EMPIRE_HOMELAND.items()}
        europe = sum(any(level1.get(k) == "1" for k in native_host.get(h, ()))
                     for h in introduced)
        other = {name: sum(bool(native_host.get(h, set()) & regions)
                           for h in introduced)
                 for name, regions in OTHER_HOMELAND.items()}
        record = {"species": butterfly, "island_code": island,
                  "level1": level1[island], "introduced_host_n": len(introduced),
                  "any_europe_native_host": int(europe > 0),
                  "europe_native_hosts_n": europe,
                  "four_metropole_native_host": int(any(metropole.values())),
                  **{f"{key}_native_host": int(n>0) for key,n in metropole.items()},
                  **{f"{key}_host_n": n for key,n in metropole.items()},
                  "France_native_host": int(other["France"] > 0)}
        cell_rows.append(record)
        grouped[island].append(record)
        for h in sorted(introduced):
            nregs = native_host.get(h, set())
            edge_rows.append({
                "species": butterfly, "island_code": island, "plant_id": h,
                "native_wgsrpd3_count": len(nregs),
                "native_europe": int(any(level1.get(k) == "1" for k in nregs)),
                **{f"native_{emp}": int(bool(nregs & regs))
                   for emp,regs in EMPIRE_HOMELAND.items()},
            })
    if not cell_rows or len(cell_rows) != expected_island_added:
        raise RuntimeError("New cells not fully decomposed into introduced plants")

    island_summary = []
    for island, rr in sorted(grouped.items()):
        row = {
            "island_code": island, "level1": level1[island],
            "new_butterfly_resource_cells": len(rr),
            "introduced_host_edges": sum(r["introduced_host_n"] for r in rr),
            "europe_native_supported_cells": sum(r["any_europe_native_host"] for r in rr),
            "four_metropole_supported_cells": sum(r["four_metropole_native_host"] for r in rr),
        }
        for emp in EMPIRE_HOMELAND:
            row[emp + "_native_supported_cells"] = sum(r[emp + "_native_host"] for r in rr)
        island_summary.append(row)
    all_cells = len(cell_rows)
    output = {
      "schema":"chocho_q2_island_introduced_host_native_source_v0.2",
      "status":"EMPIRICAL_SOURCE_ATTRIBUTION_WITHOUT_COLONIAL_EXPOSURE",
      "source_commits":q1["source_commits"],
      "corrected_island_definition":"126 strict WGSRPD3 islands including CNY Canary Is., 2026-11-02 source name audit",
      "scope":{
          "all_added_butterfly_region_cells": len(added),
          "added_strict_island_cells": all_cells,
          "islands_with_added_resource":len(grouped),
          "introduced_butterfly_plant_island_edges":len(edge_rows),
          "distinct_intro_host_taxa":len({r["plant_id"] for r in edge_rows}),
      },
      "observed":{
          "europe_native_supported_cells":sum(r["any_europe_native_host"] for r in cell_rows),
          "europe_native_supported_fraction":sum(r["any_europe_native_host"] for r in cell_rows)/all_cells,
          "any_four_metropole_supported_cells":sum(r["four_metropole_native_host"] for r in cell_rows),
          "four_metropole_source_cells":{emp:sum(r[emp+"_native_host"] for r in cell_rows) for emp in EMPIRE_HOMELAND},
          "other_homeland_source_cells":{"France":sum(r["France_native_host"] for r in cell_rows)},
          "new_cells_with_multiple_candidate_introduced_hosts":sum(r["introduced_host_n"]>1 for r in cell_rows),
      },
      "limits":[
        "A plant can be native in multiple regions; source categories overlap and sums cannot be read as mutually exclusive.",
        "Native to Britain, Spain, Portugal, or the Netherlands does not prove the metropole introduced it to a given island.",
        "Four-metropole attribution is a potential plant-source indicator, not a butterfly biogeographic origin.",
        "All 239 butterflies are the same original HOSTS–LepTraits study panel; their taxonomic and geographic sampling is not comprehensive.",
        "An introduced-host resource cell does not imply confirmed butterfly presence, colonization, or larval performance."
      ]
    }
    write_csv(a.output_dir/"island_source_summary.csv", island_summary, list(island_summary[0].keys()))
    write_csv(a.output_dir/"introduced_host_source_cells.csv", cell_rows, list(cell_rows[0].keys()))
    write_csv(a.output_dir/"introduced_host_source_edges.csv", edge_rows, list(edge_rows[0].keys()))
    (a.output_dir/"host_source_audit.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2)+"\n", encoding="utf-8")
    print(json.dumps(output,ensure_ascii=False,indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
