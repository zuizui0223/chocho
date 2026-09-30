#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import rankdata

EXPECTED_MECHANISM_SHA256 = (
    "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa"
)
EXPECTED_SPECIES_COUNT = 239
EXPECTED_SPECIES_SHA256 = (
    "c972e145e160f0c34a5fcee9f423e11c8e38072038e9765a40e79c0c5cd06157"
)
EXPECTED_REALM_LABELS = {
    "Neotropical",
    "Australian",
    "Afrotropical",
    "Madagascan",
    "Oceanian",
    "Oriental",
    "Panamanian",
    "Saharo-Arabian",
    "Nearctic",
    "Sino-Japanese",
    "Palearctic",
}


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def spearman(frame: pd.DataFrame, x: str, y: str) -> dict[str, object]:
    z = frame[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(z) < 3 or z[x].nunique() < 2 or z[y].nunique() < 2:
        return {"n": int(len(z)), "rho": None}
    rx = rankdata(z[x].to_numpy(dtype=float), method="average")
    ry = rankdata(z[y].to_numpy(dtype=float), method="average")
    rho = float(np.corrcoef(rx, ry)[0, 1])
    return {"n": int(len(z)), "rho": rho}


def family_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            family = str(row.get("Family", "")).strip()
            genus = str(row.get("Genus", "")).strip()
            species = str(row.get("Species", "")).strip()
            verbatim = str(row.get("verbatimSpecies", "")).strip()
            for name in (" ".join(x for x in (genus, species) if x), verbatim):
                if name and family:
                    out.setdefault(name, family)
    return out


def load_transport(
    occurrence_dir: Path,
    ledger_dir: Path,
    *,
    expected_shards: int,
) -> tuple[pd.DataFrame, list[dict]]:
    ledgers = []
    seen = set()
    for path in sorted(ledger_dir.rglob("*.json")):
        obj = read_json(path)
        if obj.get("schema") != "chocho_geb_occurrence_realm_shard_v0.1":
            continue
        idx = int(obj["shard_index"])
        if idx in seen:
            raise RuntimeError(f"duplicate shard ledger {idx}")
        seen.add(idx)
        ledgers.append(obj)
    if seen != set(range(expected_shards)):
        raise RuntimeError(
            f"shard ledger coverage drift: expected {expected_shards}, got {sorted(seen)}"
        )
    if any(obj["status"] != "COMPLETE" for obj in ledgers):
        raise RuntimeError("realm audit occurrence transport incomplete")

    csvs = sorted(occurrence_dir.rglob("occurrences*.csv"))
    if not csvs:
        raise RuntimeError("no occurrence shard CSVs")
    frames = [pd.read_csv(path) for path in csvs]
    occ = pd.concat(frames, ignore_index=True)
    required = {"species", "gbif_key", "latitude", "longitude"}
    if not required <= set(occ.columns):
        raise RuntimeError("occurrence schema drift")
    occ["species"] = occ["species"].astype(str)
    occ["gbif_key"] = pd.to_numeric(occ["gbif_key"], errors="raise").astype("int64")
    occ = occ.drop_duplicates(["species", "gbif_key"]).reset_index(drop=True)
    return occ, ledgers


def assign_realms(occ: pd.DataFrame, realm_geojson: Path) -> tuple[pd.DataFrame, dict]:
    realms = gpd.read_file(realm_geojson)
    if "Realm" not in realms.columns:
        raise RuntimeError(f"realm field missing: {realms.columns.tolist()}")
    labels = set(realms["Realm"].dropna().astype(str))
    if labels != EXPECTED_REALM_LABELS:
        raise RuntimeError(f"realm label drift: {sorted(labels)}")
    realms = realms[["Realm", "geometry"]].copy()
    realms = realms.to_crs("EPSG:4326")

    points = gpd.GeoDataFrame(
        occ.copy(),
        geometry=gpd.points_from_xy(occ["longitude"], occ["latitude"]),
        crs="EPSG:4326",
    )
    points["_occ_id"] = np.arange(len(points), dtype=np.int64)
    joined = gpd.sjoin(points, realms, how="left", predicate="intersects")

    records = []
    ambiguous = 0
    for occ_id, grp in joined.groupby("_occ_id", sort=True):
        vals = sorted(set(str(x) for x in grp["Realm"].dropna()))
        realm = None
        if len(vals) == 1:
            realm = vals[0]
        elif len(vals) > 1:
            ambiguous += 1
        row = points.loc[
            int(occ_id), ["species", "gbif_key", "longitude", "latitude"]
        ].to_dict()
        row["Realm"] = realm
        records.append(row)

    mapped = pd.DataFrame(records)
    diag = {
        "occurrence_rows": int(len(points)),
        "realm_mapped_rows": int(mapped["Realm"].notna().sum()),
        "realm_unmapped_or_ambiguous_rows": int(mapped["Realm"].isna().sum()),
        "ambiguous_multi_realm_boundary_rows": int(ambiguous),
        "mapped_fraction": float(mapped["Realm"].notna().mean()) if len(mapped) else None,
    }
    return mapped, diag


def summarize_species(
    mechanism: pd.DataFrame,
    occ: pd.DataFrame,
    assigned: pd.DataFrame,
    ledgers: list[dict],
) -> pd.DataFrame:
    status = {}
    for ledger in ledgers:
        for row in ledger["species"]:
            status[str(row["species"])] = str(row["status"])
    all_names = sorted(mechanism["species"].astype(str).unique())
    rows = []
    for name in all_names:
        total = int((occ["species"] == name).sum())
        sub = assigned[(assigned["species"] == name) & assigned["Realm"].notna()]
        counts = Counter(sub["Realm"].astype(str))
        mapped = int(sum(counts.values()))
        ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        primary = None
        primary_count = 0
        if ordered:
            primary_count = ordered[0][1]
            leaders = [realm for realm, n in ordered if n == primary_count]
            primary = leaders[0] if len(leaders) == 1 else "TIE"
        share = None if mapped == 0 else primary_count / mapped

        spatial = sub.copy()
        if len(spatial):
            spatial["grid_lon"] = np.floor(
                pd.to_numeric(spatial["longitude"], errors="raise") + 180.0
            ).astype(int)
            spatial["grid_lat"] = np.floor(
                pd.to_numeric(spatial["latitude"], errors="raise") + 90.0
            ).astype(int)
            spatial.loc[spatial["grid_lon"] == 360, "grid_lon"] = 359
            spatial.loc[spatial["grid_lat"] == 180, "grid_lat"] = 179
            cells = (
                spatial[["Realm", "grid_lon", "grid_lat"]]
                .drop_duplicates()
            )
            cell_counts = Counter(cells["Realm"].astype(str))
        else:
            cell_counts = Counter()
        cell_total = int(sum(cell_counts.values()))
        cell_ordered = sorted(
            cell_counts.items(), key=lambda kv: (-kv[1], kv[0])
        )
        cell_primary = None
        cell_primary_count = 0
        if cell_ordered:
            cell_primary_count = cell_ordered[0][1]
            cell_leaders = [
                realm for realm, n in cell_ordered if n == cell_primary_count
            ]
            cell_primary = (
                cell_leaders[0] if len(cell_leaders) == 1 else "TIE"
            )
        cell_share = (
            None if cell_total == 0 else cell_primary_count / cell_total
        )

        rows.append({
            "species": name,
            "transport_status": status.get(name, "MISSING_LEDGER"),
            "fetched_occurrences": total,
            "realm_mapped_occurrences": mapped,
            "realm_mapping_fraction": (mapped / total) if total else None,
            "realm_count": len(counts),
            "realm_counts_json": json.dumps(dict(ordered), sort_keys=True),
            "primary_realm": primary,
            "primary_realm_share": share,
            "realm_informative": mapped >= 30,
            "core_realm": mapped >= 30 and share is not None and share >= 0.60,
            "realm_occupied_1deg_cells": cell_total,
            "realm_cell_counts_json": json.dumps(
                dict(cell_ordered), sort_keys=True
            ),
            "cell_primary_realm": cell_primary,
            "cell_primary_realm_share": cell_share,
            "cell_core_realm": (
                mapped >= 30
                and cell_share is not None
                and cell_share >= 0.60
            ),
        })
    return pd.DataFrame(rows)


def correlations_by_group(
    frame: pd.DataFrame,
    group_col: str,
    *,
    min_interpret: int,
) -> dict[str, object]:
    out = {}
    for group, sub in frame.groupby(group_col, dropna=False):
        key = "NA" if pd.isna(group) else str(group)
        out[key] = {
            "n": int(len(sub)),
            "interpretable": bool(len(sub) >= min_interpret),
            "host_family_vs_effective": spearman(
                sub, "host_family_count", "effective_contributor_number"
            ),
            "host_family_vs_max_share": spearman(
                sub, "host_family_count", "maximum_single_host_fractional_share"
            ),
        }
    return out


def leave_one_group_out(
    frame: pd.DataFrame,
    group_col: str,
) -> dict[str, object]:
    out = {}
    for group in sorted(frame[group_col].dropna().astype(str).unique()):
        sub = frame[frame[group_col].astype(str) != group]
        out[group] = {
            "n": int(len(sub)),
            "host_family_vs_effective": spearman(
                sub, "host_family_count", "effective_contributor_number"
            ),
            "host_family_vs_max_share": spearman(
                sub, "host_family_count", "maximum_single_host_fractional_share"
            ),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--mechanism-csv", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--occurrence-dir", type=Path, required=True)
    ap.add_argument("--ledger-dir", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument("--realm-source-sha256", type=str, required=True)
    ap.add_argument("--expected-shards", type=int, default=40)
    ap.add_argument("--output-species-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = read_json(args.protocol)
    if protocol.get("schema") != "chocho_geb_occurrence_realm_protocol_v0.1":
        raise RuntimeError("unexpected realm protocol")
    if protocol.get("status") != (
        "FROZEN_POSTHOC_GENERALITY_AUDIT_BEFORE_NEW_GBIF_REALM_ACQUISITION"
    ):
        raise RuntimeError("realm protocol not frozen")

    if sha256_path(args.mechanism_csv) != EXPECTED_MECHANISM_SHA256:
        raise RuntimeError("mechanism source SHA drift")
    mechanism = pd.read_csv(args.mechanism_csv)
    names = sorted(mechanism["species"].astype(str))
    digest = hashlib.sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest()
    if len(names) != EXPECTED_SPECIES_COUNT or digest != EXPECTED_SPECIES_SHA256:
        raise RuntimeError("mechanism species panel drift")

    occ, ledgers = load_transport(
        args.occurrence_dir, args.ledger_dir, expected_shards=args.expected_shards
    )
    if not set(occ["species"]).issubset(set(names)):
        raise RuntimeError("occurrence species outside frozen panel")

    assigned, realm_diag = assign_realms(occ, args.realm_geojson)
    species_summary = summarize_species(mechanism, occ, assigned, ledgers)

    fam = family_map(args.leptraits_csv)
    mechanism["Family"] = mechanism["species"].map(fam)
    if mechanism["Family"].isna().any():
        missing = sorted(mechanism.loc[mechanism["Family"].isna(), "species"].astype(str))
        raise RuntimeError(f"unmapped butterfly families: {missing}")

    merged = mechanism.merge(species_summary, on="species", how="left", validate="one_to_one")
    adequate = merged[
        merged["host_taxonomy_lower_bound_adequate"].astype(str).str.lower().eq("true")
        & (pd.to_numeric(merged["introduced_added_units"], errors="coerce") > 0)
    ].copy()
    if len(adequate) != 191:
        raise RuntimeError(f"expected 191 expanded adequate species, got {len(adequate)}")

    informative = adequate[adequate["realm_informative"]].copy()
    record_realm_counts = informative["primary_realm"].value_counts(dropna=False)
    interpretable = record_realm_counts[record_realm_counts >= 10]
    cell_realm_counts = informative["cell_primary_realm"].value_counts(dropna=False)
    cell_interpretable = cell_realm_counts[cell_realm_counts >= 10]
    gate_pass = (
        len(informative) >= int(protocol["evaluable_gate"]["minimum_realm_informative_species"])
        and len(interpretable) >= int(protocol["evaluable_gate"]["minimum_interpretable_primary_realms"])
    )

    transport_status_counts = Counter()
    for ledger in ledgers:
        for row in ledger["species"]:
            transport_status_counts[str(row["status"])] += 1

    payload: dict[str, object] = {
        "schema": "chocho_geb_occurrence_realm_result_v0.1",
        "status": (
            "EXPLORATORY_POSTHOC_REALM_GENERALITY_AUDIT_COMPLETE"
            if gate_pass
            else "NOT_EVALUABLE_GEB_REALM_GENERALITY_AUDIT"
        ),
        "protocol_sha256": sha256_path(args.protocol),
        "realm_source_sha256": args.realm_source_sha256,
        "transport": {
            "shards": args.expected_shards,
            "status_counts": dict(sorted(transport_status_counts.items())),
            "unique_occurrence_rows": int(len(occ)),
            "species_with_occurrence_rows": int(occ["species"].nunique()),
        },
        "realm_mapping": realm_diag,
        "panel": {
            "species": 239,
            "expanded_host_taxonomy_adequate_species": 191,
            "realm_informative_expanded_adequate_species": int(len(informative)),
            "core_realm_expanded_adequate_species": int(
                adequate["core_realm"].fillna(False).sum()
            ),
            "record_weighted_primary_realm_counts_informative": {
                str(k): int(v) for k, v in record_realm_counts.items()
            },
            "record_weighted_interpretable_primary_realms_n_ge_10": {
                str(k): int(v) for k, v in interpretable.items()
            },
            "cell_weighted_primary_realm_counts_informative": {
                str(k): int(v) for k, v in cell_realm_counts.items()
            },
            "cell_weighted_interpretable_primary_realms_n_ge_10": {
                str(k): int(v) for k, v in cell_interpretable.items()
            },
        },
        "evaluable_gate": {
            "passed": bool(gate_pass),
            "minimum_realm_informative_species": 150,
            "minimum_interpretable_primary_realms": 4,
            "observed_realm_informative_species": int(len(informative)),
            "observed_interpretable_primary_realms": int(len(interpretable)),
        },
    }

    if gate_pass:
        core = adequate[adequate["core_realm"]].copy()
        payload["generality"] = {
            "overall_informative": {
                "host_family_vs_effective": spearman(
                    informative, "host_family_count", "effective_contributor_number"
                ),
                "host_family_vs_max_share": spearman(
                    informative, "host_family_count", "maximum_single_host_fractional_share"
                ),
            },
            "by_primary_realm": correlations_by_group(
                informative, "primary_realm", min_interpret=10
            ),
            "core_realm_sensitivity": {
                "n": int(len(core)),
                "by_primary_realm": correlations_by_group(
                    core, "primary_realm", min_interpret=10
                ),
            },
            "leave_one_primary_realm_out": leave_one_group_out(
                informative, "primary_realm"
            ),
            "family_x_primary_realm_counts": (
                informative.groupby(["Family", "primary_realm"])
                .size()
                .unstack(fill_value=0)
                .astype(int)
                .to_dict()
            ),
            "spatial_thinning_sensitivity": {
                "definition": (
                    "Primary realm determined by unique occupied 1-degree "
                    "longitude-latitude cells rather than raw record counts."
                ),
                "by_cell_primary_realm": correlations_by_group(
                    informative, "cell_primary_realm", min_interpret=10
                ),
                "cell_core_realm": {
                    "n": int(adequate["cell_core_realm"].fillna(False).sum()),
                    "by_cell_primary_realm": correlations_by_group(
                        adequate[adequate["cell_core_realm"]].copy(),
                        "cell_primary_realm",
                        min_interpret=10,
                    ),
                },
                "leave_one_cell_primary_realm_out": leave_one_group_out(
                    informative, "cell_primary_realm"
                ),
                "record_vs_cell_primary_realm_agreement": {
                    "n_compared": int(
                        informative[["primary_realm", "cell_primary_realm"]]
                        .dropna()
                        .shape[0]
                    ),
                    "same": int(
                        (
                            informative["primary_realm"].astype(str)
                            == informative["cell_primary_realm"].astype(str)
                        ).sum()
                    ),
                    "fraction_same": float(
                        (
                            informative["primary_realm"].astype(str)
                            == informative["cell_primary_realm"].astype(str)
                        ).mean()
                    ) if len(informative) else None,
                },
            },
        }

    args.output_species_csv.parent.mkdir(parents=True, exist_ok=True)
    species_summary.sort_values("species").to_csv(args.output_species_csv, index=False)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": payload["status"],
        "realm_informative": int(len(informative)),
        "interpretable_realms": int(len(interpretable)),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
