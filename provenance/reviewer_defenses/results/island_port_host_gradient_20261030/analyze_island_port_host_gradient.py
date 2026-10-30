#!/usr/bin/env python3
"""Descriptive within-island harbor urbanization / known-host inventory gradient.

This is NOT a matched non-port botanical inventory, a history of plant arrival,
or evidence of butterfly establishment. Prespecified exploratory specification:
provenance/reviewer_defenses/plans/island_port_within_island_host_gradient_lock_20261030.md.
Only Python standard library is required.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = Path("provenance/reviewer_defenses/results")
PORT = SOURCE / "port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
VEG = SOURCE / "coastal_port_vegetation_20261030/coastal_control_vegetation_pairs.csv"
NETWORK = SOURCE / "island_host_network_20261030/island_host_network_islands.csv"
GIT_BLOB_PINS = {
    str(PORT): "fdf73b11ddcb1644f46730123350ae52cc61cae6",
    str(VEG): "10360b39e20a720000c5ac8f9f156073c67e7c9b",
    str(NETWORK): "18cecde5bf3066f276380167f451f1fa5cc1caa1",
}
DRAW_COUNT = 2000
SEED = 20261030
FIELDS = ("xx", "xz", "zz", "xy", "zy", "yy")


def source_rows(path: Path, sha: str) -> list[dict[str, str]]:
    body = path.read_bytes()
    git_blob = hashlib.sha1(b"blob " + str(len(body)).encode() + b"\0" + body).hexdigest()
    if git_blob != sha:
        raise ValueError(f"Source Git blob differs from locked identity: {path} {git_blob}")
    return list(csv.DictReader(body.decode("utf-8").splitlines()))


def prepare(root: Path = ROOT):
    port_rows = source_rows(root / PORT, GIT_BLOB_PINS[str(PORT)])
    vegetation = source_rows(root / VEG, GIT_BLOB_PINS[str(VEG)])
    networks = source_rows(root / NETWORK, GIT_BLOB_PINS[str(NETWORK)])
    if len(port_rows) != 3630:
        raise ValueError("Port source row count differs from locked upstream source")

    ports = {}
    for row in port_rows:
        key = row["port_id"]
        if key in ports:
            raise ValueError(f"Duplicate port: {key}")
        ports[key] = row
    seen = set()
    complete = []
    classified = []
    for v in vegetation:
        if v["near_radius_km"] != "1":
            continue
        key = v["port_id"]
        if key in seen:
            raise ValueError(f"Duplicate 1-km vegetation port {key}")
        seen.add(key)
        p = ports.get(key)
        if p is None or p["query_complete"] != "True":
            continue
        complete.append(v)
        alien = float(p["contemporary_not_native_in_catalogue"])
        native = float(p["native_in_catalogue"])
        ambiguous = float(p["ambiguous_across_regions"])
        unknown = float(p["unknown_catalogue_status"])
        if any(not math.isfinite(n) or n < 0 for n in (alien, native, ambiguous, unknown)):
            raise ValueError("Invalid or missing plant-catalogue number")
        denominator = alien + native
        if denominator == 0:
            continue  # Not an observed biological absence: ratio undefined
        x = float(v["built_difference"])
        near = float(v["built_near"])
        far = float(v["built_far"])
        if not all(math.isfinite(z) for z in (x, near, far)):
            raise ValueError("Missing landcover fractions")
        if not (-1 <= x <= 1 and 0 <= near <= 1 and 0 <= far <= 1 and
                math.isclose(near - far, x, abs_tol=1e-12)):
            raise ValueError("Coastal landcover decomposition mismatch")
        classified.append({
            "port_id": key,
            "island": v["landmass_id"],
            "country": p["country"],
            "x": x,
            "near": near,
            "far": far,
            "y": alien / denominator,
            "z": math.log1p(denominator),
            "alien": alien,
            "native": native,
            "ambiguous": ambiguous,
            "unknown": unknown,
            "denominator": denominator,
        })

    groups = defaultdict(list)
    for row in classified:
        groups[row["island"]].append(row)
    eligible = {
        island: rows for island, rows in groups.items()
        if len(rows) >= 2 and len({r["x"] for r in rows}) > 1
    }
    cohort_network = {row["landmass_id"] for row in networks if row["eligible"] == "True"}
    if (len(seen), len(complete), len({v["landmass_id"] for v in complete}),
        len(classified), len(groups), len(eligible),
        sum(map(len, eligible.values()))) != (625, 210, 97, 161, 77, 24, 104):
        raise ValueError("Original pre-score cohort support differs from locked specification")

    return {
        "eligible": eligible,
        "network": cohort_network,
        "all_1km_ports": len(seen),
        "complete_ports": len(complete),
        "complete_islands": len({v["landmass_id"] for v in complete}),
        "classified_ports": len(classified),
        "classified_islands": len(groups),
    }


def centered_terms(rows: list[dict], field: str = "x") -> dict[str, float]:
    n = len(rows)
    if not n:
        raise ValueError("Empty island")
    mx = sum(r[field] for r in rows) / n
    my = sum(r["y"] for r in rows) / n
    mz = sum(r["z"] for r in rows) / n
    out = dict.fromkeys(FIELDS, 0.0)
    for row in rows:
        x, y, z = row[field] - mx, row["y"] - my, row["z"] - mz
        out["xx"] += x*x/n
        out["xz"] += x*z/n
        out["zz"] += z*z/n
        out["xy"] += x*y/n
        out["zy"] += z*y/n
        out["yy"] += y*y/n
    return out


def fit_terms(terms: list[dict[str, float]], adjusted: bool = True):
    s = {k: sum(t[k] for t in terms) for k in FIELDS}
    if s["xx"] <= 0:
        return {"ok": False, "reason": "NO_WITHIN_ISLAND_BUILT_CONTRAST"}
    if adjusted:
        det = s["xx"]*s["zz"] - s["xz"]**2
        if det <= max(1e-12, 1e-12*s["xx"]*s["zz"]):
            return {"ok": False, "reason": "EFFORT_BUILT_COLLINEARITY"}
        coefficient = (s["xy"]*s["zz"] - s["zy"]*s["xz"]) / det
        effort_coef = (s["zy"]*s["xx"] - s["xy"]*s["xz"]) / det
    else:
        coefficient = s["xy"]/s["xx"]
        effort_coef = None
    return {
        "ok": True,
        "beta_fraction_per_unit_built_contrast": coefficient,
        "percentage_points_per_10pp_built_contrast": coefficient*10,
        "effort_proxy_coefficient": effort_coef,
        "cross_products": s,
    }


def xorshift32(seed: int):
    """Portable fixed-seed RNG (only to resample *islands*, not plant records)."""
    seed &= 0xffffffff
    while True:
        seed ^= (seed << 13) & 0xffffffff
        seed ^= seed >> 17
        seed ^= (seed << 5) & 0xffffffff
        seed &= 0xffffffff
        yield seed / 2**32


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("No successful estimates")
    vals = sorted(values)
    r = fraction*(len(vals)-1)
    lo = int(math.floor(r))
    hi = int(math.ceil(r))
    return vals[lo] if lo == hi else vals[lo]*(hi-r) + vals[hi]*(r-lo)


def analyze(root: Path = ROOT, draws: int = DRAW_COUNT) -> dict:
    data = prepare(root)
    eligible = dict(sorted(data["eligible"].items(), key=lambda t: t[0]))
    islands = list(eligible)
    terms = {k: centered_terms(rows) for k, rows in eligible.items()}
    primary = fit_terms(list(terms.values()))
    raw = fit_terms(list(terms.values()), adjusted=False)
    near = fit_terms([centered_terms(rows, "near") for rows in eligible.values()])
    if not primary["ok"] or not raw["ok"] or not near["ok"]:
        raise ValueError("The fixed primary or sensitivity is unidentifiable")
    leave_one_out = {
        island: fit_terms([v for name, v in terms.items() if name != island])
        for island in islands
    }
    if any(not v["ok"] for v in leave_one_out.values()):
        raise ValueError("Leave-one-island fit is not identified")
    network_islands = [k for k in islands if k in data["network"]]
    network = fit_terms([terms[k] for k in network_islands])
    if not network["ok"]:
        raise ValueError("Network overlap fit is not identified")

    rng = xorshift32(SEED)
    boot = []
    rejected = 0
    for _ in range(draws):
        sample = [terms[islands[int(next(rng)*len(islands))]] for _ in islands]
        fit = fit_terms(sample)
        if fit["ok"]:
            boot.append(fit["percentage_points_per_10pp_built_contrast"])
        else:
            rejected += 1
    if not boot:
        raise ValueError("No successful island cluster bootstrap resamples")
    all_ports = [row for rows in eligible.values() for row in rows]
    countries = dict(sorted(Counter(r["country"] for r in all_ports).items()))
    looc = [f["percentage_points_per_10pp_built_contrast"]
            for f in leave_one_out.values()]
    return {
        "schema": "chocho_island_port_host_gradient_exploratory_20261030_v1",
        "candidate_butterflies_unchanged": 225,
        "no_nonport_host_inventory_available": True,
        "not_causal": True,
        "different_spatial_radius_and_time_windows": True,
        "input_git_blob_pins": GIT_BLOB_PINS,
        "support": {
            "harbors_with_1km_vegetation": data["all_1km_ports"],
            "plant_query_complete_harbors": data["complete_ports"],
            "plant_query_complete_islands": data["complete_islands"],
            "positive_classified_host_harbors": data["classified_ports"],
            "positive_classified_host_islands": data["classified_islands"],
            "within_island_eligible_islands": len(islands),
            "within_island_eligible_harbors": len(all_ports),
            "ports_per_eligible_island_range": [min(map(len,eligible.values())),
                                                max(map(len,eligible.values()))],
            "host_catalogue_denominator_range": [min(r["denominator"] for r in all_ports),
                                                 max(r["denominator"] for r in all_ports)],
            "countries": countries,
            "unresolved_host_catalogue_ports": sum(r["unknown"] > 0 for r in all_ports),
            "eligible_islands_shared_with_51_island_network": len(network_islands),
            "eligible_harbors_shared_with_51_island_network":
                sum(len(eligible[k]) for k in network_islands),
        },
        "primary": primary,
        "without_host_richness_effort_proxy": raw,
        "absolute_built_near_instead_of_contrast": near,
        "leave_one_island_out": {
            "min_percentage_points_per_10pp": min(looc),
            "max_percentage_points_per_10pp": max(looc),
            "positive_estimates": sum(x > 0 for x in looc),
            "negative_estimates": sum(x < 0 for x in looc),
            "number": len(looc),
        },
        "island_cluster_bootstrap": {
            "seed": SEED,
            "draws": draws,
            "successful": len(boot),
            "failed_rank": rejected,
            "percentile_95_interval_percentage_points_per_10pp":
                [percentile(boot, 0.025), percentile(boot, 0.975)],
            "median_percentage_points_per_10pp": percentile(boot, 0.5),
        },
        "network_51_island_overlap": {
            "islands": len(network_islands),
            "harbors": sum(len(eligible[k]) for k in network_islands),
            "fitted": network,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--draws", type=int, default=DRAW_COUNT)
    args = parser.parse_args()
    result = analyze(args.root, draws=args.draws)
    contents = json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    if args.output:
        dest = args.output if args.output.is_absolute() else args.root / args.output
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(contents, encoding="utf-8")
        print(f"Wrote {dest}")
        print(json.dumps(result["primary"], sort_keys=True))
    else:
        print(contents)


if __name__ == "__main__":
    main()
