#!/usr/bin/env python3
"""Source-locked native/introduced WCVP island-flora cross-check of port association.

Unit: 126 WGSRPD3 island botanical regions (and 123 HOSTS-eligible subset).
Input: committed aggregate WCVP v13 whole-flora / World Port Index data.
This tests contemporary ASSOCIATIONS, not shipborne introductions or dispersal.
Model: OLS log1p(flora count) ~ log2(1+WPI ports) + log(area_km2)
       + categorical TDWG Level1.
For adjusted models add log1p(other flora count).
HC3 intervals are descriptive, not spatially independent; nonsignificance
in the native model is NOT proof of absence or formal equivalence.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv"
SOURCE_BLOB = "53b618fa51a1a418ccc6e4363f4c866c389c90cc"
OUT = ROOT / "results/port_native_alien_cross_20261118/port_native_alien_cross_source_verified.json"


def git_blob_sha(path):
    b = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()


def fit(rows, level1_reference, outcome, additional):
    x = []
    y = []
    for r in rows:
        row = [1.0, float(r["log2ports"]), float(r["log_area"])]
        row.extend(float(r["level1"] == lev) for lev in level1_reference[1:])
        row.extend(np.log1p(float(r[v])) for v in additional)
        x.append(row)
        if outcome == "alien_minus_native_log_ratio":
            yi = np.log1p(float(r["introduced"])) - np.log1p(float(r["native"]))
        else:
            yi = np.log1p(float(r[outcome]))
        y.append(yi)
    X = np.array(x, dtype=float)
    Y = np.array(y, dtype=float)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise RuntimeError("Rank-deficient botanical region regression")
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    resid = Y - X @ beta
    gram = np.linalg.inv(X.T @ X)
    leverage = np.sum((X @ gram) * X, axis=1)
    influence = X * (resid / (1.0 - leverage))[:, None]
    # MacKinnon-White HC3 covariance; use a normal 1.96 multiplier as a
    # DESCRIPTIVE interval, not a regional-cluster confidence interval.
    cov = gram @ (influence.T @ influence) @ gram
    b = float(beta[1])
    se = float(np.sqrt(cov[1, 1]))
    return {
        "n_island_botanical_regions": len(rows),
        "response": outcome,
        "additional_flora_covariates": list(additional),
        "port_beta_per_log2_1plus_WPI_port": b,
        "HC3_standard_error": se,
        "HC3_95pct_descriptive_interval": [b - 1.96 * se, b + 1.96 * se]
    }


def main():
    actual = git_blob_sha(SOURCE)
    if actual != SOURCE_BLOB:
        raise RuntimeError(f"Source aggregate CSV changed: {actual} != {SOURCE_BLOB}")
    with SOURCE.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 126 or len({r["code"] for r in rows}) != 126:
        raise RuntimeError("Original 126 island identities changed")
    for r in rows:
        if int(r["native"]) + int(r["introduced"]) != int(r["total"]):
            raise RuntimeError(f"Whole flora count invariant failed: {r['code']}")
    host = [r for r in rows if int(r["contemporary_butterfly_resource_units"]) > 0]
    if len(host) != 123:
        raise RuntimeError("Original 123 HOSTS-eligible source subset changed")
    levels = sorted({r["level1"] for r in rows}, key=int)
    specifications = [
        ("native", ()), ("introduced", ()),
        ("alien_minus_native_log_ratio", ()),
        ("native", ("introduced",)), ("introduced", ("native",))
    ]
    results = []
    for cohort, sample in (("WCVP126", rows), ("HOSTS123", host)):
        for response, controls in specifications:
            z = fit(sample, levels, response, controls)
            z["cohort"] = cohort
            results.append(z)
    def point(cohort, resp, controls):
        return next(z["port_beta_per_log2_1plus_WPI_port"] for z in results
                    if z["cohort"] == cohort and z["response"] == resp
                    and z["additional_flora_covariates"] == list(controls))
    if abs(point("WCVP126", "introduced", ("native",)) - 0.22140536684115) > 1e-8:
        raise RuntimeError("Original 126-island port–alien conditional result not reproduced")
    if abs(point("HOSTS123", "introduced", ("native",)) - 0.23347121981392502) > 1e-8:
        raise RuntimeError("Original 123-island port–alien conditional result not reproduced")
    report = {
        "schema": "source_locked_crossed_port_native_alien_WCVP126_20261118_v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_git_blob_sha": SOURCE_BLOB,
        "estimator": "OLS with TDWG Level1 fixed effects, log area, log2(1+ports); normal-approximation HC3",
        "source_reproduced": True,
        "n_islands": 126,
        "n_HOSTS_subpanel": 123,
        "results": results,
        "interpretation": (
            "Positive port association with introduced flora; no clearly distinct "
            "native response in this specification. Direct alien-minus-native "
            "slope interval includes zero, so selective effect not established."
        ),
        "limits": [
            "Null native association is not equivalence or proof ports are irrelevant.",
            "Native-flora adjustment versus alien-flora adjustment uses different conditional estimands.",
            "Plant richness is region-level contemporary WCVP inventory and ports are current WPI proxy.",
            "HC3 does not correct spatial dependence or nonrandom botanical ascertainment.",
            "Electrical ocean resistance, shortest-water routes and directional physical currents are distinct predictors.",
            "This model does not measure delivery of seeds, larvae, or butterfly establishment."
        ]
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS_SOURCE_LOCKED",
        "WCVP126_native_port": point("WCVP126", "native", ()),
        "WCVP126_alien_port": point("WCVP126", "introduced", ()),
        "WCVP126_alien_given_native": point("WCVP126", "introduced", ("native",)),
        "WCVP126_native_given_alien": point("WCVP126", "native", ("introduced",)),
        "WCVP126_alien_minus_native_slope": point("WCVP126", "alien_minus_native_log_ratio", ())
    }, sort_keys=True))


if __name__ == "__main__":
    main()
