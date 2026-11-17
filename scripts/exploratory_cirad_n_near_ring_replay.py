#!/usr/bin/env python3
"""Offline source-locked replay: modeled N 0–10 km vs 10–25 km at WPI ports.

CIRAD 2022 agricultural supply / human excretion are modeled N-input
POTENTIAL, not observed wastewater N, terrestrial soil N or insect growth.
The response is the GBIF indexed candidate plant list, not abundance.
All raw source CSVs must come from the three pinned GitHub Actions artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

SOURCE_SHA256 = {
    "jones": "9e3d35516927741f1f176a029b9d41bb03c70c850c066d15e5fe9730c2036482",
    "worldcover": "1cb891076e13aa9d93bced1b40e6471d447771af2b8b028f134ed052d13bf071",
    "cirad": "7d934dec9bf3d62ceb1b7583ca5d84a54c57b0440147913ff27d9895026c3c0c",
}
RUNS = {"jones": 41563940689, "worldcover": 41562997326, "cirad": 41923811745}
BASE = ["log_wwp", "log_wwtp", "abs_lat", "sinlon", "coslon",
        "built_fraction_2021", "crop_fraction_2021", "log_human25", "log_crop25"]


def checked_digest(path: Path, expected: str) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    actual = h.hexdigest()
    if actual != expected:
        raise ValueError(f"Source CSV SHA256 mismatch: {path.name}: {actual}")
    return actual


def annulus_from_nested(valid10, mean10, zero10, cells10,
                        valid25, mean25, zero25, cells25):
    """Construct a masked 10–25 km annulus from original nested circle stats.

    Valid-pixel weighted sums, NOT arithmetic mean of circle means; source
    sea/no-data pixels stay missing, while genuine mapped zeros stay zeros.
    """
    numbers = [valid10, mean10, zero10, cells10, valid25, mean25, zero25, cells25]
    if any(pd.isna(x) for x in numbers):
        return dict(mean=np.nan, valid_cells=np.nan, circle_cells=np.nan,
                    coverage=np.nan, nonzero_share=np.nan, nonzero_mean=np.nan,
                    nonzero_cells=np.nan)
    n0, n1, c0, c1 = int(valid10), int(valid25), int(cells10), int(cells25)
    if min(n0, n1, c0, c1) < 0 or n0 > n1 or c0 >= c1 or n0 > c0 or n1 > c1:
        raise ValueError("Impossible nested raster coverage / count geometry")
    if n0 == 0 or n1 == 0 or n1 == n0:
        return dict(mean=np.nan, valid_cells=n1 - n0, circle_cells=c1 - c0,
                    coverage=(n1 - n0) / (c1 - c0),
                    nonzero_share=np.nan, nonzero_mean=np.nan, nonzero_cells=np.nan)
    if min(mean10, mean25) < 0 or not (0 <= zero10 <= 1 and 0 <= zero25 <= 1):
        raise ValueError("Negative N input or invalid source zero-pixel fraction")
    nr = n1 - n0
    total = n1 * float(mean25) - n0 * float(mean10)
    if total < -1e-6 * max(1, n1 * float(mean25)):
        raise ValueError("Impossible negative N total in the derived annulus")
    total = max(0., total)
    # Avoid a false inference about mapped farmland from rounded zero shares.
    nz0 = n0 * (1 - float(zero10))
    nz1 = n1 * (1 - float(zero25))
    nzr = nz1 - nz0
    if nzr < -1e-5 or nzr > nr + 1e-5:
        raise ValueError("Inconsistent nested zero fractions")
    nzr = max(0., min(float(nr), nzr))
    return dict(mean=total / nr, valid_cells=nr, circle_cells=c1 - c0,
                coverage=nr / (c1 - c0), nonzero_share=nzr / nr,
                nonzero_mean=total / nzr if nzr > 1e-7 else np.nan,
                nonzero_cells=nzr)


def add_annulus(z: pd.DataFrame) -> pd.DataFrame:
    out = z.copy()
    if not out.port_id.is_unique:
        raise ValueError("Duplicate port IDs in 2022 CIRAD source")
    for nutrient in ("human", "crop"):
        cols = [f"{nutrient}_{r}km_{kind}" for r in (10, 25)
                for kind in ("valid_cells", "mean", "zero_share_valid", "circle_cells")]
        absent = set(cols) - set(out.columns)
        if absent:
            raise ValueError(f"Original CIRAD raster feature columns missing: {sorted(absent)}")
        for col in cols:
            out[col] = pd.to_numeric(out[col], errors="coerce")
        p10, p25 = f"{nutrient}_10km_", f"{nutrient}_25km_"
        values = [annulus_from_nested(
            row[p10+"valid_cells"], row[p10+"mean"], row[p10+"zero_share_valid"], row[p10+"circle_cells"],
            row[p25+"valid_cells"], row[p25+"mean"], row[p25+"zero_share_valid"], row[p25+"circle_cells"])
            for _, row in out.iterrows()]
        for field in ("mean", "valid_cells", "circle_cells", "coverage", "nonzero_share", "nonzero_mean", "nonzero_cells"):
            out[f"{nutrient}_ring_{field}"] = [v[field] for v in values]
        out[f"{nutrient}_near_coverage"] = out[p10+"valid_cells"] / out[p10+"circle_cells"].replace(0, np.nan)
        out[f"{nutrient}_near_nonzero_cells"] = out[p10+"valid_cells"] * (1 - out[p10+"zero_share_valid"])
        out[f"{nutrient}_near_nonzero_share"] = 1 - out[p10+"zero_share_valid"]
        out[f"{nutrient}_near_nonzero_mean"] = np.where(
            out[f"{nutrient}_near_nonzero_cells"] > 1e-7,
            out[p10+"mean"] * out[p10+"valid_cells"] /
            out[f"{nutrient}_near_nonzero_cells"].replace(0, np.nan), np.nan)
        out[f"{nutrient}_delta"] = np.log1p(out[p10+"mean"]) - np.log1p(out[f"{nutrient}_ring_mean"])
        out[f"{nutrient}_nonzero_delta"] = out[f"{nutrient}_near_nonzero_share"] - out[f"{nutrient}_ring_nonzero_share"]
        out[f"{nutrient}_intensity_delta"] = np.log1p(out[f"{nutrient}_near_nonzero_mean"]) - np.log1p(out[f"{nutrient}_ring_nonzero_mean"])
    return out


def prepare_frame(source: pd.DataFrame, rasters: pd.DataFrame, *, min_coverage=0.10,
                  min_valid=8, min_nonzero=0) -> pd.DataFrame:
    d = source.merge(add_annulus(rasters), on="port_id", validate="one_to_one")
    if len(d) != 1596 or d.port_id.nunique() != 1596:
        raise ValueError("WPI1596 source ID join mismatch")
    d = d.loc[d.in_within_country_comparison & ~d.incomplete_or_missing_grid &
              (d.quality == "COMPLETE")].copy()
    d["built_fraction_2021"] = pd.to_numeric(d.built_fraction_2021, errors="coerce")
    d["crop_fraction_2021"] = pd.to_numeric(d.crop_fraction_2021, errors="coerce")
    for k in ("human", "crop"):
        for radius in ("near", "ring"):
            d = d.loc[d[f"{k}_{radius}_coverage"] >= min_coverage]
            ncol = (f"{k}_10km_valid_cells" if radius == "near" else f"{k}_ring_valid_cells")
            d = d.loc[d[ncol] >= min_valid]
            if min_nonzero:
                d = d.loc[d[f"{k}_{radius}_nonzero_cells"] >= min_nonzero]
    d["log_wwp"] = np.log1p(d.produced_25km_millionm3yr)
    d["log_wwtp"] = np.log1p(d.WWTP_count_25km)
    d["abs_lat"] = d.latitude.abs()
    d["sinlon"] = np.sin(np.deg2rad(d.longitude))
    d["coslon"] = np.cos(np.deg2rad(d.longitude))
    d["log_indexed"] = np.log1p(d.all_taxonomically_indexed)
    d["log_human25"] = np.log1p(d.human_25km_mean)
    d["log_crop25"] = np.log1p(d.crop_25km_mean)
    a, n = d.alien_candidates.astype(float), d.native_candidates.astype(float)
    if (a < 0).any() or (n < 0).any():
        raise ValueError("Negative source botanical counts")
    d["indexed_ratio"] = np.log((a + .5) / (n + .5))
    d["observed_share"] = a.div((a + n).replace(0, np.nan))
    d["log_alien"] = np.log1p(a)
    d["log_native"] = np.log1p(n)
    required = BASE + ["human_delta", "crop_delta", "indexed_ratio"]
    if min_nonzero:
        required += ["human_intensity_delta", "crop_intensity_delta",
                     "human_nonzero_delta", "crop_nonzero_delta"]
    return d.dropna(subset=required).copy()


def countries_supported(d: pd.DataFrame, *, min_ports=5):
    g = d.groupby("country").agg(n=("port_id", "size"),
        wwtp=("log_wwtp", "nunique"), crop=("crop_delta", "nunique"),
        human=("human_delta", "nunique"))
    keep = g.index[(g.n >= min_ports) & (g.wwtp >= 2) &
                   (g.crop >= 2) & (g.human >= 2)]
    return d.loc[d.country.isin(keep)].copy()


def fit(d: pd.DataFrame, response: str, *, method="ols", gradient="mean", total_proxy=False):
    if len(d) < 100 or d.country.nunique() < 8:
        return {"status": "INSUFFICIENT_SUPPORT", "ports": len(d),
                "countries": int(d.country.nunique())}
    extras = {
        "mean": ["human_delta", "crop_delta"],
        "nonzero": ["human_nonzero_delta", "crop_nonzero_delta"],
        "intensity_plus_nonzero": ["human_intensity_delta", "crop_intensity_delta",
                                   "human_nonzero_delta", "crop_nonzero_delta"],
    }[gradient]
    terms = BASE + (["log_indexed"] if total_proxy else []) + extras
    form = response + " ~ " + " + ".join(terms + ["C(country)"])
    if method == "frac_logit":
        if d[response].isna().any() or not d[response].between(0, 1).all():
            raise ValueError("Fractional logit requires defined 0–1 shares")
        model = smf.glm(form, data=d, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d.country})
    else:
        model = smf.ols(form, data=d).fit(
            cov_type="cluster", cov_kwds={"groups": d.country})
    result = {"status": "OBSERVATIONAL_NOT_CAUSAL", "ports": int(model.nobs),
              "countries": int(d.country.nunique()), "formula": form,
              "method": method, "response": response, "R2": float(model.rsquared) if method == "ols" else None}
    for term in extras:
        result[term] = {"beta": float(model.params[term]),
                        "CI95_country_cluster": list(map(float, model.conf_int().loc[term])),
                        "p_nominal_posthoc": float(model.pvalues[term])}
    return result


def study(jones: Path, worldcover: Path, cirad: Path):
    paths = {"jones": jones, "worldcover": worldcover, "cirad": cirad}
    verified = {name: checked_digest(path, SOURCE_SHA256[name]) for name, path in paths.items()}
    from exploratory_global_wwtp_jones_worldcover import read_source
    source = read_source(jones, worldcover)
    features = pd.read_csv(cirad, keep_default_na=False)
    if len(features) != 1596 or features.port_id.nunique() != 1596:
        raise ValueError("CIRAD port table does not contain the exact 1596 WPI identities")
    frames = {name: prepare_frame(source, features, min_coverage=coverage, min_nonzero=nonzero)
              for name, coverage, nonzero in
              (("main", .10, 0), ("coverage50", .50, 0), ("nonzero5", .10, 5))}
    summary = {"schema": "chocho_WPI1596_CIRAD2022_near_ring_N_gradients_offline_v1",
               "source_runs": RUNS, "verified_source_sha256": verified,
               "source_port_rows": len(source), "cohorts": {},
               "interpretation": "MODELED_N_INPUT_POTENTIAL_AND_INDEXED_HOST_COMPOSITION; NOT_SEWAGE_SOIL_N_OR_BUTTERFLY_ESTABLISHMENT"}
    for name, frame in frames.items():
        d = countries_supported(frame)
        positive = countries_supported(d.loc[d.observed_share.notna()])
        tag = {"main": "mean", "coverage50": "mean", "nonzero5": "intensity_plus_nonzero"}[name]
        models = {
          "logratio_N_intensity": fit(d, "indexed_ratio", gradient=tag, total_proxy=True),
          "alien_fraction_N_intensity": fit(positive, "observed_share", method="frac_logit", gradient=tag),
        }
        if name == "main":
            models["log_alien_negative_control"] = fit(d, "log_alien", gradient="mean")
            models["log_native_negative_control"] = fit(d, "log_native", gradient="mean")
            models["alien_fraction_N_nonzero_area"] = fit(positive, "observed_share", method="frac_logit", gradient="nonzero")
        summary["cohorts"][name] = {"n_eligible": len(d), "n_countries": int(d.country.nunique()),
                                    "n_fraction": len(positive), "n_fraction_countries": int(positive.country.nunique()),
                                    "models": models}
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--jones", required=True, type=Path)
    p.add_argument("--worldcover", required=True, type=Path)
    p.add_argument("--cirad", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    result = study(args.jones, args.worldcover, args.cirad)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("CIRAD_NEAR_RING_N_REPLAY_FINISHED", json.dumps({name: {k: v for k, v in cohort.items() if k != "models"}
          for name, cohort in result["cohorts"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
