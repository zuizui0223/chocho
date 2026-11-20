#!/usr/bin/env python3
"""Independent 239-butterfly isolation x port source audit (NOT original dd.py).
Inputs: archived accepted WCVP 126-island all-flora/WPI and frozen Q1 239 panel,
and TDWG level3 geometry (pinned at workflow download; hashes recorded).
No consumer presence, ship arrival, or causal host-rescue observations.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from pyproj import Geod
from shapely import affinity
from shapely.geometry import shape
from shapely.ops import nearest_points, unary_union

FLORA_BLOB = "53b618fa51a1a418ccc6e4363f4c866c389c90cc"
HOST_REGION_BLOB = "4aba89565354edc6c44d566523864cec66a30cbe"
URL = "https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson"
RESPONSES = {
    "plant_native": "native_plant", "plant_introduced": "introduced_plant",
    "host_native_239": "native_host_239", "host_current_239": "current_host_239"
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def blob(path):
    b = Path(path).read_bytes()
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()

def require(value, desc):
    if not value:
        raise RuntimeError(desc)

def load_geojson(path, codes):
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    acc = {}
    for ft in doc["features"]:
        code = str(ft["properties"].get("LEVEL3_COD", "")).strip()
        if code in codes and ft.get("geometry"):
            geom = shape(ft["geometry"])
            if not geom.is_valid:
                geom = geom.buffer(0)
            if not geom.is_empty:
                acc.setdefault(code, []).append(geom)
    return {code: unary_union(poly) for code, poly in acc.items()}

def geometry_features(geom_file, categories):
    geos = load_geojson(geom_file, set(categories))
    require(set(geos) == set(categories), "Missing TDWG L3 geometry: " + str(sorted(set(categories) - set(geos))))
    mainland = [geos[c] for c in sorted(categories) if categories[c] == "mainland"]
    require(len(mainland) == 237, "Expected 237 fixed mainland botanical regions")
    land = unary_union(mainland)
    # Shifts handle the Pacific date line; shortest lat/lon boundary points are
    # geodesic-converted below. This is NOT an exact ellipsoidal polygon-to-polygon
    # optimizer at high latitudes; keep it a separately disclosed proxy.
    lands = [affinity.translate(land, xoff=dx) for dx in (-360, 0, 360)]
    geod = Geod(ellps="WGS84")
    stats = {}
    for code in sorted(c for c in categories if categories[c] == "island"):
        g = geos[code]
        best = None
        for target in lands:
            a, b = nearest_points(g, target)
            lon_a = (a.x + 180) % 360 - 180
            lon_b = (b.x + 180) % 360 - 180
            _, _, metres = geod.inv(lon_a, a.y, lon_b, b.y)
            # The closest point under geographic coordinates is an approximation,
            # so store both raw and geodesic endpoint separation.
            cand = {
                "km": abs(metres) / 1000,
                "proxy_degree_distance": g.distance(target),
                "island_lon": lon_a, "island_lat": a.y,
                "mainland_lon": lon_b, "mainland_lat": b.y
            }
            if best is None or cand["km"] < best["km"]:
                best = cand
        require(best and np.isfinite(best["km"]), "No valid mainland proximity: " + code)
        stats[code] = dict(
            continent_distance_km=best["km"],
            abs_centroid_latitude=abs(g.centroid.y),
            nearest_mainland_lon=best["mainland_lon"],
            nearest_mainland_lat=best["mainland_lat"],
            proxy_degree_distance=best["proxy_degree_distance"]
        )
    return stats

def make_frame(port_path, q1_path, geo_path):
    require(blob(port_path) == FLORA_BLOB, "Whole-flora/WPI immutable source changed")
    require(blob(q1_path) == HOST_REGION_BLOB, "239 HOSTS region source changed")
    flora = pd.read_csv(port_path, dtype={"code": str, "level1": str})
    q1 = pd.read_csv(q1_path, dtype={"area_code_l3": str, "level1": str})
    require(len(flora) == 126 and flora.code.is_unique, "Expected exact original 126 botanical island source")
    require(len(q1) == 369 and q1.area_code_l3.is_unique, "Expected 369 Q1 region sources")
    require(int(flora.native.sum()) == 184997 and int(flora.introduced.sum()) == 39024,
            "Original WCVP floristic units changed")
    require(int(q1.native.sum()) == 26530 and int(q1.contemporary.sum()) == 41083,
            "Original fixed 239 native/current regional counts changed")
    categories = dict(zip(q1.area_code_l3, q1["class"]))
    require(sum(v == "island" for v in categories.values()) == 126 and
            sum(v == "mainland" for v in categories.values()) == 237,
            "Original island/mainland classification changed")
    geo = geometry_features(geo_path, categories)
    f = flora.merge(q1[["area_code_l3", "native", "contemporary", "class"]],
                    left_on="code", right_on="area_code_l3",
                    how="left", suffixes=("_plant", "_host"), validate="one_to_one")
    require(len(f) == 126 and f["class"].eq("island").all(), "Join lost fixed island cohort")
    require(f.contemporary.eq(f.contemporary_butterfly_resource_units).all(),
            "Frozen Q1 239-host count differs from plant source")
    f = f.rename(columns={"native_plant": "native_plant",
                          "introduced": "introduced_plant",
                          "native_host": "native_host_239",
                          "contemporary": "current_host_239"})
    for key in ("continent_distance_km", "abs_centroid_latitude",
                "nearest_mainland_lon", "nearest_mainland_lat", "proxy_degree_distance"):
        f[key] = f.code.map({c: v[key] for c, v in geo.items()})
    f["introduced_only_host_239"] = f.current_host_239 - f.native_host_239
    require((f.introduced_only_host_239 >= 0).all(), "Non-nested potential host resources")
    f["introduced_only_share_239"] = np.where(
        f.current_host_239 > 0, f.introduced_only_host_239 / f.current_host_239, np.nan)
    f["log_area_calc"] = np.log(f.area_km2)
    f["log_port_calc"] = np.log2(1 + f.wpi_port_count)
    f["log_distance_calc"] = np.log1p(f.continent_distance_km)
    require(np.isfinite(f[["log_area_calc","log_port_calc","log_distance_calc",
                           "abs_centroid_latitude"]]).all().all(), "Bad main source covariate")
    f["level1"] = f.level1.astype(str)
    return f

def design(data, binary_ports=False, without_lat=False):
    x = pd.DataFrame(index=data.index)
    x["intercept"] = 1.
    v = (data.wpi_port_count > 0).astype(float) if binary_ports else data.log_port_calc
    for name, val in [
        ("ports", v), ("distance", data.log_distance_calc),
        ("area", data.log_area_calc),
        ("abs_latitude", data.abs_centroid_latitude)
    ]:
        if name == "abs_latitude" and without_lat:
            continue
        mu = float(val.mean()); sd = float(val.std(ddof=1))
        require(sd > 0, "Degenerate regressor: " + name)
        x[name] = (val - mu) / sd
    levels = sorted(data.level1.unique(), key=lambda z: int(z))
    for lv in levels[1:]:
        x["level1_" + lv] = (data.level1 == lv).astype(float)
    return x

def fit_one(data, y, *, binary_ports=False, without_lat=False):
    x = design(data, binary_ports=binary_ports, without_lat=without_lat)
    model = sm.OLS(y, x).fit(cov_type="HC3")
    require(np.linalg.matrix_rank(x.to_numpy()) == x.shape[1], "Rank deficient full design")
    return model, x

def extract(m):
    out = {}
    for name in ("ports", "distance"):
        lo, hi = m.conf_int().loc[name]
        out[name] = {
            "beta": float(m.params[name]),
            "HC3_SE": float(m.bse[name]),
            "HC3_p": float(m.pvalues[name]),
            "HC3_95": [float(lo), float(hi)]
        }
    out["n"] = int(m.nobs)
    out["r2"] = float(m.rsquared)
    return out

def y_standardized(data, field):
    values = np.log1p(data[field].astype(float))
    return (values - values.mean()) / values.std(ddof=1)

def eval_suite(f, *, binary_ports=False, without_lat=False, min_area=None):
    if min_area is not None:
        f = f.loc[f.area_km2 < min_area].copy()
    result = {}
    for key, col in RESPONSES.items():
        m, x = fit_one(f, y_standardized(f, col),
                       binary_ports=binary_ports, without_lat=without_lat)
        result[key] = extract(m)
    for key, a, b in [
        ("plant_introduced_minus_native_z", "introduced_plant", "native_plant"),
        ("host_current_minus_native_z", "current_host_239", "native_host_239")
    ]:
        m, x = fit_one(f, y_standardized(f, a) - y_standardized(f, b),
                       binary_ports=binary_ports, without_lat=without_lat)
        result[key] = extract(m)
    active = f.loc[f.current_host_239 > 0].copy()
    # Log-odds of added-vs-native coverage with a Jeffreys-style continuity
    # correction, fitted at one island region per row (no 239-fold pseudo-N).
    odds = np.log((active.introduced_only_host_239 + 0.5) /
                  (active.native_host_239 + 0.5))
    m, x = fit_one(active, odds, binary_ports=binary_ports, without_lat=without_lat)
    result["introduced_only_share_logodds_239"] = extract(m)
    result["counts"] = {
        "islands": len(f), "resource_active_islands": len(active),
        "hosts_239_native_units": int(f.native_host_239.sum()),
        "hosts_239_current_units": int(f.current_host_239.sum()),
        "corr_port_count_vs_island_area": float(f.wpi_port_count.corr(f.area_km2)),
        "corr_logport_vs_logarea": float(f.log_port_calc.corr(f.log_area_calc))
    }
    return result

def bootstrap_slope(f, draws=999, seed=20261120):
    x = design(f).to_numpy()
    ys = {
        k: y_standardized(f, col).to_numpy() for k, col in RESPONSES.items()
    }
    ys["plant_introduced_minus_native_z"] = ys["plant_introduced"] - ys["plant_native"]
    ys["host_current_minus_native_z"] = ys["host_current_239"] - ys["host_native_239"]
    idx_by_level = [np.flatnonzero(f.level1.to_numpy() == lv)
                    for lv in sorted(f.level1.unique(), key=int)]
    rng = np.random.default_rng(seed)
    fits = {k: {"ports": [], "distance": []} for k in ys}
    valid = 0
    for _ in range(draws):
        idx = np.concatenate([rng.choice(z, size=len(z), replace=True) for z in idx_by_level])
        xx = x[idx]
        if np.linalg.matrix_rank(xx) < xx.shape[1]:
            continue
        valid += 1
        for key, y in ys.items():
            coef = np.linalg.lstsq(xx, y[idx], rcond=None)[0]
            fits[key]["ports"].append(float(coef[1]))
            fits[key]["distance"].append(float(coef[2]))
    require(valid >= draws * .70, "Too few full rank bootstrap draws")
    return {
        "requested_draws": draws, "valid_draws": valid,
        "scope": "Within-biogeographic-Level1 island bootstrap; no spatial independence proof",
        "effects": {
            key: {
                term: {
                    "q025_q975": [float(v) for v in np.quantile(vals, [.025, .975])],
                    "positive_fraction": float(np.mean(np.asarray(vals) > 0))
                }
                for term, vals in axes.items()
            }
            for key, axes in fits.items()
        }
    }

def conservative_geographic_cluster(f):
    """Descriptive nine-Level1 cluster standard errors, using t with G-1 df.

    Nine broad regions are a limited number of independent units: this is
    an explicit sensitivity, NOT a fully identified spatial error model.
    """
    from scipy.stats import t as t_dist
    categories = [
        ("plant_native","native_plant"),
        ("plant_introduced","introduced_plant"),
        ("host_native_239","native_host_239"),
        ("host_current_239","current_host_239"),
    ]
    ys = {key: y_standardized(f, field) for key, field in categories}
    ys["plant_introduced_minus_native_z"] = ys["plant_introduced"] - ys["plant_native"]
    ys["host_current_minus_native_z"] = ys["host_current_239"] - ys["host_native_239"]
    out = {}
    for cohort, data, subys in [
        ("all_126",f,ys),
        ("resource_active_123",f.loc[f.current_host_239>0].copy(),None)
    ]:
        if subys is None:
            subys = {"introduced_only_share_logodds_239":np.log(
                (data.introduced_only_host_239+0.5)/(data.native_host_239+0.5))}
        xx = design(data)
        g = int(data.level1.nunique())
        df = g-1
        require(g==9, "Unexpected Level1 cluster cohort")
        out[cohort] = {}
        for key, val in subys.items():
            m = sm.OLS(val,xx).fit(cov_type="cluster",
                cov_kwds={"groups": data.level1, "use_correction":True})
            effects={}
            for term in ("ports","distance"):
                b,se = float(m.params[term]),float(m.bse[term])
                critical = t_dist.ppf(.975,df)
                effects[term]={
                    "beta":b,"cluster9_SE":se,
                    "t8_two_sided_p":float(2*t_dist.sf(abs(b/se),df)),
                    "t8_95":[float(b-critical*se),float(b+critical*se)]
                }
            out[cohort][key]=effects
    return {
        "method":"TDWG Level1 9-region CR1 grouped sandwich, t8 sensitivity; only nine clusters, not geographic proof",
        "models":out
    }

def marginal_predictions(f):
    # In-range counterfactual predictions at common port, area, latitude and
    # Level1 compositions: only plotted x varies; not historical counterfactuals.
    x = design(f)
    native, _ = fit_one(f, np.log1p(f.native_host_239), without_lat=False)
    curr, _ = fit_one(f, np.log1p(f.current_host_239), without_lat=False)
    mu, sd = f.log_distance_calc.mean(), f.log_distance_calc.std(ddof=1)
    positions = np.geomspace(max(.01, float(f.continent_distance_km.min())),
                             max(100., float(f.continent_distance_km.max())), 70)
    positions = np.unique(np.concatenate([positions, [100, 1000, 4000]]))
    rows = []
    for dist in positions:
        temp = x.copy()
        temp["distance"] = (math.log1p(dist) - mu) / sd
        vals = [max(0., math.expm1(float(model.predict(temp).mean()))) for model in (native, curr)]
        rows.append({"distance_km": float(dist),
                     "pred_native_host_239": vals[0], "pred_current_host_239": vals[1]})
    return rows

def leave_one_island_out(f):
    """Check that the 126 vs unavailable 125 island mismatch is not one leverage point."""
    v = []
    for i in f.index:
        d = f.drop(index=i)
        xx = design(d)
        change = y_standardized(d,"current_host_239")-y_standardized(d,"native_host_239")
        m = sm.OLS(change,xx).fit(cov_type="HC3")
        v.append({
            "omitted_island":str(f.loc[i,"code"]),
            "slope_difference_distance_beta":float(m.params["distance"]),
            "HC3_p":float(m.pvalues["distance"])
        })
    return {
        "n_single_island_exclusions":len(v),
        "beta_min":min(z["slope_difference_distance_beta"] for z in v),
        "beta_max":max(z["slope_difference_distance_beta"] for z in v),
        "n_positive":sum(z["slope_difference_distance_beta"]>0 for z in v),
        "n_p_under_0p05":sum(z["HC3_p"]<.05 for z in v),
        "per_island":v,
        "warning":"This does not identify the specific island excluded in the inaccessible 125-island scratchpad."
    }

def adjusted_share_predictions(f):
    """Region-level, equal-island marginal fitted shares at requested distances.

    These are NOT the ratio of separately back-transformed native/current
    fitted host counts. Other recorded covariates are kept at each island's
    measured values, varying the isolation distance in the prediction only.
    """
    from scipy.special import expit
    active = f.loc[f.current_host_239>0].copy()
    x = design(active)
    logdist = active.log_distance_calc
    odds = np.log((active.introduced_only_host_239+.5)/
                  (active.native_host_239+.5))
    model = sm.OLS(odds,x).fit(cov_type="HC3")
    ret = {}
    for km in (100,1000,4000):
        temp=x.copy()
        temp["distance"]=(math.log1p(km)-float(logdist.mean()))/float(logdist.std(ddof=1))
        y=np.asarray(model.predict(temp),dtype=float)
        ret[str(km)]={
            "equal_island_mean_share":float(expit(y).mean()),
            "equal_island_median_share":float(np.median(expit(y))),
            "units":"fraction of covered butterfly-region potential resource combinations newly supported by introduced hosts"
        }
    return ret

def make_figure(f, pred, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.1))
    ax.scatter(f.continent_distance_km+1, f.native_host_239,
               s=16, alpha=.18, marker="o", color="#1b6f83")
    ax.scatter(f.continent_distance_km+1, f.current_host_239,
               s=16, alpha=.18, marker="s", color="#c07c2f")
    x = [p["distance_km"] + 1 for p in pred]
    ax.plot(x,[p["pred_native_host_239"] for p in pred],
            lw=2.4,color="#1b6f83",label="Native-host potential (239 butterflies)")
    ax.plot(x,[p["pred_current_host_239"] for p in pred],
            lw=2.4,color="#c07c2f",label="Native + introduced host potential")
    ax.set_xscale("log")
    ax.set_xlabel("Approximate nearest mainland boundary distance (km; +1 axis)")
    ax.set_ylabel("Butterflies with a mapped known host (of fixed 239; not present species)")
    ax.set_title("Source-pinned 239-butterfly opportunity vs island isolation")
    ax.grid(alpha=.16)
    ax.legend(frameon=False,loc="best")
    fig.tight_layout()
    for ext in (".svg",".png"):
        fig.savefig(out/("isolation_239_host_opportunity"+ext),dpi=220)
    plt.close(fig)

def main():
    pa = argparse.ArgumentParser()
    pa.add_argument("--port-flora",type=Path,required=True)
    pa.add_argument("--region-hosts",type=Path,required=True)
    pa.add_argument("--tdwg-geojson",type=Path,required=True)
    pa.add_argument("--out-dir",type=Path,required=True)
    a = pa.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    f = make_frame(a.port_flora,a.region_hosts,a.tdwg_geojson)
    report = {
        "schema":"chocho_fixed239_isolation_port_independent_test_20261120_v1",
        "status":"SOURCE_RECALCULATED_NOT_CLAUDE_DD_EXACT_REPLICATION",
        "population":"Original frozen 239 HOSTS butterflies, 126 WGSRPD3 islands (not original 2517/125)",
        "inputs":{"flora_git_blob":blob(a.port_flora),
                  "hosts_region_git_blob":blob(a.region_hosts),
                  "tdwg_level3_geojson_sha256":sha(a.tdwg_geojson),"tdwg_url":URL},
        "isolation_method":"Minimum WGS84 geodesic distance among nearest points of WGS84 lon/lat mainland and island WGSRPD3 boundary geometries, with ±360 longitude wrap. Pair chosen in degree-space: approximate geodesic nearest, not exact ellipsoid optimization.",
        "model":"OLS log1p(response) z-transformed per source response, standardized log2(1+ports), log1p(isolation km), log area, abs centroid latitude, TDWG Level1 fixed effects; HC3. Paired native/alien and native/current differences of standardized responses.",
        "models":eval_suite(f),
        "sensitivity": {
            "area_under_50000_km2":eval_suite(f,min_area=50000),
            "port_binary":eval_suite(f,binary_ports=True),
            "without_latitude":eval_suite(f,without_lat=True)
        },
        "bootstrap":bootstrap_slope(f),
        "marginal_adjusted_share_at_100_1000_4000km":adjusted_share_predictions(f),
        "leave_one_island_out":leave_one_island_out(f),
        "conservative_level1_cluster":conservative_geographic_cluster(f),
        "interpretation_guardrails":[
            "Geographic vegetation categories and modern port counts are not actual shipborne introduction, measured currents or host quality.",
            "Not 2517 original butterflies or original 125 islands; raw dd.py/dd_data.csv never inspected.",
            "Counts are potential available-host butterfly identities, NOT observed butterfly presence or a colonization time series.",
            "Standardized responses change the interpretation of slope differences; compare matched outcome contrasts explicitly.",
            "Non-significant contemporary slope does not prove zero slope.",
            "The nearest-boundary geographic surrogate can differ from the unavailable original mainland-isolation algorithm.",
            "Partial correlation of ports and region area and spatial autocorrelation remain risks.",
        ]
    }
    holdouts = []
    for lv in sorted(f.level1.unique(),key=int):
        pop = f.loc[f.level1 != lv].copy()
        try:
            mm = eval_suite(pop)
            holdouts.append({"excluded_level1":lv,"n":len(pop),
                             "native_host_distance_beta":mm["host_native_239"]["distance"]["beta"],
                             "current_host_distance_beta":mm["host_current_239"]["distance"]["beta"],
                             "host_slope_change_distance_beta":mm["host_current_minus_native_z"]["distance"]["beta"],
                             "introduced_plant_port_beta":mm["plant_introduced"]["ports"]["beta"]})
        except Exception as ex:
            holdouts.append({"excluded_level1":lv,"error":str(ex)})
    report["sensitivity"]["leave_one_biogeographic_group_out"]=holdouts
    pred = marginal_predictions(f)
    f.to_csv(a.out_dir/"fixed239_source_islands_with_isolation.csv",index=False)
    with (a.out_dir/"fixed239_marginal_predictions.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(pred[0]))
        writer.writeheader();writer.writerows(pred)
    make_figure(f,pred,a.out_dir)
    (a.out_dir/"SOURCE239_ISOLATION_PORT_REANALYSIS.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ISOLATION_PORT_239_SOURCE_RECALCULATED")
    print(json.dumps({"status":report["status"],"n":len(f),
                      "models":report["models"],
                      "bootstrap":report["bootstrap"]},ensure_ascii=False))
if __name__ == "__main__":
    main()
