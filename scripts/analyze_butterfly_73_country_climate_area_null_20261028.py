#!/usr/bin/env python3
"""Area-opportunity geographic reference for frozen 73 butterfly positive cases.

Does not use new butterfly outcomes, pseudoabsences, or larval use. The control is
WorldClim climate distribution among valid 10' LAND PIXEL CENTERS inside the
pinned Natural Earth 1:50m receiving-country polygons, not sampling effort.
The focal-species native-HOST-PLANT climate rectangle is held fixed in both the
butterfly observation and the geographic null. This tests geographic sorting,
not interaction-mediated niche release or physiological niche expansion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rasterio.features import rasterize
from shapely.geometry import mapping

from analyze_butterfly_73_paired_point_climate_20261028 import (
    SOURCE_CSV_SHA256, read_cases,
)
from audit_butterfly_73_worldclim_point_niche_20261027 import load_climate_from_zip
from build_country_wgsrpd_crosswalk_20261020 import (
    NATURAL_EARTH_SHA256, dissolve_countries,
)

SCHEMA = "chocho_73_country_land_climate_area_null_20261028_v01"
BOUNDS_SCHEMA = "chocho_73_native_host_worldclim_reference_bounds_20261028_v01"
EXPECTED_ORIGINAL_JSON_SHA256 = "5d3ca4137821e60aaca7e1bf965c509b1cddee76c643fb00ee15df44ae775937"
EXPECTED_WORLDCLIM_ZIP_SHA256 = "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
MIN_VALID_COUNTRY_PIXELS = 10
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASES = ROOT / "data/exploratory/butterfly_73_worldclim_case_counts_20261027.csv"
DEFAULT_BOUNDS = ROOT / "data/exploratory/butterfly_73_native_host_worldclim_reference_bounds_20261028.json"
DEFAULT_CROSSWALK = ROOT / "data/exploratory/gbif_country_20261020/country_wgsrpd_summary.json"


def file_hash(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_bounds(path, cases):
    d = json.loads(Path(path).read_text(encoding="utf8"))
    origin = d.get("origin", {})
    if d.get("schema") != BOUNDS_SCHEMA or \
            origin.get("original_source_json_sha256") != EXPECTED_ORIGINAL_JSON_SHA256 or \
            origin.get("raster_zip_sha256") != EXPECTED_WORLDCLIM_ZIP_SHA256:
        raise ValueError("Native-host WorldClim reference lacks frozen source identity")
    ref = d["references"]
    if set(ref) != {r["species"] for r in cases} or len(ref) != 20:
        raise ValueError("All 20 original species' native-host references required")
    for species, sp in ref.items():
        for which in ("area_weighted_5_95", "area_weighted_10_90", "equal_pixel_5_95"):
            for climate in ("bio1", "bio12"):
                lo, hi = sp["bounds"][which][climate]
                if not all(math.isfinite(z) for z in (lo, hi)) or not lo < hi:
                    raise ValueError("Non-finite or inverted native host climate bounds: " + species)
    return d


def load_countries(path, crosswalk, country_codes):
    """Dissolve source-pinned, unprojected Natural Earth 50m map units by ISO_A2_EH."""
    if file_hash(path) != NATURAL_EARTH_SHA256:
        raise ValueError("Natural Earth 50m source did not match 2026-10-20 pinned crosswalk")
    source = json.loads(Path(path).read_text(encoding="utf8"))
    summary = json.loads(Path(crosswalk).read_text(encoding="utf8"))
    if summary.get("schema") != "country_wgsrpd_area_crosswalk_v1" or \
            summary["source"]["natural_earth_sha256"] != NATURAL_EARTH_SHA256:
        raise ValueError("Mismatched audited country crosswalk")
    country_index = {r["country"]: r for r in summary["countries"]}
    # These 73 source-positive countries were already restricted to audited
    # crosswalk-eligible map units; fail rather than silently reassign borders.
    if not set(country_codes) <= set(country_index):
        raise ValueError("A focal country is absent from the original audited crosswalk")
    all_geo, audit = dissolve_countries(source["features"], lambda g:g)
    if not set(country_codes) <= set(all_geo):
        raise ValueError("A focal ISO country has no dissolved Natural Earth polygon")
    return {c: all_geo[c]["geometry"] for c in sorted(country_codes)}, \
        {c: country_index[c] for c in sorted(country_codes)}, audit


def area_fraction_outside(country_polygons, climate, cases, reference, *, all_touched=False):
    """Country×species WorldClim outside share, with latitude area weights.

    Country polygons are externally sourced; maps do not contain actual butterfly
    point positions or author-only plant range geometries.
    """
    shape = climate["valid"].shape
    all_ids = sorted(country_polygons)
    raster = rasterize(
        [(mapping(country_polygons[c]), i+1) for i,c in enumerate(all_ids)],
        out_shape=shape, transform=climate["transform"], fill=0,
        all_touched=all_touched, dtype="uint16",
    )
    # Country code overlap on a pixel is assigned once by rasterize; count it
    # as an approximate raster geography, NOT a fractional boundary cell.
    row_lats = climate["transform"].f + (
        np.arange(shape[0], dtype=float) + .5
    ) * climate["transform"].e
    if np.any(row_lats <= -90) or np.any(row_lats >= 90):
        raise ValueError("Unexpected worldclim grid row center latitude")
    area_weights = np.cos(np.deg2rad(row_lats))
    if np.any(area_weights <= 0):
        raise ValueError("Nonpositive land pixel area weight")
    out = {}
    required = sorted({(r["species"],r["country"]) for r in cases})
    for country, cid in zip(all_ids, range(1,len(all_ids)+1)):
        valid = (raster == cid) & climate["valid"]
        total = int(np.count_nonzero(valid))
        species = [sp for sp,c in required if c==country]
        if total < MIN_VALID_COUNTRY_PIXELS:
            for sp in species:
                out[sp,country] = {"status":"INSUFFICIENT_10MIN_LAND_PIXEL_CENTERS",
                                   "country_climate_pixels":total}
            continue
        rr,cc=np.nonzero(valid)
        denom = float(np.sum(area_weights[rr]))
        t = climate["bio1"][rr,cc]
        p = climate["bio12"][rr,cc]
        for sp in species:
            result={"status":"GEOGRAPHIC_AREA_REFERENCE_VALID",
                    "country_climate_pixels":total,
                    "country_total_latitude_area_weight":denom}
            for endpoint,bounds_key in (("area95","area_weighted_5_95"),
                                        ("area90","area_weighted_10_90"),
                                        ("equal95","equal_pixel_5_95")):
                b = reference[sp]["bounds"][bounds_key]
                outside = ((t < b["bio1"][0]) | (t > b["bio1"][1]) |
                           (p < b["bio12"][0]) | (p > b["bio12"][1]))
                result["expected_outside_"+endpoint] = float(
                    np.sum(area_weights[rr][outside]) / denom
                )
            out[sp,country] = result
    return out


def paired_compare(cases, area, *, endpoint="area95", minimum_records=1,
                   strict=False, min_country_pixels=MIN_VALID_COUNTRY_PIXELS):
    """Species-equal difference-in-differences, descriptive only; no p values."""
    if endpoint not in ("area95","area90","equal95") or minimum_records<1:
        raise ValueError("Unknown option")
    ncol = "strict_n" if strict else "n"
    ocol = ("strict_outside_" if strict else "outside_") + endpoint
    groups=defaultdict(lambda:defaultdict(list))
    rejected=Counter()
    for r in cases:
        ref=area[r["species"],r["country"]]
        if ref["status"]!="GEOGRAPHIC_AREA_REFERENCE_VALID" or \
                ref["country_climate_pixels"]<min_country_pixels:
            rejected["unscorable_country_land_pixels"] += 1
            continue
        if r[ncol]<minimum_records:
            rejected["insufficient_butterfly_points"] += 1
            continue
        bg = ref["expected_outside_"+endpoint]
        observed = r[ocol]/r[ncol]
        groups[r["species"]][r["source_class"]].append({
            "country":r["country"],"n":r[ncol],
            "observed":observed,"expected":bg,"residual":observed-bg,
        })
    paired=[]
    for species,cl in sorted(groups.items()):
        if "novel" not in cl or "analogue" not in cl:
            continue
        avg=lambda cat,key:sum(r[key] for r in cl[cat])/len(cl[cat])
        obs = avg("novel","observed")-avg("analogue","observed")
        exp = avg("novel","expected")-avg("analogue","expected")
        adj = obs-exp
        paired.append({"species":species,"n_novel_countries":len(cl["novel"]),
            "n_analogue_countries":len(cl["analogue"]),
            "n_novel_points":sum(x["n"] for x in cl["novel"]),
            "n_analogue_points":sum(x["n"] for x in cl["analogue"]),
            "observed_class_difference":obs,
            "land_available_class_difference":exp,
            "geography_adjusted_difference":adj})
    def agg(field):
        return (sum(p[field] for p in paired)/len(paired) if paired else None)
    def counts(field):
        a=[p[field] for p in paired]
        return {"positive":sum(x>1e-10 for x in a),
                "negative":sum(x < -1e-10 for x in a),
                "zero":sum(abs(x)<=1e-10 for x in a)}
    return {
        "endpoint":endpoint,"minimum_point_records_per_country_case":minimum_records,
        "strict_reported_uncertainty_le1000m":strict,
        "minimum_country_land_pixel_centers":min_country_pixels,
        "cases_used_before_pairing":sum(len(x) for v in groups.values() for x in v.values()),
        "cases_not_used":dict(rejected),
        "paired_species":len(paired),
        "paired_species_countries_each_side_ge2":sum(
            x["n_novel_countries"]>=2 and x["n_analogue_countries"]>=2 for x in paired),
        "mean_observed_country_class_difference":agg("observed_class_difference"),
        "mean_geographic_country_class_difference":agg("land_available_class_difference"),
        "mean_observed_minus_country_area_reference":agg("geography_adjusted_difference"),
        "geography_adjusted_species_directions":counts("geography_adjusted_difference"),
        "species_detail":paired,
        "p_values":None,
        "inference_gate":"Country land-climate availability only. Not GBIF observation-effort matched, not random sites of true butterfly dispersal, not a host/non-host local control, not demonstrated caterpillar use or establishment. Multiple source-positive cases and within-species country contrasts are posthoc descriptive.",
    }


def build_report(cases, bounds, nat_path, country_summary, climate, *, min_pixels=MIN_VALID_COUNTRY_PIXELS):
    countries, summaries, audit = load_countries(
        nat_path,country_summary,{r["country"] for r in cases})
    ref=bounds["references"]
    area=area_fraction_outside(countries,climate,cases,ref)
    alias=area_fraction_outside(countries,climate,cases,ref,all_touched=True)
    sensitivity=[paired_compare(cases,area,minimum_records=n) for n in (3,5,10)]
    sensitivity += [paired_compare(cases,area,endpoint="area90"),
                    paired_compare(cases,area,endpoint="equal95"),
                    paired_compare(cases,area,strict=True),
                    paired_compare(cases,alias),
                    paired_compare(cases,area,min_country_pixels=30)]
    flags={c:x["mapping_flags"] for c,x in summaries.items() if x.get("mapping_flags")}
    records=[]
    for r in cases:
        q=area[r["species"],r["country"]]
        records.append({"species":r["species"],"country":r["country"],
                        "original_source_class":r["source_class"],
                        "observed_GBIF_indexed_records":r["n"],
                        "observed_outside_native_host_5_95":r["outside_area95"],
                        "country_10min_land_climate_pixel_centers":q["country_climate_pixels"],
                        "country_status":q["status"],
                        "land_available_outside_native_host_5_95":q.get("expected_outside_area95"),
                        "land_available_outside_native_host_10_90":q.get("expected_outside_area90"),
                        "source_iso_territory_warnings":flags.get(r["country"],[])})
    out={
        "schema":SCHEMA,
        "status":"GEOGRAPHIC_AREA_REFERENCE_COMPUTED_NOT_EFFORT_MATCHED",
        "original_observation_cases":len(cases),
        "original_GBIF_point_indices":sum(r["n"] for r in cases),
        "source":{
            "input_73_case_csv_sha256":SOURCE_CSV_SHA256,
            "source_73_reference_bounds_sha256":file_hash(DEFAULT_BOUNDS),
            "source_73_reference_original_full_json_sha256":EXPECTED_ORIGINAL_JSON_SHA256,
            "natural_earth_50m_exact_sha256":NATURAL_EARTH_SHA256,
            "natural_earth_country_crosswalk_sha256":file_hash(country_summary),
            "worldclim_1970_2000_zip_sha256":EXPECTED_WORLDCLIM_ZIP_SHA256,
            "pixel_definition":"Source WorldClim 10-arc-minute valid land-climate grid; country polygon pixel CENTER only, cos(latitude) area weights; all_touched sensitivity separate"},
        "country_geography":{"source_countries":len(countries),
                             "source_map_unit_exclusions":len(audit["excluded_map_units"]),
                             "source_country_mapping_warnings":flags,
                             "country_ids_at_10min_grid_nonzero":sum(any(x["country_climate_pixels"]>0 for k,x in area.items() if k[1]==c) for c in countries)},
        "main":paired_compare(cases,area,min_country_pixels=min_pixels),
        "sensitivity":sensitivity,
        "case_summary":records,
        "result_interpretation":"A geographic control against country-wide annual temperature/precipitation supply. Residual does not show host-enabled niche shifts, observer effort bias, caterpillar survival, plant use, temporal colonization or population change."
    }
    if out["original_observation_cases"]!=73 or out["original_GBIF_point_indices"]!=6076:
        raise ValueError("Original cohort changed")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases",type=Path,default=DEFAULT_CASES)
    p.add_argument("--bounds",type=Path,default=DEFAULT_BOUNDS)
    p.add_argument("--country-crosswalk",type=Path,default=DEFAULT_CROSSWALK)
    p.add_argument("--natural-earth",type=Path,required=True)
    p.add_argument("--worldclim-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    cases,_=read_cases(args.cases)
    if file_hash(args.worldclim_zip)!=EXPECTED_WORLDCLIM_ZIP_SHA256:
        raise ValueError("WorldClim 2.1 archive SHA mismatched original butterfly point source")
    bounds=load_bounds(args.bounds,cases)
    import tempfile
    with tempfile.TemporaryDirectory(prefix="chocho-area-null-") as temp:
        climate,_=load_climate_from_zip(args.worldclim_zip,Path(temp))
        out=build_report(cases,bounds,args.natural_earth,args.country_crosswalk,climate)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"cases":out["original_observation_cases"],
        "source_countries":out["country_geography"]["source_countries"],
        "scorable_paired_species":out["main"]["paired_species"],
        "observed_contrast":out["main"]["mean_observed_country_class_difference"],
        "geographic_contrast":out["main"]["mean_geographic_country_class_difference"],
        "residual_contrast":out["main"]["mean_observed_minus_country_area_reference"]},sort_keys=True))

if __name__ == "__main__":
    main()
