#!/usr/bin/env python3
"""Audit recovered PRIVATE chocho original species×landmass×year Parquet.

Usage:
  python scripts/audit_chocho_historical_parquet_20261120.py --private-zip SOURCE.zip --out REPORT.json

Requires pandas and pyarrow when reading the private original Parquet. Do NOT
publish SOURCE.zip or attribute these aggregates to the unrelated 224 GBIF DOI.
Use --self-test without those optional dependencies to verify ID-zero semantics.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from hashlib import sha256
import csv
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

PINNED_PARQUET_SHA256 = "8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
PARQUET_NAME = "work/barriers/global_original_species_island_year.parquet"
PAIRS_NAME = "outputs/island_host_network_pairs.csv"
ISLANDS_NAME = "outputs/island_host_network_islands.csv"
FIVE = {
    "Colias croceus": (5137612, ("169", "276")),
    "Colias hyale": (5137706, ("169", "276")),
    "Vanessa cardui": (4299368, ("183", "821")),
    "Pieris rapae": (1920496, ("260", "951")),
    "Nymphalis antiopa": (5130682, ("276", "715")),
}
REQUIRED = ("speciesKey", "landmass_id", "year", "records",
            "precise1km_records", "landmass_area_km2")


def check_identifier_area(ids: list[str], areas: list[float]):
    """Numeric landmass ID=0 is allowed; negatives, null and area<=0 are not."""
    if len(ids) != len(areas):
        raise ValueError("Misaligned landmass id and area arrays")
    for index, (land, area) in enumerate(zip(ids, areas)):
        if land is None or str(land).strip() == "" or str(land).lower() == "nan":
            raise ValueError(f"Null source landmass at row {index}")
        if not str(land).isdigit():
            raise ValueError(f"Non-numeric source landmass at row {index}")
        if not math.isfinite(float(area)) or float(area) <= 0:
            raise ValueError(f"Non-positive source landmass area at row {index}")
    return {str(x) for x in ids}


def replay(frame, islands, pairs):
    """Source-audited full-panel checks, then STRICT explicit island-roster select."""
    import pandas as pd
    missing = set(REQUIRED) - set(frame.columns)
    if missing:
        raise ValueError(f"Original Parquet lacks required fields: {sorted(missing)}")
    if frame[list(REQUIRED)].isna().any().any():
        raise ValueError("Null value in original Parquet required columns")
    frame = frame.loc[:, REQUIRED].copy()
    frame["landmass_id"] = frame["landmass_id"].astype(str)
    frame["speciesKey"] = frame["speciesKey"].astype("int64")
    frame["year"] = frame["year"].astype("int64")
    check_identifier_area(frame["landmass_id"].tolist(), frame["landmass_area_km2"].tolist())
    profile = (
        len(frame), int(frame.speciesKey.nunique()), int(frame.landmass_id.nunique()),
        int(frame.year.min()), int(frame.year.max()), int(frame.records.sum()))
    expected = (83573, 224, 3384, 2001, 2024, 14425536)
    if profile != expected:
        raise ValueError(f"Historical cohort drift: {profile} vs {expected}")
    if frame.duplicated(["speciesKey", "landmass_id", "year"]).any():
        raise ValueError("Duplicate original species×landmass×year key")
    if (frame.records <= 0).any() or not frame.records.eq(frame.precise1km_records).all():
        raise ValueError("Original positive-record and precision source pin drift")
    if frame.groupby("landmass_id")["landmass_area_km2"].nunique().gt(1).any():
        raise ValueError("Inconsistent area within original landmass key")
    zero = frame[frame.landmass_id.eq("0")]
    if len(zero) != 2636 or int(zero.records.sum()) != 4386686:
        raise ValueError("Landmass ID zero retained original records unexpectedly changed")
    if zero.speciesKey.nunique() != 126:
        raise ValueError("Original landmass 0 species roster drift")
    if abs(float(zero.landmass_area_km2.iloc[0]) - 50654050.6945) > .0001:
        raise ValueError("Expected largest mainland source key does not match source")
    if frame.landmass_area_km2.max() != zero.landmass_area_km2.iloc[0]:
        raise ValueError("Original ID=0 must remain largest-area region")
    cohort = {str(r["landmass_id"]) for r in islands
              if str(r["eligible"]).strip().lower() == "true"}
    if len(cohort) != 51 or "0" in cohort:
        raise ValueError("51-island source roster drift or illegal ID0 island selection")
    selected = frame.loc[frame.landmass_id.isin(cohort)]
    if len(selected) != 9739 or selected.speciesKey.nunique() != 129:
        raise ValueError("Source-selected original-51 observational grain drift")
    years = defaultdict(list)
    for sp, lid, year in frame[["speciesKey", "landmass_id", "year"]].itertuples(index=False, name=None):
        years[(int(sp), str(lid))].append(int(year))
    nnear = nnew = nordered = 0
    species = frame.speciesKey.drop_duplicates().tolist()
    for pair in pairs:
        if float(pair["port_distance_km"]) >= 200:
            continue
        nnear += 1
        ia, ib = str(pair["island_a"]), str(pair["island_b"])
        new = ordered = 0
        for species_key in species:
            a = sorted(years.get((int(species_key), ia), ()))
            b = sorted(years.get((int(species_key), ib), ()))
            if not (a and b and 2011 <= a[0] <= 2022 and 2011 <= b[0] <= 2022):
                continue
            new += 1
            if a[0] < b[0]:
                early, late = a, b
            elif b[0] < a[0]:
                early, late = b, a
            else:
                early, late = (), ()
            ordered += int(bool(early) and any(early[0] < y < late[0] and y <= early[0]+5 for y in early))
        if new != int(pair["shared_new_species"]):
            raise ValueError(f"Original near-island first-index event count drift {ia},{ib}")
        nnew += new
        nordered += ordered
    if (len(pairs), nnear, nnew, nordered) != (1275, 14, 41, 21):
        raise ValueError("Original near-200-km temporal replay drift")
    focal = []
    for name, (species_key, endids) in FIVE.items():
        yr = {i: sorted(years.get((species_key, i), ())) for i in endids}
        if not all(yr.values()):
            raise ValueError(f"Historical focal year endpoint missing for {name}")
        first = {k: v[0] for k, v in yr.items()}
        a, b = list(yr.values())
        early, late = (a,b) if a[0] < b[0] else (b,a) if b[0] < a[0] else ([],[])
        positive = bool(early) and any(early[0] < y < late[0] and y <= early[0]+5 for y in early)
        if not positive:
            raise ValueError(f"Original preselected year order not replayed for {name}")
        focal.append({"butterfly": name, "first_year_by_landmass": first,
                      "indexed_prior_re_record_within_5y": bool(positive)})
    return {
        "schema": "chocho_original_private_history_full_row_audit_20261120_v1",
        "status": "PRIVATE_SOURCE_REPLAYED_NOT_PUBLICATION_ARCHIVE",
        "input_population": {"rows": profile[0], "butterflies": profile[1],
                             "landmasses": profile[2], "year_range": [profile[3],profile[4]],
                             "record_sum": profile[5]},
        "landmass_id_zero": {
            "rows": len(zero), "species": int(zero.speciesKey.nunique()),
            "record_sum": int(zero.records.sum()),
            "fraction_of_all_records": float(zero.records.sum()/frame.records.sum()),
            "landmass_area_km2": float(zero.landmass_area_km2.iloc[0]),
            "is_valid_populated_id_not_missing_sentinel": True,
            "in_published_51_island_subset": False},
        "precision_field": "precise1km_records equals records on all 83573 original source rows; this is a pinned input property, not new evidence about raw original coordinates",
        "original_51": {"islands":len(cohort), "rows":len(selected), "butterflies_observed":int(selected.speciesKey.nunique())},
        "near_island_temporal_replay": {"pairs_lt_200km":nnear,
              "joint_initial_index_2011_2022":nnew, "earlier_island_reindexed_before_later_initial_index":nordered,
              "original_published_pair_counts_match":True},
        "focal_preselected_five":focal,
        "gbif_original_publisher_identifiers":"NOT RECOVERED: neither gbifID nor datasetKey in source Parquet",
        "scientific_warning":"Indexed sequence is not shipborne transport, butterfly migration, breeding or true colonization; nonindependent island pair records",
    }


def check_private_archive(src: Path):
    import pandas as pd
    with ZipFile(src) as z:
        required = [PARQUET_NAME, ISLANDS_NAME, PAIRS_NAME]
        if any(path not in z.namelist() for path in required):
            raise ValueError("Original ZIP lacks required fixed paths")
        data = z.read(PARQUET_NAME)
        if sha256(data).hexdigest() != PINNED_PARQUET_SHA256:
            raise ValueError("Original private Parquet fingerprint mismatch")
        frame = pd.read_parquet(io.BytesIO(data), engine="pyarrow")
        islands = list(csv.DictReader(io.StringIO(z.read(ISLANDS_NAME).decode("utf8-sig"))))
        pairs = list(csv.DictReader(io.StringIO(z.read(PAIRS_NAME).decode("utf8-sig"))))
    out = replay(frame, islands, pairs)
    out["private_parquet_sha256"] = sha256(data).hexdigest()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--private-zip",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    if a.self_test:
        assert check_identifier_area(["0","1","2"],[50654050.6945,100.0,25.0]) == {"0","1","2"}
        try:
            check_identifier_area(["-1"],[4.0])
            raise AssertionError("negative ID should have failed")
        except ValueError:
            pass
        print("ORIGINAL_LANDMASS_ZERO_IS_ALLOWED_NOT_MISSING_SELF_TEST_PASS")
        return
    if not a.private_zip:
        ap.error("Private ZIP required for complete source replay, use --self-test for public test")
    report=check_private_archive(a.private_zip)
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    main()
