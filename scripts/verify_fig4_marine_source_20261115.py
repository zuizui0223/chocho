#!/usr/bin/env python3
"""Verify main GEB Figure 4 against original frozen marine prediction CSVs.

Figure 4 is intentionally a *subset* of the original source-locked Figure 5:
a readable comparison of positive held-island predictive associations and
negative strict both-endpoint geographic transfer. It is not a causal figure.
"""
from __future__ import annotations

import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "figures/integrated_island_port_butterfly_v03/FIGURE4_MARINE_SOURCE_RECEIPT.json"
NS = "{http://www.w3.org/2000/svg}"


def digests(path):
    raw=path.read_bytes()
    return {
        "sha256":hashlib.sha256(raw).hexdigest(),
        "git_blob":hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest(),
    }


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as file:
        return list(csv.DictReader(file))


def exactly(rs,match):
    found=[r for r in rs if all(r.get(k)==v for k,v in match.items())]
    if len(found)!=1:
        raise ValueError(f"Original marine source row missing or nonunique: {match}")
    return found[0]


def main():
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert receipt["schema"]=="chocho_GEB_v03_main_fig4_marine_20261115"
    original={}
    for k,src in receipt["sources"].items():
        f=ROOT/src["path"]
        d=digests(f)
        assert all(d[n]==src[n] for n in ("sha256","git_blob")), (k,d,src)
        original[k]=rows(f)

    v={
        "ocean":float(exactly(original["uniform"],{
            "comparison":"base_vs_uniform","target":"new_butterfly_jaccard"
        })["relative_MSE_reduction_percent"]),
        "host":float(exactly(original["joint"],{
            "comparison":"C_to_HC","target":"new_butterfly_jaccard"
        })["relative_MSE_reduction_percent"]),
        "aiswind":float(exactly(original["joint"],{
            "comparison":"H_to_HC","target":"new_butterfly_jaccard"
        })["relative_MSE_reduction_percent"]),
        "strict_ais":float(exactly(original["strict"],{
            "panel":"later_butterflies_retro_hosts",
            "contrast":"BASE_to_AIS","pair_class":"all_both_plants_indexed",
            "unique_source_pairs":"990"
        })["MSE_gain_percent"]),
        "strict_aiswind":float(exactly(original["strict"],{
            "panel":"later_butterflies_retro_hosts",
            "contrast":"BASE_to_BOTH","pair_class":"all_both_plants_indexed",
            "unique_source_pairs":"990"
        })["MSE_gain_percent"])
    }
    expect={
        "ocean":receipt["figure_panels"]["A"]["values"]["uniform_ocean_beyond_distance_MSE_reduction_pct"],
        "host":receipt["figure_panels"]["A"]["values"]["host_assemblages_beyond_AIS_wind_MSE_reduction_pct"],
        "aiswind":receipt["figure_panels"]["A"]["values"]["AIS_wind_beyond_host_assemblages_MSE_reduction_pct"],
        "strict_ais":receipt["figure_panels"]["B"]["values"]["AIS_only_MSE_reduction_pct"],
        "strict_aiswind":receipt["figure_panels"]["B"]["values"]["AIS_and_wind_MSE_reduction_pct"],
    }
    assert all(abs(v[k]-expect[k])<1e-9 for k in v),(v,expect)
    assert v["ocean"]>0 and v["host"]>0 and v["aiswind"]>0 and v["strict_ais"]>0 and v["strict_aiswind"]<0
    svg_path=ROOT/receipt["svg"]
    svg=ET.parse(svg_path).getroot()
    assert svg.tag==NS+"svg" and svg.attrib.get("viewBox")=="0 0 1100 780"
    assert "movement" in "".join(svg.itertext()).lower() or "transport" in "".join(svg.itertext()).lower()
    by_id={e.attrib["id"]:e for e in svg.iter() if e.attrib.get("id")}
    pairs={
        "fig4-ocean":v["ocean"],
        "fig4-host":v["host"],
        "fig4-aiswind":v["aiswind"],
        "fig4-strict-ais":v["strict_ais"],
        "fig4-strict-aiswind":v["strict_aiswind"],
    }
    for id_,value in pairs.items():
        expected=("+" if value>0 else "")+f"{value:.2f}%"
        assert "".join(by_id[id_].itertext())==expected,(id_,expected)
    assert receipt["figure_panels"]["A"]["different_reference_models_no_addition"] is True
    assert receipt["figure_panels"]["B"]["n_source_pairs"]==990
    assert receipt["figure_panels"]["A"]["n_source_landmasses"]==51

    print(json.dumps({
        "status":"SOURCE_LOCKED_FIGURE_4_PASS",
        "source_CSVs_verified":3,
        "SVG_source_labels_verified":len(pairs),
        "source_values":v,
        "source_restrictions":"Different holdouts and models, not added; no established transport"
    },sort_keys=True))


if __name__=="__main__":
    main()
