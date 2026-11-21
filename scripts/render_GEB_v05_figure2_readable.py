#!/usr/bin/env python3
"""Re-layout source Figure 2 without estimating new ecological quantities.

A-C use the verified numeric receipt. D retains the archived fitted paths
under an affine SVG transform: no new model, interpolation or extrapolation.
The numeric receipt remains validated against the original frozen analyses.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)
FIGDIR = Path("figures/integrated_island_port_butterfly_v03")
RECEIPT = FIGDIR / "FIGURE_V04_SOURCE_NUMBERS.json"
ARCHIVE = FIGDIR / "Fig2_ports_native_alien_floras_isolation_v05.svg"
ARCHIVE_SHA256 = "3cabb25b60227a90aede73273d36f933fab6e33edad0d34fbafd33833312fb82"
OUTPUT = FIGDIR / "Fig2_ports_native_alien_floras_readable_v06.svg"
WIDTH, HEIGHT = 1100, 1120
INK, MUTED, GRID = "#172c40", "#4f6473", "#dae3e9"
INTRO, NATIVE, CONTRAST = "#bd673c", "#2e6876", "#7e6b87"


def node(parent, tag, text=None, **attrs):
    e = ET.SubElement(parent, f"{{{NS}}}{tag}", {
        k.replace("_", "-"): str(v) for k, v in attrs.items()
    })
    if text is not None:
        e.text = str(text)
    return e


def label(parent, x, y, value, size=19, color=INK, weight=400, anchor="start", **attrs):
    return node(parent, "text", value, x=x, y=y, font_size=size,
                font_family="Arial, Helvetica, sans-serif", fill=color,
                font_weight=weight, text_anchor=anchor, **attrs)


def line(parent, x1, y1, x2, y2, color=GRID, width=1.2, **attrs):
    return node(parent, "line", x1=x1, y1=y1, x2=x2, y2=y2,
                stroke=color, stroke_width=width, **attrs)


def panel(parent, ident, x, y, w, h, title):
    g = node(parent, "g", id=f"panel-{ident}", role="group", aria_label=title)
    node(g, "rect", x=x, y=y, width=w, height=h, rx=12, fill="#f6f8fa")
    label(g, x+20, y+34, ident, size=27, weight=700)
    label(g, x+58, y+34, title, size=24, weight=700)
    return g


def interval(parent, ident, y, estimate, limits, scale, color):
    g = node(parent, "g", id=ident, data_estimate=repr(estimate),
             data_lower=repr(limits[0]), data_upper=repr(limits[1]))
    lo, hi, mid = scale(limits[0]), scale(limits[1]), scale(estimate)
    line(g, lo, y, hi, y, color, 3)
    for x in (lo, hi):
        line(g, x, y-8, x, y+8, color, 2)
    node(g, "circle", cx=mid, cy=y, r=6, fill=color)


def build(root):
    source = root / ARCHIVE
    if hashlib.sha256(source.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("Archived Figure 2 changed; recover the pinned source before rendering")
    archived = ET.parse(source).getroot()
    receipt = json.loads((root / RECEIPT).read_text())
    f = receipt["figures"]["fig2"]
    isolation_source = json.loads((root / f["isolation_panel"]["source_receipt"]).read_text())
    curves = [e for e in archived.iter(f"{{{NS}}}path")
              if e.get("stroke") in (NATIVE, INTRO)]
    if len(curves) != 2:
        raise ValueError("Expected the two original fixed239 fitted paths")
    if f["isolation_panel"]["cohort"] != {
            "islands": 126, "active_islands": 123, "fixed_butterflies": 239}:
        raise ValueError("Source isolation cohort changed")

    svg = ET.Element(f"{{{NS}}}svg", {
        "width": str(WIDTH), "height": str(HEIGHT),
        "viewBox": f"0 0 {WIDTH} {HEIGHT}", "role": "img",
        "aria-labelledby": "title desc",
    })
    node(svg, "title", "Ports, island floras and butterfly host opportunities", id="title")
    node(svg, "desc",
         "Four source-defined panels. A: port coefficients and HC3 intervals for "
         "126 island floras, including the uncertain introduced-minus-native contrast. "
         "B: native and introduced flora median Jaccard over 7,875 pairs of those islands. "
         "C: relative MSE reduction from undirected marine electrical resistance in a "
         "separate 104-island cohort; both-islands-unseen intervals include zero. "
         "D: Mainland isolation and candidate butterfly food. Fixed panel of 239 butterfly "
         "species across 126 islands; native-host coverage is nested within contemporary coverage. "
         "The unchanged fitted count curves use log1p distance. A separate share model "
         "uses 123 resource-active islands. Its nine-Level1 p is 0.014; the paired "
         "count-slope contrast p is 0.068 and its interval includes zero. "
         "Opportunities are mapped known-host availability, not butterfly occurrence.", id="desc")
    node(svg, "metadata", json.dumps({
        "renderer": "scripts/render_GEB_v05_figure2_readable.py",
        "receipt": str(RECEIPT), "archived_source": str(ARCHIVE),
        "archived_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "display_width_inches": 6.6, "minimum_font_px": 19,
        "minimum_font_pt_at_display_width": 19*6.6*72/WIDTH,
        "isolation_curve_rule": "unchanged archived SVG paths, affine layout transform only",
    }, sort_keys=True))
    node(svg, "rect", x=0, y=0, width=WIDTH, height=HEIGHT, fill="white")
    label(svg, 24, 38, "2  |  Ports, island floras and butterfly host opportunities", size=29, weight=700)
    label(svg, 24, 68, "Source cohorts and model outcomes are identified separately in each panel.", color=MUTED)

    a = panel(svg, "A", 24, 92, 518, 388, "Ports and plant richness")
    label(a, 44, 158, "126 islands · area and Level1 adjusted", color=MUTED)
    xa = lambda value: 76 + (value+.2)/.8*422
    for value in (-.2, 0, .2, .4, .6):
        x = xa(value)
        for y0, y1 in ((208, 239), (275, 306), (342, 373)):
            line(a, x, y0, x, y1, color=MUTED if value == 0 else GRID,
                 stroke_dasharray="5 4" if value == 0 else "none")
        label(a, x, 405, f"{value:.1f}", color=MUTED, anchor="middle")
    for key, name, color, y in [
        ("introduced", "Introduced flora", INTRO, 222),
        ("native", "Native flora", NATIVE, 289),
        ("introduced_minus_native", "Introduced − native", CONTRAST, 356),
    ]:
        item = f["port_coefficient"][key]
        label(a, 44, y-22, name, weight=600)
        label(a, 516, y-22, f"{item['estimate']:+.4f}", color=color, weight=700, anchor="end")
        interval(a, "port-"+key, y, item["estimate"], item["descriptive_HC3_95"], xa, color)
    label(a, 283, 437, "Port coefficient per log₂(1 + ports)", anchor="middle")
    label(a, 44, 465, "Whiskers: descriptive HC3 95% intervals", color=MUTED)

    b = panel(svg, "B", 558, 92, 518, 388, "Island floristic overlap")
    label(b, 578, 158, "126 islands · 7,875 island–island pairs", color=MUTED)
    xb = lambda value: 598 + value/.06*418
    for value in (0, .02, .04, .06):
        x = xb(value)
        line(b, x, 192, x, 378)
        label(b, x, 405, f"{value:.2f}", color=MUTED, anchor="middle")
    for key, name, color, y in [
        ("introduced_median", "Introduced", INTRO, 238),
        ("native_median", "Native", NATIVE, 335),
    ]:
        value = f["whole_flora_jaccard"][key]
        label(b, 598, y-23, name, weight=600)
        node(b, "rect", id="jaccard-"+key, x=598, y=y-6,
             width=xb(value)-598, height=28, rx=4, fill=color, data_value=repr(value))
        label(b, xb(value)+12, y+16, f"{value:.6f}", color=color, weight=700)
    label(b, 807, 437, "Median Jaccard", anchor="middle")
    ratio = f["whole_flora_jaccard"]["introduced_median"]/f["whole_flora_jaccard"]["native_median"]
    label(b, 578, 465, f"Introduced ≈{ratio:.2f}× native median", color=MUTED)

    c = panel(svg, "C", 24, 500, 518, 585, "Marine-circuit prediction")
    label(c, 44, 566, "104 islands · 5,355 available pairs", color=MUTED)
    xc = lambda value: 172 + (value+2)/8.5*260
    for value in (-2, 0, 2, 4, 6):
        x = xc(value)
        for y0, y1 in ((617, 684), (736, 804)):
            line(c, x, y0, x, y1, color=MUTED if value == 0 else GRID,
                 stroke_dasharray="5 4" if value == 0 else "none")
        label(c, x, 834, (f"{value:+d}" if value else "0")+"%", color=MUTED, anchor="middle")
    for source_key, group_title, y0 in [
        ("source_original_one_island_heldout", "One island held out", 600),
        ("strict_both_endpoints_unseen", "Both islands unseen", 719),
    ]:
        model = f["marine_circuit"][source_key]
        label(c, 44, y0, group_title, weight=700)
        for key, name, color, y in [
            ("introduced", "Introduced", INTRO, y0+38),
            ("native", "Native", NATIVE, y0+76),
        ]:
            value = model[key+"_MSE_gain_percent"]
            label(c, 44, y+6, name, color=color)
            interval(c, source_key+"-"+key, y, value,
                     model[key+"_descriptive_interval"], xc, color)
            label(c, 520, y+6, f"{value:+.3f}%", color=color, weight=700, anchor="end")
    label(c, 282, 870, "Relative MSE reduction", anchor="middle")
    label(c, 44, 914, "Both unseen: 1,089 pairs in five folds", color=MUTED)
    label(c, 44, 944, "Whiskers: descriptive source intervals", color=MUTED)
    label(c, 44, 974, "Strict intervals include zero.", weight=600)
    label(c, 44, 1024, "Undirected ocean resistance;", color=MUTED)
    label(c, 44, 1053, "ocean currents were not measured.", color=MUTED)

    d = panel(svg, "D", 558, 500, 518, 585, "Isolation and host opportunity")
    label(d, 578, 566, "126 islands · fixed 239 butterflies", color=MUTED)
    for color, name, y in [(INTRO, "Native + introduced hosts", 596),
                           (NATIVE, "Native hosts", 624)]:
        line(d, 580, y-5, 611, y-5, color, 4)
        label(d, 621, y, name, color=color)
    label(d, 613, 655, "Butterflies with ≥1 mapped known host", size=19)
    # The archived axes map 108..705 to the original log1p distance and
    # y=1176−1.9*n to fitted counts. Only display coordinates change.
    left, right, top, bottom = 622, 1040, 680, 827
    sx, sy = (right-left)/597, (bottom-top)/152
    tx, ty = left-sx*108, bottom-sy*1176
    for value in (0, 40, 80):
        y = bottom-(bottom-top)*value/80
        line(d, left, y, right, y)
        label(d, left-12, y+6, value, color=MUTED, anchor="end")
    # Copy exact archived tick coordinates, preserving log1p(km) spacing.
    for km, old_x in [(1, 109.19), (10, 236.86), (100, 402.90),
                      (1000, 574.67), (4000, 678.43)]:
        x = sx*old_x+tx
        line(d, x, top, x, bottom)
        label(d, x, 855, f"{km:,}", color=MUTED, anchor="middle")
    geometry = node(d, "g", id="source-isolation-count-curves",
                    transform=f"matrix({sx:.12f} 0 0 {sy:.12f} {tx:.12f} {ty:.12f})",
                    data_source=str(ARCHIVE))
    for e in curves:
        geometry.append(copy.deepcopy(e))
    points = [e for e in archived.iter(f"{{{NS}}}circle")
              if float(e.get("cy", 0)) > 1000]
    if len(points) != 6:
        raise ValueError("Expected all six original 100/1000/4000-km fitted anchors")
    for e in points:
        geometry.append(copy.deepcopy(e))
    label(d, 831, 883, "Mainland distance (km; log1p scale)", anchor="middle")
    label(d, 578, 920, "Introduced-only share · 123 islands", weight=700)
    ip = f["isolation_panel"]
    for km, x in [(100, 641), (1000, 810), (4000, 984)]:
        label(d, x, 949, f"{km:,} km", color=MUTED, anchor="middle")
        label(d, x, 981, f"{ip['adjusted_introduced_only_shares'][f'km{km}']:.0%}",
              color=INTRO, size=25, weight=700, anchor="middle")
    share_p = ip["introduced_only_share_nine_level1_t8_p"]
    contrast_p = ip["paired_shift_nine_level1_t8_p"]
    ci = isolation_source["models"]["paired_native_current_shift"]["distance_nine_cluster_95"]
    ci_low = f"{ci[0]:+.3f}".replace("-", "−")
    label(d, 578, 1011, f"Nine-Level1 p = {share_p:.3f} (share model)", color=MUTED)
    label(d, 578, 1041, f"Count-slope contrast: p = {contrast_p:.3f}", weight=700)
    label(d, 578, 1069, f"Nine-Level1 95% interval {ci_low} to {ci[1]:+.3f}", color=MUTED)
    label(svg, 24, 1111, "Mapped known-host opportunities; source populations, controls and uncertainty are detailed in the caption.", color=MUTED)
    ET.indent(svg, space="  ")
    return ET.tostring(svg, encoding="unicode")+"\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--check", action="store_true", help="Verify committed SVG matches its source-driven renderer")
    args = ap.parse_args()
    root = args.repo_root.resolve()
    result = build(root)
    target = root / OUTPUT
    if args.check:
        if not target.exists() or target.read_text() != result:
            raise SystemExit("Figure 2 differs from the source-driven readable layout")
        print("GEB_FIGURE2_READABLE_SOURCE_LAYOUT_PASS")
    else:
        target.write_text(result)
        print(target)


if __name__ == "__main__":
    main()
