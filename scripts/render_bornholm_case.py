"""Render Box 1 from GSHHG intermediate level 1, WPI 2019 and SKIB31 sources."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import shapefile
from shapely.geometry import shape, box

root = Path(__file__).resolve().parents[1]
out = root / "figures/integrated_island_port_butterfly_v03"
parser = argparse.ArgumentParser()
parser.add_argument("sources", type=Path, help="Directory with GSHHS_shp, ports_2019.geojson and multiyear_calls")
src = parser.parse_args().sources
ports = json.loads((src / "ports_2019.geojson").read_text(encoding="utf8"))["features"]
coords = {p["properties"]["PORT_NAME"]: p["geometry"]["coordinates"] for p in ports
          if p["properties"]["COUNTRY"] == "DK"}
trips = pd.read_csv(src / "multiyear_calls/skib31_2020_23.csv", sep=";")
rows = trips[trips.iloc[:, 0] == "Køge-Rønne"].sort_values("TID")
assert list(rows.INDHOLD.astype(int)) == [382, 395, 364, 363]
assert list(rows.TID) == [2020, 2021, 2022, 2023]
coast = shapefile.Reader(str(src / "GSHHS_shp/i/GSHHS_i_L1.shp"))
fig = plt.figure(figsize=(10, 5.6), facecolor="white")
ax = fig.add_axes([.025, .12, .55, .74])
detail = fig.add_axes([.64, .43, .32, .43])
bounds = (11.1, 54.60, 17.25, 57.35)
geoms = []
for record in coast.iterShapeRecords(bbox=bounds):
    geom = shape(record.shape.__geo_interface__).intersection(box(*bounds))
    geoms.append((record.record["id"], geom))
def draw(a, limits):
    x0,y0,x1,y1 = limits
    a.set(xlim=(x0,x1), ylim=(y0,y1), aspect=1/np.cos(np.deg2rad((y1+y0)/2)),
          xticks=[], yticks=[], facecolor="#e9f3f6")
    for ident, geom in geoms:
        clipped = geom.intersection(box(*limits))
        for poly in getattr(clipped, "geoms", [clipped]):
            if poly.geom_type == "Polygon" and not poly.is_empty:
                x,y = poly.exterior.xy
                a.fill(x,y, color="#f7edcc" if ident in ("98","262","454") else "#f9faf8",
                       ec="#86989a", lw=.5)
    for spine in a.spines.values():
        spine.set_color("#a7bbbf")
draw(ax,bounds)
draw(detail,(14.60,54.94,15.23,55.33))
route = np.array([coords["KOGE"],coords["RONNE"]])
ax.plot(*route.T,color="#70408c",lw=1.8,marker="s",ms=4)
for x,y,label in [(11.24,55.78,"Zealand\nEarlier record"),(16.0,57.10,"Öland"),(14.77,54.82,"Bornholm"),
                  (11.72,55.23,"Køge"),(12.8,56.35,"SWEDEN")]:
    ax.text(x,y,label,fontsize=10,color="#263e43")
ax.text(12.78,55.01,"Ferry service",color="#70408c",fontsize=9)
ax.annotate("Nearest earlier island\n172 km from Bornholm", (16.59,56.72), (14.40,56.57), fontsize=9,color="#263e43",arrowprops=dict(arrowstyle="-",color="#476169"))
ax.add_patch(Rectangle((14.60,54.94),.63,.39,fill=False,ec="#476169",lw=.8))
for port,label,dy in [("HASLE","Hasle",.025),("RONNE","Rønne",-.015)]:
    x,y = coords[port]
    detail.scatter(x,y,c="#00765b",s=31,zorder=3)
    detail.text(x+.045,y+dy,label,fontsize=10,weight="bold")
fig.text(.03,.93,"A   Earlier butterfly records and ferry connection",weight="bold",fontsize=12)
fig.text(.61,.93,"B   Host records on Bornholm",weight="bold",fontsize=12)
fig.text(.61,.36,"Host plants recorded in 2001–2010",fontsize=11,weight="bold")
fig.text(.61,.25,"Both ports:  Malus domestica, Pyrus communis,\n                    Salix viminalis",fontsize=10,linespacing=1.4)
fig.text(.61,.16,"Hasle:  Populus nigra\nRønne:  Salix alba",fontsize=10,linespacing=1.4)
fig.text(.03,.045,"Bornholm, 2024    •    Three records of Nymphalis polychloros in its first indexed year here since 2001.",fontsize=11)
for ext in ("png","svg"):
    fig.savefig(out / ("Box1_bornholm."+ext),dpi=220)
svg = out / "Box1_bornholm.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf8").splitlines()) + "\n", encoding="utf8")
