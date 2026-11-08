#!/usr/bin/env python3
"""Figure 3: port-associated total alien flora and distinct host-plant richness.

Sources are hashed original observational model results. Coefficients refer to
different log(1+plant species) responses. No CI is fabricated and this chart
must not be read as ship-mediated introductions or a causal mediation model.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
SOURCES=[
    ("results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json",
     "8178c6f22be868a8695ecad6d4c539c715d7eed6"),
    ("results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json",
     "b3c710384cd73c2b11bc21e7c2ea8f227be79e46"),
]


def sha_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def read_models():
    sources=[]
    for path,expected in SOURCES:
        raw=(ROOT/path).read_bytes()
        if sha_blob(raw)!=expected:
            raise ValueError(f"Original source blob has changed: {path}")
        sources.append(json.loads(raw))
    audit,conversion=sources
    models=[
        ("Entire introduced vascular flora",
         audit["models"]["all_flora_with_gbif_total_effort"],
         "regression_beta_per_log2_port_unit", "permutation", "p_two_sided_centered"),
        ("Distinct introduced larval-host taxa",
         audit["models"]["host_taxa_with_effort"],
         "regression_beta_per_log2_port_unit", "permutation", "p_two_sided_centered"),
        ("Larval-host taxa, controlling\nall introduced flora",
         conversion["model_results"]["M1eff_port_host_taxa_given_alien_flora"],
         "coefficient", "within_level1_and_area_half_label_permutation", "p_two_centered"),
    ]
    out=[]
    for label,raw,key,permk,pkey in models:
        assert raw["n"] == 98
        v=float(raw[key])
        pv=float(raw[permk][pkey])
        assert int(raw[permk]["n_permutations"] if permk=="permutation" else
                   raw[permk]["permutations"]) == 4999
        out.append({"label":label,"coefficient":v,"descriptive_p":pv,"n":raw["n"]})
    expected=[0.24447852550694815,0.16390595789001258,-0.041526655753407464]
    assert all(np.isclose(x["coefficient"],a,atol=1e-12) for x,a in zip(out,expected))
    return out


def render(rows,dest):
    dest.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":11,"axes.titlesize":14,
                         "pdf.fonttype":42,"svg.fonttype":"none"})
    fig,ax=plt.subplots(figsize=(12.2,5.5))
    fig.suptitle("Port-indexed islands contain broader introduced floras, including larval hosts",
                 fontsize=15.5,weight="bold",y=.975)
    y=[2,1,0]
    for v,item in zip(y,rows):
        x=item["coefficient"]
        ax.plot([0,x],[v,v],linewidth=5,color="#64839d" if x>0 else "#a0a5ad",
                solid_capstyle="round")
        ax.scatter([x],[v],s=115,zorder=3,color="#24435b" if x>0 else "#555a64")
        ax.text(x+(0.013 if x>=0 else -0.013),v,
                f"β {x:+.3f}   (descriptive p={item['descriptive_p']:.4f})",
                fontsize=11,va="center",ha="left" if x>=0 else "right")
    ax.axvline(0,color="black",lw=1)
    ax.set_yticks(y,[x["label"] for x in rows],fontsize=11)
    ax.set_xlim(-0.235,0.55)
    ax.set_ylim(-.62,2.55)
    ax.set_xlabel("Port-count coefficient on log(1 + distinct plant species) response")
    ax.grid(axis="x",alpha=.18)
    ax.set_axisbelow(True)
    fig.text(.5,.16,
             "Same 98 insular WGSRPD3 regions; exposure = log₂(1 + indexed WPI ports). "
             "Models adjust area, broad biogeography, native flora/resources and GBIF botanical effort.",
             ha="center",fontsize=10)
    fig.text(.5,.105,
             "The conditional host model also adjusts total introduced flora. "
             "No coefficient implies selective carriage by ships, biological plant mediation or butterfly establishment.",
             ha="center",fontsize=9.6)
    fig.text(.5,.052,
             "Descriptive within-stratum permutation p values; intervals are not shown because "
             "no suitable coefficient-confidence bounds were independently verified.",
             ha="center",fontsize=9)
    fig.subplots_adjust(left=.35,right=.95,top=.81,bottom=.34)
    checks={}
    for ext in ("png","svg","pdf"):
        path=dest/f"fig3_port_flora_host_coefficients.{ext}"
        fig.savefig(path,bbox_inches="tight",dpi=180)
        checks[ext]=hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    return checks


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,default=ROOT/"_jbi_fig3")
    args=parser.parse_args()
    values=read_models()
    checks=render(values,args.out_dir)
    receipt={"schema":"chocho_JBI_source_locked_port_flora_fig3_v1",
             "status":"SOURCE_VERIFIED_FIGURE_RENDERED",
             "input_blobs":[{"path":p,"sha":h} for p,h in SOURCES],
             "estimates":values,
             "files_sha256":checks,
             "not_measured":["ship carriage","causal mediation","established butterfly feeding"]}
    (args.out_dir/"fig3_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
