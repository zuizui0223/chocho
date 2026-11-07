"""Render a Japanese comparison table and interval plot for a validated forecast."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["temporal", "spatial"])
mode = parser.parse_args().mode
out = Path(__file__).resolve().parents[2] / "outputs"
validation = json.loads((out / f"route_terrain_{mode}_validation.json").read_text())
assert validation["status"] == "PASS_PAIRED_PREDICTION_REPLAY" and validation["mode"] == mode
data = pd.read_csv(out / f"route_terrain_{mode}_comparison.csv")
names = {"all155": "全155種", "area_supported_islands": "25島", "islands_without_GB_Ireland": "英国・アイルランドを除く23島", "at_least100km": "学習地点から100 km以上"}
cohorts = list(names)[:3] if mode == "temporal" else list(names)
refs = {"M5_host_composition": "食草・気候・土地利用", "B5_sea_terrain": "海・セル地形等を含むモデル"}
assert set(zip(data.cohort, data.reference)) == {(c,r) for c in cohorts for r in refs}
assert len(data) == len(cohorts)*2 and data["mode"].eq(mode).all()
label = {"temporal": "時間外", "spatial": "空間外"}[mode]
plt.rcParams["font.family"] = "Yu Gothic"
fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
lines = [f"# 経路の標高・起伏とチョウの記録予測（{label}検証）", "",
         "候補225種のうち適格155種を評価した。地形指標は予測時点以前の最近接記録と対象地点を結ぶ測地線上の地表高から算出した。指標の欠測は専用のフラグで保持した。", "",
         "|対象|種数|評価行|比較対象|MSE改善率|MSE差の95%区間|", "|---|---:|---:|---|---:|---|"]
for ax, (ref, title) in zip(axes, refs.items()):
    subset = data[data.reference.eq(ref)].set_index("cohort").loc[cohorts]
    for y, (cohort, row) in enumerate(subset.iterrows()):
        ax.hlines(y, row.interval_low, row.interval_high, color="#285d78")
        ax.plot(row.mse_difference, y, "o", color="#164662")
        lines.append(f"|{names[cohort]}|{row.species}|{row.rows:,}|{title}|{row.improvement_pct:+.2f}%|[{row.interval_low:+.6f}, {row.interval_high:+.6f}]|")
    ax.axvline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_title(title + "との比較", fontsize=11)
    ax.set_yticks(range(len(cohorts)), [names[c] for c in cohorts])
    ax.tick_params(axis="x", labelsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.spines[["right", "top"]].set_visible(False)
axes[0].invert_yaxis()
fig.suptitle(f"経路の標高・起伏を含むモデルの予測誤差（{label}検証）")
fig.supxlabel("平均二乗誤差の差（負の値で改善、横軸の尺度は図ごとに異なる）", fontsize=10)
fig.tight_layout()
stem = out / f"route_terrain_{mode}_ja"
fig.savefig(stem.with_suffix(".png"), dpi=160, bbox_inches="tight")
lines.extend(["", "MSE差は経路地形を含むモデルの値から各比較対象の値を引いたもの。区間は固定した予測を種と空間ブロックの両方で2,000回再標本化して算出。島の行は各検証の予測のうち、対象島に属する行を集計した。", "", f"![予測誤差の差]({stem.name}.png)"])
stem.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
