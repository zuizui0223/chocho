#!/usr/bin/env python3
"""Render fixed old/new distance summaries without rerunning either analysis."""
from pathlib import Path
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/chocho-latitude-matplotlib')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/exploratory/island_latitude_null_20261030'
SCENARIOS = ['primary', 'source_records_ge3', 'late_records_ge2']
LABELS = ['Primary', 'Early source records ≥3', 'Later target records ≥2']


def main():
    x = json.loads((BASE / 'latitude_null_diagnostic.json').read_bytes())
    rows = [x['analyses'][name] for name in SCENARIOS]
    old = np.array([r['old_longitude_only_ratio'] for r in rows])
    new = np.array([r['observed_to_null_geometric_mean_one_plus_gap_ratio'] for r in rows])
    variable = np.array([r['support']['variable_positive_cells'] for r in rows])
    forced = np.array([r['support']['forced_positive_cells'] for r in rows])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.labelcolor': '#253443', 'text.color': '#253443',
                         'xtick.color': '#465563', 'ytick.color': '#253443',
                         'axes.edgecolor': '#98a3ab', 'svg.fonttype': 'none'})
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 6.0), gridspec_kw={'width_ratios': [1.45, 1]})
    fig.subplots_adjust(left=.20, right=.975, bottom=.33, top=.71, wspace=.31)
    fig.suptitle('Recorded-source proximity after latitude conditioning', x=.04, y=.965,
                 ha='left', fontsize=17, fontweight='bold')
    fig.text(.04, .898, 'Same observed pairs and coast-gap scores; a more geographically comparable conditional reference',
             fontsize=11, color='#516170')
    old_color, new_color, grey = '#b07b27', '#205c82', '#d4dce2'
    y = np.arange(3)[::-1]
    for yi, a, b in zip(y, old, new):
        left.plot([a,b], [yi,yi], color='#aeb8c0', lw=1.7, zorder=1)
        left.scatter([a], [yi], facecolors='white', edgecolors=old_color, s=90, linewidths=1.7, zorder=3)
        left.scatter([b], [yi], color=new_color, marker='D', s=62, zorder=3)
        left.annotate(f'{a:.3f}', (a,yi), xytext=(0,-21), textcoords='offset points', ha='center', color='#866026', fontsize=11)
        left.annotate(f'{b:.3f}', (b,yi), xytext=(0,13), textcoords='offset points', ha='center', color=new_color, fontsize=11, fontweight='bold')
    left.axvline(1, color='#6a7780', lw=1.1, linestyle='--')
    left.text(1, 2.6, 'Reference = 1', ha='center', va='bottom', fontsize=10, color='#62717b')
    left.set_xlim(0, 1.1)
    left.set_ylim(-.55, 2.65)
    left.set_xticks([0,.25,.5,.75,1])
    left.set_yticks(y, LABELS)
    left.tick_params(axis='y', length=0, pad=12)
    left.set_xlabel('Observed / conditional reference\nGeometric scale of (1 + coast gap)', labelpad=11)
    left.grid(axis='x', color='#edf0f2', linewidth=.8)
    left.set_axisbelow(True)
    left.legend(handles=[Line2D([],[],color=old_color,marker='o',markerfacecolor='white',lw=0,label='Longitude only'),
                         Line2D([],[],color=new_color,marker='D',lw=0,label='Longitude + latitude zone')],
                loc='lower left', bbox_to_anchor=(-.01,1.08), frameon=False, fontsize=10, handletextpad=.5)
    right.barh(y, variable, color=new_color, height=.48)
    right.barh(y, forced, left=variable, color=grey, edgecolor='#a9b5bd', linewidth=.6, height=.48)
    for yi, a, b in zip(y, variable, forced):
        right.text(a/2,yi,str(a),va='center',ha='center',color='white',fontweight='bold')
        right.text(a+b/2,yi,str(b),va='center',ha='center',color='#364654')
    right.set_xlim(0,70)
    right.set_ylim(left.get_ylim())
    right.set_yticks(y, ['']*3)
    right.tick_params(axis='y',length=0)
    right.set_xticks([0,20,40,60])
    right.set_xlabel('Later-positive species–island pairs\nAfter latitude conditioning',labelpad=11)
    right.legend(handles=[Patch(facecolor=new_color,label='Can vary'),
                          Patch(facecolor=grey,edgecolor='#a9b5bd',label='Fixed by margins')],
                 loc='lower left',bbox_to_anchor=(0,1.08),frameon=False,fontsize=10)
    fig.text(.04,.10, 'Primary: 36 variable positive pairs across 26 species and 11 islands.\n'
             '95.9% of the remaining net distance contrast is concentrated in the tropical 90°W–0° block.',
             fontsize=10.5,color='#334957',linespacing=1.5)
    fig.text(.04,.025,'Early records: 2010–2017; later records: 2018–2025. Exploratory reanalysis; latitude zones do not fully measure climate.\n'
             'Ratios are conditional descriptive summaries; no ecological-effect confidence interval or causal port/host effect is shown.',
             fontsize=9,color='#657580',linespacing=1.4)
    for extension in ('png','svg'):
        fig.savefig(BASE / ('latitude_conditioned_distance_comparison_20261030.' + extension),dpi=170,facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    main()
