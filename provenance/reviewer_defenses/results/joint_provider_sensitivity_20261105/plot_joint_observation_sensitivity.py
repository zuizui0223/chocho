from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
O=Path(__file__).resolve().parents[2]/'outputs'
d=pd.read_csv(O/'joint_observation_sensitivity_uncertainty.csv')
labels={'temporal':'Main: 10 recording days','temporal_effort20':'20 recording days','temporal_exclude_top_provider':'Largest provider excluded'}
fig,axes=plt.subplots(1,2,figsize=(10,4.6))
for ax,(model,ref,title,color) in zip(axes,[('B5_sea_terrain','M5_host_composition','Natural geography added','#245980'),('J1_natural_transport','B5_sea_terrain','Ports added to geography','#a06320')]):
    s=d[(d.model==model)&(d.reference==ref)]
    assert len(s)==3
    for y,row in enumerate(s.itertuples()):
        ax.plot([row.interval_low,row.interval_high],[y,y],color=color,lw=2)
        ax.plot(row.mse_difference,y,'o',color=color,ms=7)
    ax.set_yticks(range(3),[labels[x] for x in s.variant]);ax.set_ylim(2.5,-.5)
    ax.axvline(0,color='#555555',ls='--',lw=1);ax.grid(axis='x',alpha=.2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3));ax.set_title(title)
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('MSE difference (added - reference)\nNegative = lower error')
fig.suptitle('Observation sensitivity: all 155 modeled species',fontsize=14)
fig.subplots_adjust(left=.25,right=.98,top=.82,bottom=.27,wspace=1.1)
fig.text(.03,.05,'95% percentile ranges: 2,000 species x spatial-block resamples of fixed predictions.\nCompare models within each condition only. Historical geography remains fixed after provider removal.',fontsize=9)
fig.savefig(O/'joint_observation_sensitivity.png',dpi=180)
