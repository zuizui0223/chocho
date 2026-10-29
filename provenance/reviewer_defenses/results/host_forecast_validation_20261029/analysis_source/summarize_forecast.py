"""Independent metric reconciliation and crossed-cluster descriptive uncertainty.

Intervals condition on the fitted predictions, not model refitting or new climates.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

R = Path(__file__).resolve().parents[2]
O = R / 'outputs'
W = R / 'work/temporal/forecast'
rng = np.random.default_rng(20261029)
comparisons = [('M5_host_composition', 'M3_both'),
               ('M5_host_composition', 'M4_host_resources'),
               ('M4_host_resources', 'M3_both'),
               ('M5_host_composition', 'M0_background'),
               ('M5_host_composition', 'persistence')]
rows, checks = [], []
for mode in ['temporal', 'species', 'spatial', 'spatial_buffer100', 'temporal_effort20', 'temporal_radius250', 'temporal_exclude_top_provider', 'temporal_precision100']:
    source_mode = 'spatial' if mode == 'spatial_buffer100' else mode
    receipt_path = O / f'forecast_{source_mode}_results.json'
    if not receipt_path.exists():
        continue
    receipt = json.loads(receipt_path.read_text())
    if receipt['status'] != 'COMPLETE':
        continue
    p = pd.read_parquet(W / f'predictions_{source_mode}.parquet')
    if mode == 'spatial_buffer100':
        p = p[p.distance_to_training_cell_km >= 100].copy()
        assert len(p)>0
        metrics={}
        for name in receipt['results']:
            error=p[name]-p.target
            metrics[name]={'rows':len(p),'species':int(p.species.nunique()),'species_equal_mse':float(error.pow(2).groupby(p.species).mean().mean()),'species_equal_mae':float(error.abs().groupby(p.species).mean().mean())}
        receipt={'status':'COMPLETE','mode':mode,'results':metrics,'cells':int(p.cell.nunique()),'interpretation':'Post-hoc geographic sensitivity: test cells >=100km from ANY training cell for their held spatial fold; same fitted models; not a separate refit'}
        (O/'forecast_spatial_buffer100_results.json').write_text(json.dumps(receipt,indent=2))
    assert not p.duplicated(['species', 'cell']).any()
    if mode == 'temporal':
        temporal = receipt
    for model, expected in receipt['results'].items():
        error = p[model] - p.target
        for name, values in [('mse', error ** 2), ('mae', abs(error))]:
            actual = values.groupby(p.species).mean().mean()
            assert np.isclose(actual, expected[f'species_equal_{name}'], atol=1e-12)
            if mode != 'spatial_buffer100':
                checks.append(dict(mode=mode, model=model, metric=name, value=float(actual)))
    species = sorted(p.species.unique())
    blocks = sorted(p.spatial_block.unique())
    count = pd.crosstab(p.species, p.spatial_block).reindex(index=species, columns=blocks, fill_value=0).to_numpy()
    # Independent resampling of both crossed dimensions preserves paired models.
    sw = rng.multinomial(len(species), np.full(len(species), 1/len(species)), size=2000)
    bw = rng.multinomial(len(blocks), np.full(len(blocks), 1/len(blocks)), size=2000)
    denominator = bw @ count.T
    for added, reference in comparisons:
        for metric in ['mse', 'mae']:
            ea, er = p[added] - p.target, p[reference] - p.target
            delta = ea**2 - er**2 if metric == 'mse' else abs(ea)-abs(er)
            q = p[['species', 'spatial_block']].assign(delta=delta)
            sums = q.groupby(['species', 'spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=species, columns=blocks, fill_value=0).to_numpy()
            means = np.divide(bw @ sums.T, denominator, out=np.zeros_like(denominator, dtype=float), where=denominator>0)
            weights = sw * (denominator>0)
            samples = (means*weights).sum(axis=1)/weights.sum(axis=1)
            byspecies = delta.groupby(p.species).mean()
            low, high = np.quantile(samples, [.025,.975])
            rows.append(dict(mode=mode, added=added, reference=reference, metric=metric,
                difference=float(byspecies.mean()), ci_low=float(low), ci_high=float(high),
                species_improved=int((byspecies<0).sum()), species=len(species), blocks=len(blocks), rows=len(p)))

pd.DataFrame(rows).to_csv(O/'forecast_paired_comparisons.csv', index=False)
(O/'forecast_metric_validation.json').write_text(json.dumps(dict(status='PASS', checks=checks,
    bootstrap='2000 crossed species x 5-degree spatial-block resamples; equal species weighting; paired fixed predictions; no model refitting; descriptive intervals',
    warning='Does not capture training, source observation, temporal replication or scenario uncertainty'), indent=2))
labels = ['Persistence','Background','+ Climate','+ Land cover','+ Climate & land','+ Host resources','+ Host composition']
fig, axes = plt.subplots(1,2,figsize=(11,5.5))
for ax, metric, title in zip(axes, ['mse','mae'], ['Mean squared error','Mean absolute error']):
    vals = [x[f'species_equal_{metric}'] for x in temporal['results'].values()]
    ax.barh(np.arange(len(vals)), vals, color=['#8b9299']+['#4c7899']*5+['#b56a35'])
    ax.set_yticks(np.arange(len(vals)), labels)
    ax.invert_yaxis()
    ax.set_xlabel(title+' (lower is better)')
    ax.set_xlim(0,max(vals)*1.23)
    for i,v in enumerate(vals): ax.text(v+max(vals)*.015,i,f'{v:.4f}',va='center',fontsize=9)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('Temporal transfer to 2019–2022: mixed evidence of predictive improvement',fontsize=13)
fig.text(.05,.025,'155 species; 208,482 species–cell pairs; species weighted equally.\nOutcome: standardized recording probability, not confirmed occupancy. Fixed-prediction uncertainty in companion CSV.',fontsize=9)
fig.tight_layout(rect=[0,.10,1,.94])
fig.savefig(O/'forecast_temporal_comparison.png',dpi=180)
plt.close(fig)
contrasts = pd.DataFrame(rows).query("mode == 'temporal'")
comparison_labels = ['All host information vs climate + land',
                     'Host composition vs host resources',
                     'Host resources vs climate + land',
                     'Full model vs background',
                     'Full model vs persistence']
fig, axes = plt.subplots(1,2,figsize=(13,5.5),sharey=True)
for ax, metric in zip(axes,['mse','mae']):
    a = contrasts[contrasts.metric==metric].reset_index(drop=True)
    for i,row in a.iterrows():
        ax.plot([row.ci_low,row.ci_high],[i,i],color='#4c7899',lw=2)
        ax.scatter(row.difference,i,color='#b56a35',s=36,zorder=3)
    ax.axvline(0,color='#777777',ls='--',lw=1)
    ax.set_yticks(range(5),comparison_labels)
    ax.set_xlabel(f'Change in {metric.upper()} (negative = improvement)')
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines[['top','right']].set_visible(False)
axes[0].invert_yaxis()
fig.suptitle('Host resource information helps; extra composition benefit remains uncertain',fontsize=13)
fig.text(.05,.025,'2019–2022 temporal test; 155 species, 100 spatial blocks. 95% intervals: 2,000 crossed species × block resamples.\nConditional on fixed fitted predictions; excludes model-fitting, observation and future-scenario uncertainty.',fontsize=9)
fig.tight_layout(rect=[0,.12,1,.94])
fig.savefig(O/'forecast_temporal_paired_uncertainty.png',dpi=180)
print(json.dumps(dict(status='PASS', metric_checks=len(checks), comparisons=len(rows))))
import subprocess,sys
subprocess.run([sys.executable,str(R/'work/temporal/reconcile_validation_coverage.py')],check=True)
subprocess.run([sys.executable,str(R/'work/temporal/summarize_initially_unrecorded.py')],check=True)
subprocess.run([sys.executable,str(R/'work/temporal/report_forecast_validation.py')],check=True)
