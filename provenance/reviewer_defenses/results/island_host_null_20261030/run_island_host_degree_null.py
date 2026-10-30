from pathlib import Path
from collections import Counter
import hashlib,json,random,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'island_host_degree_null_lock.json',W/'all225_named_host_edges_local.parquet',
       O/'island_shared_host_match_species_pairs.csv',O/'port_buffer_landmass_audit.csv',
       O/'port_alien_host_broad_all_ports.csv',W/'port_named_host_regional_status.parquet',O/'forecast_candidate_disposition.csv']
lock=json.loads(files[0].read_text());rng=random.Random(lock['seed'])
family_mode='--family' in sys.argv
prefix='island_host_family_null' if family_mode else 'island_host_degree_null'
if family_mode:
    files[0]=O/'island_host_family_null_lock.json'
    files.append(R/'work/remote_sensing/host_source/extracted/inputs/original/insect_host_accepted.csv')
    lock=json.loads(files[0].read_text())
edges=pd.read_parquet(files[1])[['species','host_name']].drop_duplicates()
c=pd.read_csv(files[6]);assert len(c)==225 and set(edges.species)<=set(c.species)
spid={s:i for i,s in enumerate(c.species)};hostid={s:i for i,s in enumerate(sorted(set(edges.host_name)))}
edge=[(spid[r.species],hostid[r.host_name]) for r in edges.itertuples()];original=set(edge);current=set(edge)
degree_s=Counter(s for s,h in edge);degree_h=Counter(h for s,h in edge)
family_by_host={}; family_groups={}
if family_mode:
    fam=pd.read_csv(files[-1],usecols=['accepted_name','family']).dropna().drop_duplicates()
    counts=fam.groupby('accepted_name').family.nunique()
    lookup=fam[fam.accepted_name.isin(counts[counts.eq(1)].index)].set_index('accepted_name').family.to_dict()
    family_by_host={j:lookup.get(name) for name,j in hostid.items()}
    for idx,(s,h) in enumerate(edge):
        family=family_by_host[h]
        if family is not None:family_groups.setdefault(family,[]).append(idx)
    initial_family_degree=Counter((s,family_by_host[h]) for s,h in edge)
geo=pd.read_csv(files[3],dtype={'single_landmass_id':str});ports=pd.read_csv(files[4],keep_default_na=False)
complete=set(ports.loc[ports.query_complete.astype(str).eq('True'),'port_id'])
geo=geo[geo.status.eq('SINGLE_LANDMASS')&geo.islands_intersecting_5km_buffer.eq(1)&geo.port_id.isin(complete)]
h=pd.read_parquet(files[5]);h=h.merge(geo[['port_id','single_landmass_id']],on='port_id',validate='many_to_one')
islandbits={}
for lid,g in h.groupby('single_landmass_id'):
    states=g.groupby('host_name').regional_catalogue_status.agg(set)
    islandbits[lid]=sum(1<<hostid[name] for name,st in states.items() if st=={'contemporary_not_native_in_catalogue'} and name in hostid)
d=pd.read_csv(files[2],dtype={'island_a':str,'island_b':str});d=d[d.new_in_both].copy()
species=[spid[x] for x in d.species];shared=[islandbits[a]&islandbits[b] for a,b in zip(d.island_a,d.island_b)]
frequency=Counter(species);weights=np.array([1/frequency[s]/len(frequency) for s in species])
def score():
    diet=[0]*225
    for s,h in edge:diet[s]|=1<<h
    matched=np.array([bool(diet[s]&bits) for s,bits in zip(species,shared)])
    return matched,float(matched.mean()),float(weights@matched)
observed,obs,obs_equal=score()
np.testing.assert_array_equal(observed,d.shared_matching_host)
accepted=0
def swaps(n):
    global accepted
    for _ in range(n):
        i=rng.randrange(len(edge));a,x=edge[i]
        if family_mode:
            family=family_by_host[x]
            if family is None:continue
            j=rng.choice(family_groups[family])
        else:j=rng.randrange(len(edge))
        b,y=edge[j]
        if a==b or x==y or (a,y) in current or (b,x) in current:continue
        current.remove((a,x));current.remove((b,y));current.add((a,y));current.add((b,x))
        edge[i]=(a,y);edge[j]=(b,x);accepted+=1
swaps(lock['burn_attempts_per_edge']*len(edge))
draw=[]
for k in range(lock['draws']):
    swaps(lock['between_draw_attempts_per_edge']*len(edge))
    assert len(current)==len(edge)==len(set(edge))
    assert Counter(s for s,h in edge)==degree_s and Counter(h for s,h in edge)==degree_h
    if family_mode:assert Counter((s,family_by_host[h]) for s,h in edge)==initial_family_degree
    _,raw,equal=score();draw.append(dict(draw=k+1,fraction_matched=raw,species_equal_fraction=equal,original_edge_retention=len(current&original)/len(edge)))
sim=pd.DataFrame(draw);sim.to_csv(O/(prefix+'_draws.csv'),index=False)
audit=dict(status='DEGREE_PRESERVING_STRUCTURAL_NULL_NOT_CAUSAL',edges=len(edge),candidate_species=225,
           species_with_edges=len(degree_s),hosts=len(degree_h),both_new_species_pairs=len(d),both_new_species=len(frequency),
           degree_checks_all_draws=True,accepted_swaps=accepted,draws=lock['draws'],
           observed_fraction=obs,observed_species_equal_fraction=obs_equal,
           null_mean_fraction=float(sim.fraction_matched.mean()),null_fraction_quantiles=sim.fraction_matched.quantile([.025,.5,.975]).to_dict(),
           null_mean_species_equal=float(sim.species_equal_fraction.mean()),null_species_equal_quantiles=sim.species_equal_fraction.quantile([.025,.5,.975]).to_dict(),
           draws_at_least_observed=int(sim.fraction_matched.ge(obs).sum()),
           draws_at_least_observed_species_equal=int(sim.species_equal_fraction.ge(obs_equal).sum()),
           lag1_fraction_autocorrelation=float(sim.fraction_matched.autocorr()),
           limits=lock['limits'],input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
if family_mode:
    audit['family_degree_checks_all_draws']=True
    audit['unresolved_family_hosts_frozen']=sum(v is None for v in family_by_host.values())
    audit['unresolved_family_edges_frozen']=sum(family_by_host[h] is None for s,h in original)
(O/(prefix+'_audit.json')).write_text(json.dumps(audit,indent=2),encoding='utf-8')
fig,axs=plt.subplots(1,2,figsize=(9,4.5))
for ax,col,ob,title in [(axs[0],'fraction_matched',obs,'Species-pair weighted'),(axs[1],'species_equal_fraction',obs_equal,'Species weighted equally')]:
    ax.hist(sim[col]*100,bins=15,color='#83b8b6',edgecolor='white');ax.axvline(ob*100,color='#b84b36',lw=2,label=f'Observed {ob*100:.1f}%')
    ax.set_title(title);ax.set_xlabel('Shared-new records with a matching shared host (%)');ax.set_ylabel('Randomized networks');ax.legend(fontsize=8)
fig.suptitle('Host matching: '+('family- and degree-preserving swaps' if family_mode else 'degree-preserving network swaps'))
fig.text(.02,.015,'199 MCMC draws; fixed island records. Null preserves degrees, not phylogeny or geography. No causal interpretation.',fontsize=8)
if family_mode:fig.texts[-1].set_text('199 MCMC draws; family counts preserved, unresolved-family edges frozen. Geography not preserved; no causal interpretation.')
fig.tight_layout(rect=[0,.05,1,.94]);fig.savefig(O/(prefix+'.png'),dpi=180);plt.close(fig)
print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))
