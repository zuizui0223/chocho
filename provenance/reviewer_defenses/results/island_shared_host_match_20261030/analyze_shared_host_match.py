from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'island_shared_host_match_lock.json',O/'island_host_network_pairs.csv',O/'forecast_candidate_disposition.csv',
       O/'port_buffer_landmass_audit.csv',O/'port_alien_host_broad_all_ports.csv',W/'port_named_host_regional_status.parquet',
       W/'all225_named_host_edges_local.parquet',O/'global_original_species_island_history.csv']
lock=json.loads(files[0].read_text());pairs=pd.read_csv(files[1],dtype={'island_a':str,'island_b':str})
c=pd.read_csv(files[2])[['species','taxon_key']];assert len(c)==225
names=c.set_index('taxon_key').species.to_dict()
geo=pd.read_csv(files[3],dtype={'single_landmass_id':str});ports=pd.read_csv(files[4],keep_default_na=False)
complete=set(ports.loc[ports.query_complete.astype(str).eq('True'),'port_id'])
geo=geo[geo.status.eq('SINGLE_LANDMASS')&geo.islands_intersecting_5km_buffer.eq(1)&geo.port_id.isin(complete)]
h=pd.read_parquet(files[5]);h=h[h.port_id.isin(geo.port_id)].copy();h['island']=h.port_id.map(geo.set_index('port_id').single_landmass_id)
hostsets={}
for lid,g in h.groupby('island'):
    states=g.groupby('host_name').regional_catalogue_status.agg(set)
    hostsets[lid]=set(states[states.map(lambda x:x=={'contemporary_not_native_in_catalogue'})].index)
e=pd.read_parquet(files[6]);diet={sp:set(g.host_name) for sp,g in e.groupby('species')}
b=pd.read_csv(files[7],dtype={'landmass_id':str});b=b[b.first_record_year.between(2011,2022)]
new={lid:set(g.speciesKey) for lid,g in b.groupby('landmass_id')}
rows=[];plantrows=[]
for pair in pairs.itertuples():
    a=new.get(pair.island_a,set());z=new.get(pair.island_b,set())
    shared=hostsets[pair.island_a]&hostsets[pair.island_b]
    for key in sorted(a|z):
        sp=names[key];matched=shared&diet.get(sp,set());both=key in a and key in z
        rows.append(dict(species=sp,taxon_key=key,island_a=pair.island_a,island_b=pair.island_b,
                         new_in_both=both,shared_matching_host=bool(matched),matching_host_count=len(matched),
                         matching_hosts=';'.join(sorted(matched))))
        if both:
            for plant in matched:plantrows.append(dict(species=sp,island_a=pair.island_a,island_b=pair.island_b,host_name=plant))
d=pd.DataFrame(rows);s=d[d.new_in_both]
d.to_csv(O/'island_shared_host_match_species_pairs.csv',index=False)
plant=pd.DataFrame(plantrows);plant.to_csv(O/'island_shared_host_match_plant_links.csv',index=False)
top=plant.groupby('host_name').agg(butterfly_species=('species','nunique'),species_pair_matches=('species','size')).reset_index()
top=top.merge(plant.drop_duplicates(['host_name','island_a','island_b']).groupby('host_name').size().rename('island_pairs'),on='host_name')
top=top.sort_values(['butterfly_species','island_pairs','host_name'],ascending=[False,False,True])
top.to_csv(O/'island_shared_host_match_top_plants.csv',index=False)
summary=[]
for r in c.itertuples():
    q=d[d.species.eq(r.species)];ss=q[q.new_in_both]
    valid=len(q)>=5 and q.new_in_both.nunique()==2 and q.shared_matching_host.nunique()==2
    summary.append(dict(species=r.species,taxon_key=r.taxon_key,union_new_pairs=len(q),both_new_pairs=len(ss),
                        matched_both_new_pairs=int(ss.shared_matching_host.sum()),
                        fraction_matched=float(ss.shared_matching_host.mean()) if len(ss) else np.nan,
                        phi=float(np.corrcoef(q.new_in_both.astype(int),q.shared_matching_host.astype(int))[0,1]) if valid else np.nan,
                        disposition='CORRELATION_EVALUABLE' if valid else 'NO_NEW_RECORD_IN_ELIGIBLE_ISLANDS' if not len(q) else 'INSUFFICIENT_VARIATION_OR_PAIRS'))
out=pd.DataFrame(summary);out.to_csv(O/'island_shared_host_match_all225.csv',index=False)
audit=dict(status='CATALOGUE_HOST_COMPATIBILITY_NOT_MECHANISM',union_species_pairs=len(d),shared_new_species_pairs=len(s),
           shared_new_species=s.species.nunique(),matched_shared_new_species_pairs=int(s.shared_matching_host.sum()),
           matched_shared_new_species=s.loc[s.shared_matching_host,'species'].nunique(),
           fraction_matched=float(s.shared_matching_host.mean()),species_equal_fraction_matched=float(out.fraction_matched.mean()),
           evaluable_species=int(out.phi.notna().sum()),mean_species_phi=float(out.phi.mean()),median_species_phi=float(out.phi.median()),
           limits=lock['limits'],input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
(O/'island_shared_host_match_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}));print(top.head(10).to_string(index=False))
