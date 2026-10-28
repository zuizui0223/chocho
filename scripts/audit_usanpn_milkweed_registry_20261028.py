#!/usr/bin/env python3
"""Conservatively resolve USA-NPN milkweed taxon/Leaves phenophase IDs only.

This stage never inspects an observation outcome. Missing or ambiguous taxon/phase
IDs are SOURCE_UNRESOLVED rather than evidence that milkweeds have no leaves.
"""
from __future__ import annotations
import argparse, hashlib, json, re, time
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASE = 'https://services.usanpn.org/npn_portal/'
SPECIES = ('Asclepias syriaca', 'Asclepias incarnata', 'Asclepias tuberosa',
           'Asclepias speciosa', 'Asclepias fascicularis', 'Asclepias curassavica')
YEARS = tuple(range(2018, 2026))
MAX_BYTES = 8_000_000


def fetch_json(url, timeout=20):
    req = Request(url, headers={'User-Agent':'chocho-usanpn-leaf-observation/1.0',
                               'Accept':'application/json'})
    with urlopen(req, timeout=timeout) as r:
        if r.status!=200:
            raise ValueError(f'HTTP {r.status}')
        raw=r.read(MAX_BYTES+1)
        if len(raw)>MAX_BYTES:raise ValueError('bounded response exceeded 8 MB')
        if 'json' not in r.headers.get('Content-Type','').lower():
            raise ValueError('service returned non-JSON content')
    return json.loads(raw),hashlib.sha256(raw).hexdigest()


def candidates(row):
    if not isinstance(row,dict):return set()
    fields=['scientific_name','scientificName','scientificname','latin_name',
            'latinName','genus_species','name','full_name','taxon_name']
    ans=set()
    for f in fields:
        v=row.get(f)
        if isinstance(v,str):ans.add(v.strip())
    gen=row.get('genus')
    ep=row.get('species')
    if isinstance(gen,str) and isinstance(ep,str):
        ans.add(gen.strip()+' '+ep.strip())
    # NPN may use genus_name/species_name as its canonical fields.
    gen=row.get('genus_name')
    ep=row.get('species_name')
    if isinstance(gen,str) and isinstance(ep,str):
        ans.add(gen.strip()+' '+ep.strip())
    return ans


def registry_match(registry, scientific):
    if not isinstance(registry,list):
        raise ValueError('Expected species registry list')
    matches=[]
    for rec in registry:
        if not isinstance(rec,dict):continue
        if scientific in candidates(rec):
            ids=[rec.get(k) for k in ('species_id','speciesId','id') if rec.get(k) is not None]
            if len({str(i) for i in ids})!=1:
                raise ValueError('Ambiguous NPN species identifier')
            ident=int(ids[0])
            if ident<=0:raise ValueError('Invalid NPN species id')
            matches.append(ident)
    matches=sorted(set(matches))
    if len(matches)>1: raise ValueError('Multiple exact NPN accepted species records')
    return matches[0] if matches else None


def leaves_id(body, source_id):
    """Require one exact Leaves phenophase, rejecting initial-growth and colored leaves."""
    if not isinstance(body,list):raise ValueError('Expected phenophase response list')
    p=[]
    for entry in body:
        if not isinstance(entry,dict):continue
        if entry.get('species_id') not in (source_id,str(source_id)):
            continue
        obj=entry.get('phenophases')
        if not isinstance(obj,list):raise ValueError('Expected species phenophases list')
        for ph in obj:
            if not isinstance(ph,dict):continue
            if str(ph.get('phenophase_name','')).strip()=='Leaves' and str(ph.get('phenophase_category','')).strip()=='Leaves':
                p.append(ph.get('phenophase_id'))
    ids=sorted({int(x) for x in p if x is not None})
    if len(ids)>1:raise ValueError('Multiple exact Leaves phenophase IDs for species/date')
    return ids[0] if ids else None


def discover(registry, phaser):
    found={}
    for name in SPECIES:
        taxid=registry_match(registry,name)
        phase={}
        if taxid is not None:
            for year in YEARS:
                phase[str(year)]=leaves_id(phaser(taxid,year),taxid)
        found[name]={'species_id':taxid,'leaves_pheno_id_by_year':phase,
                     'registry_resolved':taxid is not None,
                     'years_leaves_protocol':sum(v is not None for v in phase.values())}
    return found


def run(path):
    errors=[]
    response_hashes={}
    try:
        registry,sha=fetch_json(BASE+'species/getSpecies.json',timeout=35)
        response_hashes['species_registry']=sha
        def get_phases(taxon,year):
            data,digest=fetch_json(BASE+'phenophases/getPhenophasesForSpecies.json?'+urlencode({
                'species_id[1]':taxon,'date':f'{year}-06-15'}))
            response_hashes[f'{taxon}:{year}']=digest
            return data
        records=discover(registry,get_phases)
        code='REGISTRY_AND_EXACT_LIVE_LEAVES_PROTOCOL_RESOLVED' if any(x['years_leaves_protocol']>0 for x in records.values()) else 'NO_ELIGIBLE_MILKWEED_LEAVES_PROTOCOL'
    except Exception as e:
        code='REGISTRY_OR_SERVICE_UNAVAILABLE'
        errors.append(f'{type(e).__name__}: {str(e)[:350]}')
        records={}
    output={'schema':'chocho_usanpn_milkweed_live_leaves_registry_v01',
            'status':code,'source':'https://services.usanpn.org/npn_portal/',
            'protocol':'USANPN_LIVE_LEAF_OBSERVATION_PROTOCOL_20261028',
            'rnpn_source':'https://github.com/usa-npn/rnpn/blob/master/R/npn_data_download.R',
            'calendar_years':list(YEARS),'source_replies_sha256':response_hashes,
            'taxa':records,'errors':errors,
            'checked_phenophase_description':'Leaves: one or more live, fully unfolded leaves',
            'no_observation_outcomes_examined':True,'executed_utc':datetime.now(timezone.utc).isoformat()}
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':code,'taxa':{x:{'species_id':v['species_id'],'years_leaves_protocol':v['years_leaves_protocol']} for x,v in records.items()},'errors':errors},sort_keys=True))
    if code!='REGISTRY_AND_EXACT_LIVE_LEAVES_PROTOCOL_RESOLVED':raise SystemExit(2)

if __name__=='__main__':
    from pathlib import Path
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    run(ap.parse_args().output)
