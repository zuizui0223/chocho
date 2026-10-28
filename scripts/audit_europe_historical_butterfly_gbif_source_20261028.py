#!/usr/bin/env python3
"""Source-only 39 European candidate butterflies × three GBIF historic windows.

Never infer colonization, occupancy, climate chasing, host adoption, or managed
host presence from the count facets. No raw occurrence rows are obtained.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import json
import math
import time

SCHEMA = 'chocho_europe_39butterflies_1970_2024_gbif_source_feasibility_20261028_v01'
MANIFEST_SCHEMA = 'chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01'
MANIFEST_SHA256 = '18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2'
GBIF = 'https://api.gbif.org/v1/occurrence/search'
ORIGINAL_ROSTER = (
 'Aglais urticae', 'Anthocharis cardamines', 'Apatura ilia',
 'Aphantopus hyperantus', 'Aporia crataegi', 'Araschnia levana',
 'Argynnis aglaja', 'Argynnis paphia', 'Boloria napaea',
 'Brenthis daphne', 'Callophrys rubi', 'Carcharodus alceae',
 'Carterocephalus palaemon', 'Celastrina argiolus',
 'Coenonympha pamphilus', 'Colias croceus', 'Colias hyale',
 'Cupido minimus', 'Erebia medusa', 'Erynnis tages',
 'Gonepteryx rhamni', 'Hesperia comma', 'Hipparchia semele',
 'Iphiclides podalirius', 'Issoria lathonia', 'Lasiommata megera',
 'Leptidea sinapis', 'Lycaena dispar', 'Lycaena phlaeas',
 'Maniola jurtina', 'Nymphalis antiopa', 'Papilio machaon',
 'Pararge aegeria', 'Pieris brassicae', 'Pieris napi',
 'Pieris rapae', 'Polyommatus icarus', 'Vanessa atalanta', 'Vanessa cardui',
)
COUNTRIES = ('GB', 'DE', 'NL', 'FR', 'SE', 'ES', 'IT', 'PL')
WINDOWS = {'pre_1970_1990': '1970,1990',
           'mid_1991_2005': '1991,2005',
           'post_2010_2024': '2010,2024'}
EARLY_COUNTRY_MIN = 20
LATE_COUNTRY_MIN = 50
REQUIRED_BOTH_ERA_SPECIES = 25
REQUIRED_PAIRED_COUNTRIES = 5
MAX_BODY = 300_000
MAX_RETRIES = 3
MAX_RATE = 2.0  # starts per second, shared between all parallel workers


def sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def read_panel(path: Path) -> tuple[list[dict], str]:
    raw = path.read_bytes()
    digest = sha(raw)
    if digest != MANIFEST_SHA256:
        raise ValueError(f'Original 239 GBIF source manifest drifted: {digest}')
    parsed = json.loads(raw)
    if (parsed.get('schema') != MANIFEST_SCHEMA or
        parsed.get('species_requested') != 239 or
        len(parsed.get('species_receipts', [])) != 239):
        raise ValueError('Original 239 source roster wrong schema or sample size')
    by_name = {r['species']: r for r in parsed['species_receipts']}
    if len(by_name) != 239:
        raise ValueError('Original roster has name collisions')
    if len(ORIGINAL_ROSTER) != 39 or len(set(ORIGINAL_ROSTER)) != 39:
        raise ValueError('European cohort not frozen at 39 names')
    panel = []
    taxon_ids = []
    for name in ORIGINAL_ROSTER:
        item = by_name.get(name)
        if item is None:
            raise ValueError(f'Fixed European name not in original 239 at all: {name}')
        if item.get('status') == 'TAXON_UNRESOLVED':
            # Keep the prospectively selected taxon in the complete 39-species
            # denominator as UNSCORABLE. Never live-rematch or substitute it.
            panel.append({'species':name,'taxon_key':None,
                          'taxon_identity':None,
                          'original_retrieval_status':'TAXON_UNRESOLVED'})
            continue
        taxon = item.get('taxon') or {}
        key = taxon.get('key')
        if type(key) is not int or key <= 0:
            raise ValueError(f'Invalid original accepted taxon key for {name}')
        taxon_ids.append(key)
        panel.append({'species': name, 'taxon_key': key,
                      'taxon_identity': taxon.get('match_type'),
                      'original_retrieval_status': item.get('status')})
    if len(set(taxon_ids)) != len(taxon_ids):
        raise ValueError('Original 39 butterflies contain colliding GBIF accepted keys')
    return panel, digest


def request_url(key: int, years: str) -> str:
    if type(key) is not int or key <= 0 or years not in WINDOWS.values():
        raise ValueError('Original target key and year window not fixed')
    items = [
        ('taxonKey', str(key)), ('hasCoordinate', 'true'),
        ('hasGeospatialIssue', 'false'), ('occurrenceStatus', 'PRESENT'),
        ('year', years), ('limit', '0'), ('facet', 'country'),
        ('facet', 'basisOfRecord'), ('facetLimit', '250'),
    ]
    return GBIF + '?' + urlencode(items)


def validate_facets(data: object) -> dict:
    if not isinstance(data, dict):
        raise ValueError('GBIF returned non-object')
    n = data.get('count')
    if type(n) is not int or n < 0:
        raise ValueError('GBIF count missing/invalid')
    if data.get('limit') != 0 or data.get('results') != []:
        raise ValueError('GBIF must return zero occurrence bodies, facet only')
    facets = data.get('facets')
    if not isinstance(facets, list):
        raise ValueError('GBIF facet payload missing')
    if not facets and n == 0:
        return {'count': 0, 'countries': {}, 'bases': {}, 'missing_country': 0}
    if len(facets) != 2:
        raise ValueError('Both country and basis facets are required')
    result = {}
    for f in facets:
        if not isinstance(f, dict):
            raise ValueError('Malformed facet')
        field = str(f.get('field', '')).upper()
        if field not in ('COUNTRY', 'BASIS_OF_RECORD') or field in result:
            raise ValueError('Unexpected or duplicate facet identity')
        counts = f.get('counts')
        if not isinstance(counts, list) or len(counts) >= 250:
            raise ValueError('Facet truncated at 250 or malformed')
        decoded = {}
        for entry in counts:
            name = entry.get('name') if isinstance(entry, dict) else None
            value = entry.get('count') if isinstance(entry, dict) else None
            if not isinstance(name, str) or not name or type(value) is not int or value < 0:
                raise ValueError('Malformed source facet category')
            if name in decoded:
                raise ValueError('Duplicate source facet category')
            decoded[name] = value
        if sum(decoded.values()) > n:
            raise ValueError('Facet total exceeds source count')
        result[field] = decoded
    if set(result) != {'COUNTRY','BASIS_OF_RECORD'}:
        raise ValueError('Missing mandatory country/basis source facets')
    country = result['COUNTRY']
    if any(len(k)!=2 or k != k.upper() for k in country):
        raise ValueError('Invalid ISO country codes in GBIF country facet')
    return {'count':n, 'countries':country,
            'bases':result['BASIS_OF_RECORD'],
            'missing_country':n-sum(country.values())}


class RateLimiter:
    def __init__(self, rate: float=MAX_RATE):
        self._lock=Lock()
        self._next=0.0
        self.interval=1/rate

    def wait(self):
        while True:
            with self._lock:
                now=time.monotonic()
                rem=self._next-now
                if rem<=0:
                    self._next=now+self.interval
                    return
            time.sleep(rem)


def fetch_one(species: str, key:int, window:str, limiter:RateLimiter) -> dict:
    url = request_url(key,WINDOWS[window])
    audit={'species':species,'original_taxon_key':key,
           'window':window,'year_filter':WINDOWS[window],
           'source_query_url':url,'retrieval_status':'FAILED'}
    for attempt in range(1,MAX_RETRIES+1):
        limiter.wait()
        try:
            req=Request(url,headers={'User-Agent':'chocho-europe-historical-facet-pilot/0.1',
                                     'Accept':'application/json'})
            with urlopen(req,timeout=30) as response:
                if response.status!=200:raise ValueError('GBIF response not HTTP 200')
                raw=response.read(MAX_BODY+1)
            if len(raw)>MAX_BODY:raise ValueError('GBIF facet response exceeded payload limit')
            data=validate_facets(json.loads(raw))
            audit.update({'retrieval_status':'COMPLETE',
                          'body_sha256':sha(raw),'body_bytes':len(raw),
                          'attempts':attempt,**data})
            return audit
        except (OSError,ValueError,UnicodeError,TimeoutError) as exc:
            audit['failure_class']=type(exc).__name__
            audit['failure_description']=str(exc)[:220]
            audit['attempts']=attempt
            if attempt<MAX_RETRIES:
                time.sleep(attempt)
    return audit


def summarize(panel: list[dict], results: list[dict], *, completed:bool) -> dict:
    lookup = {(x['species'],x['window']):x for x in results
              if x['retrieval_status']=='COMPLETE'}
    eligible = []
    country_tally=Counter()
    support=defaultdict(Counter)
    for item in panel:
        sp=item['species']
        if (sp,'pre_1970_1990') not in lookup or (sp,'post_2010_2024') not in lookup:
            continue
        early=lookup[sp,'pre_1970_1990']['countries']
        later=lookup[sp,'post_2010_2024']['countries']
        mid=lookup.get((sp,'mid_1991_2005'),{}).get('countries',{})
        qualifying=[]
        for country in COUNTRIES:
            a=early.get(country,0)
            b=later.get(country,0)
            m=mid.get(country,0)
            support[country]['pre_count'] += a
            support[country]['mid_count'] += m
            support[country]['post_count'] += b
            if a>=EARLY_COUNTRY_MIN and b>=LATE_COUNTRY_MIN:
                qualifying.append(country)
                country_tally[country]+=1
        if qualifying:
            eligible.append({'species':sp,
                             'qualifying_source_countries':qualifying,
                             'qualifying_country_count':len(qualifying)})
    species_count=len(eligible)
    n_countries=sum(country_tally[c]>0 for c in COUNTRIES)
    qualified=completed and species_count>=REQUIRED_BOTH_ERA_SPECIES and n_countries>=REQUIRED_PAIRED_COUNTRIES
    basis_info={}
    for window in WINDOWS:
        entries=[x for x in results if x['retrieval_status']=='COMPLETE' and x['window']==window]
        basis=Counter()
        for x in entries:
            basis.update(x['bases'])
        reported=sum(x['count'] for x in entries)
        basis_info[window]={'completed_species_queries':len(entries),
                            'worldwide_index_records_for_fixed_panel':reported,
                            'human_observation_facet_records':basis.get('HUMAN_OBSERVATION',0),
                            'preserved_specimen_facet_records':basis.get('PRESERVED_SPECIMEN',0),
                            'other_basis_records':reported-basis.get('HUMAN_OBSERVATION',0)-basis.get('PRESERVED_SPECIMEN',0)}
    return {'historical_source_eligibility_gate':(
                'PASS_SOURCE_COUNTS_ONLY_NOT_OCCUPANCY' if qualified else
                'FAIL_SOURCE_INFORMATION_OR_INCOMPLETE'),
            'eligible_original_species_count':species_count,
            'countries_with_eligible_species':n_countries,
            'eligible_source_species':eligible,
            'paired_species_by_country':dict(sorted(country_tally.items())),
            'total_indexed_records_by_country_ERA_not_abundance':{
                 c:dict(support[c]) for c in COUNTRIES},
            'basis_source_reporting_diagnostics':basis_info}


def run(manifest:Path, output:Path, max_workers:int=3) -> dict:
    panel,digest=read_panel(manifest)
    requests=[(p['species'],p['taxon_key'],window) for p in panel
              if type(p['taxon_key']) is int for window in WINDOWS]
    limiter=RateLimiter(MAX_RATE)
    results=[]
    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        tasks={exe.submit(fetch_one,sp,key,win,limiter):(sp,win) for sp,key,win in requests}
        for future in as_completed(tasks):
            x=future.result()
            results.append(x)
            if len(results)%15==0:
                print(json.dumps({'source_queries_complete_or_failed':len(results),
                                  'source_queries_total':len(requests),
                                  'valid_completed':sum(v['retrieval_status']=='COMPLETE' for v in results)},sort_keys=True),flush=True)
    results.sort(key=lambda x:(x['species'],x['window']))
    complete=sum(x['retrieval_status']=='COMPLETE' for x in results)==len(requests)
    out={'schema':SCHEMA,'status':'COMPLETE_SOURCE_FACET_AVAILABILITY' if complete else 'PARTIAL_SOURCE_UNSCORABLE',
         'source_manifest_sha256':digest,
         'GBIF_original_239_species_count':239,
         'fixed_europe_candidate_species_count':len(panel),
         'resolved_original_taxa':sum(type(x['taxon_key']) is int for x in panel),
         'original_taxon_unresolved_names':[x['species'] for x in panel if x['taxon_key'] is None],
         'fixed_eligible_countries':list(COUNTRIES),
         'windows':WINDOWS,
         'fixed_candidate_panel':panel,
         'queries_expected':len(requests),
         'queries_complete':sum(x['retrieval_status']=='COMPLETE' for x in results),
         'individual_butterfly_occurrence_records_downloaded':0,
         'survey_effort_or_absences_estimated':False,
         'source_result_note':'GBIF country and basisOfRecord facets only; no biological range changes or host interactions derived.',
         'country_source_results':results}
    out['gate']=summarize(panel,results,completed=complete)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'queries_complete':out['queries_complete'],
                      'source_gate':out['gate']['historical_source_eligibility_gate'],
                      'eligible_species':out['gate']['eligible_original_species_count'],
                      'paired_countries':out['gate']['countries_with_eligible_species'],
                      'country_paired_species':out['gate']['paired_species_by_country']},sort_keys=True),flush=True)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-workers',type=int,default=3)
    a=p.parse_args()
    if not 1<=a.max_workers<=4:raise ValueError('At most four source API workers')
    out=run(a.manifest,a.output,a.max_workers)
    raise SystemExit(0 if out['status']=='COMPLETE_SOURCE_FACET_AVAILABILITY' else 2)
