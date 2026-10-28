"""Inspectable original GBIF facet-source receipt invariants; no inferred colonization."""
from pathlib import Path
from collections import Counter
import json,hashlib

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'provenance/archive/host_hindcast/results/europe_39_historical_source_facets_20261028.json'
SHA='3f0a819752612211158dc9f0031181254ae35f14fa67fdbfd993333b7e0ba467'


def test_complete_frozen_counts_only():
    b=RESULT.read_bytes()
    assert hashlib.sha256(b).hexdigest()==SHA
    r=json.loads(b)
    assert r['status']=='COMPLETE_SOURCE_FACET_AVAILABILITY'
    assert r['queries_expected']==r['queries_complete']==114
    assert r['fixed_europe_candidate_species_count']==39
    assert r['resolved_original_taxa']==38
    assert r['original_taxon_unresolved_names']==['Lasiommata megera']
    assert r['individual_butterfly_occurrence_records_downloaded']==0
    assert r['survey_effort_or_absences_estimated'] is False
    q=r['gate']
    assert q['historical_source_eligibility_gate']=='PASS_SOURCE_COUNTS_ONLY_NOT_OCCUPANCY'
    assert q['eligible_original_species_count']==37
    assert q['countries_with_eligible_species']==8
    assert q['paired_species_by_country']=={'DE':32,'ES':33,'FR':37,'GB':28,'IT':18,'NL':24,'PL':18,'SE':29}


def test_never_silently_interpret_source_count_zeros_as_occupancy():
    r=json.loads(RESULT.read_text())
    original_39={p['species'] for p in r['fixed_candidate_panel']}
    original_38={p['species'] for p in r['fixed_candidate_panel'] if p['taxon_key'] is not None}
    assert len(original_39)==39 and len(original_38)==38
    pairs=set()
    for x in r['country_source_results']:
        assert x['retrieval_status']=='COMPLETE'
        assert x['species'] in original_38
        assert x['count']>=sum(x['countries'].values())
        assert x['count']>=sum(x['bases'].values())
        assert len(x['body_sha256'])==64
        assert 'limit=0' in x['source_query_url']
        pair=(x['species'],x['window'])
        assert pair not in pairs
        pairs.add(pair)
    assert len(pairs)==114
    source=RESULT.read_text()
    for forbidden in ('"decimalLatitude":','"decimalLongitude":','"occurrenceId":',
                      '"gbifID":','"observation_date":'):
        assert forbidden not in source


def test_source_record_reporting_not_butterfly_demography():
    r=json.loads(RESULT.read_text())
    c=r['gate']['total_indexed_records_by_country_ERA_not_abundance']
    pre=sum(v['pre_count'] for v in c.values())
    post=sum(v['post_count'] for v in c.values())
    assert pre==1230875 and post==20656024
    assert c['GB']['pre_count']==1130062 and c['GB']['post_count']==11635336
    assert .90 < c['GB']['pre_count']/pre < .93
    assert .55 < c['GB']['post_count']/post < .58
