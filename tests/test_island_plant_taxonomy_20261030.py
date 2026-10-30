"""Scientific gate, fixed source selection and bounded public-query provenance."""
import copy
import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/bridge_island_plant_taxonomy_20261030.py'


def module():
    assert SCRIPT.exists(), 'plant taxonomy bridge implementation is missing'
    spec = importlib.util.spec_from_file_location('plant_taxonomy', SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def response(name='Acacia alba', key='AAA', status='ACCEPTED'):
    r = {'usage': {'key': key, 'canonicalName': name, 'name': name + ' Author', 'rank': 'SPECIES', 'status': status},
         'classification': [{'key': 'P', 'name': 'Plantae', 'rank': 'KINGDOM'},
                            {'key': 'F', 'name': 'Fabaceae', 'rank': 'FAMILY'}],
         'diagnostics': {'matchType': 'EXACT', 'confidence': 100, 'note': '', 'alternatives': []},
         'synonym': status == 'SYNONYM'}
    if status == 'SYNONYM':
        r['acceptedUsage'] = {'key': 'BBB', 'canonicalName': 'Acacia bella', 'rank': 'SPECIES', 'status': 'ACCEPTED'}
    return r


def test_exact_accepted_species():
    a = module().assess_match('Acacia alba', 'Fabaceae', response())
    assert a['conservative_join_eligible'] is True
    assert a['current_accepted_key'] == 'AAA'
    assert a['current_accepted_status'] == 'ACCEPTED'


def test_exact_synonym_and_independently_verified_missing_status():
    m = module()
    r = response(status='SYNONYM')
    assert m.assess_match('Acacia alba', 'Fabaceae', r)['conservative_join_eligible']
    del r['acceptedUsage']['status']
    assert 'accepted_status_not_verified' in m.assess_match('Acacia alba', 'Fabaceae', r)['rejection_reasons']
    v = response(name='Acacia bella', key='BBB')
    a = m.assess_match('Acacia alba', 'Fabaceae', r, accepted_verification=v)
    assert a['conservative_join_eligible']
    v['usage']['key'] = 'OTHER'
    assert not m.assess_match('Acacia alba', 'Fabaceae', r, accepted_verification=v)['conservative_join_eligible']


@pytest.mark.parametrize('mutation,reason', [
    ('nonexact', 'not_exact_match'), ('higher', 'not_species_rank'),
    ('kingdom', 'outside_or_unverified_Plantae'), ('family', 'source_family_not_verified'),
    ('accepted_status', 'accepted_status_not_verified'), ('numeric_key', 'accepted_key_not_nonempty_string'),
    ('canonical', 'submitted_canonical_name_not_exact'), ('ambiguity', 'diagnostics_ambiguity')])
def test_rejects_unsafe_matches(mutation, reason):
    r = response()
    if mutation == 'nonexact': r['diagnostics']['matchType'] = 'FUZZY'
    elif mutation == 'higher': r['usage']['rank'] = 'GENUS'
    elif mutation == 'kingdom': r['classification'][0]['name'] = 'Animalia'
    elif mutation == 'family': r['classification'][1]['name'] = 'Poaceae'
    elif mutation == 'accepted_status': r['usage']['status'] = 'PROVISIONALLY_ACCEPTED'
    elif mutation == 'numeric_key': r['usage']['key'] = 123
    elif mutation == 'canonical': r['usage']['canonicalName'] = 'Acacia bella'
    elif mutation == 'ambiguity': r['diagnostics']['note'] = 'Ambiguous name'
    a = module().assess_match('Acacia alba', 'Fabaceae', r)
    assert not a['conservative_join_eligible']
    assert reason in a['rejection_reasons']


def test_competing_exact_concept_including_ambiguous_synonym_is_rejected():
    r = response()
    alt = response(key='CCC', status='SYNONYM')
    alt['usage']['status'] = 'AMBIGUOUS_SYNONYM'
    r['diagnostics']['alternatives'] = [alt]
    a = module().assess_match('Acacia alba', 'Fabaceae', r)
    assert not a['conservative_join_eligible']
    assert a['competing_exact_keys'] == 'BBB'


def test_fixed_sources_keep_all_links_and_infraspecific_names_unqueried():
    m = module()
    rows, edges, pins = m.load_sources(m.ROOT)
    assert len(rows) == 1062
    assert sum(r['query_eligible'] == 'True' for r in rows) == 1040
    assert len(edges) == 1440
    assert len({r['original_name'] for r in edges}) == 93
    assert [r['wcvp_accepted_id'] for r in rows] == sorted(r['wcvp_accepted_id'] for r in rows)
    assert len(pins) == 3


def test_hash_drift_aborts_before_selection(tmp_path):
    m = module()
    p = tmp_path / next(iter(m.SOURCES))
    p.parent.mkdir(parents=True)
    p.write_text('changed')
    with pytest.raises(ValueError, match='Source hash drift'):
        m.load_sources(tmp_path)


def test_duplicate_accepted_concepts_preserve_source_rows_and_inventory():
    m = module()
    rows = [{'wcvp_accepted_id': x, 'query_name': n, 'current_accepted_key': 'AAA',
             'conservative_join_eligible': True} for x, n in [('1', 'Acacia alba'), ('2', 'Acacia bella')]]
    collisions = m.annotate_collisions(rows)
    assert len(rows) == 2
    assert all(r['conservative_join_eligible'] for r in rows)
    assert collisions == [{'current_accepted_key': 'AAA', 'source_count': 2,
                           'source_ids': '1|2', 'source_query_names': 'Acacia alba|Acacia bella'}]
    assert rows[0]['accepted_concept_source_count'] == 2


def test_request_uses_official_rank_strict_and_explicit_namespace():
    u = module().request_url('Acacia alba', 'Fabaceae')
    p = parse_qs(urlparse(u).query)
    assert p == {'checklistKey': ['7ddf754f-d193-4cc9-b351-99906754a03b'],
                 'scientificName': ['Acacia alba'], 'kingdom': ['Plantae'], 'family': ['Fabaceae'],
                 'taxonRank': ['SPECIES'], 'strict': ['true'], 'verbose': ['true']}


def test_failed_request_preserves_three_attempts_and_global_spacing(tmp_path):
    m = module()
    class Clock:
        t = 100.0
        def monotonic(self): return self.t
        def sleep(self, seconds): self.t += seconds
    clock = Clock()
    starts = []
    def fail(req, timeout):
        starts.append(clock.t)
        raise OSError('offline')
    client = m.Client(tmp_path, opener=fail, monotonic=clock.monotonic, sleep=clock.sleep)
    r, receipt = client.fetch('1', m.request_url('Acacia alba', 'Fabaceae'))
    assert r is None
    assert receipt['query_status'] == 'FAILED'
    assert len(receipt['attempts']) == 3
    assert all(b-a >= 1.05 for a,b in zip(starts, starts[1:]))
    assert (tmp_path / '1.receipt.json').exists()
    # A cached exhausted failure may not silently gain more than three attempts.
    _, again = client.fetch('1', m.request_url('Acacia alba', 'Fabaceae'))
    assert again == receipt
    assert len(starts) == 3


def test_raw_response_hash_and_cache_contract_enforced(tmp_path):
    m = module()
    class Handle:
        status = 200
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def read(self, size): return json.dumps(response()).encode()
    client = m.Client(tmp_path, opener=lambda req, timeout: Handle(), sleep=lambda _: None)
    _, receipt = client.fetch('1', m.request_url('Acacia alba', 'Fabaceae'))
    assert receipt['query_status'] == 'SUCCESS'
    (tmp_path / '1.json').write_text('{}')
    with pytest.raises(ValueError, match='Cached response hash drift'):
        client.fetch('1', m.request_url('Acacia alba', 'Fabaceae'))


def test_archive_is_deterministic_and_retains_request_receipts(tmp_path):
    m = module()
    cache = tmp_path / 'cache'; cache.mkdir()
    (cache / '1.json').write_text('{}')
    (cache / '1.receipt.json').write_text('{"url":"test"}')
    a, b = tmp_path / 'a.tar.gz', tmp_path / 'b.tar.gz'
    m.archive_cache(cache, a); m.archive_cache(cache, b)
    assert a.read_bytes() == b.read_bytes()
    import tarfile
    with tarfile.open(a) as tf: assert tf.getnames() == ['1.json', '1.receipt.json']


@pytest.mark.parametrize('accepted_rank', ['SUBSPECIES', 'VARIETY'])
def test_distinct_exact_species_homonym_to_infraspecific_concept_is_ambiguous(accepted_rank):
    # Real v2 Cardamine pratensis/Hosackia crassifolia responses expose this:
    # the alternative usage is SPECIES, although its accepted target is lower.
    r = response()
    alt = response(key='CCC', status='SYNONYM')
    alt['acceptedUsage']['rank'] = accepted_rank
    r['diagnostics']['alternatives'] = [alt]
    a = module().assess_match('Acacia alba', 'Fabaceae', r)
    assert not a['conservative_join_eligible']
    assert a['competing_exact_keys'] == 'BBB'
