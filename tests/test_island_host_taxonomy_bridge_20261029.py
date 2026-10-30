"""Conservative name-to-current-taxonomy eligibility, never numeric-key coercion."""
import importlib.util
from pathlib import Path
import copy
import json
import pytest


def module():
    path = Path(__file__).parents[1] / 'scripts/bridge_island_host_taxonomy_20261029.py'
    assert path.exists(), 'Taxonomy bridge implementation must exist'
    spec = importlib.util.spec_from_file_location('taxonomy_bridge', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def exact(name='Danaus plexippus', key='6CB8P'):
    return {'usage': {'key': key, 'canonicalName': name, 'rank': 'SPECIES', 'status': 'ACCEPTED'},
            'classification': [{'key': '5G9', 'name': 'Papilionoidea', 'rank': 'SUPERFAMILY'}],
            'diagnostics': {'matchType': 'EXACT', 'confidence': 99}, 'synonym': False}


def test_exact_string_key_survives_and_synonym_targets_accepted_species():
    m = module()
    assert m.assess_match('Danaus plexippus', exact())['provisional_eligible']
    response = exact()
    response['synonym'] = True
    response['usage']['status'] = 'SYNONYM'
    response['acceptedUsage'] = {'key': 'A1', 'canonicalName': 'Danaus example', 'rank': 'SPECIES'}
    result = m.assess_match('Danaus plexippus', response)
    assert result['current_accepted_key'] == 'A1'
    assert result['provisional_eligible']


def test_fuzzy_higher_rank_and_outside_clade_are_rejected():
    m = module()
    for field, value in [('matchType', 'FUZZY'), ('matchType', 'HIGHERRANK')]:
        response = exact(); response['diagnostics'][field] = value
        assert not m.assess_match('Danaus plexippus', response)['provisional_eligible']
    response = exact(); response['usage']['rank'] = 'GENUS'
    assert not m.assess_match('Danaus plexippus', response)['provisional_eligible']
    response = exact(); response['classification'] = []
    assert not m.assess_match('Danaus plexippus', response)['provisional_eligible']


def test_competing_exact_accepted_species_is_ambiguous_but_misapplied_is_flagged():
    m = module()
    response = exact()
    response['diagnostics']['alternatives'] = [exact(key='DIFFERENT')]
    assert not m.assess_match('Danaus plexippus', response)['provisional_eligible']
    response['diagnostics']['alternatives'][0]['usage']['status'] = 'MISAPPLIED'
    result = m.assess_match('Danaus plexippus', response)
    assert result['provisional_eligible']
    assert result['alternative_count'] == 1


def test_both_old_numeric_and_current_accepted_collisions_are_quarantined():
    m = module()
    rows = [{'original_name': 'A a', 'original_numeric_key': '1', 'current_accepted_key': 'X', 'provisional_eligible': True},
            {'original_name': 'B b', 'original_numeric_key': '1', 'current_accepted_key': 'Y', 'provisional_eligible': True},
            {'original_name': 'C c', 'original_numeric_key': '2', 'current_accepted_key': 'Z', 'provisional_eligible': True},
            {'original_name': 'D d', 'original_numeric_key': '3', 'current_accepted_key': 'Z', 'provisional_eligible': True}]
    result = m.quarantine_collisions(copy.deepcopy(rows))
    assert all(not r['conservative_join_eligible'] for r in result)
    assert result[0]['original_numeric_collision']
    assert result[2]['current_accepted_key_collision']


def test_expanded_overlap_refuses_partial_acquisition(tmp_path):
    # A RUNNING snapshot must not be reported as the fixed96 completed panel.
    folder = tmp_path / 'data/exploratory/island_expansion_20261029'
    folder.mkdir(parents=True)
    (folder / 'expanded_island_panel_manifest.json').write_text(json.dumps({'status': 'RUNNING'}))
    with pytest.raises(ValueError, match='must be FINISHED'):
        module().audit_expanded(tmp_path, tmp_path / 'output')
