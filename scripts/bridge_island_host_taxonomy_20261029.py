#!/usr/bin/env python3
"""Conservative original-239 name -> current COL XR taxonomy bridge.

No occurrences are queried. No regional host status is assigned to islands.
Raw official responses, retrieval times and request contracts are archived.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
CHECKLIST = '7ddf754f-d193-4cc9-b351-99906754a03b'
SUPERFAMILY = '5G9'
ENDPOINT = 'https://api.gbif.org/v2/species/match'
DOCS = 'https://techdocs.gbif.org/openapi/checklistbank.json'
SOURCE = 'data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json'
SOURCE_SHA = '18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2'
PILOT = 'data/exploratory/island_port_20261029/island_butterfly_pilot_species_key_facets.csv'
LAST_START = 0.0


def sha(body):
    return hashlib.sha256(body).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def assess_match(name, response):
    usage = response.get('usage') or {}
    accepted = response.get('acceptedUsage') or usage
    diagnostics = response.get('diagnostics') or {}
    synonym = response.get('synonym') is True
    reasons = []
    canonical = ' '.join(str(usage.get('canonicalName', '')).split())
    if diagnostics.get('matchType') != 'EXACT': reasons.append('not_exact_match')
    if canonical != ' '.join(name.split()): reasons.append('submitted_canonical_name_not_exact')
    if usage.get('rank') != 'SPECIES' or accepted.get('rank') != 'SPECIES': reasons.append('not_species_rank')
    status = usage.get('status', '')
    if status not in ('ACCEPTED', 'SYNONYM'): reasons.append('unsupported_usage_status')
    if synonym and not response.get('acceptedUsage'): reasons.append('synonym_without_accepted_usage')
    if (status == 'SYNONYM') != synonym: reasons.append('inconsistent_synonym_status')
    if not isinstance(accepted.get('key'), str) or not accepted.get('key'): reasons.append('accepted_key_not_nonempty_string')
    classification = response.get('classification') or []
    clade = any(str(x.get('key')) == SUPERFAMILY and x.get('rank') == 'SUPERFAMILY'
                and x.get('name') == 'Papilionoidea' for x in classification)
    if not clade: reasons.append('outside_or_unverified_Papilionoidea')
    if 'ambig' in str(diagnostics.get('note', '')).lower(): reasons.append('diagnostics_ambiguity')
    competing = set()
    alternatives = diagnostics.get('alternatives') or []
    for alt in alternatives:
        au = alt.get('usage') or {}
        aa = alt.get('acceptedUsage') or au
        ad = alt.get('diagnostics') or {}
        if (ad.get('matchType') == 'EXACT' and au.get('canonicalName') == canonical
                and au.get('rank') == 'SPECIES'
                and au.get('status') in ('ACCEPTED', 'SYNONYM', 'PROVISIONALLY_ACCEPTED')
                and aa.get('rank') == 'SPECIES' and aa.get('key') != accepted.get('key')):
            competing.add(str(aa.get('key', '')))
    if competing: reasons.append('competing_exact_current_species_concepts')
    return {'current_usage_key': str(usage.get('key', '')),
            'current_usage_canonical_name': canonical,
            'current_usage_full_name': usage.get('name', ''),
            'current_usage_rank': usage.get('rank', ''), 'current_usage_status': status,
            'current_accepted_key': str(accepted.get('key', '')),
            'current_accepted_canonical_name': accepted.get('canonicalName', ''),
            'current_accepted_full_name': accepted.get('name', ''),
            'current_accepted_rank': accepted.get('rank', ''),
            'synonym': synonym, 'match_type': diagnostics.get('matchType', ''),
            'confidence': diagnostics.get('confidence'), 'diagnostic_note': diagnostics.get('note', ''),
            'alternative_count': len(alternatives), 'competing_exact_keys': '|'.join(sorted(competing)),
            'Papilionoidea_verified': clade, 'rejection_reasons': '|'.join(reasons),
            'provisional_eligible': not reasons,
            'classification_json': json.dumps(classification, separators=(',', ':'))}


def quarantine_collisions(rows):
    old, new = defaultdict(list), defaultdict(list)
    for row in rows:
        if row['original_numeric_key']: old[row['original_numeric_key']].append(row['original_name'])
        if row['current_accepted_key']: new[row['current_accepted_key']].append(row['original_name'])
    for row in rows:
        oc = old.get(row['original_numeric_key'], [])
        nc = new.get(row['current_accepted_key'], [])
        row['original_numeric_collision'] = len(oc) > 1
        row['current_accepted_key_collision'] = len(nc) > 1
        row['original_numeric_collision_names'] = '|'.join(oc) if len(oc) > 1 else ''
        row['current_accepted_collision_names'] = '|'.join(nc) if len(nc) > 1 else ''
        row['conservative_join_eligible'] = bool(row['provisional_eligible'] and len(oc) <= 1 and len(nc) <= 1)
        row['old_numeric_occurrence_join_eligible'] = bool(row['conservative_join_eligible'] and row['original_numeric_key'])
    return rows


def cached_get(cache, label, url):
    global LAST_START
    body_path, receipt_path = cache / (label + '.json'), cache / (label + '.receipt.json')
    if body_path.exists() and receipt_path.exists():
        body = body_path.read_bytes(); receipt = json.loads(receipt_path.read_text())
        if receipt['url'] != url or receipt['body_sha256'] != sha(body):
            raise ValueError('Cached request contract or body changed')
        return json.loads(body), receipt
    attempts = []
    for attempt in range(3):
        time.sleep(max(0, 1.05 - (time.monotonic() - LAST_START)))
        LAST_START = time.monotonic()
        started = utc()
        try:
            with urlopen(Request(url, headers={'User-Agent': 'chocho-taxonomy-bridge/1.0',
                                               'Accept': 'application/json'}), timeout=45) as handle:
                body = handle.read(8_000_001)
                if len(body) > 8_000_000: raise ValueError('Response size exceeds bound')
                status = handle.status
            parsed = json.loads(body)
            receipt = {'url': url, 'started_utc': started, 'finished_utc': utc(),
                       'http_status': status, 'body_sha256': sha(body), 'body_bytes': len(body),
                       'previous_attempts': attempts}
            body_path.write_bytes(body)
            receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
            return parsed, receipt
        except Exception as exc:
            attempts.append({'started_utc': started, 'error_type': type(exc).__name__, 'error': str(exc)[:250]})
            if attempt == 2: raise
            time.sleep(2 ** (attempt + 1))


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def audit_existing_host_evidence(root, output):
    """Expose exact existing evidence grains without inventing island joins."""
    base = 'provenance/reviewer_defenses/results/'
    paths = {
        'regional_candidates': base + 'gbif_country_20261020/gbif_all26_source_host_country_inventory_20261023.json',
        'plant_gadm': base + 'gbif_country_20261020/gbif_all26_host_gadm_source_20261023.json',
        'point_proximity': base + 'butterfly_host_point_proximity_summary_20261026.json',
        'juvenile_annotations': base + 'gbif_country_20261020/gbif_juvenile_host_annotation_coverage_20261023.json'}
    evidence = {k: json.loads((root / p).read_text()) for k,p in paths.items()}
    inventory, botanical, points, juvenile = (evidence[k] for k in paths)
    examples = []
    for cell in inventory['all_48_positive_country_cells']:
        if (cell['butterfly'], cell['country']) in (('Aguna asander', 'CU'), ('Anteos clorinde', 'CU')):
            examples.append({'butterfly': cell['butterfly'], 'country': cell['country'],
                             'source_candidate_plants': cell['candidate_accepted_plants']})
    result = {
        'schema': 'chocho_original_host_evidence_grain_audit_20261029_v1',
        'sources': {k: {'path': p, 'sha256': sha((root / p).read_bytes()),
                        'schema': evidence[k]['schema'], 'status': evidence[k]['status']}
                    for k,p in paths.items()},
        'regional_source': {'added_positive_butterfly_country_cells': inventory['frozen_added_host_cells'],
            'candidate_butterfly_country_plant_links': inventory['candidate_butterfly_country_plant_links'],
            'distinct_WCVP_plant_IDs': inventory['distinct_source_plant_ids'],
            'geographic_grain': 'WGSRPD3 introduced host-plant ranges and derived country flags, not individual island plant presence',
            'specific_WCVP_ID_examples': examples},
        'existing_botanical_observation_evidence': {
            'indexed_plant_country_records': botanical['indexed_plant_source_records'],
            'plant_country_pairs': botanical['plant_country_source_pairs'],
            'cases_with_GADM1_butterfly_plant_overlap': botanical['cases_with_butterfly_GADM1_overlap_with_indexed_host'],
            'geographic_grain': 'GADM1/country facets for the original26 selected positive country cases; no keyed link to current24/96 GSHHG islands',
            'source_original_zip_sha256': botanical['source']['source_original_zip_sha256']},
        'existing_point_proximity_evidence': {'workflow_run': points['run'],
            'derived_artifact': points['original_derived_artifact'], 'cohort': points['cohort'],
            'within5km_3years_butterfly_record_indices': points['summary']['near5'],
            'strict_uncertainty_within5km_3years_indices': points['summary']['strict5'],
            'validated_host_use_literal_annotations': points['summary']['validated_larval_host_use_literal_annotated_indices'],
            'grain': 'Derived point-proximity aggregates in original26 country cases. Raw point/record identities are not present in this public summary; it does not resolve current island membership.'},
        'existing_juvenile_annotation_evidence': {
            'records_checked': juvenile['original_juvenile_stage_records'],
            'records_with_exact_candidate_plant_binomial': juvenile['records_with_exact_accepted_candidate_binomial_in_text'],
            'records_with_candidate_binomial_and_host_action_word': juvenile['records_with_candidate_plant_name_and_host_word_not_verified'],
            'interpretation': 'Missing host annotations are missing evidence, not no feeding or no recruitment.'},
        'conclusion': 'Taxonomic compatibility can connect original butterfly labels to current recording keys. Existing regional, GADM1 and point-proximity evidence must not be assigned to sampled islands without original coordinates/geometry, timing and plant evidence checks. Current taxonomy matches do not establish actual local host availability or host use.',
        'coordinate_availability_audit': {
            'point_analysis_script': 'scripts/analyze_butterfly_host_point_proximity_20261026.py',
            'point_analysis_script_sha256': sha((root / 'scripts/analyze_butterfly_host_point_proximity_20261026.py').read_bytes()),
            'source_script_policy': 'Unredacted coordinates stay in GitHub Actions memory ONLY; exported acquisition no_raw_coordinates_preserved=true.',
            'checked_local_scope': ['repository data/ and provenance/', 'workspace data-file extensions zip/parquet/feather/jsonl/ndjson',
                '/tmp/chocho-island-gbif-cache', '/tmp/chocho-island-expansion-gbif-cache', '/tmp/chocho-island-sources'],
            'observed': 'No stored decimalLatitude/decimalLongitude fields or point-level host/butterfly input files found in the checked data/provenance/task caches. Current island occurrence responses have results=[] and facets only. GSHHG/WPI files contain geography/ports, not botanical occurrence points.',
            'direct_existing_point_to96_island_join_currently_executable': False,
            'next_required_input': 'Original coordinate-level botanical/butterfly records, or a separately authorized fresh point acquisition using the archived exact taxon-country source queries. Country/GADM facets alone cannot recover point membership.'},
        'current_panel_island_host_availability_join_validated': False}
    (output / 'existing_host_evidence_grain_audit.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result


def audit_expanded(root, output):
    """Reuse cached-name bridge only after fixed96 occurrence acquisition ends."""
    base = 'data/exploratory/island_expansion_20261029/'
    manifest_path = root / (base + 'expanded_island_panel_manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if manifest['status'] != 'FINISHED':
        raise ValueError('Expanded96 panel must be FINISHED before bridge overlap audit')
    if manifest['taxonomy']['checklistKey'] != CHECKLIST or manifest['taxonomy']['taxonKey'] != SUPERFAMILY:
        raise ValueError('Expanded96 taxonomy namespace differs')
    if manifest['summary']['requested_islands'] != 96:
        raise ValueError('Expanded source is not the fixed96 panel')
    bridge_path = output / 'original239_to_colxr_bridge.csv'
    bridge_receipt = json.loads((output / 'taxonomy_bridge_manifest.json').read_text())
    if sha(bridge_path.read_bytes()) != bridge_receipt['output_sha256']['original239_to_colxr_bridge.csv']:
        raise ValueError('Cached bridge CSV differs from source-verified bridge receipt')
    with bridge_path.open() as handle: bridge = list(csv.DictReader(handle))
    if len(bridge) != 239: raise ValueError('Cached bridge does not retain original239')
    facets_path = root / (base + 'expanded_island_panel_species_key_facets.csv')
    with facets_path.open() as handle: facets = list(csv.DictReader(handle))
    by_key = defaultdict(list)
    for row in bridge:
        if row['current_accepted_key']: by_key[row['current_accepted_key']].append(row)
    overlap, by_window = [], defaultdict(set)
    eligible_names, numeric_names = set(), set()
    for facet in facets:
        candidates = by_key.get(facet['facet_key'], [])
        eligible = [r['original_name'] for r in candidates if r['conservative_join_eligible'] == 'True']
        old_numeric = [r['original_name'] for r in candidates if r['old_numeric_occurrence_join_eligible'] == 'True']
        eligible_names.update(eligible); numeric_names.update(old_numeric)
        if eligible: by_window[facet['island_key'], facet['window']].add(facet['facet_key'])
        overlap.append({**facet, 'original_panel_candidate_names': '|'.join(r['original_name'] for r in candidates),
                        'conservative_eligible_original_names': '|'.join(eligible),
                        'old_numeric_occurrence_join_eligible_names': '|'.join(old_numeric),
                        'conservative_eligible': bool(eligible), 'old_numeric_occurrence_join_eligible': bool(old_numeric)})
    islands = {r['island_key'] for r in manifest['results']}
    total = sum(int(r['record_count']) for r in facets)
    assigned = sum(int(r['record_count']) for r in overlap if r['conservative_eligible'])
    keys = {r['facet_key'] for r in overlap if r['conservative_eligible']}
    csv_path = output / 'expanded96_taxonomy_overlap.csv'
    write_csv(csv_path, overlap)
    result = {'schema': 'chocho_fixed96_cached_taxonomy_bridge_overlap_20261029_v1',
        'status': 'FINISHED_OCCURRENCE_PANEL_TAXONOMIC_OVERLAP_NOT_ISLAND_HOST_EVIDENCE',
        'source_manifest': {'path': str(manifest_path.relative_to(root)), 'sha256': sha(manifest_path.read_bytes())},
        'source_species_facets': {'path': str(facets_path.relative_to(root)), 'sha256': sha(facets_path.read_bytes())},
        'source_bridge_sha256': sha(bridge_path.read_bytes()),
        'occurrence_completeness': manifest['summary'],
        'summary': {'eligible_original_host_panel_names_recorded': len(eligible_names),
            'eligible_current_species_keys_recorded': len(keys), 'eligible_original_names': sorted(eligible_names),
            'old_numeric_occurrence_join_eligible_names_recorded': len(numeric_names),
            'originally_unresolved_names_now_mapped_and_recorded': sorted(eligible_names - numeric_names),
            'all_current_species_keys_recorded': len({r['facet_key'] for r in facets}),
            'species_assigned_record_denominator': total, 'eligible_original_panel_record_numerator': assigned,
            'eligible_original_panel_fraction_of_species_assigned_records': assigned / total if total else None,
            'islands_with_any_eligible_original_panel_record': len({r['island_key'] for r in overlap if r['conservative_eligible']}),
            'islands_atleast3_eligible_original_panel_species_each_window': sum(
                all(len(by_window[i,w]) >= 3 for w in ('2010_2017', '2018_2025')) for i in islands)},
        'caveats': ['Distinct from the original24 pilot overlap; expanded96 includes reused pilot windows and is not an independent validation cohort.',
            'Eligible taxa are original HOSTS/WCVP host-panel labels, not measured contemporary hosts on an island.',
            'Species-assigned record fraction excludes records without SPECIES_KEY. Zero indexed records are not biological absence.',
            'Old numeric eligibility remains false for original unresolved keys and either namespace collision. No original counts were merged into current counts.'],
        'output_csv_sha256': sha(csv_path.read_bytes())}
    (output / 'expanded96_taxonomy_overlap_manifest.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result['summary'], sort_keys=True), flush=True)
    return result


def run(root, cache, output):
    cache.mkdir(parents=True, exist_ok=True); output.mkdir(parents=True, exist_ok=True)
    source_bytes = (root / SOURCE).read_bytes()
    if sha(source_bytes) != SOURCE_SHA: raise ValueError('Original239 manifest SHA changed')
    originals = json.loads(source_bytes)['species_receipts']
    if len(originals) != 239 or len({r['species'] for r in originals}) != 239:
        raise ValueError('Original239 name panel changed')
    docs, docs_receipt = cached_get(cache, 'official_species_api_openapi', DOCS)
    contract = docs['paths']['/v2/species/match']['get']
    parameters = {r['name'] for r in contract['parameters']}
    if not {'scientificName', 'checklistKey', 'strict', 'verbose', 'taxonRank'} <= parameters:
        raise ValueError('Official API contract changed')
    metadata, metadata_receipt = cached_get(cache, 'colxr_matching_service_metadata',
                                           ENDPOINT + '/metadata?' + urlencode({'checklistKey': CHECKLIST}))
    if metadata.get('mainIndex', {}).get('datasetKey') != CHECKLIST:
        raise ValueError('Matching service metadata does not identify requested COL XR checklist')
    rows, receipts = [], []
    for index, original in enumerate(originals):
        name = original['species']
        params = {'scientificName': name, 'checklistKey': CHECKLIST,
                  'taxonRank': 'SPECIES', 'strict': 'true', 'verbose': 'true'}
        url = ENDPOINT + '?' + urlencode(params)
        label = 'match_' + sha(name.encode())[:20]
        response, receipt = cached_get(cache, label, url)
        assessed = assess_match(name, response)
        taxon = original.get('taxon') or {}
        rows.append({'original_name': name, 'original_source_status': original['status'],
                     'original_numeric_key': str(taxon.get('key') or ''),
                     'original_canonical_name': taxon.get('canonical', ''),
                     'original_taxon_status': taxon.get('taxon_status', ''),
                     **assessed, 'response_sha256': receipt['body_sha256'],
                     'retrieved_utc': receipt['finished_utc'], 'archive_raw_member': label + '.json'})
        receipts.append({'original_name': name, 'archive_raw_member': label + '.json', **receipt})
        if (index + 1) % 25 == 0: print(json.dumps({'taxonomy_names_completed': index + 1}), flush=True)
    quarantine_collisions(rows)
    expected_old_collision = {r['original_name'] for r in rows if r['original_numeric_collision']}
    if expected_old_collision != {'Oeneis bore', 'Oeneis polixenes'}:
        raise ValueError('Original numeric collision audit changed')
    pilot_path = root / PILOT
    with pilot_path.open() as handle: pilot = list(csv.DictReader(handle))
    pilot_manifest_path = root / 'data/exploratory/island_port_20261029/island_butterfly_pilot_manifest.json'
    pilot_manifest = json.loads(pilot_manifest_path.read_text())
    if pilot_manifest['taxonomy']['checklistKey'] != CHECKLIST or pilot_manifest['taxonomy']['taxonKey'] != SUPERFAMILY:
        raise ValueError('Pilot taxonomy namespace differs')
    matched = defaultdict(list)
    for row in rows:
        if row['current_accepted_key']: matched[row['current_accepted_key']].append(row)
    overlap = []
    for facet in pilot:
        matches = matched.get(facet['facet_key'], [])
        eligible = [r['original_name'] for r in matches if r['conservative_join_eligible']]
        overlap.append({**facet, 'original_panel_candidate_names': '|'.join(r['original_name'] for r in matches),
                        'conservative_eligible_original_names': '|'.join(eligible),
                        'conservative_eligible': bool(eligible)})
    eligible_keys = {r['current_accepted_key'] for r in rows if r['conservative_join_eligible']}
    pilot_keys = {r['facet_key'] for r in pilot}
    detected = sorted({name for r in overlap for name in r['conservative_eligible_original_names'].split('|') if name})
    species_by_window = defaultdict(set)
    for row in overlap:
        if row['conservative_eligible']:
            species_by_window[row['island_key'], row['window']].add(row['facet_key'])
    pilot_islands = {r['island_key'] for r in pilot_manifest['results']}
    summary = {'original_names': len(rows), 'match_type_counts': dict(Counter(r['match_type'] for r in rows)),
               'provisionally_eligible_names': sum(r['provisional_eligible'] for r in rows),
               'conservative_join_eligible_names': sum(r['conservative_join_eligible'] for r in rows),
               'exact_synonym_names': sum(r['match_type'] == 'EXACT' and r['synonym'] for r in rows),
               'original_numeric_collision_names': sorted(expected_old_collision),
               'current_accepted_collision_groups': {k: [r['original_name'] for r in v] for k,v in matched.items() if len(v)>1},
               'rejection_reason_counts': dict(Counter(reason for r in rows for reason in r['rejection_reasons'].split('|') if reason)),
               'pilot_unique_species_keys': len(pilot_keys), 'pilot_total_species_bin_records': sum(int(r['record_count']) for r in pilot),
               'pilot_overlap_unique_conservative_keys': len(eligible_keys & pilot_keys),
               'pilot_overlap_original_names': detected,
               'pilot_overlap_islands': len({r['island_key'] for r in overlap if r['conservative_eligible']}),
               'pilot_overlap_species_island_window_bins': sum(bool(r['conservative_eligible']) for r in overlap),
               'pilot_islands_atleast3_eligible_panel_species_each_window': sum(
                   all(len(species_by_window[i, w]) >= 3 for w in ('2010_2017', '2018_2025'))
                   for i in pilot_islands),
               'pilot_noneligible_primary_candidate_names': sorted({r['original_name'] for r in rows
                   if not r['conservative_join_eligible'] and r['current_accepted_key'] in pilot_keys}),
               'pilot_overlap_species_bin_records': sum(int(r['record_count']) for r in overlap if r['conservative_eligible']),
               'actual_island_host_availability_rows_supplied_by_taxonomy_bridge': 0}
    write_csv(output / 'original239_to_colxr_bridge.csv', rows)
    write_csv(output / 'original24_pilot_taxonomy_overlap.csv', overlap)
    # Normalize archive metadata so a fully cached replay is byte-reproducible.
    with (output / 'official_taxonomy_raw_responses.tar.gz').open('wb') as archive_file:
        with gzip.GzipFile(fileobj=archive_file, mode='wb', filename='', mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w') as archive:
                labels = ['official_species_api_openapi', 'colxr_matching_service_metadata']
                labels.extend(Path(r['archive_raw_member']).stem for r in receipts)
                paths = [cache / (label + suffix) for label in labels
                         for suffix in ('.json', '.receipt.json')]
                for path in sorted(paths):
                    raw = path.read_bytes(); info = tarfile.TarInfo(path.name)
                    info.size = len(raw); info.mtime = 0; info.mode = 0o644
                    archive.addfile(info, io.BytesIO(raw))
    report = {'schema': 'chocho_original239_current_COLXR_conservative_taxonomy_bridge_20261029_v1',
        'status': 'TAXONOMY_BRIDGE_ACQUIRED_HOST_ISLAND_AVAILABILITY_UNVALIDATED',
        'original_manifest': {'path': SOURCE, 'sha256': SOURCE_SHA},
        'pilot_source': {'path': PILOT, 'sha256': sha(pilot_path.read_bytes()),
                         'manifest_sha256': sha(pilot_manifest_path.read_bytes()),
                         'recording_completeness_summary': pilot_manifest['summary']},
        'taxonomy_namespace': {'checklistKey': CHECKLIST, 'superfamilyKey': SUPERFAMILY, 'keys_are_strings': True},
        'official_api_contract': {'docs_url': DOCS, 'docs_sha256': docs_receipt['body_sha256'],
                                  'matching_metadata': metadata, 'metadata_receipt': metadata_receipt,
                                  'query_parameters': {'scientificName': '<original panel name>', 'checklistKey': CHECKLIST,
                                      'taxonRank': 'SPECIES', 'strict': 'true', 'verbose': 'true'}},
        'rules': ['Require EXACT diagnostics and unchanged submitted canonical name, SPECIES usage and accepted target, ACCEPTED or SYNONYM status, verified Papilionoidea5G9 classification.',
                  'Reject diagnostics ambiguity or competing exact ACCEPTED/SYNONYM/PROVISIONALLY_ACCEPTED concepts; MISAPPLIED alternatives retained in raw diagnostics but are not accepted biological concepts.',
                  'Quarantine every original numeric or current accepted-key collision. No numeric-to-string conversion is treated as a namespace bridge.',
                  'A newly resolved original name without an original numeric key remains ineligible for old numeric occurrence joins even when its current name mapping qualifies.',
                  'Conservative eligibility allows taxonomic record joins only; it is not proof that biological host evidence transfers between every synonym or concept.'],
        'host_evidence_audit': {'meaning': 'Original HOSTS larval associations plus WCVP native/introduced plant ranges reconstruct potential WGSRPD3 regional resource opportunity, not actual host plants on a sampled island.',
            'available_original_inputs': 'Public original239 name/GBIF audit manifest and regional/country derived summaries. Author-only HOSTS/WCVP/TDWG/CHELSA source root is not supplied by this bridge.',
            'island_grain_host_status': 'UNVALIDATED: no independent island plant presence, contemporaneous introduction date, larval feeding or recruitment evidence is created by taxonomy matching.',
            'no_country_host_status_assigned_to_islands': True},
        'summary': summary, 'query_receipts': receipts,
        'output_sha256': {n: sha((output/n).read_bytes()) for n in ('original239_to_colxr_bridge.csv', 'original24_pilot_taxonomy_overlap.csv', 'official_taxonomy_raw_responses.tar.gz')}}
    (output / 'taxonomy_bridge_manifest.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    audit_existing_host_evidence(root, output)
    print(json.dumps(summary, sort_keys=True), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--cache-dir', type=Path, default=Path('/tmp/chocho-island-taxonomy-cache'))
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'data/exploratory/island_taxonomy_20261029')
    parser.add_argument('--expanded-overlap-only', action='store_true', help='Read completed fixed96 panel and existing bridge; no network calls')
    args = parser.parse_args()
    if args.expanded_overlap_only:
        audit_expanded(args.root, args.output_dir)
        audit_existing_host_evidence(args.root, args.output_dir)
    else:
        run(args.root, args.cache_dir, args.output_dir)
