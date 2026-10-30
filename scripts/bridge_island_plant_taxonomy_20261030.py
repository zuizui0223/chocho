#!/usr/bin/env python3
"""Resolve the fixed early93 host catalogue conservatively in current COL XR.

Only public taxonomy JSON is acquired. Source links remain intact when several
WCVP IDs resolve to one concept; occurrence joins must deduplicate per butterfly.
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
from urllib.error import HTTPError

ROOT = Path(__file__).resolve().parents[1]
CHECKLIST = '7ddf754f-d193-4cc9-b351-99906754a03b'
ENDPOINT = 'https://api.gbif.org/v2/species/match'
OUTPUT = 'data/exploratory/island_dated_host_20261030'
BASE = 'data/exploratory/island_host_roster_20261030/'
SOURCES = {
    BASE + 'plant_taxon_query_candidates.csv': 'a1cc87b20fbd19415fa81ca4269b8414b0c32167bb6033ca0425354c89c1906f',
    BASE + 'original230_host_edges.csv': 'e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22',
    BASE + 'baseline93_host_coverage.csv': 'f702ce482ca6b6411bb4530aa8b829400300cb5ef2464e367eef4e5a365ce9aa',
}
MAX_BYTES = 8_000_000


def sha(body):
    return hashlib.sha256(body).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def load_sources(root):
    tables, pins = [], {}
    for path, expected in SOURCES.items():
        body = (root / path).read_bytes()
        if sha(body) != expected:
            raise ValueError('Source hash drift: ' + path)
        pins[path] = expected
        tables.append(list(csv.DictReader(io.StringIO(body.decode('utf-8')))))
    plants, all_edges, early = tables
    names = {r['original_name'] for r in early}
    edges = [r for r in all_edges if r['original_name'] in names]
    ids = {r['wcvp_accepted_id'] for r in edges}
    rows = sorted((r for r in plants if r['wcvp_accepted_id'] in ids), key=lambda r: r['wcvp_accepted_id'])
    if (len(names), len(edges), len(rows), sum(r['query_eligible'] == 'True' for r in rows)) != (93, 1440, 1062, 1040):
        raise ValueError('Fixed early93 catalogue counts changed')
    if len({r['wcvp_accepted_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate source plant IDs')
    return rows, edges, pins


def request_url(name, family):
    return ENDPOINT + '?' + urlencode({'checklistKey': CHECKLIST, 'scientificName': name,
        'kingdom': 'Plantae', 'family': family, 'taxonRank': 'SPECIES',
        'strict': 'true', 'verbose': 'true'})


def canonical(value):
    return ' '.join(str(value or '').split())


def assess_match(name, family, response, accepted_verification=None):
    usage = response.get('usage') or {}
    accepted = response.get('acceptedUsage') or usage
    diagnostics = response.get('diagnostics') or {}
    classification = response.get('classification') or []
    alternatives = diagnostics.get('alternatives') or []
    synonym = response.get('synonym') is True
    reasons = []
    if diagnostics.get('matchType') != 'EXACT': reasons.append('not_exact_match')
    if canonical(usage.get('canonicalName')) != canonical(name): reasons.append('submitted_canonical_name_not_exact')
    if usage.get('rank') != 'SPECIES' or accepted.get('rank') != 'SPECIES': reasons.append('not_species_rank')
    status = usage.get('status', '')
    if status not in ('ACCEPTED', 'SYNONYM'): reasons.append('unsupported_usage_status')
    if (status == 'SYNONYM') != synonym: reasons.append('inconsistent_synonym_status')
    if status == 'SYNONYM' and not response.get('acceptedUsage'): reasons.append('synonym_without_accepted_usage')
    key = accepted.get('key', '')
    if not isinstance(key, str) or not key.strip(): reasons.append('accepted_key_not_nonempty_string')
    kingdom = any(x.get('rank') == 'KINGDOM' and x.get('name') == 'Plantae' and x.get('key') == 'P' for x in classification)
    verified_family = bool(family) and any(x.get('rank') == 'FAMILY' and x.get('name') == family for x in classification)
    if not kingdom: reasons.append('outside_or_unverified_Plantae')
    if not verified_family: reasons.append('source_family_not_verified')
    if 'ambig' in json.dumps({k: v for k,v in diagnostics.items() if k != 'alternatives'}).lower():
        reasons.append('diagnostics_ambiguity')
    competing = set()
    for alt in alternatives:
        au = alt.get('usage') or {}
        aa = alt.get('acceptedUsage') or au
        ad = alt.get('diagnostics') or {}
        # Every distinct exact species concept competes, including ambiguous
        # synonym alternatives. Their low ranking is not disambiguating evidence.
        if (ad.get('matchType') == 'EXACT' and canonical(au.get('canonicalName')) == canonical(name)
                and au.get('rank') == 'SPECIES'
                and aa.get('key') and aa.get('key') != key):
            competing.add(str(aa['key']))
    if competing: reasons.append('competing_exact_current_species_concepts')
    accepted_status = accepted.get('status', '')
    accepted_status_source = 'primary_response' if accepted_status else 'unverified'
    if accepted_status != 'ACCEPTED' and not accepted_status and accepted_verification:
        checked = assess_match(canonical(accepted.get('canonicalName')), family, accepted_verification)
        if (checked['conservative_join_eligible'] and checked['current_usage_status'] == 'ACCEPTED'
                and checked['current_usage_key'] == key and checked['current_accepted_key'] == key):
            accepted_status = 'ACCEPTED'
            accepted_status_source = 'independent_accepted_name_match'
    if accepted_status != 'ACCEPTED': reasons.append('accepted_status_not_verified')
    return {'current_usage_key': str(usage.get('key', '')),
        'current_usage_name': canonical(usage.get('canonicalName')),
        'current_usage_full_name': usage.get('name', ''), 'current_usage_rank': usage.get('rank', ''),
        'current_usage_status': status, 'current_accepted_key': str(key),
        'current_accepted_name': canonical(accepted.get('canonicalName')),
        'current_accepted_full_name': accepted.get('name', ''), 'current_accepted_rank': accepted.get('rank', ''),
        'current_accepted_status': accepted_status, 'accepted_status_evidence': accepted_status_source,
        'synonym': synonym, 'match_type': diagnostics.get('matchType', ''),
        'confidence': diagnostics.get('confidence', ''), 'diagnostic_note': diagnostics.get('note', ''),
        'alternative_count': len(alternatives), 'competing_exact_keys': '|'.join(sorted(competing)),
        'Plantae_verified': kingdom, 'source_family_verified': verified_family,
        'conservative_join_eligible': not reasons, 'rejection_reasons': '|'.join(reasons),
        'classification_json': json.dumps(classification, sort_keys=True, separators=(',', ':')),
        'diagnostics_json': json.dumps(diagnostics, sort_keys=True, separators=(',', ':'))}


class Client:
    """One sequential client, global start spacing including verification/retries."""
    def __init__(self, cache, opener=urlopen, monotonic=time.monotonic, sleep=time.sleep):
        self.cache = cache
        self.cache.mkdir(parents=True, exist_ok=True)
        self.opener, self.monotonic, self.sleep = opener, monotonic, sleep
        self.last_start = None

    def fetch(self, label, url):
        body_path = self.cache / (label + '.json')
        receipt_path = self.cache / (label + '.receipt.json')
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt['url'] != url: raise ValueError('Cached request contract drift')
            if receipt['query_status'] == 'FAILED': return None, receipt
            body = body_path.read_bytes()
            if sha(body) != receipt['response_sha256']: raise ValueError('Cached response hash drift')
            return json.loads(body), receipt
        attempts = []
        for attempt in range(1, 4):
            if self.last_start is not None:
                self.sleep(max(0, 1.05 - (self.monotonic() - self.last_start)))
            self.last_start = self.monotonic()
            started = utc()
            try:
                with self.opener(Request(url, headers={'User-Agent': 'chocho-conservative-plant-taxonomy/1.0',
                                                      'Accept': 'application/json'}), timeout=45) as handle:
                    body = handle.read(MAX_BYTES + 1)
                    status = handle.status
                if len(body) > MAX_BYTES: raise ValueError('Response exceeds 8MB bound')
                # Preserve successful HTTP bytes even if JSON/schema validation fails.
                attempt_body = self.cache / f'{label}.attempt-{attempt}.json'
                attempt_body.write_bytes(body)
                parsed = json.loads(body)
                if not isinstance(parsed, dict): raise ValueError('Taxonomy response is not an object')
                attempts.append({'attempt': attempt, 'started_utc': started, 'finished_utc': utc(),
                                 'http_status': status, 'response_sha256': sha(body), 'response_bytes': len(body)})
                body_path.write_bytes(body)
                receipt = {'url': url, 'method': 'GET', 'query_status': 'SUCCESS',
                           'response_sha256': sha(body), 'attempts': attempts}
                receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
                return parsed, receipt
            except Exception as exc:
                error = {'attempt': attempt, 'started_utc': started, 'finished_utc': utc(),
                         'error_type': type(exc).__name__, 'error': str(exc)[:500]}
                if isinstance(exc, HTTPError):
                    error_body = exc.read(MAX_BYTES + 1)
                    error_path = self.cache / f'{label}.attempt-{attempt}.error-body'
                    error_path.write_bytes(error_body)
                    error.update(http_status=exc.code, response_sha256=sha(error_body), response_bytes=len(error_body))
                attempts.append(error)
                if attempt < 3: self.sleep(2 ** attempt)
        receipt = {'url': url, 'method': 'GET', 'query_status': 'FAILED', 'response_sha256': '', 'attempts': attempts}
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
        return None, receipt


def annotate_collisions(rows):
    by_key = defaultdict(list)
    for row in rows:
        if row.get('conservative_join_eligible'):
            by_key[row['current_accepted_key']].append(row)
    collisions = []
    for key, group in sorted(by_key.items()):
        ids = '|'.join(r['wcvp_accepted_id'] for r in group)
        names = '|'.join(r['query_name'] for r in group)
        for row in group:
            row['accepted_concept_source_count'] = len(group)
            row['accepted_concept_source_ids'] = ids
        if len(group) > 1:
            collisions.append({'current_accepted_key': key, 'source_count': len(group),
                               'source_ids': ids, 'source_query_names': names})
    for row in rows:
        row.setdefault('accepted_concept_source_count', 0)
        row.setdefault('accepted_concept_source_ids', '')
    return collisions


def archive_cache(cache, destination):
    with destination.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w') as tar:
                for path in sorted(cache.iterdir()):
                    if not path.is_file(): continue
                    body = path.read_bytes()
                    info = tarfile.TarInfo(path.name)
                    info.size, info.mtime, info.mode = len(body), 0, 0o644
                    tar.addfile(info, io.BytesIO(body))


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='', encoding='utf-8') as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def acquire(root, cache, output):
    plants, edges, pins = load_sources(root)
    output.mkdir(parents=True, exist_ok=True)
    client = Client(cache)
    acquired, verifications = {}, {}
    started = utc()
    eligible = [r for r in plants if r['query_eligible'] == 'True']
    print(json.dumps({'event': 'START', 'primary_names': len(eligible), 'started_utc': started}), flush=True)
    for i, plant in enumerate(eligible, 1):
        pid = plant['wcvp_accepted_id']
        acquired[pid] = client.fetch('primary-' + pid, request_url(plant['query_name'], plant['family']))
        if i % 25 == 0 or i == len(eligible):
            counts = Counter(rec['query_status'] for _,rec in acquired.values())
            progress = {'event': 'PRIMARY_PROGRESS', 'completed': i, 'expected': len(eligible),
                        'status_counts': dict(counts), 'utc': utc()}
            (output / 'plant_taxonomy_progress.json').write_text(json.dumps(progress, indent=2) + '\n')
            print(json.dumps(progress), flush=True)
    # v2 acceptedUsage omits status for many synonym targets. Verify those
    # targets through independent exact accepted-name matches in the SAME namespace.
    requests = {}
    for plant in eligible:
        response, _ = acquired[plant['wcvp_accepted_id']]
        if not response: continue
        accepted = response.get('acceptedUsage') or {}
        if accepted and not accepted.get('status') and accepted.get('canonicalName') and accepted.get('key'):
            requests[(str(accepted['key']), plant['family'])] = canonical(accepted['canonicalName'])
    print(json.dumps({'event': 'ACCEPTED_VERIFICATION_START', 'unique_requests': len(requests)}), flush=True)
    for i, ((key, family), name) in enumerate(sorted(requests.items()), 1):
        label = 'accepted-' + sha((key + '\t' + family).encode())[:24]
        verifications[key, family] = client.fetch(label, request_url(name, family))
        if i % 25 == 0 or i == len(requests):
            print(json.dumps({'event': 'VERIFICATION_PROGRESS', 'completed': i, 'expected': len(requests), 'utc': utc()}), flush=True)
    rows = []
    empty = assess_match('', '', {})
    for plant in plants:
        pid = plant['wcvp_accepted_id']
        row = {'wcvp_accepted_id': pid, 'query_name': plant['query_name'],
            'source_accepted_name': plant['accepted_name'], 'source_family': plant['family'],
            'source_name_gate': plant['source_name_gate'], 'source_name_query_eligible': plant['query_eligible'] == 'True',
            'source_ambiguity_flags': plant['ambiguity_flags'], 'query_status': 'UNQUERIED_SOURCE_NAME',
            'response_sha256': '', 'request_url': '', 'attempt_count': 0,
            'accepted_verification_query_status': '', 'accepted_verification_response_sha256': '',
            'accepted_verification_request_url': '', **empty}
        row['rejection_reasons'] = 'unqueried_infraspecific_or_unresolved_source_name'
        if pid in acquired:
            response, receipt = acquired[pid]
            row.update(query_status=receipt['query_status'], response_sha256=receipt['response_sha256'],
                       request_url=receipt['url'], attempt_count=len(receipt['attempts']))
            if response is not None:
                accepted = response.get('acceptedUsage') or {}
                verify, vr = verifications.get((str(accepted.get('key', '')), plant['family']), (None, None))
                row.update(assess_match(plant['query_name'], plant['family'], response, verify))
                if vr:
                    row.update(accepted_verification_query_status=vr['query_status'],
                        accepted_verification_response_sha256=vr['response_sha256'], accepted_verification_request_url=vr['url'])
            else:
                row['rejection_reasons'] = 'taxonomy_request_failed'
        rows.append(row)
    collisions = annotate_collisions(rows)
    bridge_path = output / 'plant_taxonomy_bridge.csv'
    collision_path = output / 'plant_taxonomy_collisions.csv'
    write_csv(bridge_path, rows)
    write_csv(collision_path, collisions, ['current_accepted_key', 'source_count', 'source_ids', 'source_query_names'])
    archive_path = output / 'plant_taxonomy_raw_responses.tar.gz'
    archive_cache(cache, archive_path)
    good = [r for r in rows if r['conservative_join_eligible']]
    reason_counts = Counter(reason for row in rows for reason in row['rejection_reasons'].split('|') if reason)
    manifest = {'schema': 'chocho_early93_current_plant_taxonomy_20261030_v1', 'status': 'FINISHED',
        'started_utc': started, 'finished_utc': utc(), 'source_sha256': pins,
        'taxonomy': {'endpoint': ENDPOINT, 'checklistKey': CHECKLIST, 'kingdom': 'Plantae', 'taxonRank': 'SPECIES',
            'strict': True, 'verbose': True, 'family': 'exact source family per request'},
        'request_bounds': {'attempts_per_request_max': 3, 'global_minimum_request_start_interval_seconds': 1.05,
                           'timeout_seconds': 45, 'response_bytes_max': MAX_BYTES},
        'summary': {'source_rows': len(rows), 'source_query_eligible_rows': len(eligible),
            'early_butterfly_count': 93, 'early_catalogue_edges': len(edges),
            'queried_rows': len(acquired), 'successful_rows': sum(r['query_status'] == 'SUCCESS' for r in rows),
            'failed_rows': sum(r['query_status'] == 'FAILED' for r in rows),
            'unqueried_source_name_rows': sum(r['query_status'] == 'UNQUERIED_SOURCE_NAME' for r in rows),
            'conservative_join_eligible_rows': len(good), 'unique_eligible_accepted_concepts': len({r['current_accepted_key'] for r in good}),
            'duplicate_accepted_concept_groups': len(collisions),
            'source_rows_in_duplicate_accepted_concept_groups': sum(c['source_count'] for c in collisions),
            'accepted_verification_requests': len(verifications),
            'accepted_verification_failures': sum(rec['query_status'] == 'FAILED' for _,rec in verifications.values()),
            'primary_http_attempts': sum(len(rec['attempts']) for _,rec in acquired.values()),
            'accepted_verification_http_attempts': sum(len(rec['attempts']) for _,rec in verifications.values()),
            'rejection_reason_counts': dict(sorted(reason_counts.items()))},
        'api_schema_adaptation': 'If acceptedUsage omits status, a separate exact accepted-canonical-name match must return usage.status=ACCEPTED and the same string key, species rank, Plantae, source family and no competing exact concept. Missing verification remains unresolved.',
        'duplicate_concept_policy': 'All source IDs and edges retained; deduplicate eligible accepted concepts per butterfly when joining occurrence facets. Collisions do not themselves make an unambiguous synonym invalid.',
        'missing_state_policy': 'Unqueried source names, failed requests and rejected matches are unresolved taxonomy, never biological zero.',
        'interpretation': 'Catalogue-compatible taxonomy only; it does not establish island host availability, introduced status or feeding.',
        'output_sha256': {p.name: sha(p.read_bytes()) for p in (bridge_path, collision_path, archive_path)}}
    (output / 'plant_taxonomy_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'event': 'FINISHED', 'summary': manifest['summary']}), flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--cache', type=Path, default=Path('/tmp/chocho-plant-taxonomy-20261030'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    acquire(args.root, args.cache, args.output or args.root / OUTPUT)


if __name__ == '__main__':
    main()
