#!/usr/bin/env python3
"""Audit archived species scores and historical metadata without exposing records."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
import zipfile

import numpy as np
import pandas as pd

import recover_GEB_v05_original155_private_inputs_20261122 as recovery
import verify_GEB_v05_original155_refit_20261122 as check

PREDICTION = 'work/temporal/forecast/predictions_species.parquet'
SCORES = 'outputs/forecast_species_results.json'
NATIVE_SCORE_SHA256 = '3eaaeeb8b77088314bf42051c5a51c07895c17615bf5f27fc3a33424f25e125a'
SCRIPT = 'scripts/temporal/fit_forecast.py'
RUNTIME_FILES = {
    'requirements.txt', 'requirements.lock', 'requirements-dev.txt', 'requirements.in',
    'environment.yml', 'environment.yaml', 'conda-lock.yml', 'conda-lock.yaml',
    'pyproject.toml', 'poetry.lock', 'uv.lock', 'Pipfile', 'Pipfile.lock',
    'pip-freeze.txt', 'pip_freeze.txt', 'runtime.txt', '.python-version',
    'sessionInfo.txt', 'runtime_versions.json', 'environment.json',
}
VERSION = re.compile(r'[0-9]{1,4}(?:\.[0-9]{1,4}){1,3}(?:(?:a|b|rc)[0-9]{1,4})?')
WRITER = re.compile(r'(parquet-cpp-arrow|parquet-mr|fastparquet) version (' + VERSION.pattern + ')')


def digest(value):
    return hashlib.sha256(value).hexdigest()


def archive_metadata(archive, members, records):
    panel = recovery.find_zip_member(members, recovery.PANEL)
    prefix = panel.filename[:-len(recovery.PANEL)]
    entries = {name: [] for name in (SCRIPT, SCORES, PREDICTION, recovery.PANEL, recovery.LOCK)}
    environment = []
    unbound = 0
    for member in members:
        name = recovery.safe_zip_member(member)
        if member.is_dir() or not name.startswith(prefix):
            continue
        rel = name[len(prefix):]
        base = PurePosixPath(rel).name
        if rel not in entries and base not in RUNTIME_FILES:
            continue
        try:
            record = recovery.find_manifest_record(records, rel)
        except recovery.RecoveryError:
            unbound += 1
            continue
        h = hashlib.sha256()
        with archive.open(member) as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        if member.file_size != record['bytes'] or h.hexdigest() != record['sha256']:
            raise recovery.RecoveryError('Diagnostic member failed manifest verification')
        item = dict(bytes=member.file_size, sha256=h.hexdigest(),
                    zip_local_datetime=list(member.date_time))
        if rel in (SCRIPT, SCORES) and member.file_size <= 2 * 1024 * 1024:
            content = archive.read(member)
            item['crlf_to_lf_sha256'] = digest(content.replace(b'\r\n',b'\n'))
            item['crlf_pair_count'] = content.count(b'\r\n')
            if rel == SCORES:
                historical = json.loads(content)
                item['header_fields_present'] = {key:key in historical for key in
                                                 ('sensitivity','input_panel','rarefaction_n')}
        if rel in entries:
            entries[rel].append(item)
        else:
            environment.append(dict(allowlisted_basename=base, **item))
    source = entries[SCRIPT]
    scores = entries[SCORES]
    return dict(selected_members=entries, runtime_file_evidence=environment,
                unbound_candidate_count=unbound,
                private_script_matches_public_sha256=source[0]['sha256'] == recovery.SCRIPT_SHA256 if source else None,
                private_scores_match_frozen_sha256=scores[0]['sha256'] == check.FROZEN['species'] if scores else None,
                private_scores_match_manifest_native_sha256=(
                    scores[0]['sha256'] == NATIVE_SCORE_SHA256 if scores else None),
                private_script_crlf_normalized_matches_public=(
                    source[0]['crlf_to_lf_sha256']==recovery.SCRIPT_SHA256
                    if source and 'crlf_to_lf_sha256' in source[0] else None),
                private_scores_crlf_normalized_match_frozen=(
                    scores[0]['crlf_to_lf_sha256']==check.FROZEN['species']
                    if scores and 'crlf_to_lf_sha256' in scores[0] else None),
                zip_datetime_interpretation='ZIP local timestamp, timezone unknown. No execution timestamp inferred.')


def parquet_metadata(path):
    import pyarrow.parquet as pq
    f = pq.ParquetFile(path)
    schema = f.schema_arrow
    metadata = json.loads((schema.metadata or {}).get(b'pandas', b'{}'))
    ver = metadata.get('pandas_version', '')
    writer = f.metadata.created_by or ''
    dtypes = [str(field.type) for field in schema]
    return dict(rows=f.metadata.num_rows, columns=len(schema), row_groups=f.metadata.num_row_groups,
                arrow_schema_sha256=digest(schema.remove_metadata().serialize().to_pybytes()),
                ordered_dtype_sha256=digest(json.dumps(dtypes).encode()),
                pandas_version=ver if isinstance(ver,str) and VERSION.fullmatch(ver) else None,
                created_by=writer if WRITER.fullmatch(writer) else None,
                created_by_sha256=digest(writer.encode()))


def inputs(root):
    original = root / 'original'
    panel, preds = original / recovery.PANEL, original / PREDICTION
    receipt = json.loads((root / 'restore_receipt.json').read_text())
    if recovery.sha256_file(panel) != recovery.PANEL_SHA256 or recovery.sha256_file(preds) != receipt['verified_source_files'][PREDICTION]['sha256']:
        raise recovery.RecoveryError('Diagnostic inputs differ from restored sources')
    return panel, preds


def transition_evidence(panel, predictions):
    cols = [*check.KEYS, 'recording_probability10']
    before = panel.loc[panel.period.eq(1), cols]
    after = panel.loc[panel.period.eq(2), cols]
    pairs = before.merge(after, on=check.KEYS, validate='one_to_one',
                         suffixes=('_before', '_after')).dropna()
    ordered = pd.MultiIndex.from_frame(pairs[check.KEYS])
    existing = pd.MultiIndex.from_frame(predictions[check.KEYS])
    check.require(ordered.is_unique and existing.is_unique, 'duplicate_transition_keys')
    same = (len(existing) == len(ordered) and ordered.difference(existing).empty)
    evidence = dict(reconstructed_rows=len(pairs), archived_rows=len(predictions),
                    same_key_set=same, same_key_order=ordered.equals(existing))
    if same:
        old = predictions.set_index(check.KEYS).reindex(ordered)
        evidence.update(target_exact=np.array_equal(pairs.recording_probability10_after, old.target),
                        baseline_exact=np.array_equal(pairs.recording_probability10_before,
                                                      old.baseline_probability))
    return evidence


def rescore(frame, frozen):
    check.keyed(frame)
    check.require(len(frame) == 208482 and frame.species.nunique() == 155, 'species_support')
    audit = check.Audit()
    for model in check.MODELS:
        check.require(np.isfinite(frame[[model,'target','baseline_probability']]).all().all(),
                      'nonfinite_archived_predictions')
        metrics, _ = check.scores(frame, model)
        check.metrics_compare(audit, model, frozen[model], metrics)
    return audit


def lineage(root, snapshot, output):
    root.mkdir(parents=True, exist_ok=False)
    if recovery.sha256_file(snapshot / 'analysis_source/fit_forecast.py') != recovery.SCRIPT_SHA256:
        raise recovery.RecoveryError('Public original fit source hash mismatch')
    release = recovery.gh_json(f'repos/{recovery.REPO}/releases/{recovery.RELEASE_ID}')
    assets = {x['id']: x for x in release['assets']}
    check.require(release.get('draft') is True, 'private_draft_required')
    for asset, size, sha in [recovery.RECEIPT, *recovery.PARTS]:
        item = assets.get(asset, {})
        check.require(item.get('size') == size and item.get('digest') == 'sha256:'+sha, 'asset_binding')
    asset, size, sha = recovery.RECEIPT
    receipt = root / 'manifest.json'
    with receipt.open('xb') as out:
        recovery.stream_asset(asset,out,size,sha,hashlib.sha256(),root)
    records = json.loads(receipt.read_text())['files']
    check.require(len(records)==3032,'source_manifest_count')
    archive_path, combined = root / 'source.zip', hashlib.sha256()
    selected = {}
    try:
        with archive_path.open('xb') as out:
            for asset,size,sha in recovery.PARTS:
                recovery.stream_asset(asset,out,size,sha,combined,root)
        check.require(combined.hexdigest()==recovery.ARCHIVE_SHA256,'combined_zip_hash')
        with zipfile.ZipFile(archive_path) as archive:
            members = archive.infolist()
            names = [recovery.safe_zip_member(m) for m in members]
            check.require(len(names)==len(set(names)), 'duplicate_normalized_members')
            archive_evidence = archive_metadata(archive,members,records)
            for rel in (recovery.PANEL, PREDICTION):
                record = recovery.find_manifest_record(records, rel)
                member = recovery.find_zip_member(members, rel)
                selected[rel] = recovery.extract_verified_member(archive,member,root/'original'/rel,record)
    finally:
        archive_path.unlink(missing_ok=True)
    recovery.write_json(root/'restore_receipt.json', dict(verified_source_files=selected))
    panel, preds = inputs(root)
    frozen = check.load_frozen(snapshot, 'species')
    predictions = pd.read_parquet(preds)
    identity = transition_evidence(pd.read_parquet(panel, columns=[
        *check.KEYS, 'period', 'recording_probability10']), predictions)
    audit = rescore(predictions, frozen)
    report = dict(schema='geb_species_archived_prediction_lineage_v1', fresh_fit_executed=False,
                  archive_sha256=recovery.ARCHIVE_SHA256, manifest_sha256=recovery.RECEIPT[2],
                  panel_sha256=recovery.PANEL_SHA256, frozen_species_scores_sha256=check.FROZEN['species'],
                  archived_species_predictions_sha256=recovery.sha256_file(preds), archive=archive_evidence,
                  test_transition_alignment=identity, panel_metadata=parquet_metadata(panel), predictions_metadata=parquet_metadata(preds),
                  absolute_tolerance=check.ATOL, relative_tolerance=check.RTOL,
                  archived_rows_match_frozen_scores=audit.passed, comparisons=audit.comparisons)
    recovery.write_json(output, report)
    return audit.passed


def main():
    parser = argparse.ArgumentParser()
    for key in ['root','snapshot','output']:
        parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args()
    ok = lineage(a.root.resolve(), a.snapshot.resolve(), a.output)
    print('DIAGNOSTIC_COMPARISON_MATCH' if ok else 'DIAGNOSTIC_COMPARISON_MISMATCH')
    return 0 if ok else 1


if __name__=='__main__':
    try:
        sys.exit(main())
    except (recovery.RecoveryError,check.Invalid) as exc:
        sys.exit(str(exc))
    except Exception as exc:
        sys.exit('DIAGNOSTIC_ERROR_'+type(exc).__name__)
