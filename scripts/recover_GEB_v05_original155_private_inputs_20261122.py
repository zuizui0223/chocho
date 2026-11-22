#!/usr/bin/env python3
"""Restore verified original155 inputs privately and attest a fresh original fit.

Only fixed aggregate receipts leave the runner. Recovery assets are data, never
executed code. Existing joblib bytes are scanned statically, never unpickled.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile
import zlib

REPO = 'zuizui0223/chocho'
RELEASE_ID = 453667492
SOURCE_COMMIT = '42f8fcc3a388e96456c76180d7e56d75dd2ed48c'
ARCHIVE_SHA256 = '9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f'
PANEL = 'work/temporal/forecast/species_cell_period_panel.parquet'
PANEL_SHA256 = '5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1'
LOCK = 'outputs/forecast_model_lock.json'
LOCK_SHA256 = '0c7fc34d7858d6f62fe31a914188aee63e56ab1d41872abe40c314f6c5ed8b1e'
SCRIPT_SHA256 = 'a34ad4cd1e4f40460c4f6bf6ed0b22ba3885b17d20ac17f74bd6564d0618ce6d'
RECEIPT = (775476029, 569496, 'b8eff912a43373424e5f31ee8714c89d672f6d44c4011bbd81514e8a9afc3aa2')
PARTS = [
    (775475943, 1500000000, '4fca3d7d1d20c431455cea62ffd40d7074008921a5500586a305aafa8e20812b'),
    (775475944, 1500000000, '23ff9c59c23ffbb8f0fb34367c0871e4661489ca7bf5b4d55b581b41f4f047c3'),
    (775475941, 1500000000, 'b6f7c09f29c0dd0117a72ed49a931b85b0aff254235df90392b9f069e1446bbb'),
    (775475945, 1052304553, 'f0e0eccfbe24ac1e9d285b0331ea159c6fc3da7e7e0d0f5549b9e2874dce6ce7'),
]
OLD_GATE_PATHS = [PANEL,
    'work/temporal/forecast/species_recording_days.parquet',
    'work/temporal/forecast/target_group_recording_days.parquet',
    'work/temporal/forecast/candidate_cells.parquet',
    'work/temporal/forecast/cell_environment_period.parquet',
    'work/temporal/forecast/cell_climate_grid.parquet',
    'work/temporal/forecast/predictions_temporal.parquet',
    'work/temporal/forecast/predictions_spatial.parquet',
    'work/temporal/species_year_cells_10km.parquet',
    'work/temporal/host_region_composition_metrics.parquet',
    'work/temporal/host_region_family_composition.parquet',
    'work/temporal/occurrence_parts',
    'work/remote_sensing/host_region_counts.csv',
]
MODEL_SETTINGS = {'n_estimators': 160, 'max_depth': 12, 'min_samples_leaf': 40,
                  'max_features': 0.7, 'random_state': 20261029, 'n_jobs': 4}


class RecoveryError(ValueError):
    """An approved constant diagnostic, never a third-party exception message."""


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def append_verified_part(stream, output, expected_sha256, expected_bytes, combined):
    """Stream a part; callers must discard the joined file on any failure."""
    digest, size = hashlib.sha256(), 0
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        size += len(block)
        if size > expected_bytes:
            raise RecoveryError('Source part exceeds its registered size')
        digest.update(block)
        combined.update(block)
        output.write(block)
    if size != expected_bytes or digest.hexdigest() != expected_sha256:
        raise RecoveryError('Source part byte count or SHA256 mismatch')
    return size


def safe_zip_member(member):
    name = member.filename
    path = PurePosixPath(name)
    if (not name or '\\' in name or path.is_absolute() or '..' in path.parts
            or path == PurePosixPath('.') or ':' in path.parts[0]):
        raise RecoveryError('Unsafe archive member path')
    if stat.S_ISLNK(member.external_attr >> 16):
        raise RecoveryError('Archive symlink is not an input file')
    return path.as_posix()


def find_manifest_record(records, relative_path):
    matches = [entry for entry in records
               if isinstance(entry, dict) and isinstance(entry.get('path'), str)
               and (entry['path'].replace('\\', '/') == relative_path
                    or entry['path'].replace('\\', '/').endswith('/' + relative_path))]
    if len(matches) != 1:
        raise RecoveryError('Source manifest entry is missing or ambiguous')
    entry = matches[0]
    if (not isinstance(entry.get('bytes'), int) or entry['bytes'] < 0
            or not re.fullmatch('[0-9a-f]{64}', str(entry.get('sha256', '')))):
        raise RecoveryError('Source manifest size or digest is invalid')
    return entry


def find_zip_member(members, relative_path, required=True):
    matches = [member for member in members
               if not member.is_dir() and
               (member.filename == relative_path or member.filename.endswith('/' + relative_path))]
    if not matches and not required:
        return None
    if len(matches) != 1:
        raise RecoveryError('Archive input is missing or ambiguous')
    return matches[0]


def extract_verified_member(archive, member, destination, record):
    safe_zip_member(member)
    if member.file_size != record['bytes']:
        raise RecoveryError('Archive member size differs from source manifest')
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with archive.open(member) as src, destination.open('wb') as out:
            append_verified_part(src, out, record['sha256'], record['bytes'], hashlib.sha256())
    except Exception:
        destination.unlink(missing_ok=True)
        raise
    return {'sha256': record['sha256'], 'bytes': record['bytes']}


def static_sklearn_version(payload):
    """Read only a pickle string opcode adjacent to the known version key."""
    if payload.startswith(b'\x78'):
        decoder = zlib.decompressobj()
        payload = decoder.decompress(payload, 512 * 1024 * 1024)
        if not decoder.eof:
            raise RecoveryError('Compressed model metadata exceeds the inspection bound')
    versions = set()
    for match in re.finditer(b'_sklearn_version', payload):
        tail = payload[match.end():match.end() + 80]
        version = re.search(rb'\x8c([\x03-\x20])(\d+\.\d+\.\d+(?:rc\d+)?)(?=\x94|q|r)', tail)
        if version and version.group(1)[0] == len(version.group(2)):
            versions.add(version.group(2).decode('ascii'))
    if len(versions) > 1:
        raise RecoveryError('Recovered model contains conflicting library versions')
    return next(iter(versions)) if versions else None


def gh_json(endpoint):
    result = subprocess.run(['gh', 'api', endpoint], capture_output=True, check=False)
    if result.returncode:
        status = re.search(r'HTTP (\d{3})', result.stderr.decode(errors='replace'))
        raise RecoveryError('Authenticated metadata read failed; HTTP '
                         + (status.group(1) if status else 'unknown'))
    return json.loads(result.stdout)


def stream_asset(asset_id, output, size, digest, combined, error_dir):
    with tempfile.TemporaryFile(dir=error_dir) as err:
        process = subprocess.Popen(['gh', 'api', f'repos/{REPO}/releases/assets/{asset_id}',
                                    '-H', 'Accept: application/octet-stream'],
                                   stdout=subprocess.PIPE, stderr=err)
        try:
            copied = append_verified_part(process.stdout, output, digest, size, combined)
            process.stdout.close()
            if process.wait() != 0:
                raise RecoveryError('Authenticated asset download failed')
            return copied
        except Exception:
            process.kill()
            process.wait()
            raise


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def restore(root, mode, fit_script, output):
    import pandas as pd
    import pyarrow.parquet as pq

    root.mkdir(parents=True, exist_ok=False)
    original, fit = root / 'original', root / 'refit'
    original.mkdir()
    fit.mkdir()
    release = gh_json(f'repos/{REPO}/releases/{RELEASE_ID}')
    if release.get('draft') is not True:
        raise RecoveryError('Source release is no longer the expected private draft')
    assets = {asset['id']: asset for asset in release['assets']}
    for asset_id, size, digest in [RECEIPT] + PARTS:
        if (assets.get(asset_id, {}).get('size') != size
                or assets[asset_id].get('digest') != 'sha256:' + digest):
            raise RecoveryError('Registered source asset size or digest changed')
    receipt_path = root / 'source_receipt.json'
    with receipt_path.open('wb') as out:
        stream_asset(RECEIPT[0], out, RECEIPT[1], RECEIPT[2], hashlib.sha256(), root)
    records = json.loads(receipt_path.read_text())['files']
    if len(records) != 3032:
        raise RecoveryError('Original source manifest file count changed')
    archive_path, combined = root / 'original.zip', hashlib.sha256()
    try:
        with archive_path.open('wb') as out:
            for number, (asset_id, size, digest) in enumerate(PARTS, 1):
                stream_asset(asset_id, out, size, digest, combined, root)
                print(f'VERIFIED_ARCHIVE_PART {number}/4 bytes={size}', flush=True)
        if combined.hexdigest() != ARCHIVE_SHA256:
            raise RecoveryError('Complete restored ZIP SHA256 mismatch')
        verified, optional_csvs = {}, []
        prediction = f'work/temporal/forecast/predictions_{mode}.parquet'
        model_metadata = 'work/temporal/forecast/model_M0_background.joblib'
        with zipfile.ZipFile(archive_path) as archive:
            members = archive.infolist()
            names = [safe_zip_member(member) for member in members]
            if len(names) != len(set(names)):
                raise RecoveryError('Archive contains duplicate normalized paths')
            availability = {path: (any('/' + path + '/' in '/' + n for n in names)
                             if path.endswith('occurrence_parts')
                             else find_zip_member(members, path, required=False) is not None)
                            for path in OLD_GATE_PATHS}
            selected = [p for p in OLD_GATE_PATHS if p.endswith('.parquet') and availability[p]]
            selected += [LOCK, prediction]
            if find_zip_member(members, model_metadata, required=False):
                selected.append(model_metadata)
            for model in ['M0_background', 'M1_climate', 'M2_landcover', 'M3_both',
                          'M4_host_resources', 'M5_host_composition']:
                path = f'outputs/forecast_{mode}_{model}_species_metrics.csv'
                if find_zip_member(members, path, required=False):
                    selected.append(path)
                    optional_csvs.append(path)
            for path in dict.fromkeys(selected):
                member = find_zip_member(members, path)
                record = find_manifest_record(records, path)
                verified[path] = extract_verified_member(archive, member, original / path, record)
            total_uncompressed = sum(m.file_size for m in members if not m.is_dir())
            file_count = sum(not m.is_dir() for m in members)
    except Exception:
        archive_path.unlink(missing_ok=True)
        raise
    # The draft release remains the persistent source; avoid retaining a second 5.55GB copy.
    archive_path.unlink()
    if verified[PANEL]['sha256'] != PANEL_SHA256 or verified[LOCK]['sha256'] != LOCK_SHA256:
        raise RecoveryError('Original panel or pre-fit model lock changed')
    if sha256_file(fit_script) != SCRIPT_SHA256:
        raise RecoveryError('Frozen original fit script checksum mismatch')
    lock = json.loads((original / LOCK).read_text())
    if lock.get('status') != 'LOCKED_BEFORE_MODEL_FIT' or lock.get('random_forest') != MODEL_SETTINGS:
        raise RecoveryError('Original fixed-model settings changed')
    columns = ['speciesKey', 'species', 'cell', 'period', 'recording_probability10']
    panel = pd.read_parquet(original / PANEL, columns=columns)
    if len(panel) != 842964 or panel.species.nunique() != 155:
        raise RecoveryError('Original source panel cohort changed')
    if panel.duplicated(['speciesKey', 'cell', 'period']).any():
        raise RecoveryError('Source panel has duplicate species-cell-period keys')
    transitions = {}
    for period, name, expected in [(0, 'training', 193703), (1, 'test', 208482)]:
        before = panel.loc[panel.period.eq(period), ['speciesKey', 'species', 'cell', 'recording_probability10']]
        after = panel.loc[panel.period.eq(period+1), ['speciesKey', 'cell', 'recording_probability10']]
        pairs = before.merge(after, on=['speciesKey', 'cell'], validate='one_to_one', suffixes=('_before', '_after'))
        pairs = pairs.dropna(subset=['recording_probability10_before', 'recording_probability10_after'])
        if len(pairs) != expected or pairs.species.nunique() != 155:
            raise RecoveryError('Original evaluable transition cohort changed')
        transitions[name] = {'rows': len(pairs), 'species': int(pairs.species.nunique())}
    recording_path = original / 'work/temporal/forecast/species_recording_days.parquet'
    dataset_key_retained = recording_path.exists() and 'datasetKey' in pq.read_schema(recording_path).names
    if not dataset_key_retained:
        raise RecoveryError('Recording-day source no longer retains datasetKey attribution')
    metadata_path = original / model_metadata
    source_version = static_sklearn_version(metadata_path.read_bytes()) if metadata_path.exists() else None
    replay_version = source_version or '1.8.0'
    (root / 'model_requirement.txt').write_text('scikit-learn==' + replay_version + '\n')
    for path in [PANEL, LOCK]:
        (fit / path).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original / path, fit / path)
    dest_script = fit / 'scripts/temporal/fit_forecast.py'
    dest_script.parent.mkdir(parents=True)
    shutil.copyfile(fit_script, dest_script)
    restored = {
        'schema': 'chocho_original155_verified_private_source_restore_v1',
        'status': 'ORIGINAL_PANEL_AND_LOCK_RESTORED_AND_COHORT_VERIFIED;FIT_PENDING',
        'source_commit': SOURCE_COMMIT, 'release_id': RELEASE_ID, 'release_remains_draft': True,
        'source_receipt_sha256': RECEIPT[2], 'archive_sha256': ARCHIVE_SHA256,
        'archive_bytes': sum(p[1] for p in PARTS), 'source_manifest_files': len(records),
        'archive_files': file_count, 'archive_uncompressed_bytes': total_uncompressed,
        'verified_archive_parts': 4, 'verified_source_files': verified,
        'original_missing_gate_paths_present_in_complete_archive': availability,
        'panel_rows': len(panel), 'panel_species': int(panel.species.nunique()),
        'transitions': transitions, 'recording_days_datasetKey_retained': dataset_key_retained,
        'model_lock_sha256': LOCK_SHA256, 'source_fit_script_sha256': SCRIPT_SHA256,
        'model_settings': MODEL_SETTINGS, 'mode_to_refit': mode,
        'original_prediction_sha256': verified[prediction]['sha256'],
        'original_direct_species_metric_csvs': len(optional_csvs),
        'source_sklearn_version_from_static_pickle_string': source_version,
        'requested_replay_sklearn_version': replay_version,
        'original_complete_runtime_lock_available': False,
        'model_refit_completed': False, 'no_original_records_exported': True,
        'reviewer_redistribution_and_anonymity_cleared': False,
        'scope': 'Complete original archive byte identity and selected input payload hashes; no claim of all-file re-execution or upstream raw-to-panel rebuild.',
    }
    write_json(root / 'restore_receipt.json', restored)
    write_json(output, restored)
    print('ORIGINAL155_PRIVATE_RESTORE_PASS', mode, '842964 rows; 155 species; training193703/test208482', flush=True)


def execute_fit(root, mode, output):
    fit = root / 'refit'
    script = fit / 'scripts/temporal/fit_forecast.py'
    work = fit / 'work/temporal/forecast'
    outputs = fit / 'outputs'
    if sorted(p.name for p in outputs.iterdir()) != ['forecast_model_lock.json']:
        raise RecoveryError('Fit output directory is not fresh')
    if sorted(p.name for p in work.iterdir()) != ['species_cell_period_panel.parquet']:
        raise RecoveryError('Fit work directory contains pre-existing fitted outputs')
    if sha256_file(script) != SCRIPT_SHA256 or sha256_file(work / 'species_cell_period_panel.parquet') != PANEL_SHA256:
        raise RecoveryError('Prepared model source or input changed before execution')
    started = datetime.now(timezone.utc).isoformat()
    with (root / 'private_model_execution.log').open('w') as log:
        result = subprocess.run([sys.executable, str(script), mode], cwd=fit,
                                stdout=log, stderr=subprocess.STDOUT, check=False)
    completed = datetime.now(timezone.utc).isoformat()
    receipt = {'run_id': os.environ.get('GITHUB_RUN_ID', 'local'), 'mode': mode,
               'source_script_sha256': SCRIPT_SHA256, 'panel_sha256': PANEL_SHA256,
               'model_lock_sha256': LOCK_SHA256, 'started_at': started, 'completed_at': completed,
               'returncode': result.returncode, 'preexisting_outputs': False}
    write_json(output, receipt)
    for line in (root / 'private_model_execution.log').read_text().splitlines():
        if line.startswith((mode + ' ', 'COMPLETE ', 'M0_background ', 'M1_climate ',
                            'M2_landcover ', 'M3_both ', 'M4_host_resources ', 'M5_host_composition ')):
            print(line, flush=True)
    if result.returncode:
        raise RecoveryError('Original model execution failed; no successful refit claimed')
    print('ORIGINAL155_FRESH_MODEL_EXECUTION_COMPLETE', mode, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['restore', 'fit'], required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--mode', choices=['temporal', 'spatial', 'species'], required=True)
    parser.add_argument('--fit-script', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'restore':
        if args.fit_script is None:
            parser.error('--fit-script required for restore')
        restore(args.root.resolve(), args.mode, args.fit_script.resolve(), args.output)
    else:
        execute_fit(args.root.resolve(), args.mode, args.output)


def cli():
    try:
        main()
        return 0
    except RecoveryError as error:
        print(str(error), file=sys.stderr)
        return 1
    except Exception as error:
        # Unexpected library errors can contain raw data values or private paths.
        print('ORIGINAL155_RECOVERY_OR_EXECUTION_FAILED: ' + type(error).__name__, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(cli())
