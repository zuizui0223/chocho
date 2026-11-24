#!/usr/bin/env python3
"""Audit scikit-learn version strings in the 30 ORIGINAL private joblib objects.

Never unpickle or execute original model bytes. Fetch the immutable draft release,
check complete ZIP+per-object SHA256, statically parse a bounded version string,
and export aggregate version metadata only. No source file paths or observations
are included in the output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile

import recover_GEB_v05_original155_private_inputs_20261122 as source

MODELS = (
    'M0_background', 'M1_climate', 'M2_landcover', 'M3_both',
    'M4_host_resources', 'M5_host_composition'
)
CONDITIONS = ('primary', 'effort20', 'radius250', 'exclude_top_provider', 'precision100')


class VersionAuditError(ValueError):
    pass


def require(condition, reason):
    if not condition:
        raise VersionAuditError(reason)


def model_relative_path(condition, model):
    require(condition in CONDITIONS and model in MODELS, 'invalid_expected_model')
    suffix = '' if condition == 'primary' else '_' + condition
    return 'work/temporal/forecast/model_' + model + suffix + '.joblib'


def read_version_safely(data):
    # The original helper decompresses with a 512-MB bound and uses pickle
    # literal-string inspection ONLY; there is no pickle.load/joblib.load.
    version = source.static_sklearn_version(data)
    if version is not None:
        require(isinstance(version, str) and len(version) <= 24, 'invalid_version_string')
    return version


def run(output_path):
    release = source.gh_json(f'repos/{source.REPO}/releases/{source.RELEASE_ID}')
    require(release.get('draft') is True, 'source_release_not_draft')
    assets = {item['id']: item for item in release['assets']}
    for asset_id, size, digest in [source.RECEIPT, *source.PARTS]:
        asset = assets.get(asset_id, {})
        require(asset.get('size') == size and asset.get('digest') == 'sha256:' + digest,
                'source_asset_identity_changed')

    with tempfile.TemporaryDirectory(prefix='original155-version-audit-') as directory:
        root = Path(directory)
        manifest_path = root / 'manifest.json'
        with manifest_path.open('wb') as out:
            source.stream_asset(
                source.RECEIPT[0], out, source.RECEIPT[1], source.RECEIPT[2],
                hashlib.sha256(), root
            )
        entries = json.loads(manifest_path.read_text())['files']
        require(len(entries) == 3032, 'source_manifest_record_count')

        combined = hashlib.sha256()
        archive_path = root / 'source.zip'
        with archive_path.open('wb') as out:
            for asset_id, size, digest in source.PARTS:
                source.stream_asset(asset_id, out, size, digest, combined, root)
        require(combined.hexdigest() == source.ARCHIVE_SHA256, 'original_zip_digest')

        items = []
        with zipfile.ZipFile(archive_path) as zipped:
            members = zipped.infolist()
            names = [source.safe_zip_member(m) for m in members]
            require(len(names) == len(set(names)), 'archive_duplicate_names')

            for condition in CONDITIONS:
                for model in MODELS:
                    rel = model_relative_path(condition, model)
                    receipt = source.find_manifest_record(entries, rel)
                    member = source.find_zip_member(members, rel)
                    require(0 < receipt['bytes'] <= 100_000_000, 'model_size_guard')
                    target = root / 'single_model.joblib'
                    try:
                        source.extract_verified_member(zipped, member, target, receipt)
                        # Historical compressed pickle bytes are NEVER deserialized.
                        version = read_version_safely(target.read_bytes())
                    finally:
                        target.unlink(missing_ok=True)
                    items.append({
                        'condition': condition,
                        'model': model,
                        'source_sha256': receipt['sha256'],
                        'source_bytes': receipt['bytes'],
                        'static_sklearn_version': version,
                        'version_readable': version is not None,
                    })
        require(len(items) == 30, 'incorrect_model_count')

    summary = {}
    for condition in CONDITIONS:
        group = [i for i in items if i['condition'] == condition]
        counts = {}
        for x in group:
            key = x['static_sklearn_version'] or 'UNREADABLE'
            counts[key] = counts.get(key, 0) + 1
        summary[condition] = counts

    result = {
        'schema': 'GEB_original155_model_pickle_static_version_audit_v1',
        'source_commit': source.SOURCE_COMMIT,
        'source_archive_sha256': source.ARCHIVE_SHA256,
        'source_manifest_sha256': source.RECEIPT[2],
        'model_lock_sha256': source.LOCK_SHA256,
        'status': 'SOURCE_30_MODEL_OBJECTS_HASH_VERIFIED_STATIC_VERSION_METADATA_ONLY',
        'models_verified': len(items),
        'condition_version_counts': summary,
        'per_model': items,
        'no_unpickling_or_historical_model_execution': True,
        'no_observations_or_local_private_paths_exported': True,
        'source_release_remained_draft': True,
        'interpretation': (
            'Static version strings describe the saved object metadata; they do not '
            'prove which code revision generated the archived species predictions or '
            'establish historical fitting order. A matching version is not evidence '
            'of numerical reproduction.'
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n',
                           encoding='utf-8')
    print('STATIC_JOBLIB_VERSION_AUDIT_PASS', len(items), json.dumps(summary, sort_keys=True),
          flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.output)
    except (VersionAuditError, source.RecoveryError) as exc:
        # source.RecoveryError labels contain no third-party source rows or paths.
        raise SystemExit('STATIC_VERSION_AUDIT_FAILED_' + str(exc))
    except Exception as exc:
        raise SystemExit('STATIC_VERSION_AUDIT_ERROR_' + type(exc).__name__)
