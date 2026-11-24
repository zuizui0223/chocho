#!/usr/bin/env python3
"""Bounded original155 partial panel reconstruction from independently saved upstream inputs.

No panel row, source observation, source path, or generated panel is published.
The missing original host-count CSV is NEVER imputed from the final panel.
"""
import argparse
import hashlib
import json
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.neighbors import BallTree
from pyproj import Transformer

import recover_GEB_v05_original155_private_inputs_20261122 as recovery

SOURCE_COMMIT = recovery.SOURCE_COMMIT
ORIGINAL_BUILD_SCRIPT_SHA = '99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6'
PUBLISHED_SELECTION_SHA = 'ae0a874b1d743657fdd95a1e653d7ca0939fc41a56b947f8a077be7785c4332d'
PANEL = recovery.PANEL
SOURCE_FILES = (
    'work/temporal/forecast/candidate_cells.parquet',
    'work/temporal/forecast/cell_environment_period.parquet',
    'work/temporal/forecast/cell_climate_grid.parquet',
    'work/temporal/forecast/species_recording_days.parquet',
    'work/temporal/forecast/target_group_recording_days.parquet',
    'work/temporal/species_year_cells_10km.parquet',
    'work/temporal/host_region_composition_metrics.parquet',
    'work/temporal/host_region_family_composition.parquet',
    recovery.LOCK,
    PANEL,
)
MISSING = (
    'outputs/forecast_candidate_disposition.csv',
    'outputs/temporal_host_family_support.csv',
    'work/remote_sensing/host_region_counts.csv',
)


class RebuildError(Exception):
    pass


def require(value, label):
    if not value:
        raise RebuildError(label)


def restore_upstream(root):
    release = recovery.gh_json(f'repos/{recovery.REPO}/releases/{recovery.RELEASE_ID}')
    require(release.get('draft') is True, 'private_release_is_not_draft')
    by_id = {asset['id']: asset for asset in release['assets']}
    for aid, nbytes, digest in [recovery.RECEIPT, *recovery.PARTS]:
        current = by_id.get(aid, {})
        require(current.get('size') == nbytes and current.get('digest') == 'sha256:' + digest,
                'private_source_asset_identity')
    manifest = root / 'manifest.json'
    with manifest.open('wb') as out:
        recovery.stream_asset(
            recovery.RECEIPT[0], out, recovery.RECEIPT[1], recovery.RECEIPT[2],
            hashlib.sha256(), root)
    records = json.loads(manifest.read_text())['files']
    require(len(records) == 3032, 'source_manifest_record_count')
    for missing in MISSING:
        try:
            recovery.find_manifest_record(records, missing)
        except recovery.RecoveryError:
            pass
        else:
            raise RebuildError('previously_missing_original_input_now_present')
    archive = root / 'source.zip'
    hasher = hashlib.sha256()
    try:
        with archive.open('wb') as out:
            for aid, nbytes, sha in recovery.PARTS:
                recovery.stream_asset(aid, out, nbytes, sha, hasher, root)
        require(hasher.hexdigest() == recovery.ARCHIVE_SHA256, 'source_archive_sha')
        checked = {}
        with zipfile.ZipFile(archive) as z:
            members = z.infolist()
            names = [recovery.safe_zip_member(m) for m in members]
            require(len(names) == len(set(names)), 'duplicate_source_member_names')
            for path in SOURCE_FILES:
                entry = recovery.find_manifest_record(records, path)
                item = recovery.find_zip_member(members, path)
                checked[path] = recovery.extract_verified_member(
                    z, item, root / 'upstream' / path, entry)
        return checked
    finally:
        archive.unlink(missing_ok=True)
        manifest.unlink(missing_ok=True)


def group_support(df, only_positive=True):
    if only_positive:
        df = df[pd.to_numeric(df.share, errors='coerce').gt(0)]
    groups = df.groupby('family').agg(
        butterfly_species=('species', 'nunique'),
        regions=('region', 'nunique'))
    return sorted(groups.index[groups.butterfly_species.ge(10) & groups.regions.ge(10)].tolist())


def compare_column(left, right, name):
    la, rb = left.to_numpy(), right.to_numpy()
    require(len(la) == len(rb), 'column_length_mismatch')
    if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
        u, v = np.asarray(la, dtype='float64'), np.asarray(rb, dtype='float64')
        finite = np.isfinite(u) & np.isfinite(v)
        equal = (np.isnan(u) & np.isnan(v)) | (
            finite & np.isclose(u, v, rtol=1e-6, atol=1e-8))
        err = np.abs(u[finite] - v[finite])
        return {'matched': int(equal.sum()),
                'mismatched': int((~equal).sum()),
                'max_absolute_difference': float(err.max()) if len(err) else None}
    eq = left.astype('string').fillna('<MISSING>').eq(
        right.astype('string').fillna('<MISSING>'))
    return {'matched': int(eq.sum()), 'mismatched': int((~eq).sum())}


def reproduce(root, snapshot):
    upstream = root / 'upstream'
    snap = Path(snapshot)
    require(recovery.sha256_file(snap / 'analysis_source/build_forecast_panel.py') ==
            ORIGINAL_BUILD_SCRIPT_SHA, 'original_panel_builder_script_changed')
    selection_path = snap / 'forecast_all_225_disposition.csv'
    require(recovery.sha256_file(selection_path) == PUBLISHED_SELECTION_SHA,
            'original_published_selection_digest')
    selection = pd.read_csv(selection_path)
    eligible = selection[selection.stage_disposition.str.startswith('BASELINE_ELIGIBLE')]
    require(len(selection) == 225 and len(eligible) == 155 and
            eligible.taxon_key.nunique() == 155, 'eligible_species_cohort')

    cells = pd.read_parquet(upstream / SOURCE_FILES[0])
    env = pd.read_parquet(upstream / SOURCE_FILES[1])
    base = env[env.period.eq(0)].set_index('cell')
    cells = cells[(cells.center_region_matches == 1) &
                  (~cells.region_boundary_flag) &
                  cells.cell.map(base.unclassified).le(.05)].copy()
    grid = pd.read_parquet(upstream / SOURCE_FILES[2])
    cells = cells.merge(grid[['cell', 'climate_id', 'climate_lat', 'climate_lon']],
                        on='cell', validate='one_to_one')
    cells['spatial_block'] = (np.floor(cells.climate_lat / 5).astype(int).astype(str) +
                              '_' + np.floor(cells.climate_lon / 5).astype(int).astype(str))

    past = pd.read_parquet(upstream / SOURCE_FILES[5],
                           columns=['speciesKey', 'year', 'cell_x_10km', 'cell_y_10km'])
    past = past[past.year.le(2014)].drop_duplicates(
        ['speciesKey', 'cell_x_10km', 'cell_y_10km'])
    inv = Transformer.from_crs(6933, 4326, always_xy=True)
    targets = np.radians(cells[['latitude', 'longitude']].values)
    pieces = []
    for sp in eligible.itertuples(index=False):
        g = past[past.speciesKey.eq(sp.taxon_key)]
        require(len(g) > 0, 'baseline_species_has_no_location')
        lon, lat = inv.transform((g.cell_x_10km.values + .5) * 10000,
                                 (g.cell_y_10km.values + .5) * 10000)
        tree = BallTree(np.radians(np.column_stack([lat, lon])), metric='haversine')
        distances = tree.query(targets, k=1)[0].ravel() * 6371.0088
        chosen = cells[distances <= 500].copy()
        chosen['distance_km'] = distances[distances <= 500]
        chosen['speciesKey'] = sp.taxon_key
        chosen['species'] = sp.species
        pieces.append(chosen)
    domain = pd.concat(pieces, ignore_index=True)
    require(domain.speciesKey.nunique() == 155, 'domain_species')
    require(len(domain) == 280988, 'domain_species_cells')

    days = pd.read_parquet(upstream / SOURCE_FILES[3])
    days['period'] = ((days.year - 2011) // 4).astype(int)
    effort = pd.read_parquet(upstream / SOURCE_FILES[4])
    effort['period'] = ((effort.year - 2011) // 4).astype(int)
    denominator = effort.groupby(['cell', 'period']).size().rename('days').reset_index()
    numerator = days.groupby(['speciesKey', 'cell', 'period']).size().rename('success').reset_index()
    quarters = effort.assign(quarter=(effort.month.astype(int) - 1) // 3).groupby(
        ['cell', 'period', 'quarter']).size().unstack(fill_value=0).reindex(
        columns=range(4), fill_value=0)
    quarters = quarters.div(quarters.sum(axis=1), axis=0)
    quarters.columns = ['season_q' + str(i) for i in range(4)]
    quarters = quarters.reset_index()
    provider = effort.groupby(['cell', 'period', 'datasetKey']).size().rename('n').reset_index()
    dominance = (provider.groupby(['cell', 'period']).n.max() /
                 provider.groupby(['cell', 'period']).n.sum()).rename(
                 'provider_dominance').reset_index()

    env = (env.merge(denominator, on=['cell', 'period'], how='left', validate='one_to_one')
              .merge(quarters, on=['cell', 'period'], how='left', validate='one_to_one')
              .merge(dominance, on=['cell', 'period'], how='left', validate='one_to_one'))

    host_metrics = pd.read_parquet(upstream / SOURCE_FILES[6]).drop(
        columns=['composition_hash'])
    family = pd.read_parquet(upstream / SOURCE_FILES[7])
    independent_families = group_support(family, only_positive=True)
    alternate_families = group_support(family, only_positive=False)
    frozen_panel = pd.read_parquet(upstream / PANEL)
    frozen_family_cols = sorted(c[len('host_family_'):] for c in frozen_panel.columns
                                if c.startswith('host_family_'))
    chosen = independent_families
    if not chosen:
        raise RebuildError('no_host_family_candidates')
    fp = family[family.family.isin(chosen)].pivot(
        index=['species', 'region'], columns='family', values='share').fillna(0)
    fp.columns = ['host_family_' + str(c) for c in fp.columns]
    host_metrics = host_metrics.merge(
        fp.reset_index(), on=['species', 'region'], how='left', validate='one_to_one'
    ).fillna({c: 0 for c in fp.columns})
    domain = domain.merge(host_metrics, on=['species', 'region'], how='left',
                          validate='many_to_one')

    parts = []
    for period in (0, 1, 2):
        p = domain.copy()
        p['period'] = period
        p = p.merge(env, on=['cell', 'period'], how='left', validate='many_to_one')
        p = p.merge(numerator, on=['speciesKey', 'cell', 'period'], how='left',
                    validate='one_to_one')
        p['success'] = p.success.fillna(0)
        p['days'] = p.days.fillna(0)
        require(bool((p.success <= p.days).all()), 'success_exceeds_provider_days')
        p['recording_probability10'] = np.nan
        supported = p.days.ge(10)
        p.loc[supported, 'recording_probability10'] = hypergeom.sf(
            0, p.loc[supported, 'days'].astype(int),
            p.loc[supported, 'success'].astype(int), 10)
        parts.append(p)
    computed = pd.concat(parts, ignore_index=True)
    require(len(computed) == 842964, 'derived_panel_row_count')

    key = ['speciesKey', 'cell', 'period']
    require(not computed.duplicated(key).any() and
            not frozen_panel.duplicated(key).any(), 'duplicated_model_panel_keys')
    computed = computed.sort_values(key).reset_index(drop=True)
    frozen_panel = frozen_panel.sort_values(key).reset_index(drop=True)
    require(computed[key].equals(frozen_panel[key]), 'model_panel_key_mismatch')

    exclude = {'native_host_count', 'introduced_host_count', 'introduced_share'}
    common = [x for x in computed.columns if x in frozen_panel.columns and x not in exclude]
    stats = {x: compare_column(frozen_panel[x], computed[x], x) for x in common}
    unmatched = {k: v['mismatched'] for k, v in stats.items() if v['mismatched']}
    categories = {
        'space': ['distance_km', 'spatial_block', 'region', 'latitude', 'longitude'],
        'effort_and_response': ['days', 'success', 'recording_probability10', 'provider_dominance',
                                'season_q0', 'season_q1', 'season_q2', 'season_q3'],
        'original_upstream_host_metrics': list(host_metrics.columns),
    }
    return {
        'status': ('PARTIAL_PANEL_REPLAY_EXACT_ON_COMPARABLE_COLUMNS'
                   if not unmatched else 'PARTIAL_PANEL_REPLAY_WITH_COLUMN_DIFFERENCES'),
        'original_source_selected_species': len(eligible),
        'original_selected_source': 'frozen public 225 candidate disposition, stage labels only',
        'current_panel_rows': len(frozen_panel),
        'reconstructed_keyed_rows': len(computed),
        'species_by_cell_pairs': len(domain),
        'source_family_positive_selection_count': len(independent_families),
        'source_family_all_row_selection_count': len(alternate_families),
        'frozen_family_count': len(frozen_family_cols),
        'positive_family_set_matches_original': independent_families == frozen_family_cols,
        'unfiltered_family_set_matches_original': alternate_families == frozen_family_cols,
        'source_column_count_independently_computed': len(common),
        'columns_different_from_original': unmatched,
        'column_match_counts': {
            c: {n: stats[n] for n in fields if n in stats}
            for c, fields in categories.items()
        },
        'unreconstructed_original_upstream_columns': sorted(exclude),
        'raw_GBIF_occurrences_replayed': False,
        'daily_record_intermediates_reused': True,
        'original_host_count_csv_not_found': True,
        'copying_final_host_count_values_to_new_panel_prohibited': True,
        'complete_independent_raw_to_panel_replay': False,
        'no_source_rows_or_private_paths_exported': True,
        'method_boundary': ('Recomputed original cohort/domain, distance, provider-day '
                            'counts, target hypergeometric probabilities, environmental '
                            'joins and host-composition fields from hash-verified derived '
                            'source inputs. Three original host-count columns require the '
                            'missing independent source and are deliberately omitted.'),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        'schema': 'GEB_original155_partial_model_panel_upstream_replay_v1',
        'state': 'STARTED',
        'original_source_archive_sha256': recovery.ARCHIVE_SHA256,
        'original_panel_sha256': recovery.PANEL_SHA256,
        'expected_original_build_script_sha256': ORIGINAL_BUILD_SCRIPT_SHA,
        'full_raw_GBIF_to_panel_rebuilt': False,
    }
    try:
        with tempfile.TemporaryDirectory(prefix='original155-partial-replay-') as tmp:
            root = Path(tmp)
            verified = restore_upstream(root)
            require(verified[PANEL]['sha256'] == recovery.PANEL_SHA256,
                    'original_panel_hash_changed')
            require(verified[recovery.LOCK]['sha256'] == recovery.LOCK_SHA256,
                    'original_lock_hash_changed')
            result['source_inputs_sha256'] = {k: v['sha256'] for k, v in verified.items()}
            result['replay'] = reproduce(root, args.snapshot)
            result['state'] = 'COMPLETED'
    except RebuildError as exc:
        result['state'] = 'INCOMPLETE_WITH_SAFE_FAILURE'
        result['failure'] = str(exc)
    except recovery.RecoveryError:
        result['state'] = 'INCOMPLETE_WITH_SAFE_FAILURE'
        result['failure'] = 'private_source_integrity_or_read_failure'
    except Exception as exc:
        result['state'] = 'INCOMPLETE_WITH_SAFE_FAILURE'
        result['failure'] = 'source_rebuild_exception_' + type(exc).__name__
    finally:
        args.output.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    print('ORIGINAL155_PARTIAL_PANEL_REBUILD_' + result['state'], flush=True)
    if result.get('replay'):
        print('PARTIAL_REPLAY_SUMMARY', json.dumps({
            k: v for k, v in result['replay'].items()
            if k in ('status','current_panel_rows','source_column_count_independently_computed',
                     'columns_different_from_original','positive_family_set_matches_original')
        }), flush=True)
    return 0 if result['state'] == 'COMPLETED' else 1


if __name__ == '__main__':
    sys.exit(main())
