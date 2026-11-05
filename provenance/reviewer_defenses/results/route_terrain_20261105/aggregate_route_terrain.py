"""Aggregate the locked DSM route design only after every tile has a receipt."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[2]
W = R / 'work/barriers'
O = R / 'outputs'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    lock_hash = digest(O / 'route_terrain_design_lock.json')
    design = json.loads((O / 'route_terrain_design_audit.json').read_text())
    assert design['lock_sha256'] == lock_hash
    assert design['source_mapping_sha256'] == digest(W / 'sea_route_mapping.parquet')
    tiles = pd.read_csv(O / 'route_terrain_required_tiles.csv')
    acquisition = O / 'route_terrain_acquisition.json'
    if not acquisition.exists():
        raise RuntimeError('Full acquisition pass not finished; no partial route features written')
    receipt = json.loads(acquisition.read_text())
    assert receipt['requested_tiles'] == len(tiles) and receipt['partial_limit'] == 0
    frames = []
    statuses = {}
    for t in tiles.itertuples(index=False):
        p = W / 'route_terrain_values' / (t.tile + '.json')
        d = json.loads(p.read_text())
        assert d['lock_sha256'] == lock_hash and d['points'] == t.route_land_sample_points
        statuses[d['status']] = statuses.get(d['status'], 0) + 1
        if d['status'] == 'ERROR':
            raise RuntimeError('Resolve failed retrieval before feature construction: ' + t.tile)
        assert d['status'] in ['COMPLETE', 'NO_TILE']
        data = p.with_suffix('.parquet')
        if 'data_sha256' in d:
            assert digest(data) == d['data_sha256']
        q = pd.read_parquet(data)
        assert len(q) == d['points']
        coords = q[['route_index', 'longitude', 'latitude']]
        actual = hashlib.sha256(pd.util.hash_pandas_object(coords, index=False).to_numpy().tobytes()).hexdigest()
        assert actual == d['point_hash']
        assert np.isfinite(q.elevation_m).sum() == d['valid_points']
        if d['status'] == 'NO_TILE':
            assert q.elevation_m.isna().all()
        frames.append(q[['route_index', 'elevation_m']])
    values = pd.concat(frames, ignore_index=True)
    support = pd.read_parquet(W / 'route_terrain_sampling_support.parquet')
    assert support.route_index.is_unique
    counts = values.groupby('route_index').size().reindex(support.route_index, fill_value=0)
    np.testing.assert_array_equal(counts.to_numpy(), support.land_samples.to_numpy())
    assert len(values) == design['planned_land_samples']
    grouped = values.groupby('route_index').elevation_m
    quantiles = grouped.quantile([.1, .9]).unstack().rename(columns={.1:'route_elevation_p10_m', .9:'route_elevation_p90_m'})
    support = support.join(quantiles, on='route_index')
    support['valid_land_samples'] = support.route_index.map(grouped.count()).fillna(0).astype(int)
    support['valid_land_fraction'] = support.valid_land_samples / support.land_samples.replace(0, np.nan)
    support['route_terrain_unknown'] = (~(support.valid_land_samples.ge(3) & support.valid_land_fraction.ge(.8))).astype(int)
    support['route_elevation_spread_m'] = support.route_elevation_p90_m - support.route_elevation_p10_m
    for raw, feature in [('route_elevation_p90_m','log_route_elevation_p90'), ('route_elevation_spread_m','log_route_elevation_spread')]:
        support[feature] = np.where(support.route_terrain_unknown.eq(0), np.log1p(support[raw].clip(lower=0)), 0)
    routes = pd.read_parquet(W / 'route_terrain_routes.parquet')[['route_key','route_index']]
    assert routes.route_key.is_unique and routes.route_index.is_unique
    metrics = routes.merge(support, on='route_index', validate='one_to_one')
    mapping = pd.read_parquet(W / 'sea_route_mapping.parquet')[['speciesKey','period','cell','route_key']]
    features = mapping.merge(metrics, on='route_key', validate='many_to_one')
    assert len(features) == len(mapping) and not features.duplicated(['speciesKey','period','cell']).any()
    assert np.isfinite(features[['log_route_elevation_p90','log_route_elevation_spread','route_terrain_unknown']]).all().all()
    metrics.to_parquet(W / 'route_terrain_metrics.parquet', index=False)
    features.to_parquet(W / 'route_terrain_features.parquet', index=False)
    audit = dict(status='FEATURES_VALIDATED_NOT_MODEL_RESULT', routes=len(metrics), feature_rows=len(features), supported_routes=int(metrics.route_terrain_unknown.eq(0).sum()), tile_statuses=statuses, lock_sha256=lock_hash, source_mapping_sha256=design['source_mapping_sha256'], feature_sha256=digest(W / 'route_terrain_features.parquet'), limits='Native DSM sampled every <=2km on geometric routes; not actual flight, least-cost resistance, historical elevation or maximum mountain barrier.')
    (O / 'route_terrain_feature_audit.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit))

if __name__ == '__main__':
    main()
