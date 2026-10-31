"""Retain named hosts and conservative regional status; no arrival inference."""
from pathlib import Path
import json
import pandas as pd
import geopandas as gpd
import shapely

R = Path(__file__).resolve().parents[2]
W, O = R/'work/barriers', R/'outputs'
S = R/'work/remote_sensing/host_source/extracted/inputs'
names = pd.read_csv(S/'original/insect_host_accepted.csv', dtype=str)[['accepted_plant_name_id', 'accepted_name']].drop_duplicates()
mapping = dict(names.itertuples(index=False, name=None))
for butterfly in json.loads((S/'bce/accepted_host_edges.json').read_text())['per_butterfly']:
    for item in butterfly['details']:
        if item.get('match_status') == 'ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY':
            mapping[str(item['accepted_plant_name_id'])] = item['accepted_name']

def pairs(a, b):
    q = pd.concat([pd.read_csv(S/a, dtype=str), pd.read_csv(S/b, dtype=str)])
    return set(zip(q.accepted_plant_name_id.map(mapping), q.area_code_l3))

native = pairs('original/native_extant_nondoubtful_wgsrpd3.csv', 'bce/additional_native_wgsrpd3.csv')
present = pairs('original/extant_nondoubtful_wgsrpd3.csv', 'bce/additional_contemporary_wgsrpd3.csv')
q = pd.read_parquet(W/'host_flora_name_join_audit.parquet')
q = q[q.join_status.eq('SPECIES_FACET_COMPATIBLE')]
assert not q.accepted_gbif_species_key.duplicated().any()
keynames = dict(zip(q.accepted_gbif_species_key.astype(int), q.host_name))
regions = gpd.read_file(S/'level3.geojson').to_crs(4326)
rows, audits = [], []
for f in sorted((W/'port_flora_inventory').glob('*.json')):
    try:
        a = json.loads(f.read_text(encoding='utf-8'))
    except json.JSONDecodeError:
        continue  # acquisition may be writing this file
    if a['status'] != 'COMPLETE_FACETS':
        continue
    geom = shapely.from_wkt(a['geometry'])
    ix = regions.sindex.query(geom, predicate='intersects')
    codes = sorted(set(regions.iloc[ix].LEVEL3_COD.dropna()))
    audits.append(dict(port_id=a['port_id'], candidate_regions=';'.join(codes), region_count=len(codes)))
    for item in a['facets']['SPECIES_KEY']:
        key = int(item['name'])
        if key not in keynames:
            continue
        name = keynames[key]
        states = ['native_in_catalogue' if (name,c) in native else
                  'contemporary_not_native_in_catalogue' if (name,c) in present else
                  'unknown_catalogue_status' for c in codes]
        # A buffer may cross regional boundaries: only unanimous status is usable.
        status = states[0] if states and len(set(states)) == 1 else (
            'ambiguous_across_regions' if states else 'no_region_intersection')
        rows.append(dict(port_id=a['port_id'], host_name=name, gbif_species_key=key,
                         plant_records=item['count'], candidate_regions=';'.join(codes),
                         regional_catalogue_status=status, source_response_sha256=a['response_sha256']))
out = pd.DataFrame(rows)
assert not out.duplicated(['port_id','gbif_species_key']).any()
assert out.plant_records.gt(0).all()
out.to_parquet(W/'port_named_host_regional_status.parquet', index=False)
pd.DataFrame(audits).to_csv(O/'port_flora_regional_support.csv', index=False)
summary = out.groupby('regional_catalogue_status').agg(port_host_pairs=('host_name','size'),
    ports=('port_id','nunique'), host_species=('gbif_species_key','nunique'),
    plant_records=('plant_records','sum')).reset_index()
summary.to_csv(O/'port_named_host_status_summary.csv', index=False)
receipt = dict(status='PARTIAL_INVENTORY_REGIONAL_STATUS_NOT_CONFIRMED_INTRODUCTION',
    complete_ports_read=len(audits), ports_with_hosts=int(out.port_id.nunique()),
    host_species=int(out.gbif_species_key.nunique()), port_host_pairs=len(out),
    status_counts=out.regional_catalogue_status.value_counts().to_dict(),
    limits=['regional contemporary-minus-native is a non-native candidate, not independent introduced-status evidence',
            'regional status does not establish introduction through a port or introduction year',
            'pooled 2001-2024 species facets do not date individual hosts',
            'region determined from full query polygon, not port centre; multi-region discordance retained',
            'port buffers may include multiple islands; no island assignment of pooled records',
            'snapshot is partial and acquisition order is geographically structured',
            'record counts measure recording, not plant abundance; overlapping ports may share records'])
(O/'port_named_host_status_audit.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
print(json.dumps(receipt), flush=True)
