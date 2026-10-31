"""Connect a literature transport case to frozen local data, without imputing routes."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs'
INPUTS = {
    'hosts': OUT / 'port_alien_top_hosts_previous_relationship.csv',
    'forecast': OUT / 'port_alien687_previous_forecast_overlap.csv',
    'ports': OUT / 'port_alien_host_broad_all_ports.csv',
    'edges': ROOT / 'work/barriers/all225_named_host_edges_local.parquet',
    'named': ROOT / 'work/barriers/port_named_host_regional_status.parquet',
    'dates': ROOT / 'work/barriers/host_butterfly_first_record_pairs.parquet',
    'catalogue': ROOT / 'work/barriers/ports_2019.geojson',
}
hosts = pd.read_csv(INPUTS['hosts'])
forecast = pd.read_csv(INPUTS['forecast'])
ports = pd.read_csv(INPUTS['ports'], keep_default_na=False)
edges = pd.read_parquet(INPUTS['edges'])
named = pd.read_parquet(INPUTS['named'])
dates = pd.read_parquet(INPUTS['dates'])
catalogue = json.loads(INPUTS['catalogue'].read_text(encoding='utf-8'))
nelson = [x['properties'] for x in catalogue['features']
          if x['properties']['PORT_NAME'] == 'NELSON' and x['properties']['COUNTRY'] == 'NZ']
assert len(nelson) == 1
port_id = nelson[0]['OBJECTID']
species = 'Pieris brassicae'
selected_hosts = ['Tropaeolum majus', 'Brassica oleracea']
selected = hosts[hosts.host_name.isin(selected_hosts)].copy()
assert len(selected) == 2
selected['host_edge_to_pieris_brassicae'] = selected.host_name.map(
    lambda name: bool(((edges.species == species) & (edges.host_name == name)).any()))
assert selected.host_edge_to_pieris_brassicae.all()
selected.to_csv(OUT / 'nelson_transport_case_host_links.csv', index=False)
f = forecast[forecast.species == species].iloc[0]
p = ports[ports.port_id == port_id]
assert len(p) == 1
n = named[named.port_id == port_id]
date_pairs = dates[(dates.species == species) & (dates.Region == 'New Zealand') & dates.host_name.isin(selected_hosts)]
date_pairs[['species', 'Region', 'host_name', 'FirstRecord_orig_butterfly',
            'FirstRecord_orig_plant', 'ordering', 'Source_butterfly', 'Source_plant']].to_csv(
    OUT / 'nelson_transport_case_country_dates.csv', index=False)
audit = {
    'status': 'LITERATURE_CASE_LINKED_NOT_ROUTE_INFERRED',
    'port': nelson[0],
    'port_query_complete': str(p.iloc[0].query_complete),
    'local_named_host_rows_at_port': len(n),
    'forecast_mse_without_hosts': float(f.mse_without_hosts),
    'forecast_mse_with_hosts': float(f.mse_with_hosts),
    'forecast_relative_mse_improvement_pct': 100 * (1 - float(f.mse_with_hosts) / float(f.mse_without_hosts)),
    'country_date_pair_rows': len(date_pairs),
    'pathway_values_in_screened_pairs': {
        col: dates[col].fillna('MISSING').value_counts().to_dict()
        for col in ['Pathway_butterfly', 'Pathway_plant']},
    'sources': {
        'field_study': 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236791',
        'agency': 'https://www.doc.govt.nz/nature/pests-and-threats/animal-pests-and-threats/great-white-butterfly/'},
    'source_accessed': '2026-10-30',
    'claim_limits': [
        'Shipping-container introduction is probable, not an observed original interception.',
        'Local host use is documented; global port host counts are not Nelson observations.',
        'Eradication in 2016 prevents treating this as uninterrupted expansion.',
        'Country first records do not date host establishment in Nelson.',
        'No island-to-island relay or joint transport is demonstrated.',
        'Case selected for documented pathway relevance, not positive forecast performance.'],
    'input_sha256': {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in INPUTS.items()},
}
(OUT / 'nelson_transport_case_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(audit, ensure_ascii=True, indent=2))
