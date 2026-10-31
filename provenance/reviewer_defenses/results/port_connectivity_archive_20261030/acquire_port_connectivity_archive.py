from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,time
import requests,pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs';D=W/'unctad_plsciq_archive';D.mkdir(exist_ok=True)
endpoint='https://api.db.nomics.world/v22/series/UNCTAD/PLSCIQ'
docs=[];receipts=[];offset=0;meta=None
while True:
    file=D/f'offset_{offset}.json'
    if file.exists():raw=file.read_bytes()
    else:
        res=requests.get(endpoint,params={'limit':500,'offset':offset,'observations':1},timeout=60)
        res.raise_for_status();raw=res.content;file.write_bytes(raw);time.sleep(1)
    j=json.loads(raw);dataset=j['dataset'];batch=j['series']['docs'];total=j['series']['num_found']
    if meta is None:meta=dataset
    assert dataset['dir_hash']==meta['dir_hash'] and total==meta['nb_series']
    assert j['series']['offset']==offset and batch
    docs.extend(batch);receipts.append(dict(offset=offset,rows=len(batch),sha256=hashlib.sha256(raw).hexdigest()))
    offset+=len(batch)
    if offset>=total:break
assert len(docs)==total and len({d['series_code'] for d in docs})==total
rows=[]
for d in docs:
    assert d['dimensions']['frequency']=='Q'
    label=meta['dimensions_values_labels']['port'][d['dimensions']['port']]
    country,port=label.split(', ',1)
    assert len(d['period'])==len(d['value'])
    for period,value in zip(d['period'],d['value']):
        rows.append(dict(series_code=d['series_code'],port_code=d['dimensions']['port'],country_label=country,
                         port_label=port,period=period,year=int(period[:4]),quarter=int(period[-1]),
                         value=pd.to_numeric(value,errors='coerce')))
p=pd.DataFrame(rows);assert not p.duplicated(['series_code','period']).any()
p.to_parquet(W/'unctad_port_connectivity_2006_2021.parquet',index=False)
p[['series_code','port_code','country_label','port_label']].drop_duplicates().to_csv(O/'port_connectivity_archive_catalogue.csv',index=False)
audit=dict(status='ARCHIVE_ACQUIRED_NOT_JOINED_OR_FITTED',primary_provider='UNCTAD',retrieval_provider='DBnomics',
           primary_page='https://unctadstat.unctad.org/datacentre/dataviewer/US.PLSCI',primary_access='web fetch returned403',
           archive_endpoint=endpoint,provider_updated_at=meta['updated_at'],archive_indexed_at=meta['indexed_at'],
           archive_dir_hash=meta['dir_hash'],retrieved_at=datetime.now(timezone.utc).isoformat(),
           series=len(docs),rows=len(p),nonmissing=int(p.value.notna().sum()),first_period=p.period.min(),last_period=p.period.max(),
           measure=meta['dimensions_values_labels']['measure'],pages=receipts,
           license_status='Primary usage terms not yet directly verified; raw archive kept local, not redistributed',
           limits=['older2021 archived edition; do not mix index base with current series',
                   'container liner connectivity, not cargo volume, passenger traffic or actual port-pair route',
                   'NA remains missing, not zero; port identity requires country-aware validation'])
(O/'port_connectivity_archive_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit))
