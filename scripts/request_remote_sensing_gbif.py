"""Submit the frozen 225-taxon 2021 source predicate; never print credentials."""
import os,json,base64,time,urllib.request
from pathlib import Path
out=Path('artifacts');out.mkdir(exist_ok=True)
predicate=json.loads(Path('data/exploratory/remote_sensing_20261029/predicate.json').read_text())
body={'creator':os.environ['GBIF_USERNAME'],'notificationAddresses':[],'sendNotification':False,'format':'SIMPLE_CSV','predicate':predicate}
token=base64.b64encode((os.environ['GBIF_USERNAME']+':'+os.environ['GBIF_PASSWORD']).encode()).decode()
req=urllib.request.Request('https://api.gbif.org/v1/occurrence/download/request',data=json.dumps(body).encode(),headers={'Content-Type':'application/json','Authorization':'Basic '+token})
with urllib.request.urlopen(req,timeout=120) as response:key=response.read().decode().strip().strip('"')
assert key and '/' not in key
(out/'download_key.txt').write_text(key+'\n')
print('Download key:',key,flush=True)
for _ in range(120):
    with urllib.request.urlopen('https://api.gbif.org/v1/occurrence/download/'+key,timeout=60) as response:meta=json.load(response)
    receipt={k:meta.get(k) for k in ['key','status','doi','downloadLink','totalRecords','size','created','modified','format']}
    (out/'remote_sensing_gbif_download.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('Status:',receipt['status'],flush=True)
    if receipt['status']=='SUCCEEDED':break
    if receipt['status'] in ('KILLED','CANCELLED','FAILED'):raise RuntimeError(receipt['status'])
    time.sleep(20)
else:raise RuntimeError('Download pending after 40 minutes; resume from recorded key, do not resubmit')
