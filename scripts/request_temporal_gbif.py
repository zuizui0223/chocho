"""Submit the frozen 225-taxon 2001-2024 source predicate; never print credentials."""
import os,json,base64,time,urllib.request
from pathlib import Path
out=Path('artifacts');out.mkdir(exist_ok=True)
predicate=json.loads(Path('data/exploratory/temporal_20261029/predicate.json').read_text())
body={'creator':os.environ['GBIF_USERNAME'],'notificationAddresses':[],'sendNotification':False,'format':'SIMPLE_CSV','predicate':predicate}
token=base64.b64encode((os.environ['GBIF_USERNAME']+':'+os.environ['GBIF_PASSWORD']).encode()).decode()
req=urllib.request.Request('https://api.gbif.org/v1/occurrence/download/request',data=json.dumps(body).encode(),headers={'Content-Type':'application/json','Authorization':'Basic '+token})
with urllib.request.urlopen(req,timeout=120) as response:key=response.read().decode().strip().strip('"')
assert key and '/' not in key
(out/'download_key.txt').write_text(key+'\n')
print('Download key:',key,flush=True)
receipt={'key':key,'status':'REQUEST_ACCEPTED','source_predicate':'data/exploratory/temporal_20261029/predicate.json'}
(out/'temporal_gbif_download.json').write_text(json.dumps(receipt,indent=2)+'\n')
