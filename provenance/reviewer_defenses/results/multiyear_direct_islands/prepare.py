"""Assemble directed arrivals by year. Local raw downloads stay outside the repository."""
from pathlib import Path
from collections import Counter
from zipfile import ZipFile,is_zipfile
import hashlib, re, json, unicodedata, csv
import requests, pandas as pd, numpy as np, shapefile
from lxml import etree
from shapely.geometry import Point,shape
from shapely.prepared import prep
from concurrent.futures import ThreadPoolExecutor

out=Path(__file__).resolve().parent
w=Path('work/barriers')
cache=w/'multiyear_calls'
cache.mkdir(exist_ok=True)
sources=[]
def fetch(name,url):
    f=cache/name
    if not f.exists():
        r=requests.get(url,timeout=60);r.raise_for_status();f.write_bytes(r.content)
    sources.append(dict(file=name,url=url,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    return f
def norm(v):
    return re.sub('[^A-Z0-9]','',unicodedata.normalize('NFKD',str(v)).encode('ascii','ignore').decode().upper())
pin='0fb7229b64476966c67cb4b90f000d554154591c'
loc=pd.read_csv(fetch('locodes.csv','https://raw.githubusercontent.com/datasets/un-locode/'+pin+'/data/code-list.csv'),dtype=str,keep_default_na=False)
loc=loc[loc.Function.str.startswith('1')].copy()
loc['code']=loc.Country+loc.Location
loc['normal']=loc.NameWoDiacritics.map(norm)
codes=loc.drop_duplicates('code').set_index('code').Coordinates.to_dict()
unique=loc.drop_duplicates(['code','normal'])
unique=unique[~unique.normal.duplicated(keep=False)].set_index('normal').code.to_dict()
known=set(pd.read_csv('outputs/island_port_early_inventory_status.csv').landmass_id)-{15,26,456}
geoms={}
shp=shapefile.Reader(str(w/'GSHHS_shp/h/GSHHS_h_L1.shp'))
for rec in shp.iterRecords(fields=['sibling_id']):
    i=int(rec[0])
    if i in known:geoms.setdefault(i,[]).append(shape(shp.shape(rec.oid).__geo_interface__))
polys=[(i,g.bounds,prep(g)) for i,gs in geoms.items() for g in gs]
inventory=pd.read_csv('outputs/island_port_early_inventory_status.csv').set_index('port_id')
coords=[]
for f in json.loads((w/'ports_2019.geojson').read_text(encoding='utf-8'))['features']:
    prop=f['properties']
    if prop['OBJECTID'] in inventory.index:
        coords.append((prop['COUNTRY'],inventory.at[prop['OBJECTID'],'landmass_id'],*f['geometry']['coordinates']))
coords=pd.DataFrame(coords,columns=['country','id','lon','lat'])
from functools import lru_cache
@lru_cache(None)
def island(code):
    c=re.fullmatch(r'(\d{2})(\d{2})([NS])\s+(\d{3})(\d{2})([EW])',codes.get(code,'').strip())
    if not c:return None
    a,b,h,x,y,k=c.groups()
    point=Point((int(x)+int(y)/60)*(1 if k=='E' else -1),(int(a)+int(b)/60)*(1 if h=='N' else -1))
    hits={i for i,(l,d,r,u),g in polys if l<=point.x<=r and d<=point.y<=u and g.covers(point)}
    adjacent=coords[coords.country.eq(code[:2])]
    dx=np.radians(adjacent.lon-point.x);dy=np.radians(adjacent.lat-point.y)
    arc=np.sin(dy/2)**2+np.cos(np.radians(adjacent.lat))*np.cos(np.radians(point.y))*np.sin(dx/2)**2
    hits.update(adjacent.loc[12742.0176*np.arcsin(np.sqrt(arc.clip(0,1)))<=5,'id'])
    return hits.pop() if len(hits)==1 else None
edges=Counter()
def add(source,yr,origin,destination,a,b):
    if a in known and b in known and a!=b:
        edges[source,yr,origin,destination,a,b]+=1
# Balearic completed berths; code and arrival date identify the inbound call.
months=[f'{y}{m:02}' for y in [2020,2021] for m in range(1,13)]+['202201']
paths=[fetch('atraques_'+m+'.csv','https://posidoniaweb.portsdebalears.com/docs/files/atraques_'+m+'.csv') for m in months]
calls=[]
for f in paths:
    with f.open(encoding='utf-8-sig',newline='') as fp:
        for r in csv.DictReader(fp,delimiter=';'):
            calls.append([r['Código puerto'],r['Locode origen'],r['Fecha/Hora llegada'],r['Escala']])
d=pd.DataFrame(calls,columns=['port','origin','arrival','call']).dropna()
d['arrival']=pd.to_datetime(d.arrival,dayfirst=True,format='mixed',errors='raise')
d=d.sort_values('arrival').drop_duplicates(['port','call'])
for r in d.itertuples():
    if r.port in ('P','A') and r.arrival.year in (2020,2021):
        add('Balearic',r.arrival.year,r.origin,{'P':'ESPMI','A':'ESALD'}[r.port],island(r.origin),152)
# Taiwan's international commercial-port archive covers the main island.
# The home port field is vessel registry, so retain island-level destination only.
u='https://data.gov.tw/api/v2/rest/dataset/99641'
items=requests.get(u,timeout=30).json()['result']['distribution']
ts=[]
for x in items:
    m=re.match(r'(109|110)\D+(\d+)',x['resourceDescription'])
    if m:ts.append((f'tw{int(m[1])+1911}{int(m[2]):02}.csv',x['resourceDownloadUrl']))
assert len(set(n for n,_ in ts))==24
def twrows(f):
    if not is_zipfile(f):
        for enc in ['utf-8-sig','cp950']:
            try:return pd.read_csv(f,encoding=enc,dtype=str).values.tolist()
            except UnicodeDecodeError:pass
        data=pd.read_csv(f,encoding='cp950',encoding_errors='replace',dtype=str)
        assert data.columns[11]=='出發港' and data.columns[13]=='實際進港日'
        return data.values.tolist()
    with ZipFile(f) as z:xml=etree.fromstring(z.read('content.xml'))
    data=[]
    for row in xml.xpath('//*[local-name()="table-row"]')[1:]:
        vals=[]
        for c in row:
            key='{urn:oasis:names:tc:opendocument:xmlns:table:1.0}number-columns-repeated'
            vals.extend([''.join(c.itertext()).strip()]*min(18,int(c.get(key,1))))
            if len(vals)>=17:break
        if len(vals)>=14:data.append(vals)
    return data
seen=set()
for name,url in sorted(ts):
    for vals in twrows(fetch(name,url)):
        date=str(vals[13]).replace('/','').replace('-','').strip()
        if not re.fullmatch('202[01][0-9]{4}',date):continue
        eng=re.split(r'[‧·\u4e00-\u9fff]',str(vals[11]))[0].strip()
        code=unique.get(norm(eng),'')
        key=tuple(str(vals[i]) for i in [0,1,11,13,14] if len(vals)>i)
        if key in seen:continue
        seen.add(key)
        add('Taiwan',int(date[:4]),code,'Taiwan international ports',island(code),44)
# US domestic AIS-derived port-pair records, year-specific origin and arrival,
# Hawaii port names checked within state. Aggregate route entries, not vessel counts.
us=w/'freight2018_2022.xlsx'
if not us.exists():us=fetch(us.name,'https://ndownloader.figshare.com/files/46387069')
meta=requests.get('https://api.figshare.com/v2/articles/25846111/versions/1',timeout=30).json()
assert hashlib.md5(us.read_bytes()).hexdigest()==meta['files'][0]['computed_md5']
sources.append(dict(file=us.name,url='https://doi.org/10.6084/m9.figshare.25846111.v1',sha256=hashlib.sha256(us.read_bytes()).hexdigest()))
haw={'HONOLULU':242,'PEARL HARBOR':242,'BARBERS POINT':242,'HILO':81,'KAWAIHAE':81,'KAHULUI':214,'LAHAINA':214,'NAWILIWILI':257,'PORT ALLEN':257}
uscache=cache/'us_2020_21_state_port_pairs.csv'
if not uscache.exists():
    counts=Counter()
    with ZipFile(us) as z:
        shared=[''.join(s.itertext()) for s in etree.fromstring(z.read('xl/sharedStrings.xml'))]
        with z.open('xl/worksheets/sheet1.xml') as fp:
            for _,row in etree.iterparse(fp,events=('end',),tag='{*}row'):
                vals={}
                for c in row[:6]:
                    col=re.sub('[0-9]','',c.get('r'));v=c.findtext('{*}v','')
                    vals[col]=shared[int(v)] if c.get('t')=='s' else v
                if vals.get('A') in ['2020','2021']:counts[tuple(vals.get(i,'') for i in 'ACDEF')]+=1
                row.clear()
                while row.getprevious() is not None:del row.getparent()[0]
    pd.DataFrame([list(k)+[v] for k,v in counts.items()],columns=['year','origin','ostate','dest','dstate','n']).to_csv(uscache,index=False)
for r in pd.read_csv(uscache).itertuples():
    if r.ostate==r.dstate=='Hawaii' and r.origin in haw and r.dest in haw and haw[r.origin]!=haw[r.dest]:
        edges['US',r.year,r.origin,r.dest,haw[r.origin],haw[r.dest]]+=r.n
e=pd.DataFrame([list(k)+[v] for k,v in edges.items()],columns=['source','year','origin','destination','a','b','source_rows'])
e.sort_values(['source','year','a','b','origin']).to_csv(out/'verified_links.csv',index=False)
(out/'source_receipts.json').write_text(json.dumps(sources,indent=2)+'\n')
print(e.groupby(['source','year']).agg(origin_islands=('a','nunique'),destinations=('b','nunique'),rows=('source_rows','sum')))
