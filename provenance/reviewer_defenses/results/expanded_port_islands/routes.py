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
portmap=pd.read_csv('outputs/port_buffer_landmass_audit.csv').query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1").copy()
portmap['landmass_id']=portmap.single_landmass_id.astype(int)
known=set(portmap.landmass_id)-{15,26,456}
geoms={}
shp=shapefile.Reader(str(w/'GSHHS_shp/h/GSHHS_h_L1.shp'))
for rec in shp.iterRecords(fields=['sibling_id']):
    i=int(rec[0])
    if i in known:geoms.setdefault(i,[]).append(shape(shp.shape(rec.oid).__geo_interface__))
polys=[(i,g.bounds,prep(g)) for i,gs in geoms.items() for g in gs]
inventory=portmap.set_index('port_id')
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
def twrows(f):
    if not is_zipfile(f):
        for enc in ['utf-8-sig','cp950']:
            try:return pd.read_csv(f,encoding=enc,dtype=str).values.tolist()
            except UnicodeDecodeError:pass
        data=pd.read_csv(f,encoding='cp950',encoding_errors='replace',dtype=str)
        assert data.columns[11]=='出發港' and data.columns[13]=='實際進港日'
        return data.values.tolist()
    with ZipFile(f) as z:
        if 'content.xml' not in z.namelist():
            data=pd.read_excel(f,engine='openpyxl',dtype=str)
            assert data.columns[7]=='船種' and data.columns[11]=='出發港'
            return data.values.tolist()
        xml=etree.fromstring(z.read('content.xml'))
    data=[]
    for row in xml.xpath('//*[local-name()="table-row"]')[1:]:
        vals=[]
        for c in row:
            key='{urn:oasis:names:tc:opendocument:xmlns:table:1.0}number-columns-repeated'
            vals.extend([''.join(c.itertext()).strip()]*min(18,int(c.get(key,1))))
            if len(vals)>=17:break
        if len(vals)>=14:data.append(vals)
    return data

# Binary connectivity only; published round-trip totals are never split into
# imputed directional counts. endpoint counts track observing coverage separately.
from collections import defaultdict
from io import StringIO
links=Counter()
nodes=set()
def call(src,year,a,b,kind):
    if a in known:nodes.add((src,kind,year,a,'out'))
    if b in known:nodes.add((src,kind,year,b,'in'))
    if a in known and b in known and a!=b:links[src,kind,year,int(a),int(b)]+=1
def roundtrip(src,year,a,b):
    call(src,year,a,b,'ferry')
    call(src,year,b,a,'ferry')

api='https://data.gov.tw/api/v2/rest/dataset/99641'
items=requests.get(api,timeout=30).json()['result']['distribution']
files={}
for x in items:
    m=re.match(r'(109|110|111|112)\D+(\d+)',x['resourceDescription'])
    if m:files[f'tw{int(m[1])+1911}{int(m[2]):02d}.csv']=x['resourceDownloadUrl']
assert len(files)==48
esfiles={f'atraques_{y}{m:02}.csv' for y in range(2020,2024) for m in range(1,13)}|{'atraques_202401.csv'}
allfiles={**files,**{name:'https://posidoniaweb.portsdebalears.com/docs/files/'+name for name in esfiles}}
with ThreadPoolExecutor(4) as pool:list(pool.map(lambda kv:fetch(*kv),allfiles.items()))
seen=set()
vessel_audit=Counter()
for name in sorted(files):
    if b'<html' in (cache/name).read_bytes()[:500].lower():
        print('Unavailable monthly table:',name)
        continue
    for r in twrows(cache/name):
        if len(r)<15:continue
        date=str(r[13]).strip().replace('/','').replace('-','')
        if not re.fullmatch('202[0-3][0-9]{4}',date):continue
        typ=str(r[7]).strip()
        if typ not in ['全貨櫃船','貨櫃船']:continue
        key=tuple(str(r[j]) for j in [0,1,11,13,14])
        if key in seen:continue
        seen.add(key)
        eng=re.split(r'[‧·\u4e00-\u9fff]',str(r[11]))[0].strip()
        a=island(unique.get(norm(eng),''))
        call('Taiwan official calls',int(date[:4]),a,44,'container')
        vessel_audit['TW',int(date[:4]),typ]+=1
seen=set()
for name in sorted(esfiles):
    with (cache/name).open(encoding='utf-8-sig',newline='') as stream:
        for r in csv.DictReader(stream,delimiter=';'):
            if r['Código puerto'] not in ['A','P']:continue
            typ=r['Tipo buque'].strip().upper()
            if 'FERRY' in typ or 'FERRI' in typ:kind='ferry'
            elif 'CONTENEDOR' in typ:kind='container'
            else:continue
            dates=re.findall(r'\b20[0-9]{2}\b',r['Fecha/Hora llegada'])
            if len(dates)!=1 or not 2020<=int(dates[0])<=2023:continue
            key=(r['Código puerto'],r['Escala'])
            if not key[1] or key in seen:continue
            seen.add(key);year=int(dates[0])
            call('Balearic completed calls',year,island(r['Locode origen']),152,kind)
            vessel_audit['ES',year,typ]+=1
# Two-terminal domestic lines among islands already in the botanical frame.
danish={'Fåborg-Søby':(169,1577),'Svendborg-Ærøskøbing':(169,1577),
        'Søby-Fynshav':(1577,715),'Søby-Mommark':(1577,715),'Bøjden-Fynshav':(169,715),
        'Køge-Rønne':(98,454),'København-Rønne':(98,454),'Kastrup-Rønne':(98,454)}
u='https://api.statbank.dk/v1/data/SKIB31/CSV?lang=en&F%C3%86RGE=*&ENHED=8000&Tid=2020,2021,2022,2023'
dk=pd.read_csv(fetch('skib31_2020_23.csv',u),sep=';')
for route,(a,b) in danish.items():
    g=dk[dk.iloc[:,0].eq(route)]
    assert set(g.TID)==set(range(2020,2024))
    for r in g.itertuples():
        n=pd.to_numeric(r.INDHOLD,errors='coerce')
        if n>0:roundtrip('Statistics Denmark SKIB31',r.TID,a,b)
dk32=pd.read_csv(fetch('skib32_2020_23.csv',u.replace('SKIB31','SKIB32')),sep=';')
for route,(a,b) in {'Rødby Færgehavn-Puttgarden':(276,997),'Rønne-Mukran':(454,321),'Rønne-Fährhafen Sassnitz':(454,321)}.items():
    g=dk32[dk32.iloc[:,0].eq(route)]
    assert set(g.TID)==set(range(2020,2024))
    for r in g.itertuples():
        if pd.to_numeric(r.INDHOLD,errors='coerce')>0:roundtrip('Statistics Denmark SKIB32',r.TID,a,b)
# Annual report p.6: the operator's entire route system links Cape Breton
# (North Sydney) with Newfoundland (Port aux Basques or Argentia).
for end,n in [(2021,1498),(2022,1642),(2023,1750)]:
    roundtrip('Marine Atlantic FY ending March',end,79,21)
    vessel_audit['CA',end,'single crossings, both directions combined']+=n
# MLIT Hokkaido: Tomakomai--Oarai, observed passenger totals FY2020--FY2022.
for end,n in [(2021,97548),(2022,114142),(2023,154679)]:
    roundtrip('MLIT Hokkaido FY ending March',end,27,13)
    vessel_audit['JP',end,'Tomakomai-Oarai passengers, both directions combined']+=n
edges=pd.DataFrame([dict(source=s,kind=k,year=y,a=a,b=b) for s,k,y,a,b in sorted(links)])
coverage=pd.DataFrame(sorted(nodes),columns=['source','kind','year','island','direction'])
lng=out/'through_2023'
edges=pd.concat([edges,pd.read_csv(lng/'links.csv').drop(columns='voyages').assign(source='LNG-T3',kind='LNG')],ignore_index=True)
coverage=pd.concat([coverage,pd.read_csv(lng/'endpoint_coverage.csv').drop(columns='voyages').assign(source='LNG-T3',kind='LNG')],ignore_index=True)
edges.to_csv(out/'links.csv',index=False);coverage.to_csv(out/'endpoint_coverage.csv',index=False)
pd.DataFrame([(*k,n) for k,n in sorted(vessel_audit.items())],columns=['source','year','original_type','n']).to_csv(out/'source_types.csv',index=False)
(out/'sources.json').write_text(json.dumps(sorted(sources,key=lambda x:x['file']),indent=2))
print(edges.groupby(['source','kind']).size())
