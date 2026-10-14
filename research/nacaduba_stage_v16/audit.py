"""Post-hoc flowering-census sensitivity; NOT a butterfly feeding or colonization test.

Source: Kato et al. 2012 Table 4/5/7, pinned at research/nacaduba_stage_v15/data/
"""
from pathlib import Path
import csv, itertools, json, hashlib
ROOT=Path(__file__).resolve().parents[1]/"nacaduba_stage_v15"/"data"
MONTHS="Apr May Jun Jul Aug Sep Oct Nov Dec Jan Feb Mar".split()
PARTS="early mid late".split()
YEARS=(2008,2009,2010)
EXPECTED={
  "ramets_kato2012_table7.csv":"ed8425db6cbc6830249824a2edd41fbfc7ec2db3",
  "clones_kato2012_table4_table5.csv":"47f30f98807f27c3f47f8797731f7e05cc9eca9f"
}
def audit():
  for f,h in EXPECTED.items():
    b=(ROOT/f).read_bytes()
    assert hashlib.sha1((f"blob {len(b)}\0").encode()+b).hexdigest()==h
  def read(f):
    with (ROOT/f).open(newline="") as s:return list(csv.DictReader(s))
  clones=read("clones_kato2012_table4_table5.csv")
  data=read("ramets_kato2012_table7.csv")
  assert len(clones)==23 and len(data)==36 and sum(int(c["ramets"]) for c in clones)==39
  ns={int(r["clone"]):int(r["ramets"]) for r in clones}
  prov={int(r["clone"]):r["seed_source"] for r in clones}
  masks={c:0 for c in ns}
  ramet_masks=[]
  for r in data:
    cid=int(r["clone"]);bits=0
    assert prov[cid]==r["seed_source"]
    for yi,y in enumerate(YEARS):
      b=r[f"begin_{y}"];n=r[f"duration_{y}"]
      assert bool(b)==bool(n)
      if not b:continue
      part,mon=b.split()
      first=36*yi+3*MONTHS.index(mon)+PARTS.index(part)
      for i in range(first,first+int(n)):bits|=1<<i
    masks[cid]|=bits
    ramet_masks.append((cid,bits))
  total=0
  for v in masks.values():total|=v
  assert total.bit_count()==57
  bysource={p:0 for p in "abcdefgh"}
  for cid,bits in masks.items():bysource[prov[cid]]|=bits
  pairs=[(a,b,(bysource[a]|bysource[b]).bit_count()) for a,b in itertools.combinations(bysource,2)]
  g=frozenset(cid for cid in ns if prov[cid]=="g")
  gh=frozenset(cid for cid in ns if prov[cid] in ("g","h"))
  assert g=={146,147,148,149} and len(gh)==7
  singles=[cid for cid,n in ns.items() if n==1]
  doubles=[cid for cid,n in ns.items() if n==2]
  assert len(singles)==7 and len(doubles)==16
  def num_bins_removed(subset):
    r=0
    for c,m in masks.items():
      if c not in subset:r|=m
    return (total&~r).bit_count()
  g_unique=num_bins_removed(g)
  rare=sum(num_bins_removed(frozenset(aa+bb))>=g_unique
    for aa in itertools.combinations(singles,2) for bb in itertools.combinations(doubles,2))
  success=0
  for aa in itertools.combinations(singles,2):
    for bb in itertools.combinations(doubles,5):
      m=0
      for c in aa+bb:m|=masks[c]
      success+=(m==total)
  triples=[]
  for n in (1,2,3):
    triples=[comb for comb in itertools.combinations(masks,n)
             if _union(comb,masks)==total]
    if triples:break
  def flowering_density(subset):
    values=[sum(bool(m&(1<<t)) for cid,m in ramet_masks if cid in subset) for t in range(108)]
    return {"at_least_1":sum(z>=1 for z in values),
            "at_least_3":sum(z>=3 for z in values),
            "at_least_5":sum(z>=5 for z in values),
            "flowering_ramet_decads":sum(values)}
  metrics={name:flowering_density(s) for name,s in [
        ("all",frozenset(ns)),("g_plus_h",gh),("three_clones",frozenset([146,148,152]))]}
  return {
    "full_calendar_binyears":total.bit_count(),
    "unique_g_binyears":g_unique,
    "season_unique_g_binyears":[num_bins_removed_year(g,yi,masks,total) for yi in range(3)],
    "two_provenance_pairs_complete":[[a,b] for a,b,n in pairs if n==57],
    "three_clone_full_calendar_subsets":[list(comb) for comb in triples],
    "g_plus_h_ramets":sum(ns[c] for c in gh),
    "three_clones_ramets":sum(ns[c] for c in [146,148,152]),
    "g_matched_size_subsets_total":2520,
    "g_matched_size_subsets_as_extreme":rare,
    "seven_clone_matched_portfolios_total":91728,
    "seven_clone_matched_portfolios_complete":success,
    "concurrent_flowering_proxy":metrics,
    "butterfly_larval_success_effect":None,
    "flower_bud_mass":None,
    "biological_colonization_lag":None,
    "causal_p_value":None
  }
def _union(ids,masks):
  m=0
  for cid in ids:m|=masks[cid]
  return m
def num_bins_removed_year(subset,yi,masks,total):
  others=_union((cid for cid in masks if cid not in subset),masks)
  win=((1<<36)-1)<<(yi*36)
  return (total&~others&win).bit_count()
if __name__=="__main__":print(json.dumps(audit(),indent=2,ensure_ascii=False))
