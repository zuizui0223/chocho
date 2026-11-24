#!/usr/bin/env python3
"""Hash-verify newly added 3-input recovery bundle; export aggregate metadata only."""
import csv, hashlib, io, json, sys, tempfile, zipfile, traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import recover_GEB_v05_original155_private_inputs_20261122 as orig

ASSETS={"zip":(785479033,726857,"50e659381bcd5628cd2110361f2fe9a6ba895ffd8fb1a2a1dce1655235f8a6c6"),
        "receipt":(785479036,1302,"9c76155c4765458ef1e51d31e5593830e4dbea99455ef14720d69482fca5e36b")}
PATHS=[
"outputs/forecast_candidate_disposition.csv",
"outputs/temporal_host_family_support.csv",
"work/remote_sensing/host_region_counts.csv",
]
COLUMNS=[
("species","taxon_key","stage_disposition"),
("family","butterfly_species","regions"),
("species","region","native_host_count","introduced_host_count","introduced_share"),
]

def verify():
  meta=orig.gh_json(f"repos/{orig.REPO}/releases/{orig.RELEASE_ID}")
  assert meta.get("draft") is True
  files={a["id"]:a for a in meta["assets"]}
  for id_,size,sha in ASSETS.values():
    a=files[id_]
    assert a["size"]==size and a["digest"]=="sha256:"+sha
  with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)
    paths={}
    for kind,(id_,size,sha) in ASSETS.items():
      path=root/(kind+".bin")
      with path.open("wb") as stream:
        orig.stream_asset(id_,stream,size,sha,hashlib.sha256(),root)
      paths[kind]=path
    receipt=json.loads(paths["receipt"].read_text())
    receipt_text=json.dumps(receipt,sort_keys=True)
    summary={}
    member_count=0
    with zipfile.ZipFile(paths["zip"]) as archive:
      members=[m for m in archive.infolist() if not m.is_dir()]
      member_count=len(members)
      assert 3 <= len(members) <= 30
      allpaths=[orig.safe_zip_member(m) for m in members]
      assert len(set(allpaths))==len(allpaths)
      for i,name in enumerate(PATHS):
        found=[m for m in members if m.filename==name or m.filename.endswith("/"+name)]
        assert len(found)==1 and found[0].file_size < 20000000
        raw=archive.read(found[0])
        sha=hashlib.sha256(raw).hexdigest()
        assert sha in receipt_text and Path(name).name in receipt_text
        tab=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        assert tab and set(COLUMNS[i]).issubset(tab[0]) and all(None not in row for row in tab)
        entry={"sha256":sha,"bytes":len(raw),"rows":len(tab),"columns":list(tab[0])}
        if i==0:
          selected=[row for row in tab if row["stage_disposition"].startswith("BASELINE_ELIGIBLE")]
          assert len(tab)==225 and len(selected)==155
          assert len({row["taxon_key"] for row in selected})==155
          entry["selected_species"]=155
        if i==1:
          families=[row for row in tab if float(row["butterfly_species"])>=10 and float(row["regions"])>=10]
          assert len(families)==18
          entry["selected_families"]=18
        if i==2:
          pairs={(row["species"],row["region"]) for row in tab}
          assert len(pairs)==len(tab)
          missing=0
          for row in tab:
            assert float(row["native_host_count"])>=0 and float(row["introduced_host_count"])>=0
            if not row["introduced_share"]: missing+=1
            else: assert 0<=float(row["introduced_share"])<=1
          entry["missing_introduced_share"]=missing
          entry["unique_species_region_pairs"]=len(pairs)
        summary[name]=entry
  return {"status":"PASS","source_release_draft":True,"files":summary,
          "archive_members_non_directory":member_count,
          "additional_archive_members_not_yet_reviewed":member_count-len(PATHS),
          "zip_sha256":ASSETS["zip"][2],"receipt_sha256":ASSETS["receipt"][2],
          "source_record_rows_exported":False,
          "original_host_counts_independent_origin_not_established_by_hash_match":True}

if __name__=="__main__":
  out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True)
  result={"status":"INCOMPLETE"}
  try:result=verify()
  except Exception as e:result={"status":"FAIL","error":"verification_"+type(e).__name__,"checkpoint_line":traceback.extract_tb(e.__traceback__)[-1].lineno}
  out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
  print("EXTRA3_VERIFICATION_"+result["status"],{k:v["rows"] for k,v in result.get("files",{}).items()},flush=True)
  sys.exit(0 if result["status"]=="PASS" else 1)
