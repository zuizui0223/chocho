#!/usr/bin/env python3
"""Re-run ORIGINAL 155-species panel builder with recovered 3 source CSVs.

No reconstructed or original data rows/CSVs/parquets leave the runner.
The 2026-11-24 recovered host-count source was absent from the old manifest;
model-panel match establishes numerical reconstruction, not historic file origin.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import recover_GEB_v05_original155_private_inputs_20261122 as source
import replay_GEB_v05_original155_partial_panel_20261124 as prev
import audit_GEB155_new_extra3_20261124 as additional

EXPECTED_SHA={
 'outputs/forecast_candidate_disposition.csv':'f275181f49a7d8da99cae815717c2f2ce608eb14ca3be5609482d82f307a7003',
 'outputs/temporal_host_family_support.csv':'4644a2ed2f233052089b8880c25f22072158421f5f0396d668398899414cbb0e',
 'work/remote_sensing/host_region_counts.csv':'5c074721f4ebb3e2e413d6e4926320b97f5cb1b78f7004b0138ceeecd33d2b8d',
}
SCRIPTHASH='99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6'
class GateError(ValueError):pass
def require(ok,label):
 if not ok:raise GateError(label)
def add_three(root, execution):
 rel=source.gh_json(f'repos/{source.REPO}/releases/{source.RELEASE_ID}')
 require(rel.get('draft') is True,'source_release_draft')
 assets={x['id']:x for x in rel['assets']}
 for _,(id_,size,digest) in additional.ASSETS.items():
  obj=assets.get(id_,{})
  require(obj.get('size')==size and obj.get('digest')=='sha256:'+digest,
          'new_assets_changed')
 rootzip=root/'extras.zip';receiptpath=root/'extras.json'
 for kind, path in [('zip',rootzip),('receipt',receiptpath)]:
  id_,size,digest=additional.ASSETS[kind]
  with path.open('xb') as out:
   source.stream_asset(id_,out,size,digest,hashlib.sha256(),root)
 receipt=json.loads(receiptpath.read_text(encoding='utf-8'))
 receipt_text=json.dumps(receipt,sort_keys=True)
 verify={}
 try:
  with zipfile.ZipFile(rootzip) as z:
   members=[m for m in z.infolist() if not m.is_dir()]
   require(len(members)==4,'extra_source_member_count_drift')
   names=[source.safe_zip_member(x) for x in members]
   require(len(names)==len(set(names)),'extra_source_duplicate_paths')
   for path in additional.PATHS:
    matches=[m for m in members if m.filename==path or m.filename.endswith('/'+path)]
    require(len(matches)==1,'missing_extra_member')
    member=matches[0]
    require(member.file_size < 20000000,'extra_member_size')
    destination=execution/path
    destination.parent.mkdir(parents=True,exist_ok=True)
    with z.open(member) as src, destination.open('xb') as out:
     sha=hashlib.sha256();total=0
     for chunk in iter(lambda:src.read(1024*1024),b''):
      total+=len(chunk);sha.update(chunk);out.write(chunk)
    require(total==member.file_size and sha.hexdigest()==EXPECTED_SHA[path],
      'extra_source_sha_mismatch')
    require(sha.hexdigest() in receipt_text,'extra_receipt_lacks_sha')
    verify[path]={'sha256':sha.hexdigest(),'bytes':total}
 finally:
  rootzip.unlink(missing_ok=True);receiptpath.unlink(missing_ok=True)
 return verify
def check_parquets(original,fresh):
 key=['speciesKey','cell','period']
 require(len(original)==842964 and len(fresh)==842964,'panel_row_count')
 require(original.species.nunique()==fresh.species.nunique()==155,'panel_species_count')
 require(set(original.columns)==set(fresh.columns),'panel_column_set')
 require(not original.duplicated(key).any() and not fresh.duplicated(key).any(),
         'duplicate_panel_keys')
 a=original.sort_values(key).reset_index(drop=True)
 b=fresh.sort_values(key).reset_index(drop=True)
 require(a[key].equals(b[key]),'panel_key_disagreement')
 cols=list(original.columns)
 data={}
 for c in cols:
  left,right=a[c],b[c]
  if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
   x=np.asarray(left,dtype='float64');y=np.asarray(right,dtype='float64')
   nan=np.isnan(x) & np.isnan(y)
   finite=np.isfinite(x)&np.isfinite(y)
   with np.errstate(invalid='ignore'):
    diff=np.abs(x[finite]-y[finite])
   good=nan.copy()
   good[finite]=np.isclose(x[finite],y[finite],rtol=1e-6,atol=1e-8)
   count=int((~good).sum())
   maximum=float(diff.max()) if len(diff) else None
  else:
   count=int((~left.astype('string').fillna('<NA>').eq(
     right.astype('string').fillna('<NA>'))).sum())
   maximum=None
  data[c]={'mismatched_rows':count,'max_absolute_difference':maximum}
 return {
  'status':'ALL_COLUMNS_REPRODUCED' if all(x['mismatched_rows']==0 for x in data.values())
    else 'ORIGINAL_PANEL_VALUE_MISMATCH',
  'columns':len(data),'rows':len(a),
  'all_column_mismatch_count':sum(v['mismatched_rows'] for v in data.values()),
  'mismatched_columns':[k for k,v in data.items() if v['mismatched_rows']>0],
  'origin_credit_columns':{name:data[name] for name in (
    'native_host_count','introduced_host_count','introduced_share')},
  'numeric_rtol':1e-6,'numeric_atol':1e-8,
  'source_panel_byte_identity_not_assumed':True
 }
def run(snapshot, output):
 report={'status':'INCOMPLETE','no_source_rows_exported':True,
         'original_source_commit':source.SOURCE_COMMIT,
         'old_source_archive_sha256':source.ARCHIVE_SHA256,
         'three_missing_csv_sha256':EXPECTED_SHA,
         'historic_three_csv_source_provenance_proven':False}
 try:
  require(source.sha256_file(snapshot/'analysis_source/build_forecast_panel.py')==
    SCRIPTHASH,'original_script_digest')
  with tempfile.TemporaryDirectory(prefix='geb155-full-panel-replay-') as d:
   root=Path(d)
   verified=prev.restore_upstream(root)
   require(verified[source.PANEL]['sha256']==source.PANEL_SHA256,
           'original_target_panel_digest')
   execution=root/'execution'
   original_target=root/'upstream'/source.PANEL
   for path in prev.SOURCE_FILES:
    if path==source.PANEL:continue
    dest=execution/path
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(root/'upstream'/path,dest)
   report['extra_source_checksums']=add_three(root,execution)
   original_script=snapshot/'analysis_source/build_forecast_panel.py'
   script=execution/'scripts/temporal/build_forecast_panel.py'
   script.parent.mkdir(parents=True,exist_ok=True)
   shutil.copyfile(original_script,script)
   proc=subprocess.run([sys.executable,str(script)],cwd=execution,
                      stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=False)
   report['actual_original_builder_executed']=(proc.returncode==0)
   require(proc.returncode==0,'original_source_builder_execution')
   freshpath=execution/source.PANEL
   require(freshpath.is_file(),'new_panel_not_written')
   report['comparison']=check_parquets(pd.read_parquet(original_target),
                                        pd.read_parquet(freshpath))
   report['original_panel_sha256']=source.PANEL_SHA256
   report['fresh_parquet_sha256']=source.sha256_file(freshpath)
   report['status']='COMPLETE' if report['comparison']['status']=='ALL_COLUMNS_REPRODUCED' \
                   else 'NUMERICAL_COMPARISON_FAILED'
  report['raw_GBIF_occurrence_partitions_reread']=False
  report['recovered_source_intermediates_used']=True
  report['complete_raw_observation_to_panel_replay']=False
  report['original_host_count_provenance_independently_validated']=False
 except GateError as exc:
  report['failure']='bounded_'+str(exc)
 except source.RecoveryError:
  report['failure']='verified_archive_source_error'
 except Exception as exc:
  report['failure']='safe_'+type(exc).__name__
 finally:
  output.parent.mkdir(parents=True,exist_ok=True)
  output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
 print('GEB155_THREE_SOURCE_FULL_PANEL_REPLAY_'+report['status'],flush=True)
 if 'comparison' in report:
  print('GEB155_PANEL_AGGREGATE',json.dumps({k:v for k,v in report['comparison'].items()
    if k in ('status','rows','columns','all_column_mismatch_count','mismatched_columns')}),
    flush=True)
 return 0 if report['status']=='COMPLETE' else 1
if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--snapshot',required=True,type=Path)
 parser.add_argument('--output',required=True,type=Path)
 z=parser.parse_args();sys.exit(run(z.snapshot.resolve(),z.output.resolve()))
