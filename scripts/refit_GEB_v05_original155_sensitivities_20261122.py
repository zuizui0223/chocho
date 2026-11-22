#!/usr/bin/env python3
"""Re-fit an original 155-species temporal sensitivity in a private runner; export aggregate QA only."""
import argparse,hashlib,json,shutil,subprocess,sys,zipfile
from pathlib import Path
import pandas as pd
import recover_GEB_v05_original155_private_inputs_20261122 as rec
import verify_GEB_v05_original155_refit_20261122 as v
S={
'effort20':(122547,159084,20,'78c3b5151f3e1f9b4584aa5ae89ae757f8ceeb91fc923c9e7b37d3180cd13a06'),
'radius250':(167192,180654,10,'b6d8e6b0a064a040eed7c1f434ebd8ec494c14eceeccca6a214f9121d6c958d6'),
'exclude_top_provider':(173414,191004,10,'6a2e00321a0ef027131a90a0eced4d4bba6f629d1b8869c29c2d7eb77fa93831'),
'precision100':(98838,135532,10,'1fac06a6564958128a7191ceb61c43e51d22ffb28cedab02603fda25f56de4f2')}
class SafeError(Exception):pass
def ok(cond,label):
 if not cond:raise SafeError(label)
def source_names(k):
 return f'work/temporal/forecast/panel_{k}.parquet',f'work/temporal/forecast/predictions_temporal_{k}.parquet'
def restore(root,k,fit_script):
 panel,pred=source_names(k)
 ok(rec.sha256_file(fit_script)==rec.SCRIPT_SHA256,'fit_script_identity')
 rel=rec.gh_json(f'repos/{rec.REPO}/releases/{rec.RELEASE_ID}')
 ok(rel.get('draft') is True,'source_release_not_draft')
 assets={a['id']:a for a in rel['assets']}
 for asset,size,sha in [rec.RECEIPT,*rec.PARTS]:
  ok(assets.get(asset,{}).get('size')==size and assets[asset]['digest']=='sha256:'+sha,'release_assets_changed')
 root.mkdir(parents=True,exist_ok=False)
 manifest=root/'manifest.json'
 with manifest.open('xb') as out:
  rec.stream_asset(rec.RECEIPT[0],out,rec.RECEIPT[1],rec.RECEIPT[2],hashlib.sha256(),root)
 records=json.loads(manifest.read_text())['files']
 ok(len(records)==3032,'source_manifest_count')
 targets=[panel,pred,rec.LOCK]
 expected={p:rec.find_manifest_record(records,p) for p in targets}
 arc=root/'archive.zip'
 combined=hashlib.sha256()
 try:
  with arc.open('xb') as out:
   for asset,size,sha in rec.PARTS:rec.stream_asset(asset,out,size,sha,combined,root)
  ok(combined.hexdigest()==rec.ARCHIVE_SHA256,'complete_source_archive_sha')
  with zipfile.ZipFile(arc) as z:
   members=z.infolist(); names=[rec.safe_zip_member(m) for m in members]
   ok(len(names)==len(set(names)),'archive_duplicate_names')
   checked={}
   for p in targets:
    checked[p]=rec.extract_verified_member(z,rec.find_zip_member(members,p),root/'original'/p,expected[p])
 finally:
  arc.unlink(missing_ok=True);manifest.unlink(missing_ok=True)
 ok(checked[rec.LOCK]['sha256']==rec.LOCK_SHA256,'source_lock_sha')
 lock=json.loads((root/'original'/rec.LOCK).read_text())
 ok(lock['status']=='LOCKED_BEFORE_MODEL_FIT' and lock['random_forest']==rec.MODEL_SETTINGS,'source_lock_settings')
 d=pd.read_parquet(root/'original'/panel,columns=['speciesKey','cell','species','period','recording_probability10'])
 def ntrans(t):
  a=d[d.period==t][['speciesKey','cell','recording_probability10']]
  b=d[d.period==t+1][['speciesKey','cell','recording_probability10']]
  return len(a.merge(b,on=['speciesKey','cell'],validate='one_to_one',suffixes=('_a','_b')).dropna())
 train,test,nrare,_=S[k]
 ok(d.species.nunique()==155 and ntrans(0)==train and ntrans(1)==test,'source_panel_cohort')
 fit=root/'fit'
 for name in (panel,rec.LOCK):
  target=fit/name;target.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(root/'original'/name,target)
 script=fit/'scripts/temporal/fit_forecast.py'
 script.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(fit_script,script)
 return {'panel_sha256':checked[panel]['sha256'],'predictions_sha256':checked[pred]['sha256'],
         'panel_rows':len(d),'train_rows':train,'test_rows':test,'species':155,
         'archive_sha256':rec.ARCHIVE_SHA256,'manifest_sha256':rec.RECEIPT[2]}
def compare(root,k,source,snapshot):
 train,test,nrare,score_sha=S[k]
 original_score=snapshot/f'forecast_temporal_{k}_results.json'
 ok(rec.sha256_file(original_score)==score_sha,'source_score_sha')
 old_s=json.loads(original_score.read_text())
 new_s=json.loads((root/'fit/outputs'/original_score.name).read_text())
 for score in (old_s,new_s):
  for name,value in {'status':'COMPLETE','mode':'temporal','sensitivity':k,'rarefaction_n':nrare,
   'train_rows':train,'test_rows':test,'train_species':155,'test_species':155,
   'settings':rec.MODEL_SETTINGS,'model_lock_sha256':rec.LOCK_SHA256}.items():
   ok(score.get(name)==value,'score_header_'+name)
  ok(set(score['results'])==set(v.MODELS),'model_groups')
 _,file=source_names(k)
 old_path=root/'original'/file; new_path=root/'fit'/file
 ok(rec.sha256_file(old_path)==source['predictions_sha256'],'original_prediction_sha')
 old,new=pd.read_parquet(old_path),pd.read_parquet(new_path)
 ok(len(old)==len(new)==test and old.species.nunique()==new.species.nunique()==155,'test_row_support')
 audit=v.Audit();v.compare_predictions(audit,old,new,'temporal')
 for model in v.MODELS:
  metrics_old,by_old=v.scores(old,model)
  metrics_new,by_new=v.scores(new,model)
  v.metrics_compare(audit,'archived_rows.'+model,old_s['results'][model],metrics_old)
  v.metrics_compare(audit,'fresh_score.'+model,old_s['results'][model],new_s['results'][model])
  v.metrics_compare(audit,'fresh_rows.'+model,old_s['results'][model],metrics_new)
  if model!='persistence':
   v.compare_species(audit,'species_error.'+model,by_old,by_new)
   csv=root/'fit/outputs'/f'forecast_temporal_{k}_{model}_species_metrics.csv'
   ok(csv.is_file(),'fresh_species_csv_missing')
   v.compare_species(audit,'species_csv.'+model,by_old,pd.read_csv(csv,index_col='species'))
 def gain(s):
  a=s['results']['M3_both']['species_equal_mse'];b=s['results']['M5_host_composition']['species_equal_mse']
  return (a-b)/a*100
 return {'status':'PASS' if audit.passed else 'NUMERIC_MISMATCH',
  'comparisons':len(audit.comparisons),
  'comparison_groups_failed':sum(x['outside_tolerance']>0 for x in audit.comparisons.values()),
  'per_model_row_mismatches':{m:audit.comparisons['row.'+m]['outside_tolerance'] for m in v.MODELS},
  'archived_M5_M3_MSE_reduction_percent':gain(old_s),
  'fresh_M5_M3_MSE_reduction_percent':gain(new_s),
  'relative_tolerance':v.RTOL,'absolute_tolerance':v.ATOL,'frozen_scores_sha256':score_sha}
def main():
 p=argparse.ArgumentParser()
 p.add_argument('--sensitivity',choices=S,required=True)
 for name in ('root','snapshot','fit-script','output'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args(); root=a.root.resolve(); output=a.output.resolve();output.parent.mkdir(parents=True,exist_ok=True)
 result={'schema':'GEB_original155_four_sensitivity_private_refit_v1',
  'sensitivity':a.sensitivity,'status':'INCOMPLETE','source_commit':rec.SOURCE_COMMIT,
  'source_panel_reconstruction_not_executed':True,'no_private_rows_or_paths_exported':True}
 try:
  source=restore(root,a.sensitivity,a.fit_script.resolve());result['verified_inputs']=source
  with (root/'local_private_fit.log').open('w') as log:
   r=subprocess.run([sys.executable,str((root/'fit/scripts/temporal/fit_forecast.py').resolve()),
                     'temporal','--sensitivity',a.sensitivity],cwd=root/'fit',stdout=log,stderr=subprocess.STDOUT)
  result['fresh_script_executed']=r.returncode==0
  ok(r.returncode==0,'fresh_fit_failed')
  result['numerical_comparison']=compare(root,a.sensitivity,source,a.snapshot.resolve())
  result['status']=result['numerical_comparison']['status']
 except SafeError as e:result['failure']=str(e)
 except rec.RecoveryError:result['failure']='verified_private_source_error'
 except Exception as e:result['failure']='safe_error_'+type(e).__name__
 finally:
  output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
  if root.is_dir():shutil.rmtree(root)
 print('ORIGINAL155_SENSITIVITY_'+result['status'],a.sensitivity,flush=True)
 return 0 if result['status']=='PASS' else 1
if __name__=='__main__':sys.exit(main())
