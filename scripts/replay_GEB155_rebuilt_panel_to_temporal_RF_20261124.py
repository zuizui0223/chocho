#!/usr/bin/env python3
"""Derived-source to original155 temporal model, no private records exported."""
import argparse, json, shutil, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd
import recover_GEB_v05_original155_private_inputs_20261122 as rec
import replay_GEB_v05_original155_partial_panel_20261124 as base
import replay_GEB_v05_original155_all_columns_20261124 as panel
import verify_GEB_v05_original155_refit_20261122 as verify

ORIGINAL_PRED='work/temporal/forecast/predictions_temporal.parquet'
class Guard(ValueError):pass
def require(b,label):
 if not b:raise Guard(label)

def check_model(root,snapshot,old_sha):
 oldpath=root/'upstream'/ORIGINAL_PRED
 newpath=root/'execution'/ORIGINAL_PRED
 require(rec.sha256_file(oldpath)==old_sha,'source_prediction_sha')
 archive_score=snapshot/'forecast_temporal_results.json'
 require(rec.sha256_file(archive_score)==verify.FROZEN['temporal'],'source_score_sha')
 orig=json.loads(archive_score.read_text())['results']
 out=root/'execution/outputs/forecast_temporal_results.json'
 fresh=json.loads(out.read_text())
 for key,expected in dict(status='COMPLETE',mode='temporal',train_rows=193703,
          train_species=155,test_rows=208482,test_species=155,
          settings=rec.MODEL_SETTINGS,model_lock_sha256=rec.LOCK_SHA256).items():
  require(fresh.get(key)==expected,'fitted_header_'+key)
 audit=verify.Audit()
 old,new=pd.read_parquet(oldpath),pd.read_parquet(newpath)
 require(len(old)==len(new)==208482,'original_test_cohort')
 verify.compare_predictions(audit,old,new,'temporal')
 for model in verify.MODELS:
  origm,origsp=verify.scores(old,model)
  newm,newsp=verify.scores(new,model)
  verify.metrics_compare(audit,'original_scored.'+model,orig[model],origm)
  verify.metrics_compare(audit,'new_scored.'+model,orig[model],newm)
  verify.metrics_compare(audit,'new_receipt.'+model,orig[model],fresh['results'][model])
  if model!='persistence':
   verify.compare_species(audit,'new_species.'+model,origsp,newsp)
 def gain(a):
  x=a['M3_both']['species_equal_mse']
  return (x-a['M5_host_composition']['species_equal_mse'])/x*100
 return dict(status='PASS' if audit.passed else 'NUMERICAL_MISMATCH',
   original_M5_M3_reduction_percent=gain(orig),
   fresh_M5_M3_reduction_percent=gain(fresh['results']),
   comparison_groups=len(audit.comparisons),
   failed_comparison_groups=sum(x['outside_tolerance']>0 for x in audit.comparisons.values()),
   predictions_outside_tolerance={
     m:audit.comparisons['row.'+m]['outside_tolerance'] for m in verify.MODELS
   },
   keyed_test_rows=len(new),relative_tolerance=verify.RTOL,absolute_tolerance=verify.ATOL,
   original_predictions_sha256=old_sha,
   fresh_predictions_sha256=rec.sha256_file(newpath),
   fresh_scores_sha256=rec.sha256_file(out))

def run(snapshot,output):
 result=dict(schema='GEB_original155_rebuilt_all96_panel_to_fresh_temporal_RF_v1',
             status='INCOMPLETE',no_source_rows_or_fits_exported=True,
             source_reconstructed_from_original_derived_intermediates=True,
             raw_original_gbif_occurrence_parquets_not_recomputed=True,
             historical_three_csv_creation_lineage_unbound=True)
 try:
  originalbuild=snapshot/'analysis_source/build_forecast_panel.py'
  originalfit=snapshot/'analysis_source/fit_forecast.py'
  require(rec.sha256_file(originalbuild)==panel.SCRIPTHASH,'builder_source_SHA')
  require(rec.sha256_file(originalfit)==rec.SCRIPT_SHA256,'fit_source_SHA')
  old_sourcefiles=base.SOURCE_FILES
  base.SOURCE_FILES=tuple(old_sourcefiles)+(ORIGINAL_PRED,)
  with tempfile.TemporaryDirectory(prefix='original155-panel-through-fit-') as tmp:
   root=Path(tmp)
   verified=base.restore_upstream(root)
   require(verified[rec.PANEL]['sha256']==rec.PANEL_SHA256,'panel_reference_SHA')
   execution=root/'execution'
   for rel in old_sourcefiles:
    if rel==rec.PANEL:continue
    target=execution/rel
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(root/'upstream'/rel,target)
   result['new_source_csv_sha256']=panel.add_three(root,execution)
   bind=execution/'scripts/temporal/build_forecast_panel.py'
   bind.parent.mkdir(parents=True,exist_ok=True)
   shutil.copyfile(originalbuild,bind)
   with (root/'private_builder_log').open('w') as log:
    p=subprocess.run([sys.executable,str(bind)],cwd=execution,stdout=log,
                     stderr=subprocess.STDOUT,check=False)
   result['original_builder_executed']=p.returncode==0
   require(p.returncode==0,'original_builder_process')
   rebuilt=execution/rec.PANEL
   require(rebuilt.is_file(),'rebuilt_panel_missing')
   comparison=panel.check_parquets(pd.read_parquet(root/'upstream'/rec.PANEL),
                                   pd.read_parquet(rebuilt))
   result['reconstructed_panel']=comparison
   result['regenerated_panel_sha256']=rec.sha256_file(rebuilt)
   require(comparison['status']=='ALL_COLUMNS_REPRODUCED','panel_column_mismatch')
   fit=execution/'scripts/temporal/fit_forecast.py'
   shutil.copyfile(originalfit,fit)
   with (root/'private_RF_log').open('w') as log:
    p=subprocess.run([sys.executable,str(fit),'temporal'],cwd=execution,stdout=log,
                     stderr=subprocess.STDOUT,check=False)
   result['original_temporal_RF_executed']=p.returncode==0
   require(p.returncode==0,'original_RF_fit_process')
   result['fit_comparison']=check_model(root,snapshot,verified[ORIGINAL_PRED]['sha256'])
   result['status']=result['fit_comparison']['status']
 except (Guard,panel.GateError,base.RebuildError) as err:
  result['status']='GUARD_FAILED';result['failure']=type(err).__name__
 except rec.RecoveryError:
  result['status']='SOURCE_RESTORE_FAILED';result['failure']='original_source_integrity'
 except Exception as err:
  result['status']='INCOMPLETE';result['failure']='bounded_'+type(err).__name__
 finally:
  output.parent.mkdir(parents=True,exist_ok=True)
  output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
  print('GEB155_REBUILT_PANEL_TO_TEMPORAL_'+result['status'],flush=True)
  if 'fit_comparison' in result:
   d=result['fit_comparison']
   print('REBUILT_PANEL_TEMPORAL_AGGREGATE',json.dumps({
      x:d[x] for x in ['keyed_test_rows','comparison_groups','failed_comparison_groups',
                       'fresh_M5_M3_reduction_percent']},sort_keys=True),flush=True)
 return 0 if result['status']=='PASS' else 1

if __name__=='__main__':
 a=argparse.ArgumentParser()
 a.add_argument('--snapshot',required=True,type=Path)
 a.add_argument('--output',required=True,type=Path)
 x=a.parse_args();sys.exit(run(x.snapshot.resolve(),x.output.resolve()))
