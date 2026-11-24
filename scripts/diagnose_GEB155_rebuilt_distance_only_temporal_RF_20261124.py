#!/usr/bin/env python3
"""Derived-source to original155 temporal model, no private records exported."""
import argparse, json, shutil, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd
import numpy as np
import recover_GEB_v05_original155_private_inputs_20261122 as rec
import replay_GEB_v05_original155_partial_panel_20261124 as base
import replay_GEB_v05_original155_all_columns_20261124 as panel
import verify_GEB_v05_original155_refit_20261122 as verify

ORIGINAL_PRED='work/temporal/forecast/predictions_temporal.parquet'
class Guard(ValueError):pass
def require(b,label):
 if not b:raise Guard(label)

def diagnose_row_order(original, regenerated):
 keys=['speciesKey','cell','period']
 oldkeys=original[keys].reset_index(drop=True)
 newkeys=regenerated[keys].reset_index(drop=True)
 require(len(oldkeys)==len(newkeys),'ordered_key_support')
 matches=(oldkeys['speciesKey'].eq(newkeys['speciesKey']) &
          oldkeys['cell'].eq(newkeys['cell']) &
          oldkeys['period'].eq(newkeys['period']))
 orig_species=original.species.drop_duplicates().tolist()
 fresh_species=regenerated.species.drop_duplicates().tolist()
 keyed_a=original.set_index(keys).sort_index()
 keyed_b=regenerated.set_index(keys).sort_index()
 numeric_exact_diffs={}
 numeric_float32_diffs={}
 for name in sorted(set(original.columns)-set(keys)):
  left=keyed_a[name]
  right=keyed_b[name]
  if not (pd.api.types.is_numeric_dtype(left) and
          pd.api.types.is_numeric_dtype(right)):
   continue
  a=np.asarray(left,dtype='float64')
  b=np.asarray(right,dtype='float64')
  good=(np.isnan(a)&np.isnan(b))|(a==b)
  mismatched=int((~good).sum())
  if mismatched:numeric_exact_diffs[name]=mismatched
  a32=a.astype('float32')
  b32=b.astype('float32')
  good32=(np.isnan(a32)&np.isnan(b32))|(a32==b32)
  count32=int((~good32).sum())
  if count32:numeric_float32_diffs[name]=count32
 return {
  'ordered_key_rows':len(oldkeys),
  'rows_with_exact_original_position':int(matches.sum()),
  'rows_with_different_position':int((~matches).sum()),
  'original_vs_rebuilt_species_first_appearance_order_equal':orig_species==fresh_species,
  'species_first_appearance_positions_different':sum(
      a!=b for a,b in zip(orig_species,fresh_species)),
  'numeric_exact_differences_by_column':numeric_exact_diffs,
  'numeric_float32_differences_by_column':numeric_float32_diffs,
  'compared_by_sorted_unique_keys':True,
  'no_original_predictor_values_copied':True,
 }

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
 result=dict(schema='GEB_original155_rebuilt_panel_distance_only_counterfactual_RF_v1',
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
   archived_panel=pd.read_parquet(root/'upstream'/rec.PANEL)
   freshly_built=pd.read_parquet(rebuilt)
   comparison=panel.check_parquets(archived_panel,freshly_built)
   result['reconstructed_panel']=comparison
   result['regenerated_panel_sha256_before_row_alignment']=rec.sha256_file(rebuilt)
   require(comparison['status']=='ALL_COLUMNS_REPRODUCED','panel_column_mismatch')
   result['ordering_diagnostic']=diagnose_row_order(archived_panel,freshly_built)
   keys=['speciesKey','cell','period']
   reference_order=pd.MultiIndex.from_frame(archived_panel[keys])
   aligned=(freshly_built.set_index(keys).loc[reference_order].reset_index())
   require(aligned[keys].reset_index(drop=True).equals(
      archived_panel[keys].reset_index(drop=True)),'alignment_failed')
   aligned=aligned[freshly_built.columns]
   # COUNTERFACTUAL ONLY: replace one original-model predictor from the
   # original archived target; all other 95 columns remain independently rebuilt.
   aligned['distance_km']=archived_panel['distance_km'].to_numpy(copy=True)
   result['distance_counterfactual_original_values_reused']=True
   result['distance_counterfactual_only_one_column_replaced']='distance_km'
   result['recreated_non_distance_columns_from_source']=95
   result['counterfactual_exact_match_diagnostic']=diagnose_row_order(archived_panel,aligned)
   require(not result['counterfactual_exact_match_diagnostic']['numeric_exact_differences_by_column'],
           'non_distance_numerical_difference_after_control')
   require(not result['counterfactual_exact_match_diagnostic']['numeric_float32_differences_by_column'],
           'non_distance_float32_difference_after_control')
   aligned.to_parquet(rebuilt,index=False,compression='zstd')
   result['aligned_regenerated_panel_sha256']=rec.sha256_file(rebuilt)
   result['alignment_used_only_original_row_key_order_except_declared_distance_column']=True
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
