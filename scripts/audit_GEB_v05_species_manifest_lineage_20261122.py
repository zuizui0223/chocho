#!/usr/bin/env python3
"""Read the original155 private source MANIFEST only; export bounded hash aggregates."""
import argparse, hashlib, json, re, tempfile
from pathlib import Path, PurePosixPath
import recover_GEB_v05_original155_private_inputs_20261122 as source

SHA = re.compile(r"^[0-9a-f]{64}$")
ALLOW = {
    "species_cell_period_panel.parquet": "panel",
    "panel_effort20.parquet": "sensitivity_panel_effort20",
    "panel_radius250.parquet": "sensitivity_panel_radius250",
    "panel_exclude_top_provider.parquet": "sensitivity_panel_exclude_top_provider",
    "panel_precision100.parquet": "sensitivity_panel_precision100",
    "predictions_temporal_effort20.parquet": "sensitivity_prediction_effort20",
    "predictions_temporal_radius250.parquet": "sensitivity_prediction_radius250",
    "predictions_temporal_exclude_top_provider.parquet": "sensitivity_prediction_exclude_top_provider",
    "predictions_temporal_precision100.parquet": "sensitivity_prediction_precision100",

    "predictions_species.parquet": "species_predictions",
    "predictions_temporal.parquet": "temporal_predictions",
    "predictions_spatial.parquet": "spatial_predictions",
    "forecast_species_results.json": "species_scores",
    "forecast_model_lock.json": "model_lock",
}
EXPECTED_SHA = {
    "panel": source.PANEL_SHA256,
    "model_lock": source.LOCK_SHA256,
    "species_scores": "3eaaeeb8b77088314bf42051c5a51c07895c17615bf5f27fc3a33424f25e125a",
}
class ScanError(ValueError): pass

def category(raw):
    if not isinstance(raw, str): return None
    name = PurePosixPath(raw.replace("\\", "/")).name.lower()
    if name in ALLOW: return ALLOW[name]
    if re.fullmatch(r"fit_forecast[a-z0-9_.-]*\.py", name): return "fit_forecast_script_candidates"
    if re.fullmatch(r"model_m[0-5]_[a-z0-9_]+\.joblib", name): return "saved_model_candidates"
    if name in ("requirements.txt", "requirements.lock", "environment.yml",
                "environment.yaml", "pyproject.toml", "uv.lock",
                "poetry.lock", "pip-freeze.txt", "pip_freeze.txt"): return "runtime_manifest_candidates"
    return None

def summarize(records):
    if not isinstance(records, list): raise ScanError("bad_manifest")
    groups = {}
    invalid = 0
    for item in records:
        if not isinstance(item, dict): raise ScanError("invalid_record")
        c = category(item.get("path"))
        if c is None: continue
        value,size=item.get("sha256"),item.get("bytes")
        if not isinstance(value,str) or not SHA.fullmatch(value) or type(size) is not int or size < 0:
            invalid += 1
            continue
        key=(value,size)
        groups.setdefault(c,{}).setdefault(key,0)
        groups[c][key]+=1
    if invalid: raise ScanError("invalid_matching_record")
    summary={}
    for kind,counts in sorted(groups.items()):
        summary[kind] = {
            "entries": sum(counts.values()),
            "distinct_content_hashes": len(counts),
            "hash_groups": [{"sha256": sha,"bytes": sz,"count": n,
               "is_expected_digest": (sha == EXPECTED_SHA[kind]) if kind in EXPECTED_SHA else None}
                for (sha,sz),n in sorted(counts.items())],
        }
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    release=source.gh_json(f"repos/{source.REPO}/releases/{source.RELEASE_ID}")
    if release.get("draft") is not True: raise ScanError("draft_release_not_verified")
    assets={x["id"]:x for x in release.get("assets",[])}
    receipt_id,nbytes,expected=source.RECEIPT
    asset=assets.get(receipt_id,{})
    if asset.get("size")!=nbytes or asset.get("digest")!="sha256:"+expected:
        raise ScanError("source_manifest_asset_binding_failed")
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/"original_manifest.json"
        with path.open("wb") as stream:
            source.stream_asset(receipt_id,stream,nbytes,expected,hashlib.sha256(),Path(temp))
        data=json.loads(path.read_text())
        records=data.get("files")
        if len(records)!=3032: raise ScanError("source_manifest_count_changed")
        groups=summarize(records)
        exact_suffix_presence = {}
        for sensitivity in ('effort20','radius250','exclude_top_provider','precision100'):
            for role, relative in (
                ('panel',f'work/temporal/forecast/panel_{sensitivity}.parquet'),
                ('predictions',f'work/temporal/forecast/predictions_temporal_{sensitivity}.parquet')):
                try:
                    source.find_manifest_record(records,relative)
                    match = True
                except source.RecoveryError:
                    match = False
                exact_suffix_presence[sensitivity+'_'+role] = match

    result={
        "schema":"geb_original155_manifest_only_code_lineage_v1",
        "release_is_draft":True,
        "manifest_sha256":expected,
        "manifest_bytes":nbytes,
        "source_record_count":3032,
        "original_zip_not_downloaded":True,
        "matching_source_groups":groups,
        "expected_relative_path_suffixes":exact_suffix_presence,
        "no_private_paths_or_source_rows_exported":True,
        "limits":"Verified manifest paths and registered digests only; existence does not establish historical execution time or code version.",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MANIFEST_ONLY_CLASS_COUNTS",json.dumps({k:v["entries"] for k,v in groups.items()},sort_keys=True))
    print("MANIFEST_ONLY_SOURCE_SCRIPT_HASHES",groups.get("fit_forecast_script_candidates",{}).get("distinct_content_hashes",0))

if __name__=="__main__":
    try: main()
    except (ScanError,source.RecoveryError) as exc: raise SystemExit("MANIFEST_AUDIT_FAILED_"+str(exc))
    except Exception as exc: raise SystemExit("MANIFEST_AUDIT_ERROR_"+type(exc).__name__)
