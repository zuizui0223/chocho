"""Wait for one identified acquisition, then run locked comparisons sequentially."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
import psutil

R = Path(__file__).resolve().parents[2]
O = R / 'outputs'
parser = argparse.ArgumentParser()
parser.add_argument('--acquisition-pid', type=int, required=True)
parser.add_argument('--created', type=float, required=True)
parser.add_argument('--check-only', action='store_true')
args = parser.parse_args()

def source_alive():
    try:
        p = psutil.Process(args.acquisition_pid)
        if abs(p.create_time() - args.created) > .01:
            return False
        assert any(Path(x).name == 'acquire_route_terrain.py' for x in p.cmdline())
        return p.is_running()
    except psutil.NoSuchProcess:
        return False

def state(status, **extra):
    dest = O / 'route_terrain_pipeline_progress.json'
    temp = dest.with_suffix('.tmp')
    temp.write_text(json.dumps(dict(status=status, source_pid=args.acquisition_pid,
        source_created=args.created, **extra), indent=2))
    temp.replace(dest)
    print(status, extra, flush=True)

def run(script, *arguments):
    state('RUNNING_STAGE', script=script, arguments=arguments)
    subprocess.run([sys.executable, str(R/'work/barriers'/script), *arguments],
                   cwd=R, check=True)

alive = source_alive()
if args.check_only:
    assert alive, 'Expected acquisition process is not live'
    print('PASS: exact PID, creation time and acquisition command; no stages started')
    sys.exit(0)

try:
    assert alive, 'Start only while the identified acquisition is live'
    state('WAITING_FOR_IDENTIFIED_ACQUISITION')
    while source_alive():
        time.sleep(10)
    # A crashed acquisition without its complete-pass receipt must not trigger fitting.
    path = O/'route_terrain_acquisition.json'
    assert path.exists() and path.stat().st_mtime >= args.created
    expected = json.loads((O/'route_terrain_design_audit.json').read_text())['required_tiles']
    for attempt in range(3):
        receipt = json.loads(path.read_text())
        assert receipt['status'] == 'ACQUISITION_PASS_FINISHED_NOT_FEATURE_VALIDATION'
        assert receipt['requested_tiles'] == expected and receipt['partial_limit'] == 0
        assert sum(receipt['counts'].values()) == expected
        if receipt['counts'].get('ERROR', 0) == 0:
            break
        if attempt == 2:
            raise RuntimeError('Retrieval errors remain after two retries; no model fitted')
        # Successful tile caches are reused; only unsuccessful tiles are retrieved again.
        run('acquire_route_terrain.py')
    run('aggregate_route_terrain.py')
    for mode in ['temporal','spatial']:
        run('fit_route_terrain_forecast.py', mode)
        run('summarize_route_terrain_forecast.py', mode)
    state('COMPARISONS_COMPLETE_REPORT_AND_REVIEW_PENDING')
except Exception as exc:
    state('STOPPED_REQUIRES_INSPECTION', error_type=type(exc).__name__)
    raise
