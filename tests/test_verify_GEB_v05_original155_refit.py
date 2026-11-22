import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import verify_GEB_v05_original155_refit_20261122 as v


def predictions():
    table = pd.DataFrame(dict(speciesKey=[1,1,2,2], species=['a','a','b','b'],
                              cell=[4,5,4,5], target=[0.,.5,.2,.3],
                              baseline_probability=[0.,.4,.1,0.],
                              distance_to_training_cell_km=[99.,100.,111.,5.]))
    for model in v.MODELS:
        table[model] = table.baseline_probability if model == 'persistence' else [.1,.3,.4,.1]
    return table


class Comparisons(unittest.TestCase):
    def test_keys_align_independently_of_order(self):
        frame = predictions()
        audit = v.Audit()
        v.compare_predictions(audit, frame, frame.iloc[::-1], 'spatial')
        self.assertTrue(audit.passed)

    def test_wrong_rows_even_when_mean_error_equal(self):
        original = predictions()
        original['target'] = 0.
        original.loc[[0,1], 'baseline_probability'] = 0.
        original['persistence'] = original.baseline_probability
        changed = original.copy()
        changed.loc[[0,1], v.MODELS[-1]] = changed.loc[[1,0], v.MODELS[-1]].to_numpy()
        audit = v.Audit()
        v.compare_predictions(audit,original,changed,'temporal')
        self.assertFalse(audit.passed)
        self.assertEqual(audit.comparisons['row.'+v.MODELS[-1]]['outside_tolerance'],2)
        self.assertEqual(v.scores(original,v.MODELS[-1])[0],v.scores(changed,v.MODELS[-1])[0])

    def test_missing_duplicate_misaligned_and_nan(self):
        orig = predictions()
        for new in [orig.iloc[:-1],pd.concat([orig.iloc[:-1], orig.iloc[[1]]]),
                    orig.assign(cell=[8,9,10,11]),orig.assign(speciesKey=[1,1,np.nan,2])]:
            with self.assertRaises(v.Invalid):
                v.compare_predictions(v.Audit(),orig,new,'temporal')
        nan = orig.copy()
        nan.loc[0,'target'] = np.nan
        a = v.Audit(); v.compare_predictions(a,orig,nan,'temporal')
        self.assertFalse(a.passed)

    def test_lock_and_support(self):
        header = dict(mode='temporal',status='COMPLETE',train_rows=193703,test_rows=208482,
                      train_species=155,test_species=155,settings=v.SETTINGS,model_lock_sha256=v.LOCK)
        v.header(header,'temporal')
        for corrupt in [dict(train_rows=193704),dict(test_species=154),dict(settings={**v.SETTINGS,'n_estimators':159}),dict(status='IN_PROGRESS')]:
            with self.assertRaises(v.Invalid): v.header({**header,**corrupt},'temporal')

    def test_boundary_and_incorrect_spatial_membership(self):
        orig = predictions()
        new = orig.copy()
        new.loc[1,'distance_to_training_cell_km']=99.99999999
        a=v.Audit();v.compare_predictions(a,orig,new,'spatial')
        self.assertFalse(a.passed)
        self.assertEqual(a.comparisons['row.buffer100_membership']['outside_tolerance'],1)

    def test_per_species_error_and_tolerance(self):
        orig=predictions()
        metrics, per=v.scores(orig,'M3_both')
        self.assertAlmostEqual(metrics['species_equal_mse'],.0325)
        changed=per.copy()
        changed.loc['a','mse']+=.001
        a=v.Audit();v.compare_species(a,'group',per,changed)
        self.assertFalse(a.passed)
        a=v.Audit();a.compare('inside',[1e-3],[1e-3+1e-10])
        self.assertTrue(a.passed)

class ExecutionChecks(unittest.TestCase):
    def test_attestation_guards(self):
        receipt = dict(run_id=12, mode='spatial', source_script_sha256=v.FIT,
                       panel_sha256=v.PANEL, started_at='2026-01-01T00:00:02Z',
                       completed_at='2026-01-01T00:00:07+00:00', returncode=0,
                       preexisting_outputs=False)
        run_id, start, end = v.execution(receipt, 'spatial')
        self.assertEqual((run_id, end-start), ('12', 5))
        for change in [dict(preexisting_outputs=True), dict(returncode=False), dict(returncode=1),
                       dict(mode='species'), dict(panel_sha256='incorrect'), dict(source_script_sha256='wrong'),
                       dict(run_id='private-label'), dict(completed_at=start-1),
                       dict(started_at='2026-01-01T00:00:02')]:
            with self.subTest(change=list(change)), self.assertRaises(v.Invalid):
                v.execution({**receipt, **change}, 'spatial')

    def test_missing_outputs_or_outside_fit_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'fresh-file'
            with self.assertRaises(v.Invalid): v.fresh_file(p, 50, 75)
            p.write_bytes(b'file')
            os.utime(p, (55, 55))
            self.assertEqual(v.fresh_file(p, 50, 75), p)
            for t in [35, 85]:
                os.utime(p, (t, t))
                with self.assertRaises(v.Invalid): v.fresh_file(p, 50, 75)

    def test_cli_failure_is_sanitized_and_never_attests(self):
        secret = 'PRIVATE_machine_name_and_label'
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root/'public.json'
            args = ['verify', '--mode', 'temporal', '--fit-root', str(root/'fit'),
                    '--original-predictions', str(root/secret), '--snapshot-dir', str(root/'original'),
                    '--execution-receipt', str(root/secret), '--original-predictions-sha256', '0'*64,
                    '--output', str(output)]
            log = io.StringIO()
            with patch.object(v.sys, 'argv', args), contextlib.redirect_stdout(log):
                code = v.main()
            report = json.loads(output.read_text())
            self.assertEqual(code, 1)
            self.assertEqual(report['status'], 'FAIL')
            self.assertFalse(report['model_refit_completed'])
            for text in (log.getvalue(), output.read_text()):
                self.assertNotIn(secret, text)
                self.assertNotIn(str(root), text)


if __name__ == '__main__':
    unittest.main()
