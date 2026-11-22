import hashlib
import io
from pathlib import Path
import sys
import unittest
import zipfile
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import diagnose_GEB_v05_original155_species_20261122 as d


class Lineage(unittest.TestCase):
    def test_archive_report_uses_only_allowlisted_names_and_verified_bytes(self):
        raw = io.BytesIO()
        files = {d.recovery.PANEL:b'panel', d.SCRIPT:b'private fit',
                 'personal-name/requirements.txt': b'private contents',
                 'personal-name/unrecognized.json': b'private'}
        records = []
        with zipfile.ZipFile(raw, 'w') as z:
            for name, data in files.items():
                z.writestr('private-directory/'+name,data)
                records.append(dict(path='/original-machine/'+name,bytes=len(data),
                                    sha256=hashlib.sha256(data).hexdigest()))
        with zipfile.ZipFile(raw) as z:
            result = d.archive_metadata(z,z.infolist(),records)
            self.assertIs(result['private_script_matches_public_sha256'],False)
            self.assertIsNone(result['private_scores_match_frozen_sha256'])
            self.assertEqual(result['runtime_file_evidence'][0]['allowlisted_basename'],'requirements.txt')
            for word in ['private-directory','personal-name','private contents','original-machine','unrecognized']:
                self.assertNotIn(word,str(result))
            records[1]['sha256']='0'*64
            with self.assertRaises(d.recovery.RecoveryError):
                d.archive_metadata(z,z.infolist(),records)

    def test_key_order_distinguished_from_matching_targets_and_support(self):
        panel=pd.DataFrame(dict(speciesKey=[1,2,2,1],cell=[3,4,4,3],period=[1,1,2,2],
                                recording_probability10=[.2,.3,.8,.9]))
        preds=pd.DataFrame(dict(speciesKey=[2,1],cell=[4,3],baseline_probability=[.3,.2],
                                target=[.8,.9]))
        r=d.transition_evidence(panel,preds)
        self.assertTrue(r['same_key_set'] and r['target_exact'] and r['baseline_exact'])
        self.assertFalse(r['same_key_order'])
        self.assertFalse(d.transition_evidence(panel,preds.iloc[:1])['same_key_set'])

    def test_archived_rescore_uses_all_groups_and_detects_stale_scores(self):
        import numpy as np
        n=208482
        idx=np.arange(n)
        f=pd.DataFrame(dict(speciesKey=idx%155,species=idx%155,cell=idx//155,
                            target=(idx%13)/13,baseline_probability=(idx%7)/7))
        for model in d.check.MODELS:
            f[model]=f.baseline_probability
        metrics={model:d.check.scores(f,model)[0] for model in d.check.MODELS}
        audit=d.rescore(f,metrics)
        self.assertEqual(len(audit.comparisons),49)
        self.assertTrue(audit.passed)
        metrics[d.check.MODELS[-1]]['species_equal_mse']+=.01
        self.assertFalse(d.rescore(f,metrics).passed)

    def test_parquet_version_allowlist(self):
        self.assertTrue(d.WRITER.fullmatch('parquet-cpp-arrow version 23.0.1'))
        self.assertFalse(d.WRITER.fullmatch('parquet-cpp-arrow version 23.0.1/private/path'))
        self.assertFalse(d.VERSION.fullmatch('private text'))


if __name__=='__main__':
    unittest.main()
