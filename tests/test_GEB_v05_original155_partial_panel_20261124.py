import unittest
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))

import numpy as np
import pandas as pd
import replay_GEB_v05_original155_partial_panel_20261124 as a

class PartialSourceBoundaries(unittest.TestCase):
    def test_must_not_construct_missing_native_introduced_counts_from_target(self):
        self.assertEqual(set(a.MISSING), {
            'outputs/forecast_candidate_disposition.csv',
            'outputs/temporal_host_family_support.csv',
            'work/remote_sensing/host_region_counts.csv',
        })
        self.assertEqual(len(a.SOURCE_FILES),10)

    def test_numeric_and_null_column_comparison(self):
        left = pd.Series([np.nan, 0., 1., 0.25])
        right = pd.Series([np.nan, 0., 1., .250000001])
        self.assertEqual(a.compare_column(left,right,'prob')['mismatched'],0)
        wrong = pd.Series([np.nan, 0., 1., .2])
        self.assertEqual(a.compare_column(left,wrong,'prob')['mismatched'],1)

    def test_host_family_selection_uses_positive_observed_credit(self):
        data = []
        for family in ('A','B'):
            for i in range(12):
                data.append({'family':family,'species':str(i),
                             'region':'r'+str(i),'share':1. if family=='A' else 0.})
        frame = pd.DataFrame(data)
        self.assertEqual(a.group_support(frame, only_positive=True),['A'])
        self.assertEqual(a.group_support(frame, only_positive=False),['A','B'])

if __name__ == '__main__':
    unittest.main()
