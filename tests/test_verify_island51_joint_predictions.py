"""Catch silent errors in the independent saved-prediction auditor."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'results/island51_host_circuit_joint_20261106/verify_saved_predictions.py'
TARGETS = ('new_butterfly_jaccard', 'earlier_butterfly_jaccard')
MODELS = ('G', 'H', 'C', 'HC', 'I')
COMPARISONS = {
    'G_to_H': ('G', 'H', 'contextual'),
    'G_to_C': ('G', 'C', 'contextual'),
    'C_to_HC': ('C', 'HC', 'primary'),
    'H_to_HC': ('H', 'HC', 'primary'),
    'HC_to_I': ('HC', 'I', 'secondary_interaction'),
}


def fixtures():
    # Hand-chosen errors differ between the two held-out endpoints of a pair.
    # Averaging predictions before scoring incorrectly divides G's MSE by four.
    ids = ['1', '2', '3']
    pairs = [dict(island_a=a, island_b=b, **{t: .5 for t in TARGETS})
             for a, b in [('1', '2'), ('1', '3'), ('2', '3')]]
    predicted = {'G': [.3, .5, .7], 'H': [.4, .6, .8],
                 'C': [.2, .4, .6], 'HC': [.4, .5, .6], 'I': [.4, .5, .6]}
    square = {'G': [.04, 0., .04], 'H': [.01, .01, .09],
              'C': [.09, .01, .01], 'HC': [.01, 0., .01], 'I': [.01, 0., .01]}
    absolute = {'G': [.2, 0., .2], 'H': [.1, .1, .3],
                'C': [.3, .1, .1], 'HC': [.1, 0., .1], 'I': [.1, 0., .1]}
    rows = []
    tables = {name: [] for name in ('model_fold_losses', 'comparison_fold_losses',
                                    'model_summary', 'comparison_summary')}
    for target in TARGETS:
        for model in MODELS:
            for index, pair in enumerate(pairs):
                for heldout in (pair['island_a'], pair['island_b']):
                    k = int(heldout) - 1
                    rows.append(dict(target=target, model=model, heldout_island_id=heldout,
                                     source_pair_index=index,
                                     pair_id=f"{pair['island_a']}:{pair['island_b']}",
                                     island_a=pair['island_a'], island_b=pair['island_b'],
                                     observed=.5, predicted=predicted[model][k],
                                     squared_error=square[model][k], absolute_error=absolute[model][k]))
            for k, heldout in enumerate(ids):
                tables['model_fold_losses'].append(dict(target=target, model=model,
                    heldout_island_id=heldout, train_pairs=1, test_pairs=2,
                    MSE=square[model][k], MAE=absolute[model][k]))
            tables['model_summary'].append(dict(target=target, model=model,
                MSE=sum(square[model])/3, MAE=sum(absolute[model])/3,
                n_folds=3, n_predictions=6))
        for label, (ref, aug, role) in COMPARISONS.items():
            for k, heldout in enumerate(ids):
                r, a = square[ref][k], square[aug][k]
                tables['comparison_fold_losses'].append(dict(target=target, comparison=label,
                    reference_model=ref, augmented_model=aug, heldout_island_id=heldout,
                    train_pairs=1, test_pairs=2, reference_mse=r, augmented_mse=a,
                    delta_mse_augmented_minus_reference=a-r,
                    reference_mae=absolute[ref][k], augmented_mae=absolute[aug][k],
                    improved=a<r, worse=a>r, tied=a==r))
            r, a = sum(square[ref])/3, sum(square[aug])/3
            tables['comparison_summary'].append(dict(comparison=label, target=target,
                comparison_role=role, reference_model=ref, augmented_model=aug,
                MSE_reference=r, MSE_augmented=a,
                relative_MSE_reduction_percent=100*(r-a)/r,
                MAE_reference=sum(absolute[ref])/3, MAE_augmented=sum(absolute[aug])/3,
                absolute_MSE_gain=r-a,
                islands_improved=sum(x<y for x,y in zip(square[aug],square[ref])),
                islands_worse=sum(x>y for x,y in zip(square[aug],square[ref])),
                islands_tied=sum(x==y for x,y in zip(square[aug],square[ref]))))
    return ids, pairs, rows, tables


def receipt_fixture(directory):
    repo, results = directory / 'repo', directory / 'results'
    results.mkdir()
    source = dict(input_lock_sha256='a'*64, protocol_sha256='b'*64,
                  implementation_plan_sha256='c'*64,
                  frozen_at_utc='2026-11-06T12:10:44+00:00')
    validation = dict(source, completed_before_joint_fit_at_utc='2026-11-06T12:11:00+00:00')
    prior = dict(status='PRIOR_NATIVE_G_AND_C_FOLDS_AND_SUMMARIES_REPRODUCED')
    files = ['heldout_predictions.csv.gz', 'model_fold_losses.csv', 'comparison_fold_losses.csv',
             'model_summary.csv', 'comparison_summary.csv']
    for name in files:
        (results / name).write_bytes(('fixture '+name).encode())
    for name, data in [('source_validation.json', validation), ('prior_circuit_reproduction.json', prior)]:
        (results / name).write_text(json.dumps(data))
        files.append(name)
    implementation = {}
    for name in ['run_island51_host_circuit_joint.py', 'exploratory_island51_circuit_network_holdout.py',
                 'run_island51_observed_ocean_circuits.py', 'island51_ocean_circuit.py']:
        path = repo / 'scripts' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('# fixture '+name+'\n')
        implementation['scripts/'+name] = hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = dict(schema='chocho_fixed_island51_joint_host_circuit_run_v1',
                   status='FIVE_FIXED_JOINT_MODELS_AND_ALL_HELDOUT_PREDICTIONS_COMPLETE',
                   source_validation=validation, prior_circuit_reproduction=prior,
                   completed_at_utc='2026-11-06T12:12:00+00:00',
                   implementation_sha256=implementation,
                   outputs_sha256={name: hashlib.sha256((results/name).read_bytes()).hexdigest() for name in files})
    (results / 'joint_run_receipt.json').write_text(json.dumps(receipt))
    return repo, results, source


class IndependentPredictionAuditTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.exists(), 'Independent saved-prediction verifier is not implemented')
        spec = importlib.util.spec_from_file_location('joint_prediction_audit', SCRIPT)
        self.audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.audit)
        self.ids, self.pairs, self.rows, self.tables = fixtures()

    def test_losses_are_scored_before_endpoint_predictions_are_averaged(self):
        result = self.audit.reconstruct(self.rows, self.pairs, self.ids)
        g = next(x for x in result['model_summary']
                 if x['target'] == TARGETS[0] and x['model'] == 'G')
        self.assertAlmostEqual(g['MSE'], 2/75, places=14)
        self.assertNotAlmostEqual(g['MSE'], 1/150, places=14)
        check = self.audit.compare_tables(result, self.tables)
        self.assertEqual(check['model_fold_losses']['checked_rows'], 30)
        self.assertEqual(check['comparison_summary']['checked_rows'], 10)
        interaction = next(x for x in result['comparison_summary'] if x['comparison'] == 'HC_to_I')
        self.assertEqual((interaction['islands_improved'], interaction['islands_worse'], interaction['islands_tied']), (0, 0, 3))

    def test_duplicate_endpoint_cannot_replace_missing_endpoint(self):
        self.rows[-1] = copy.deepcopy(self.rows[-2])
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.audit.reconstruct(self.rows, self.pairs, self.ids)

    def test_fold_must_hold_out_an_actual_pair_endpoint(self):
        self.rows[0]['heldout_island_id'] = '3'
        with self.assertRaisesRegex(ValueError, 'endpoint'):
            self.audit.reconstruct(self.rows, self.pairs, self.ids)

    def test_self_consistent_modified_observation_is_rejected_against_source(self):
        self.rows[0].update(observed=.6, squared_error=.09, absolute_error=.3)
        with self.assertRaisesRegex(ValueError, 'observed'):
            self.audit.reconstruct(self.rows, self.pairs, self.ids)

    def test_individual_loss_is_recomputed_instead_of_trusted(self):
        self.rows[0]['squared_error'] = .09
        with self.assertRaisesRegex(ValueError, 'squared_error'):
            self.audit.reconstruct(self.rows, self.pairs, self.ids)

    def test_pair_index_and_canonical_identifier_cannot_disagree(self):
        self.rows[0]['pair_id'] = '1:3'
        with self.assertRaisesRegex(ValueError, 'pair_id'):
            self.audit.reconstruct(self.rows, self.pairs, self.ids)

    def test_wrong_percentage_denominator_is_detected(self):
        result = self.audit.reconstruct(self.rows, self.pairs, self.ids)
        self.tables['comparison_summary'][0]['relative_MSE_reduction_percent'] = -27.2727272727
        with self.assertRaisesRegex(ValueError, 'relative_MSE_reduction_percent'):
            self.audit.compare_tables(result, self.tables)

    def test_omitted_adverse_fold_is_detected(self):
        result = self.audit.reconstruct(self.rows, self.pairs, self.ids)
        self.tables['comparison_fold_losses'].pop(1)
        with self.assertRaisesRegex(ValueError, 'coverage'):
            self.audit.compare_tables(result, self.tables)

    def test_tie_is_not_counted_as_worse(self):
        result = self.audit.reconstruct(self.rows, self.pairs, self.ids)
        row = next(x for x in self.tables['comparison_fold_losses'] if x['comparison'] == 'HC_to_I')
        row['worse'] = True
        with self.assertRaisesRegex(ValueError, 'worse'):
            self.audit.compare_tables(result, self.tables)

    def test_receipt_hashes_bind_outputs_to_the_recorded_implementation(self):
        self.assertTrue(hasattr(self.audit, 'verify_receipt_links'), 'Receipt hash-link audit is missing')
        with tempfile.TemporaryDirectory() as tmp:
            repo, results, source = receipt_fixture(Path(tmp))
            check = self.audit.verify_receipt_links(results, source, repo=repo)
            self.assertEqual(check['output_hashes_verified'], 7)
            self.assertEqual(check['implementation_hashes_verified'], 4)
            self.assertTrue(check['freeze_precedes_preflight_and_completion'])

    def test_modified_predictions_cannot_retain_an_old_receipt(self):
        self.assertTrue(hasattr(self.audit, 'verify_receipt_links'), 'Receipt hash-link audit is missing')
        with tempfile.TemporaryDirectory() as tmp:
            repo, results, source = receipt_fixture(Path(tmp))
            (results/'heldout_predictions.csv.gz').write_bytes(b'changed predictions')
            with self.assertRaisesRegex(ValueError, 'output.*hash'):
                self.audit.verify_receipt_links(results, source, repo=repo)

    def test_changed_fit_code_cannot_claim_the_previous_recorded_implementation(self):
        self.assertTrue(hasattr(self.audit, 'verify_receipt_links'), 'Receipt hash-link audit is missing')
        with tempfile.TemporaryDirectory() as tmp:
            repo, results, source = receipt_fixture(Path(tmp))
            (repo/'scripts/run_island51_host_circuit_joint.py').write_text('# changed\n')
            with self.assertRaisesRegex(ValueError, 'implementation.*hash'):
                self.audit.verify_receipt_links(results, source, repo=repo)

    def test_preflight_before_the_frozen_plan_is_rejected(self):
        self.assertTrue(hasattr(self.audit, 'verify_receipt_links'), 'Receipt hash-link audit is missing')
        with tempfile.TemporaryDirectory() as tmp:
            repo, results, source = receipt_fixture(Path(tmp))
            validation = json.loads((results/'source_validation.json').read_text())
            validation['completed_before_joint_fit_at_utc'] = '2026-11-06T12:00:00+00:00'
            (results/'source_validation.json').write_text(json.dumps(validation))
            receipt = json.loads((results/'joint_run_receipt.json').read_text())
            receipt['source_validation'] = validation
            receipt['outputs_sha256']['source_validation.json'] = hashlib.sha256((results/'source_validation.json').read_bytes()).hexdigest()
            (results/'joint_run_receipt.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'freeze.*preflight'):
                self.audit.verify_receipt_links(results, source, repo=repo)


if __name__ == '__main__':
    unittest.main()
