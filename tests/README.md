# Test map

The offline test suite protects both the scientific results and the submission/release package.

## Scientific-analysis tests

- resource-envelope reconstruction and host-resource utilities;
- anthropogenic resource expansion;
- host-contribution decomposition;
- within-family specialization hierarchy;
- climate quality gates, cross-fit and independent hypothesis test.

## Manuscript and release tests

- `test_butterfly_specialization_manuscript_bundle.py` — frozen values, claim boundaries, evidence order, data-scale communication and independent-paper wording.
- `test_butterfly_specialization_submission_bundle.py` — blinded/submission package checks.
- `test_reproducibility_snapshot_integrity.py` — pinned frozen inputs and exact hashes.
- `test_release_metadata.py` — v1 version/tag consistency.
- `test_repository_layout.py` — independent package namespace and repository organization.

Run the full offline suite with:

```bash
python -m pip install -e ".[test]"
pytest -q
```
