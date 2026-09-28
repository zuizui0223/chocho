# Test suite

The offline test suite protects both scientific values and the submission package.

## Scientific/reconstruction tests

Tests for host-resource reconstruction, anthropogenic expansion, portfolio architecture, within-family hierarchy and the independent climate analysis.

## Manuscript/release tests

- `test_butterfly_specialization_manuscript_bundle.py` — claim values, inference boundaries, narrative/figure order and manuscript metadata.
- `test_butterfly_specialization_submission_bundle.py` — blinded/submission package checks.
- `test_reproducibility_snapshot_integrity.py` — frozen source identities and hashes.
- `test_release_metadata.py` — v1 release identifiers.
- `test_repository_layout.py` — repository structure and independent package namespace.

Run all tests with:

```bash
python -m pip install -e ".[test]"
pytest -q
```
