# Migration provenance

Source repository: zuizui0223/TTF

Frozen source commit: `1a112334cca2f2ef5e234c3ae1fc1a80b8266956`

Migration rule:

1. Copy all files whose path contains `butterfly` or `lepidoptera` under manuscript, scripts, tests, docs, benchmarks, results, and workflows.
2. Copy `src/` initially as a conservative dependency closure.
3. Copy Python project/test configuration when present.
4. Establish a passing paper-specific test suite in this repository.
5. Only then prune unrelated source modules and remove the migrated ecology surface from TTF.

This two-stage rule prevents a repository split from silently changing the frozen analysis.
