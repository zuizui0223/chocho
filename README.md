# chocho

Reproducible research repository for the butterfly specialization ecology manuscript targeted to *Global Ecology and Biogeography*.

## Provenance

This repository is being split from `zuizui0223/TTF` so that the ecological paper has a small, paper-specific reproducibility surface while TTF remains a methods repository.

The initial migration is pinned to TTF commit:

`1a112334cca2f2ef5e234c3ae1fc1a80b8266956`

The import intentionally copies the butterfly/Lepidoptera manuscript, scripts, tests, frozen/exploratory receipts, workflows, and the current `src/` implementation needed to preserve dependency closure. After import, unused TTF modules will be removed only after the paper-specific test suite establishes the dependency boundary.

## Migration status

Do not cite this repository yet. The migration snapshot must first pass the butterfly/Lepidoptera tests, dependency cleanup, license decision, citation metadata, and a frozen release tag.
