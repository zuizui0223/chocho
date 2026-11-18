# One-paper source-receipt consistency check

Run `python3 scripts/verify_integrated_positive_evidence.py` (Python 3.10+ stdlib only).

This gate re-derives **arithmetic and numerical relationships in the frozen source receipts** rather than refitting raw HOSTS, rWCVPdata, GBIF or Circuitscape. It safeguards source-typed claims across 239 known-host butterflies; 126 island/237 mainland botanical regions (+6 excluded mixed); 98/123/126 port-flora cohorts; 104 flora coastal terminals; 51 butterfly island groups; 155 known-butterfly forecast species; independently listed 41 living/18 endemic island butterflies; and 224 global GBIF butterfly-record taxa. 18 source files are SHA-256 fingerprinted in the workflow artifact.

**Not equivalent to a full reproducible computation from raw third-party data.** The original full source scripts and upstream original data/CI provenance remain in their pinned PR branches and verified GitHub Actions. A separate selective-code integration is required before retiring them.

Specific inference safeguards: an original 126-region port β+0.221 source receipt lists **native host RESOURCE count** as a covariate, not full native-flora richness; a separate 123-region full-flora regression β+0.233 uses genuine native flora richness. Their sample/control models must not be silently interchanged. The +10.70pp matched island added-share contrast comes from **15 level1×area strata**; an earlier **7 Level1-only** contrast is a different statistic. Sea shortest paths, uniform-ocean circuit resistance and observed NOAA currents are different predictors, not equivalent measured seed-delivery paths. Source-documented wild/garden larvae do not prove demographic establishment.

Run by `.github/workflows/verify-integrated-positive-evidence.yml`, whose Actions artifact lists group-level outputs and exact file SHA-256. On failure preserve the historical output and correct the source data/claim provenance rather than tweaking a cutoff to restore positivity.
