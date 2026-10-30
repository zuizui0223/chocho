# Broad multispecies landmass extension: temporal snapshot

This is a draft continuation of PR #101, separate from the93-species island panel in PR #102. Candidate225 / eligible155 butterfly species are retained. Read `barrier_first_results_ja.md` for the actual result and limits.

B2 adds history updates and landmass membership to M5. Overall temporal MSE barely changes. Island-subset estimates improve relative to M5, but crossed species/block intervals include zero; the additional landmass contrast against B1 is near zero. No positive causal island effect is established.

The2,000-resample intervals condition on fixed predictions. Excluding the two dominant sampled landmasses is post-hoc. The2019 port inventory is not used as an earlier forecast input. DEM, sea routes, shipping and independently dated arrivals are not completed by this snapshot. Spatial fits are now complete: B2 shows no supported incremental predictive improvement against M5 or B1. DEM fits and route validation remain in progress.

All225 disposition is in the parent PR #101 snapshot. Scripts are archived analysis sources, not turnkey: they require the existing local directory structure, large occurrence inputs and restricted host-derived panels. Raw observations, people, plant source tables and models are not included.

Run `python verify_snapshot.py` to verify published bytes and aggregate arithmetic only; it does not independently rerun the biological analysis. Derived source geography uses GSHHG2.3.7 (LGPL); no shoreline files redistributed.
