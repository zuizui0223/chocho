# Frozen source-matched island–mainland butterfly-resource convergence contrast

**Protocol date:** 2026-11-18, before examining island/mainland pairwise resource Jaccard deltas or restricted-null scores. **Exploratory macroecological follow-up** (the global +0.18510 versus Level1-restricted null +0.18159 result was already known). No independent data or preregistration claim.

## Question

Does the *same fixed 239-butterfly potential larval-host resource geography* become more similar between distinct botanical islands than expected from merely adding the observed number of new butterfly×region opportunities, especially relative to mainland–mainland pairs? A difference in plant species pool and present-day accepted host distributions **cannot by itself identify historical colonization, larval quality, butterfly range change or actual biotic homogenization**.

## Original sources and science units

- Freeze the exact prior PR #118 commit `bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2`, source original `provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/beta_baseline_v3.matrix.csv.gz` and `beta_baseline_v3.json`. Original source-verified WCVP v13 + HOSTS 239 species; 369 WGSRPD3 official botanical regions; **26,530** native, **14,553** introduced-added and **41,083** contemporary resource combination bits. Retain the original species identity bit order, region codes, source checksum and the explicit `native & added == 0` rule. The original matrix is source-derived and committed on the archived analysis branch, not a new butterfly observation.
- Use the existing main `results/island-mainland-resource-q1/region_metrics.csv` for original fixed **126 island and 237 mainland** botanical region codes. Of the 369 official region codes, six have no original binary island/mainland assignment; leave these **unclassified** and never guess their assignment.
- Use the original official `data/frozen/wgsrpd3_level1_lookup_v0.1.json` from archived source to define the TDWG Level1 strata, not user-inferred continents.
- Restrict the **primary response** to the 355 regions with at least one native accepted known-host butterfly resource (the exact original global common-support cohort: **62,835 unordered pairs**). Describe excluded zero-native resource regions and report them by island/mainland status; do **not** call their original absence true biological absence. Do not impute missing-source regions as native absence.
- Unit of Jaccard is a set of candidate butterfly **species having ≥1 documented host plant within that botanical region**, not host *plant taxa* and not a GBIF-indexed butterfly assemblage. Native and contemporary plant maps represent a **counterfactual nested snapshot**, not dates.

## Planned results, held fixed before observing class-specific scores

For the original native matrix N and current C=N|A, compute paired island–island (II), mainland–mainland (MM), island–mainland (IM), and **same-Level1 subgroups of II and MM**:
1. Numbers of eligible regions and unique non-independent pairs; mean native Jaccard, mean current Jaccard and `delta = mean(current-native)`, sample-pair median delta, fraction positive and number of undefined comparisons.
2. Absolute raw II–MM delta difference `D_II-D_MM`, **not** a risk ratio or a multiplication by plant richness. Also retain the original global 355-region mean delta **0.18510150276498433**.
3. For every null draw, compute exactly those original class-specific delta statistics and their II–MM contrast on the same **fixed original region pairs**. Avoid counting each source dyad as independent for significance.

## Nulls (fixed margins first, morphology and shipping not entered)

Keep **all** original N bits fixed; preserve each region's exact added butterfly count, each butterfly's total added regions, and forbid overlap with N using the original source `analyze_global_resource_gain_fixed_margins.py` constrained bipartite 2×2 switches. Use source-preserving fixed `attempts_per_draw = 5×14,553 = 72,765`, **199 sequential MCMC draws**. Draws are correlated Monte Carlo states, not independent biological replicates.

- **A (primary):** swaps only within original TDWG Level1, preserving also each butterfly's gain count within each Level1. Exact seed **41131**, reproducing the original fully source-locked global within-Level1 null mean `0.18158108553894184` (tolerance 1e-8). This is the main geographic null; refuse to interpret island contrast if that source recheck fails.
- **B (strict sensitivity):** swaps only within `TDWG Level1 × original island/mainland/unclassified`; this additionally preserves each species' count of gains in islands versus mainland within each Level1, as well as regional gain counts. Exact seed **20261118**, same 199 sequential draws and attempts. Report accepted swaps and fail/flag strata with insufficient scope, not replace strata after seeing results.
- For each cohort, output observed delta; null mean/median/2.5–97.5% *draw quantiles*; observed-minus-null median; and a descriptive high-tail fraction including an add-one correction. These are **not p-values for independent plant introductions** and do not control phylogeny/spatial autocorrelation inside Level1.
- Evaluate the II–MM contrast under both nulls rather than treating II and MM observed gains as independent samples. If an apparent island excess vanishes under source B, interpret that as differences in **which butterflies enter islands**, not an additional intra-island spatial assembly pattern.

## Release/interpretation guards

Never convert a potential resource into confirmed host use, butterfly abundance, GBIF colonization, human-beneficial introduction, or movement through ports/ocean circuits. Existing high-confidence +54.9% global resource addition, +118.9% island-relative baseline amplification, and photographed Canary-endemic *Pieris cheiranthi* using introduced *Tropaeolum majus* remain intact as different estimates. The independent, conditional GBIF 224-butterfly observed community analysis showed only +0.157% MSE improvement from added known-host composition and poor region transfer; it must not be displaced with this source-derived outcome.

Reproduce exact original regional source sums and code lists before fitting. Source work is anchored to PR #118's original archived Git commit; derived **summary-only** receipts and no raw original HOSTS/WCVP files may be added to the GEB manuscript branch or main. Preserve both positive and adverse results in supplementary evidence. Manuscript changes only after an executed source-verified CI.

**Source acknowledgments:** Original whole-global resource Jaccard analysis [PR #118](https://github.com/zuizui0223/chocho/pull/118), independent GBIF [source readout](../../provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/GBIF_GLOBAL224_EMPIRICAL_READOUT_20261101_ja.md), current fixed 239 resource [GEB evidence ledger](../../docs/synthesis/POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md).
