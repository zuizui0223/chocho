# Frozen test: is native vs introduced plant-area HOSTS evidence overlap lower than source-reporting margins predict?

**Frozen 2026-10-27 before computing the new host-pair overlap null statistics.** Source original and geographic status classification are unchanged from source-audited chocho PR #76 and reporting-bias PR #77; this work is posthoc, exploratory and independent of the GEB submission.

## Central target and risk

The previous **239 focal butterfly** original-host database yielded 156 distinct butterfly×accepted-plant pairs with HOSTS Location reports from WCVP-introduced plant regions, but only 11 of those same pairs also had a separately geographically qualified WCVP-native plant-region report. An absent native report is **not** an evolved novel host adoption. Before giving this number mechanistic meaning, ask whether the *degree of overlap of documented native and introduced source host lists* is exceptional relative to **the same source plant list and number of documented pairs per butterfly**.

The analysis is an artificial conditional **source-recording** reference, not a biological randomized null for colonization or host evolution.

## Source authority and fixed estimand

- Original 239 accepted butterfly species, exactly 140,485 HOSTS raw rows / 86,826 binomial rows, original accepted HOSTS-to-WCVP plant IDs and geography of documented original `Location`. Source original inner archive SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`, original country WCVP/TDWG 969 disjoint equal-area atoms and their SHA, classification from PR #76.
- Only original 95%-whole-country status + source-qualified Taiwan/Hawaii/West Malaysia localities. Keep all other Location statuses explicitly unresolved, do not promote them to source-native or introduced. 80% country-status threshold is the *only* alternative geography run.
- For each species `s`, let `U_s` be its unique WCVP-accepted host plant IDs in frozen `insect_host_accepted.csv` (global catalogue, accepted IDs). Let `N_s` be the **subset of those host IDs with at least one HOSTS record geographically classifiable as WCVP native-plant-area**, and `I_s` similarly for WCVP introduced-plant-area. `X_s` is the subset documented in WCVP-*unlisted*-plant-area as a negative reporting control.
- Main response `O = sum_s |N_s ∩ I_s|` distinct butterfly×plant pairs in both WCVP origin-status classes. Sanity check `sum_s|I_s|=156`, `O=11`, `# species with I>0=62`, and previous unique location triple count `3586`. Also track `sum_s|N_s ∩ X_s|`, `sum_s|X_s|` and presence recording. Do not treat HOSTS rows, multiple same-species countries, or plant pair locations as independent ecological replicates.

## Artificial conditional source references — outcome-blind

1. **Within-butterfly host-list null (main)**: fix `U_s`, `n_s=|N_s|`, `i_s=|I_s|`, and actual source introduced set identities only as fixed margins. Under an exchangeability model that labels exactly `n_s` of the `U_s` host plant IDs as native-documented uniformly, `|N_s∩I_s|∼Hypergeom(U_s,i_s,n_s)`. Convolve all exact per-species hypergeometric distributions. Report observed, exact conditional mean, central 95% reference interval, left/right *mathematical* tails (never ecological p value), number of informative species and constrained variance.
2. **Host plant family-conditioned source null**: within *each butterfly×original accepted host plant family*, fix globally documented host species count and both report margins; convolve hypergeom per species×plant family. Source uncertain host-family records must be a separate stratum or unscorable, not guessed.
3. **Generic source documentation control**: re-run (1) and (2) with WCVP-unlisted pair set `X_s` instead of introduced `I_s`, preserving native pair evidence. Under-reporting and broad published source bias can generate the same apparently-low documented native overlap for unlisted reports.
4. **Geographic-status sensitivity**: repeat all four nulls with the original predefined WCVP 80%-area country dominant threshold; do not recalibrate.
5. **Species concentration**: report number of host-list-informative butterfly species, species providing any overlap, total source-pair coverage, and taxon-genus contribution. A species with `N_s=0` contributes expectation zero and is NOT evidence of a new host adoption. Report an explicitly conditioned panel with `N_s>0, I_s>0`, which is not an unbiased sample of butterfly invasions.

## Information and effect gate

- If original 239-source identities, 156/11 margins, 3,586 triple identity or 95%-area origin labels change, fail instead of silently using another population.
- Distinguish an extreme mathematical artificial null tail from an ecological p-value; regional reporting, source host preference, and plant availability violate exchangeability. The family null partly controls host identity but may leave few permutable strata. A "compelling" result requires **at least 20 host-list-informative butterfly species** and **at least 5** expected co-documented pairs from the artificial null; otherwise do not interpret an apparent deficit/excess.
- Even with a strong conditional excess/deficit, call it **HOSTS native-versus-introduced documentation overlap** only; not new host use, host switching, successful invasive butterfly establishment, historical mechanism, selection, or ecological niche expansion.
- No raw original source butterfly×host×location triplets, per-plant accepted IDs, full WCVP table, original TDWG geometry or individual literature records in public Git or artifacts. Export cumulative source-only results and genus/species-level numerical *aggregate* reporting when privacy/licensing permits, not raw host identities.

This work tests whether the earlier "145 pairs with no native-source report" is merely what sparse native-area documentation would ordinarily produce, and whether source family and generic WCVP-unlisted controls change interpretation.
