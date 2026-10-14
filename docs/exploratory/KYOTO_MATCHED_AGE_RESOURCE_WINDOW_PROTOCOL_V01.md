# Kyoto food-demand shape: matched biological-age falsification gate

**Date:** 2026-10-14. **Classification:** POST-HOC DIAGNOSTIC ONLY, frozen before the dedicated matching script outcome; not preregistered before viewing the 2016 feeding series. **No new butterfly-fitness claim.** Work exclusively in \`theory/temporal-resource-competition-v01\`; GEB PR #38 and main unchanged.

## Why this comparison matters

The earlier analysis of Hashimoto & Ohgushi (2023) independently recovered Table 1 total late-instar dry-leaf intake (+0.77% for *Atrophaneura alcinous*, +0.41% for *Sericinus montela*) and reported a 2-day difference in median half-consumption age and 4-day difference in the calendar-relative peak. Those totals and medians are informative but cannot tell whether the difference is **within a common observed life-history window** or comes from the fact that S was observed at day 22–24 when most A larvae no longer generated matched 2-day feeding observations.

An independent experimental resource-timing hypothesis must be supported by nontrivially different normalized food-demand curves at **the same biological age**, and by replication across multiple maternal groups and host environments. The one-cohort-per-species source cannot supply the latter.

## Frozen computation before new profile estimates

- Use the original 2016 author CSV and MD5s from the prior audit; do not modify the source-data parsing, dry-mass formula or adjacent 2-day inclusion rules.
- Map dates to original founder-group start date **per species**, not absolute dates (species founder dates differ by 16 days).
- Require exact individual IDs with one valid original 2-day intake observation at each age-bin endpoint **8, 10, 12, 14, 16, 18 days**. Missing/invalid bins exclude that individual from this *descriptive* complete-window comparison; **report the exclusions** and retain the full original sample in the v01 reconstruction. This matched 6-bin window is chosen because the source-audited interval counts show near-complete support before the A sample falls to 8 at day 20 and 3 at day 22. It is a post-hoc completeness decision, not a biologically optimized response window.
- For each included larva, calculate six-bin percentages \`p_i(t)=intake_i(t)/sum_{t=8,10,...,18} intake_i(t)\`; average **each larva's share equally within its species**, rather than giving the biggest larvae more weight.
- Quantify observed curve overlap \`O=sum_t min(mean_p_A(t),mean_p_S(t))\`, total variation \`TV=1-O\` and *average date of food intake* \`mu=sum_t t*mean_p(t)\`.
- For an illustrative early/late split **fixed at the six-bin midpoint**, define early=8,10,12 and late=14,16,18. Report \`a=early_share_A\`, \`b=early_share_S\`, their difference, and whether the simple model's relative-disadvantage reversal is mathematically possible for some positive scarcity sensitivities. **No choice of the split after seeing outcome** and no biological emergence parameters inferred.
- Single-individual-delete descriptive diagnostics: recompute early-share difference with every included A larva omitted in turn, then every included S larva omitted in turn; report min/max and any sign change. These are finite-sample sensitivity bounds, **not** confidence intervals for populations or independent trials.
- Alternate coverage: 8–16 (five bins) and 10–18 (five bins), with fresh matched-complete selection. Report mean-time differences and TV only. Do not choose one after seeing which looks largest.
- Compare source per-larva cumulative totals against previous reconstruction as a guard, not as a new outcome. Any absence of complete participants, duplicates or non-positive total prevents promotion.

## Deliberate inferential limitations

The timing of the **unobserved foliage supply** is not in these CSVs. These larvae were raised in two distinct 2016 laboratory cohorts, not simultaneous competitors; previous 2012 mixed cages contain stages but no contemporaneous food-access time series and no adult-emergence outcomes. Mortality/attrition and earlier group rearing limit what counts as full lifetime food demand. The plant's phenological tissues, leaf age, maternal origin and the location of the larvae's actual seasonal oviposition are not independently sampled as resource windows.

**No stage-duration or time-shift significance test**, no 10 individual larvae as 10 independent maternal populations, no proof of resource-mediated competitive sign reversal. The 2016 totals are an already-published effect; author data reanalysis cannot be sold as new adult-demographic evidence.

## GO/NO-GO for an actual experimental model

**GO for an independent pilot** if the matched-age difference has a stable sign after leave-one-out, but only as a feasibility target; seek repeatable and ecologically important divergence in replicate maternal families and naturally comparable egg-cohort starts. Independent actual edible-leaf supply and dry-mass demand must be measured across larvae and phenological host stages.

**NO-GO for a competitive sign-reversal claim on the existing source**, regardless of overlap value, because resource timing was never experimentally manipulated and no flight-capable adult production was followed. A null or near-identical matched-window curve is additionally a reason to **deprioritize the Kyoto temporal inversion mechanism** instead of sweeping early/late boundaries to rescue an effect.

This preserves the independent GEB host-distribution submission and keeps Kyoto as an explicitly experimental ecological extension.