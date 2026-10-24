# Gain–single-host support correlation is structurally expected: completed null verdict

**24 October 2026 | source-locked, two-chain exploratory test | primary finding: NO ADDITIONAL RECORDED-HOST SPECIFIC GEOGRAPHIC ASSOCIATION.**

## Decision at a glance
The previously reported association that larger proportional potential butterfly larval-resource expansion coincides with a smaller share of single-host-supported new regions is **real in the reconstructed catalogue but is NOT exceptional under a plant-family × continent fixed-degree spatial structural null**.

The observed coefficient is almost exactly reproduced by two independently seeded 199-draw null chains, both passing all predeclared numerical quality gates. The 191-butterfly host-taxonomy-adequate subset shows the same non-exceptional outcome.

**Do not promote "resource expansion creates demographic buffering" to a supported ecological conclusion.** Neither a new ecological selection mechanism nor genuine host redundancy protecting butterfly populations was identified. The conclusion is about potential plant host resource maps, not actual caterpillars or host-quality assurance.

## Question, frozen comparison, and source provenance
The observed species-level regression from [the exploratory gain–fragility result](butterfly_geographic_gain_fragility_20261024.json) was β = −0.220718729049... across 206 butterflies with newly opened regions (a source-conditional genus-bootstrap interval −0.298 to −0.165). This coefficient describes a within-catalogue association and alone does **not** test against a plant geography null.

Before the new null simulation, the [protocol](../../provenance/reviewer_defenses/protocols/GAIN_SUPPORT_PLANT_RANGE_NULL_20261024.md) fixed both numerical seeds, original source identities, candidate fauna, plant-graph constraints, predictor/outcome definitions and inferential stop rule. The null was selected **after** seeing the old observed coefficient, so it is an explicit **posthoc falsification test**, not a registered independent confirmation.

Exact sources:
- HOSTS pinned Git commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`
- WCVP pinned Git commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`
- Frozen 239-species CSV SHA-locked resource and host metrics
- Frozen WGSRPD3 Level1 mapping source commit `52da7828aba9d461dd133c27b3bd7a4407161f54`
- Existing source builder, exact host–plant–region graph and reusable double-edge swaps in `scripts/analyze_butterfly_shared_host_structural_null.py`.

The fresh workflow reconstructed 239 butterflies, 365 geographic regions, 30,373 introduced host-plant×region edges, 26,530 native butterfly×region resources, 14,553 newly opened resources, **8,574** new regions with just one recorded introduced host, and original β values at single-species precision. Each of the 239 butterfly native/added resource counts and 206 positive-gain source-specific single-host shares matched the immutable Figure/source CSVs exactly.

## Null design
Each swap moves two introduced botanical plant–region edges within the same plant family and WGSRPD Level1 while excluding the plants' own native-range regions and any duplicate introduced edge. It preserves all focal butterfly–host identities, each introduced plant's per-continent region count, each geographic recipient's per-family introduced plant count and all native host areas. The source includes **28,670 eligible introduced edges** among 30,373. No cherry-picked plant hosts or extra organismal observations enter.

**Not conditioned upon** fine-scale climate, host quality, habitat, phenology, actual butterfly occupancy or successful reproduction. Butterfly added resource unit counts change as a mathematical result of swapping plant edges. Both chains retained 206 species with U>0 on every draw, so observed and null butterfly panel identity was identical.

The species-equal adjusted OLS response is `fraction of new regions with exactly one documented introduced host` and its predictor is `log(current/native resource regions)`; nuisance controls are baseline native resource count, recorded host count, and host-family breadth. Species cluster bootstrap was the previous analysis; **the new null statistic is the same regression coefficient recomputed at every shuffled plant graph**. Distinguish these methods.

## Main source-verified result

| Panel | Observed coefficient β | Seed 20261024 null median (2.5–97.5%) | Seed 20261025 null median (2.5–97.5%) | Evidence assessment |
|---|---:|---:|---:|---|
| **206 butterflies gaining resources** | **−0.220719** | −0.222145 (−0.227286, −0.217066) | −0.222255 (−0.227699, −0.216364) | Observed inside both null intervals |
| **191 host-taxonomy-adequate butterflies** | **−0.271129** | −0.270004 (−0.278527, −0.262430) | −0.270131 (−0.278395, −0.261019) | Observed inside both null intervals |

Observed-minus-null median β is only **+0.001426 / +0.001536** for the full 206 species and **−0.001125 / −0.000997** for the adequate-host subset. These effect departures are tiny on the coefficient scale. Descriptive serial empirical tail fractions also show no extreme observed departure (206-species left tails 0.690/0.705; adequate-species left tails 0.385/0.450), but **these are NOT independent permutation p-values**.

## Numerical qualification and exact source replay

Two independent MCMC chains used seeds 20261024 and 20261025, each with 5 × 30,373 burn proposals and 2 × 30,373 swap attempts between **199** states. Both source movement rates were >30%; burn accepted 20,270 / 20,547 swaps; lowest per-sampling-interval accepted 7,900 / 7,826. All fixed native plant exclusions and plant × continent and region × plant-family degree margins were verified against original source after burn and throughout sampling.

Approximate positive-sequence ESS for the full-panel coefficients was **119.8 / 150.4**, and for adequate-host coefficients **97.1 / 150.0**, all above the prospectively required 30. The two-chain median agreement test passed for both responses. The exploratory MCMC is therefore **numerically qualified** under the declared constraints.

[Full numeric source receipt](../../provenance/reviewer_defenses/results/butterfly_gain_support_plant_null_20261024.json) (contains each draw's slope, panel consistency diagnostics, original source identities and chain quality).
[Source workflow, successful](https://github.com/zuizui0223/chocho/actions/runs/40081217645).

## What this adds, what must be retracted and what survives

**Retract/promote-down:**
- Do not say "the strongest geographic gain is exceptionally linked to multiple-host buffering." The observed cross-species association is already reproduced when introduced host ranges are randomized subject to major source constraints.
- The prior 5,000-genus-bootstrap confidence interval excluding zero is **not** evidence against this stronger structural baseline.
- Do not infer larger larval production, stronger population stability, adaptive host portfolios or ecological homogenization of actual butterfly abundances.

**Retain as source-conditional structural facts, not novel mechanisms:**
- 206/239 species gained potential larval host geography; total geographic opportunity +54.9% in the original fixed HOSTS/WCVP representation.
- 8,574/14,553 newly opened butterfly×region cells are singly supported by one documented introduced host; in the old aggregate family×continent plant randomization, this *fraction* was itself near null expectation.
- The species-level negative gain–single-support slope is a **catalogue-architecture pattern** and its nonexceptionality is a useful falsification of the intuitive gain-versus-vulnerability mechanism.
- BCE/Clarke supplementary host-list sensitivity altered recorded candidate support intensity, but the positive slopes there were for **different response metrics**. They do not rescue this observed-original single-host-fraction null or prove local host substitution.

## Implication for the GEB single-paper strategy

The *main* submission claim should remain **partner redistribution changes butterfly potential larval-resource geography across the specialization spectrum**, with source-sensitive geographic novelty and disproportionate plant contributions as supporting ecological context. This plant-null analysis is best retained as a **transparent exploratory negative control** in a limited supplementary note if relevant to reviewer interpretation; it should **not be marketed as a second independent butterfly ecology discovery**.

Substantive ecological next data must come from source-independent local larval feeding and survival, dated host arrival, or independently validated community abundance/occupancy—none are inferred from this randomized host map. Stop further arbitrary host/country strata and do not retune null thresholds to manufacture extra signal.

## Branch and files
- [Frozen prospective-for-this-null protocol](../../provenance/reviewer_defenses/protocols/GAIN_SUPPORT_PLANT_RANGE_NULL_20261024.md)
- `scripts/analyze_butterfly_gain_support_plant_null_20261024.py`: replay exact frozen source with 2×199 draws
- `tests/test_gain_support_plant_range_null_20261024.py`: source-independent bit/margin/ESS guards
- `.github/workflows/gain-support-plant-null-20261024.yml`: CI rebuild from pinned HOSTS and WCVP
- [Raw prior gain–support exploratory analysis](../../docs/exploratory/BUTTERFLY_GAIN_AND_HOST_SUPPORT_ARCHITECTURE_20261024.md)
- GEB `main` and PR #38 remain unchanged.
