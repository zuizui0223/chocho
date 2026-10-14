# Kyoto two butterflies: dated larval intake windows, source reconstruction and scientific gate

**2026-10-14 | SOURCE-VERIFIED POST-HOC DESCRIPTIVE REANALYSIS. No randomized timing effect or viable-adult result.**

**Branch:** `theory/temporal-resource-competition-v01`, **draft PR #39**. The GEB `main` and submission PR #38 are unchanged.

## Decision first

There is **a small age-aligned difference in the temporal shape of recorded larval food intake**, while the first observed third-instar dates are largely synchronous in the mixed cages. It is **insufficient to claim that alternative host-plant supply times reverse the competitive disadvantage of the two species**. The original experiment contains no food-supply-time randomization, time-resolved available leaf biomass, adult emergence, or locally established introduced host PLANT.

More research by fitting models to the same final defoliation outcomes will not identify the hypothetical mechanism. A future direct timing manipulation is necessary. A pre-experiment biological pilot should assess whether the relative-age intake offset is repeatable across genuinely independent cohorts, source families and temperatures.

## Two independent source programmes and their exact identities

Original: Hashimoto & Ohgushi (2023), `Ecology and Evolution` 13:e10164, DOI https://doi.org/10.1002/ece3.10164. Original data: https://doi.org/10.6084/m9.figshare.23170898.

| Experimental source | Period | Original records | Unit of inference |
|---|---|---|---|
| Competition cages | Sep–Oct **2012** | 30 cages × 23 censuses = 690 records | **Cage**, not larvae or 4 plants per cage |
| Independently reared larval food demand | Jun–Jul **2016** | *A. alcinous* 98 rows; *S. montela* 105 rows | One 10-neonate founder group per species, then **10 A and 9 S** tracked individual larvae in source |

The 2016 cohort start dates are **2016-06-27** for A and **2016-06-11** for S (16 calendar days apart). Each species' ten-larva founder group is split to individually recorded third-instar cohorts six days after the original start. Naively overlaying raw June/July dates between these separate laboratory cohorts fabricates a biological temporal difference. They are compared **as days since each cohort's own hatch group started**, not on absolute calendar dates.

Exact original source file MD5s:

- 2012 `larvaldata.csv`, Figshare file `41146991`: `0e9ff924e4ae38cb78529ec540887c42`.
- 2016 `AFoodDemand.csv`, file `41147915`: `08c2bcfb9a896f7eb89bebd74153f769`.
- 2016 `SFoodDemand.csv`, file `41147924`: `174f0bd6ba056bcd49b348427da9f749`.

## 1. 2012 cages: direct stage occupancy, NOT intake or survival

Of 30 cages, **6 A-only, 6 S-only and 18 mixed** were assigned. **16/18** mixed cages had a later observation of both species in instars 3–5.

Of these 16 mixed cages:
- first observed 3–5 instar: **S earlier in 7**, **same census day in 8**, **A earlier in 1**;
- median paired first-observation lag `first_S − first_A = 0 days`.
- Among all 18 mixed cages × 23 visits = **414 cage-visit snapshots**, the 3–5 instars were recorded together at **171**, A only at 73, S only at 17, neither at 153.

The snapshot-overlap frequency (171/414) is **not** actual feeding or within-plant contact, and the 16-cage lag analysis conditions on *both species surviving until late instars*. The census grid resolves appearances at approximately 2-day intervals, not actual molts.

**Source audit:** https://github.com/zuizui0223/chocho/actions/runs/38502905650 and subsequent https://github.com/zuizui0223/chocho/actions/runs/38504162822.

## 2. 2016 intake: published method and numeric validation

The paper, Methods §2.2, states that fresh *Aristolochia debilis* leaves were replenished **every two days**, weighed before supply, and weighed when removed. Authors fitted the leaf fresh-to-dry conversion:

`ln(dry mass in g) = −1.48 + 0.90 × ln(fresh mass in g)`.

For a source row pair, the code independently converts *previous* `given` and *current* `taken` and then subtracts dry masses. Correct explicit units are:

`dry_leaf_mg = 1000 × exp(−1.48 + 0.90 × ln(fresh_leaf_g))`.

Using mg *inside* the log instead of converting **after** calibration gives a factor-of-1.995 error; this was found and corrected using the published Table 1 values before interpreting the chronology.

Only genuine adjacent 2-day single-larva pairs are included. The parser retains, reports and excludes missing/unpairable, physically negative, and `memo=redundant` intervals rather than silently imputing them. A same-day `S` founder-group stage 2/3 transition is retained as provenance; it is excluded from the separate single-larva late-intake calculation.

| Outcome | *A. alcinous* | *S. montela* |
|---|---:|---:|
| Individuals with source pupation label | 10 | 9 |
| Valid late-instar 2-day intake intervals | 70 | 76 |
| Published Table 1 instars 3–5 dry consumption per larva | 1,003.23 mg | 394.25 mg |
| Independently reconstructed with conservative pairing | **1,011.00 mg** | **395.86 mg** |
| Discrepancy from Table 1 | +0.77% | +0.41% |
| Median day since founder-group start to 50% of own recorded late-instar intake | **day 18** | **day 20** |
| Highest observed mean per tracked 2-day interval | **day 18**, 353.98 mg (n=10) | **day 22**, 105.21 mg (n=9) |

**Interpretation:** There is a small **descriptive 2-day age-aligned lag** in the median half-consumption time, and peak per-interval mean occurs four days later in S than A. This was computed from source records the original publication already collected; it is **not independent biological replication or a new experiment**. There is one source cohort per species, with 2-day interval measurement resolution. The apparently shifted peak should be interpreted in the context of late-stage survivor/observation censoring (e.g. A day 22 has only three recorded intervals), unequal overall dry intake, the stage transition and different founding calendar dates.

There are also **three physically negative matched wet-to-dry intake intervals in each species**, which are reported and omitted as invalid pairs rather than truncated to zero. The ~1% Table 1 agreement is a **computation verification gate** under author-source assumptions, not evidence of a causal or population-generalizable effect.

Source CI: https://github.com/zuizui0223/chocho/actions/runs/38504162822. The exact original-file summary, per-individual data and day-by-day masses are retained as downloadable JSON in that workflow artifact, `kyoto-two-species-calendar-source-v01`.

## 3. Bridge to the temporal competition theory, and a stringent no-go condition

The toy model in `TEMPORAL_RESOURCE_COMPETITION_SIGN_REVERSAL_V01.md` predicts which consumer receives a larger **relative competition penalty** when precisely the **same total edible supply** is scheduled early versus late. The structural mathematical result is conditional on stage-aligned demand and shortage-to-emergence response functions. The 2016 chronology suggests the **possibility** of different relative intake shapes even though both species entered third instar at similar developmental ages.

However, **do not retroactively choose whichever early/late split gives a modeled reversal**, or label the published intake curves as experiments manipulating supply. The empirical difference can be due to source family, developmental timing, censoring or tissue handling, and neither source provides an independent adult-emergence response as a function of randomized foliage supply. In particular, absolute-date differences across nonoverlapping laboratory cohorts are irrelevant to natural temporal resource competition.

Before choosing a decisive window, run a **new blinded feasibility pilot** with synchronized eggs from independent mothers and site/population blocks in which both species have daily or 2-day true per-capita edible dry intake and stage; measure a shared *A. debilis* food replenishment and contemporaneous plant quality. Compare *age-since-hatch* as well as actual seasonal calendar dates in ecologically relevant same-location cohort starts. Freeze early/late windows using **the pilot data only**, reserving separate independently randomized cages for the primary test.

The future causal experiment then crosses predeclared food timing arms (early/continuous/late at equal total actually edible leaf mass) and presence of a competitor, ideally with a same-added-larva-count conspecific control and separate consumption-matched controls if the claim is specifically *heterospecific identity* rather than total crowding. The primary outcome is **flight-capable adults per originally assigned cohort**, with live diapausing pupae separated from confirmed deaths. Cage is the randomized unit; replicate source populations/maternal families. A sign change in species-specific *relative* competitor penalties requires two separate supported contrasts, not merely noisy opposite point estimates.

**STOP / NO-GO:** If independent pilots fail to confirm an ecologically meaningful relative demand shift, **do not market Kyoto as a test of time-induced competitive sign reversal**. Prefer the already-defined quality/hidden shortage/diapause mechanism test, which is closer to the unresolved nondepletion suppression in Hashimoto & Ohgushi (2023). This is competition involving an introduced *butterfly*, on the **native** host *A. debilis*; it cannot itself validate introduced-host-*plant* globalization claims.

## Reproduce

Original MD5-pinned data are fetched from Figshare by the branch-only Actions workflow `.github/workflows/kyoto-two-species-calendar-source.yml`.

```bash
pytest -q tests/test_kyoto_calendar_demand_identifiability.py tests/test_kyoto_late_instar_intake.py
python scripts/audit_kyoto_calendar_demand_identifiability.py --receipt results/stage_original.json
python scripts/reconstruct_kyoto_late_instar_intake.py --receipt results/intake_original.json
```

**No new randomized animal experiment has been done; GEB manuscript remains separate.**
