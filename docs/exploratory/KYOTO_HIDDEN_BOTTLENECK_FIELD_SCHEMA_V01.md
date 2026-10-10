# Field-ready recording contract: hidden transient resource limitation in the Kyoto butterfly system

**2026-10-10 | PRE-EXPERIMENT PROCEDURAL DESIGN, not collected data or a new result.**

Scientific design is frozen in `KYOTO_HIDDEN_RESOURCE_BOTTLENECK_FACTORIAL_V01.json`. The working question is whether *Sericinus montela* reduces native *Atrophaneura alcinous* output through **transient edible-host shortage**, which a FINAL ordinal defoliation rating cannot detect.

## A. Randomization (only after local plant suitability and approvals are checked)

Prepare a UTF-8 CSV **with exactly these headings**, one row per *independent plant/enclosure* to randomize:

```csv
plant_id,block_id
```

No dummy plants or survival observations have been entered. Choose `block_id` from a predeclared combination of trial date, plant source/genotype batch and butterfly parental-source stratum. A complete block must contain **at least four real eligible plant/enclosures, and a multiple of four**; this is a balanced-allocation requirement, **not a sample-size or power recommendation**. Determine independent plant/enclosure sample size after a pilot estimating adult-eclosion variance and plant supply feasibility.

Use the validated script:

```bash
python scripts/randomize_kyoto_bottleneck_factorial.py \
  --inventory eligible_plants.csv \
  --allocation randomized_plants.csv \
  --receipt allocation_receipt.json \
  --seed <chosen_integer_recorded_before_randomization>
```

Each block contains equal numbers in four arms: `NATURAL_NO_COMPETITOR`, `NATURAL_COMPETITOR`, `CLAMP_NO_COMPETITOR`, `CLAMP_COMPETITOR`. Record seed, plant IDs, family/block and both factors before exposing plants. Neither the script nor its unit tests produces or predicts butterfly outcomes.

## B. Trial-assignment ledger (after randomization)

One row per randomized cage/plant, keyed by `plant_id`:

`plant_id, block_id, treatment, competitor_present, food_clamp, trial_date, native_family_id, competitor_family_id, plant_source_id, native_initial_n, competitor_initial_n, trial_length_days, resource_target_method, target_accessible_food_units, handling_protocol_id, observer_id`.

- Keep the native initial cohort size constant across treatments; record *every initial neonate*.
- Competitor density for the positive arms must be chosen in a blinded feasibility pilot, **not optimized from survival results**.
- Resource target and measurement units (`cm2` accessible edible leaf area, edible fresh mass, or calibrated combination) must be fixed from plant age/instar demand before unblinding.
- Extra food added in clamp arms is an experimental resource intervention; the sham visits must reproduce handling in natural arms without adding food.

## C. Repeated resource access, biology and handling

A separate longitudinal ledger is required, **one row per cage × prospectively fixed observation occasion**:

`plant_id, observation_datetime, elapsed_hours, native_alive_n, competitor_alive_n, native_instar_stage, competitor_instar_stage, accessible_leaf_area_cm2, fresh_leaf_area_added_cm2, edible_fresh_mass_g, fresh_mass_added_g, leaf_age_class, plant_water_status, temperature_C, relative_humidity_pct, handling_sham_done, visible_contact_count, observer_id, missingness_reason`.

Distinguish all measures:
- `accessible_leaf_area_cm2`: simultaneously consumable live leaf area, not total leaves or visually scored percent remaining;
- `fresh_leaf_area_added_cm2`: intervention amount, **not** the amount remaining after later feeding;
- `native_instar_stage`: life-stage demand varying over time, not a single calendar date proxy;
- `visible_contact_count`: direct observational candidate for interference, **not** proof of interference;
- missing observations remain missing; **do not infer accessible food = 0** from a skipped visit.

Daily/resource observation schedule and target threshold must be set after a separate pilot (without fitting the focal survival contrasts) and followed for all four arms. Keep plant source and tissue quality/age matched where possible. Compensation may alter leaf nutrients, moisture and handling; any putative "quantity" causal inference must check these sources of bias.

## D. Prospective individual outcome ledger

Key on the uniquely identified originally assigned larva, **not** only surviving larvae:

`plant_id, larva_id, species, assignment_date, original_cohort, fate, last_seen_datetime, pupation_datetime, adult_emergence_datetime, flight_capable, parasitoid_emergence, notes`.

Allowed primary `fate`: `flight_capable_adult`, `adult_nonflight`, `dead_immature`, `pupa_no_adult`, `lost_to_followup`, `other_documented`. Define how censoring will be treated *before* fitting the primary contrast; do **not** treat missing fate as success or complete-case-drop it without sensitivity. Every assigned focal larva appears once.

**Primary endpoint** at plant/enclosure level: number of viable flight-capable focal adults / *all originally randomized focal larvae*. Report raw numerators, randomized cages and attrition. Do not count each larva or final plant leaf as independently randomized units.

## E. Trial quality gates

1. Legal rearing/transport/containment of introduced *S. montela*, authorized use of host plants and no release in natural habitats.
2. Distinct cage/plant IDs and no treatment swaps after allocation; block balance and dates recorded.
3. Repeated accessible-leaf measurements cover the focal larvae's instars, including high-demand periods; terminal remaining-food ordinal scores are insufficient.
4. Food clamp **actually changes the transient resource constraint** in competitor-present arms. If not, this is an inconclusive manipulation rather than a rejection of hidden food limitation.
5. Same handling across conditions, validated plant-stage comparability, no differential loss caused by the supplement procedure.
6. Clear individual fates and original-cohort denominators. If adult eclosion not reliably observable, describe a pupation study, **not** full recruitment.
7. If the resource clamp eliminates the negative competitor effect, conclude a supply-access-moderated competition effect, **not** a unique proof of mass depletion; fresh tissue chemistry and handling remain alternatives.
8. If suppression persists in clamped food, follow the distinct randomized prior-feeding/standardized wounding design in `KYOTO_NONRESOURCE_COMPETITION_MECHANISM_V01.json`; no automatic claim about induced plant defenses.

## F. Previous published source vs new experiment

The MD5-exact 2023 Kyoto paper source has 30 cages × 23 larval observational dates and 120 **end-of-trial ordinal defoliation scores** for four different plants per cage; it has **zero repeated measurements of edible mass/area or randomized food replenishment**. The 2012 ESJ poster and 2017 plant-regrowth experiment establish temporal availability and plant response as prior hypotheses. The new experiment is an independent causal intervention, not a reanalysis of 2023 cage ratings.

**All records here are schemas only. No live animals have been allocated, observed or handled through this repository workflow. GEB PR #38 remains untouched.**

## G. Validated four-file outcome analysis contract (2026-10-10 extension)

The original longer-form ledgers above are a field collection template. For a **future real experiment**, the analysis tool expects four small, locked and consistently keyed exports in this exact order. Export *without editing treatment assignments after outcome inspection*:

```text
randomized_plants.csv
plant_id,block_id,treatment,competitor_present,food_clamp,random_seed,design_version

trial_assignments.csv
plant_id,block_id,treatment,competitor_present,food_clamp,native_initial_n,competitor_initial_n,food_floor_cm2,expected_n_resource_visits

native_fates.csv
plant_id,larva_id,species,fate

resource_visits.csv
plant_id,observation_datetime,accessible_leaf_area_cm2,fresh_leaf_area_added_cm2
```

**Before animals enter cages:** `native_initial_n` is the originally assigned focal cohort size; `food_floor_cm2` is a **positive, prospectively selected accessible-leaf-area target**, fixed within each randomization block; `expected_n_resource_visits` is a preset sampling schedule, the same for all four treatments in a block. The condition `food_clamp` is the *randomized strategy*, not a post-hoc indicator that a measured threshold happened to be achieved.

For `resource_visits.csv`, a date must be ISO 8601 with an explicit timezone offset (example structure `YYYY-MM-DDTHH:MM:SS+09:00`). Missing leaf-area observations remain **empty cells**, never replaced with 0. A missing visit is a protocol violation to audit, not permission to delete the affected cage. A biologically present leaf area of zero is numerically `0`, not blank. Where supplemented foliage is supplied, record its measured area separately from the contemporaneously accessible area.

For `native_fates.csv`, permitted `species` is `Atrophaneura alcinous`. Every originally assigned focal larva requires its own unique `larva_id` and one fate: `flight_capable_adult`, `adult_nonflight`, `dead_immature`, `pupa_no_adult`, `lost_to_followup`, or `other_documented`. Unknown outcomes remain in the randomized denominator; the tool produces **true worst-case interaction bounds** accounting for positive and negative arms of the difference-in-differences. A complete-fate dataset is preferable. Any missing cohort member is an input failure rather than a complete-case exclusion.

After valid original data exist, run:

```bash
python scripts/analyze_kyoto_food_clamp_outcomes.py \
  --allocation randomized_plants.csv \
  --assignments trial_assignments.csv \
  --fates native_fates.csv \
  --resource resource_visits.csv \
  --receipt results/kyoto_food_clamp_trial_v01.json
```

This writes per-arm cage counts, initial cohort denominators, adult success, attrition bounds, prospectively logged leaf access and intervention, and the **natural** and **clamped** competitor effects *separately*, along with their interaction. Block bootstrap intervals are withheld if there are fewer than three distinct independent source blocks; three is a computational minimum for output, **not** a statistically adequate sample-size recommendation. Report independent cage sample size and feasibility pilot precision separately. The model does not automatically declare that the resource clamp succeeded: a clamp's effect on fresh tissue quality, humidity or labor remains a potential alternative mechanism.

Code quality check: [GitHub Actions outcome-contract tests](https://github.com/zuizui0223/chocho/actions/workflows/kyoto-food-clamp-outcome-contract.yml) use **synthetic ledger fixtures only**, strictly for input, randomization and attrition integrity. They are not field observations, simulations of ecological effects or evidence of hypothesis confirmation.

**Do not calculate interaction or confidence intervals from the 2023 experiment by substituting its four terminal plant scores; it did not contain these randomized treatment arms.**
