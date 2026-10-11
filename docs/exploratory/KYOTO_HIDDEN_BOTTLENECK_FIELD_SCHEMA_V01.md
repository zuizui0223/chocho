# Field-ready recording contract: hidden transient resource limitation in the Kyoto butterfly system

**2026-10-10 | PRE-EXPERIMENT PROCEDURAL DESIGN, not collected data or a new result.**

**Primary replication clarification:** The input field `plant_id` is a legacy variable name. It MUST identify **one independently randomized cage/enclosure containing all its potted host plants**. If four plants are in a cage as in the 2023 experiment, all four belong to ONE randomized unit; never enter them as four `plant_id` rows. Use a corresponding pot-level plant inventory for within-cage measures if appropriate, keeping `plant_id` as the parent cage ID.

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
plant_id,observation_datetime,native_larvae_alive_n,native_late_instar_alive_n,accessible_leaf_area_cm2,fresh_leaf_area_added_cm2,present_leaf_age_class,added_leaf_age_class
```

**Source-grounded stage targeting:** Hashimoto & Ohgushi (2023), Table 1, measured dry-leaf consumption of **21.66 mg (instars 1–2) versus 1,003.23 mg (instars 3–5)** for `Atrophaneura`, and **1.88 versus 394.25 mg** for `Sericinus`. Those are cumulative measured consumption amounts, **not** measured standing leaf supply minima. A food clamp must be checked during 3rd–5th instars (and earlier), with the `food_floor_cm2` calibrated by a separate blinded edible leaf area/instar demand pilot. One constant floor selected at random or inferred from final leaf classes is invalid.

The analysis-ready longitudinal export now includes `native_larvae_alive_n` and `native_late_instar_alive_n` (count of **living native 3rd–5th instars**, not all assigned larvae or survivors at the end). A skipped observation is not a zero. When all native caterpillars died before late instars, later "no food shortage detected" is uninformative about earlier shortage. `present_leaf_age_class` must be `young`, `mature`, `mixed`, `none` or `unknown`; `added_leaf_age_class` must be supplied whenever fresh leaf area was added. Added young tissue and resident mature tissue can have different quality even if accessible cm² is equal. All of these are **postrandomization diagnostics**, never criteria for excluding cages from intention-to-treat comparisons.

**Preferred quality-matched sham, where biologically feasible:** Include the same quantity, source and leaf-age distribution of extra *Aristolochia* tissue in all four arms, but place it behind a food-inaccessible, airflow/space-matched perforated barrier in the natural-access arms and allow larval access in the clamp arms. Then the randomized distinction is physical access to tissue rather than simply the delivery of young foliage. Validate the barrier's effects on microclimate and movement, and record inaccessible sham foliage in the full observation ledger. `fresh_leaf_area_added_cm2` in the concise analysis ledger denotes **accessible** supplemental leaf area; inaccessible sham material is tracked separately in the full field ledger. This is a design control, not a claim that leaf age or secondary chemistry has been fully neutralized.

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

**Audit reporting:** The script outputs arm-specific living-late-instar visit counts, late-instar food-floor breaches, missing late-stage area measurements, and added foliage leaf-age classes. Comparisons based only on cages where third instars survived would condition on post-treatment survival and can manufacture apparent support for the food-shortage hypothesis. Interpret the clamp's nutritional composition and sham treatment separately. Stage-specific measures are missing in the 2023 source and cannot be reconstructed from its four final plant ratings.

**Do not calculate interaction or confidence intervals from the 2023 experiment by substituting its four terminal plant scores; it did not contain these randomized treatment arms.**

## H. Distinguish HOST PLANT ontogeny from age of added leaves (independently published 2023 control)

[Jeong et al. (2023), *Frontiers in Plant Science*, DOI 10.3389/fpls.2023.1145363](https://doi.org/10.3389/fpls.2023.1145363) tested *Sericinus montela* on ***A. contorta*** leaves from **1-, 2-, and 3-year-old whole plants**, deliberately standardizing sampled leaves to the **5th–7th leaves below the shoot tip**. Feeding trials started at instar 2, ran six days at 24°C, and replaced leaves as soon as over half the leaf had been consumed. Nonetheless, published butterfly relative growth depended on plant age × season (F=8.952); host leaf C/N had a large July-versus-September effect (F=270.727). Measured July leaf AAI means (ng/mg) were **83.80** for first-year plants, **7.61** for second-year, **2.26** for third-year plants. These are **published effects**, verified against original article XML, **not** new chocho causal analyses or necessarily transferable to Japanese `A. debilis`.

**Key trap:** The currently recorded `present_leaf_age_class`/`added_leaf_age_class` describe *the condition of the individual leaf*, not *the age of the whole donor plant*. They cannot control this known ecological axis. Amend the **real pre-experiment donor/source ledger** with these columns while preserving original randomized cage IDs:

```csv
plant_id,donor_plant_unique_id,host_species,donor_source_population,whole_donor_plant_age_years,plant_age_known,leaf_shoot_node_number,leaf_development_class,harvest_datetime,food_access_state,foliage_fresh_mass_g,edible_leaf_area_cm2,leaf_N_mass_percent,leaf_C_to_N,AAI_ng_per_mg,AAII_ng_per_mg,chemical_assay_sample_id,quality_measure_missing_reason
```

One row per **donor plant × harvested leaf batch**, nested within each **whole-enclosure `plant_id`** (which is a historical misnomer for whole cage ID). The `food_access_state` distinguishes leaf made available for eating from fresh but caged-off/inaccessible sham foliage; never treat the latter as added edible supply. All four randomized arms need *balanced sources, same whole-plant maturity where identifiable, comparable leaf nodes and harvest time*. For wild established `A. debilis` whose true plant age cannot be known, use `plant_age_known=false`, do not impute years; randomize/block on matched known accession/provenance and observed growth/phenophase instead. Chemical assays and N/AA content are optional if unavailable, but the missingness must be explicitly reported. An unknown donor age cannot be silently reinterpreted as a one-year plant.

Even if matched `A. debilis` foliage is provided throughout the high-demand 3rd–5th instars, the 2×2 contrast still identifies an intervention of **access to matched donor tissue** rather than a pure mathematical effect of edible biomass grams alone: barriers, foraging area, residual plant state, light and unmeasured defense traits are alternative pathways.

**Data gate:** Never use Korean `A. contorta` numeric AAI concentrations as values for *A. debilis*. Nor treat 20 larvae per cell or six-day specialist growth as independent wild-population fitness evidence. This study did not simultaneously test *A. alcinous* against *S. montela*, and neither a joint competition mechanism nor flight-capable adult recruitment is identified by it.

## I. Host quality can change diapause, not just survival — confirmed 2008 prior art

**Source:** Takagi & Miyashita (2008), *Host Plant Quality Influences Diapause Induction of Byasa alcinous*, *Annals of the Entomological Society of America* 101:392–396, DOI [10.1603/0013-8746(2008)101[392:HPQIDI]2.0.CO;2](https://academic.oup.com/aesa/article/101/2/392/8474). NCBI Taxonomy 85299 lists *Atrophaneura alcinous* as a homotypic combination of *Byasa alcinous*; this is the **same butterfly taxon**, not a study of a different species. The experiment used a **Chiba mountain population** and ***Aristolochia kaempferi*** (NOT the Kyoto native *A. debilis*) at 20 °C and 16:8 light:dark. Soft-versus-tough host leaf treatments changed early survival, larval duration and pupal **diapause** under controlled conditions; nearly all pupae from the tough-leaf treatment entered diapause. Human mowing altered young-versus-old host leaf abundance at the monitored site. These outcomes were **published in 2008** and cannot be claimed as a newly discovered butterfly mechanism or transferred automatically to the Kyoto source population.

**Biological consequence for the current trial:** a pupa still alive but in diapause at the end of a four-week experiment cannot be called `pupa_no_adult` and treated as dead. The same treatment can decrease *short-term adult emergence* while leaving eventual postdiapause adult emergence unresolved; early adult abundance, diapause occurrence, and adult production after winter are distinct ecological endpoints. Moreover, relative leaf toughness, tissue age and larval developmental delay must be measured independently before attributing changed pupal fate to a unique chemical signal.

**Revised accepted `native_fates.csv` codes** (each originally assigned neonate retains one terminal/current state):

| `fate` | Biological interpretation in analysis |
| --- | --- |
| `flight_capable_adult` | Confirmed successfully emerged flying adult |
| `adult_nonflight` | Emerged but confirmed unable to fly under predefined criteria |
| `dead_immature` | Confirmed mortality before/at pupation |
| `pupa_confirmed_dead` | Documented pupal death, NOT merely no adult by an arbitrary census date |
| **`diapause_alive`** | **Pupa alive and in diapause**, adult emergence pending; keep in original cohort |
| **`pupa_alive_pending_emergence`** | Living pupa, diapause status not yet resolved |
| `lost_to_followup` | Actual unobserved condition |
| `other_documented` | Other explicit condition requiring review |

**Previous `pupa_no_adult` is now REJECTED** as ambiguous. The code distinguishes confirmed deaths from still-living pupae, and counts the latter as **pending eventual adult outcomes**, contributing to uncertainty bounds but never to confirmed mortality.

Before randomization, register one **common within-block pupal/adult observation cutoff**, a predefined way to recognize diapause with uncertainty noted, and what happens to pupae requiring overwintering. Use a supplementary **per-assigned-larva follow-up ledger** in addition to the simple fate export:

```csv
plant_id,larva_id,pupation_datetime,pupa_alive_at_census,diapause_diagnosis,diagnosis_method,diagnosis_datetime,followup_census_datetime,emergence_datetime,adult_flight_test_date,adult_flight_capable,death_confirmation_datetime,last_seen_datetime,followup_unknown_reason
```

Record `diapause_diagnosis` as `diapause`, `non_diapause`, or `unresolved`, not based solely on failure to emerge by a date. Dates must include time zone. A second trial after overwintering, if permitted and logistically feasible, can determine actual emergence for pupae initially scored `diapause_alive`; analysis must retain original randomized treatment, plant/cage ID, source family and denominator.

The main 2×2 interaction at a **fixed early census date** is only a treatment effect on *adult appearance by that date*. A valid interpretation as *eventual adult reproductive success* requires sufficiently complete postdiapause follow-up. Otherwise only **lower/upper possible eventual adult success fractions**, plus separately reported percentages in diapause, are warranted. Even when the pupal fate codes are correct, a diapause shift is not automatically an adaptive evolutionary response without a fitness advantage in the relevant seasonal environment.

**Source and code boundaries:** the 2008 PDF/article did not provide original individual-level rows in this continuation, and the current code tests use synthetic data integrity fixtures only. They provide no evidence that *S. montela* induces diapause in the Kyoto butterfly population. The GEB global host-plant manuscript remains separate.

## J. Independent 2020 density×food study: record individual life-history choices and direct causes

**Published prior art, not our discovery:** Nakahara, Horita, Booton & Yamaguchi (2020), DOI [10.1111/ens.12399](https://doi.org/10.1111/ens.12399), experimentally altered food supply and **conspecific** native larval density, finding an interaction between extra moults, cannibalism, disease, and diapause. Kozuki & Takeda (2004), DOI [10.11257/jjeez.15.157](https://doi.org/10.11257/jjeez.15.157), already linked starvation, density, local host-use populations and split life cycles. Therefore the novel question is NOT that starvation or overcrowding can affect diapause; it is whether **adding heterospecific larvae instead of the same number of conspecific larvae** shifts these alternatives while source host quality and food access are controlled.

**The old 2×2 assignment tests the total effect of adding *Sericinus* relative to adding no one. It CANNOT uniquely identify heterospecific identity** because adding *Sericinus* also changes total larval density. A separate whole-cage randomized 2(resource support)×3(additions: none/extra native *Atrophaneura*/*Sericinus*) identity design would permit an equal-headcount interspecific-minus-intraspecific contrast while the **initial tagged focal native recipients** stay constant; two focal native density levels would be needed to test its moderation by native density (12 combinations). Equal headcounts do NOT equal added feeding demand: the 2023 study's consumption ratio was ≈2.6× higher for native *A. alcinous* per larva, so the food-use pathway is still a separate biological mechanism.

For **any** prospectively observed 2026+ new trial, the optional **individual-level life-history event ledger** should use the following independent event-type codes and must not infer a direct event simply because a larva disappeared:

```csv
plant_id,larva_id,observation_datetime,life_history_event,evidence_type,observed_molt_number_or_instar,cannibal_victim_id,cannibal_aggressor_species,disease_confirmed_method,observer_id,observation_complete_at_this_instar
```

Acceptable `life_history_event` includes `ecdysis_observed`, `sixth_instar_observed`, `cannibalism_directly_observed`, `disease_confirmed`, `pupation_observed`, `pupal_diapause_diagnosed`, `adult_emergence_observed`, or `other_observation`. `cannibalism_directly_observed` requires direct observation or video with recorded victim and observed attacker (if identity remains unknown, mark it unknown rather than infer from the census). `disease_confirmed` requires a recorded confirmation method; a vanished larva cannot be relabelled as disease. `sixth_instar_observed` must be based on individual-stage follow-up or verified moult evidence; another day with a fifth-instar caterpillar in the **cage-level count** is not proof of an extra ecdysis. When the event-monitoring completeness is unknown, absence of a recorded event is **not a biological zero**.

The [original 2023 cage CSV](https://doi.org/10.6084/m9.figshare.23170898) contains 690 rows from 30 independently randomized whole cages. The rows include **counts** of first through fifth instars, prepupae/pupae, cumulative pupation and two-species source densities, but **no individual larval IDs, sixth instar, individual moults, cannibalism observer identity, disease confirmation, diapause diagnosis or flight-capable adults**. See [original source-specific audit](https://github.com/zuizui0223/chocho/actions/workflows/kyoto-original-stage-demand-audit.yml). It is not legitimate to estimate any of those unmeasured individual-level processes from the 690 group-level entries.

**Data/interpretation stop:** The 2020 experiment's original individual-level raw files were not identified in indexed public source searches here. That does not prove they do not exist; no outside data were fabricated. No new butterfly life-history outcomes or dose-response effects were generated. See [2020/2004/2023 mechanism cross-study novelty decision](KYOTO_HETEROSPECIFIC_DENSITY_AND_LIFE_HISTORY_IDENTIFICATION_20261011.md).
