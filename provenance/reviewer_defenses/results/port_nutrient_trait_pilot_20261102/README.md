# Nutrient-affinity traits of candidate alien butterfly hosts at ports — exploratory pilot

**2026-11-02 · separate branch · NOT a new confirmed port mechanism or a GEB/JBI submission result**

## November 2 status update: island-matched coast data now acquired

The [strict 78-port / 35-island 5km observational source audit](STRICT_WPI78_SOURCE_AUDIT_20261102.md) is completed and reproducible. It supersedes the earlier missing-nonport-plant-data status below **for the newly reselected sample only**, not for the original nonport coordinates. The old 79-port 5km comparison was invalid for independent contrasts because 48/79 buffers overlapped; an outcome-blind 78-site geography revision removed that issue. Strict pair-local means of host Ellenberg N differences, island-weighted: early **+0.1129** (9 eligible islands; bootstrap 95% [−0.0340,+0.3020]); later **−0.0123** (22 eligible islands; [−0.3166,+0.2811]). Same 30 WPI pairs in 9 islands: late−early **−0.1622**, bootstrap [−0.5912,+0.1721]. No consistent port-specific nutrient-affinity difference is established. The positive **global species-level rho +0.341** is preserved. **No sewage or soil nitrogen exposure was measured.**

The research plan that follows remains a statement of *what additional evidence would be necessary* for any causal or JBI port gateway claim; it is not a claim that botanical geography is still completely unacquired.

## Global prevalence control — additional positive result

A [source-locked global GBIF prevalence audit](GLOBAL_GBIF_PREVALENCE_CONTROL_20261102.md) has since verified that the **raw port frequency × N correlation (rho=+0.338; n=182 accepted taxa) remains positive when country-wide geographic recording and total GBIF occurrence indexing are controlled (partial rank rho=+0.210, genus-bootstrap descriptive interval +0.065 to +0.348)**. This is stronger evidence of a *trait-specific port-recording association* than the unadjusted global correlation, but not evidence for local sewage nitrogen, plant imports, butterfly reproduction, or port-specific enrichment; the strict same-island comparisons still remain inconclusive.

## Scope and research question

The main chocho paper shows +54.9% additional potential butterfly × WGSRPD3
regional resource opportunities when contemporary introduced-host geography is
retained. The time-out-of-sample 155-butterfly evaluation has MSE improvement
of 4.12% when host information is added. Island-port analyses also retain the
positive 51-island host–butterfly assemblage correspondence (Spearman ρ=0.300).
**All of these results are preserved and unchanged.**

This small new analysis asks a different, trait-based question:

> Among candidate introduced larval host plants already indexed near ports,
> are species preferring nutrient-rich environments recorded at more ports?

It tests the proposed **nutrient-affinity mechanism** as a trait association.
It does not measure N inputs, the presence of wastewater, soil nutrient
concentrations, host quality, butterfly larval survival, port arrival, or
realized butterfly colonization.

## Exact, auditable source inputs

- Source 1: chocho pre-existing **candidate non-native host plants indexed at
  ports**, 687 host-species rows; source file:
  [`port_alien_host_broad_species.csv`](../port_alien_broad_inventory_20261030/port_alien_host_broad_species.csv),
  Git blob `9712b8e22a7e74f7e8809e7212ed6a8cfa2158fb`.
  Underlying audit: 1,596 ports with complete named-host queries; 1,150
  had ≥1 candidate non-native host indexed. These are 2001–2024 pooled
  GBIF catalogues within 5 km of port locations; a zero observation
  is not plant absence.
- Source 2: Germany-centered Ellenberg vascular-plant indicator table
  `Ellenberg_VascularPlants.csv` of `bathyscapher/EVA`, immutable repo
  commit `73b323983702a74cb9b6867bcebe00a307ea6566` and exact
  Git blob `23c27041588d62dd3c8dcff056858834cde56d0d`.
  Original indicators were developed chiefly for German/central-European
  vegetation. Importantly, Ellenberg N is better interpreted as **site
  nutrient status / fertility**, not a direct measurement of N supply or
  foliar tissue nitrogen. See Tichý et al. 2023,
  DOI `10.1111/jvs.13168`.

Join exactly by the complete species name (case-insensitive). No synonym
expansion, no genus imputation, no taxonomic imputation, no conversion
of unscored `x` into zero. Do not redistribute the external indicator table
in this repository.

## New results and scope restrictions

| Diagnostic | Observed value |
|---|---:|
| Candidate introduced host species in port catalogue | 687 |
| With exact-match nutrient values | **183 (26.6%)** |
| Spearman ρ: nutrient score N × number of indexed ports | **+0.341** |
| Nutrient-affinity group N=7–9 | 62 spp.; median **50** indexed ports |
| Low-affinity group N=1–3 | 41 spp.; median **9** indexed ports |
| 90th-percentile cap on port counts (≤129 ports) | 164 spp.; ρ≈**+0.380** |
| Delete-one-genus sensitivity | ρ **+0.327 to +0.372** |
| 127 recorded genera | 27 with ≥2 matched spp. |
| Within-genus comparison (genus has ≥2 spp. and differing N) | **24 genera; 77 spp.** |
| Within-genus fixed-effect slope, log(1+ports) per N point | **+0.105** |
| 4,999 within-genus N permutations, two-sided | **p≈0.38**, NOT supported |
| Median ports: N-matched / N-unmatched species | **28 / 12** |
| Same-species subset: nutrient N versus light L | ρ=+0.341 vs **+0.025**; Δρ≈+0.316 (182 spp.) |
| Same-species subset: nutrient N versus moisture F | ρ=+0.395 vs **+0.043**; Δρ≈+0.352 (155 spp.) |
| Same-species subset: nutrient N versus salinity S | ρ=+0.344 vs **−0.023**; Δρ≈+0.366 (182 spp.) |

The overall positive rank correlation is robust to removal of very frequent
species and to leaving out any one genus, but the small eligible within-genus
sample does not independently resolve a nutrient-affinity effect. The Germany-
centric matching is selective: matched species have more indexed ports even
before any trait model. A genus-cluster resampling interval for the overall
correlation was positive in exploratory calculations, but it does not correct
selection, geographic range, observation effort, or causal confounding.

The positive nutrient pattern is stronger than the similarly calculated
light, moisture, or salinity trait–port associations on their respective
identical-species subsets. This makes the *trait-specific descriptive pattern*
more interesting, but does not substitute for measured nitrogen exposure,
plant range size, matched non-port coasts, or independently verified host use.

**Conclusion:** a positive association in **port-indexed plant catalogues**,
not a demonstrated nutrient-enrichment mechanism. The existing positive GEB,
island and host-forecast results are not contradicted.

## Reproduce

```bash
curl -fLsS \
  'https://raw.githubusercontent.com/bathyscapher/EVA/73b323983702a74cb9b6867bcebe00a307ea6566/EllenbergsValuesAutomaton/Ellenberg_VascularPlants.csv' \
  -o /tmp/ellenberg_germany_2020.tsv

python scripts/exploratory_port_nutrient_ellenberg.py \
  --ellenberg /tmp/ellenberg_germany_2020.tsv \
  --bootstrap 2000 \
  --permutations 4999 \
  --output /tmp/port_nutrient_trait_pilot.json
```

Read and run `tests/test_exploratory_port_nutrient_ellenberg.py` to validate
strict name matching, site-indicator parsing, tied ranks, and genus-restricted
label permutations without any network request. This branch neither modifies
nor supersedes the original GEB frozen manuscript.

## Competing hypotheses and next decisive observational design

The current positive pattern might reflect (a) nutrient-affinity traits,
(b) globally common ruderal plants entering many ports, (c) botanical and
geographic relatedness, (d) stronger GBIF observation effort near accessible
ports, or (e) actual nitrogen-rich habitat. Without non-port sampling, these
cannot be separated.

**The next paper-critical dataset is still the existing 36 islands / 79
same-WPI-ID ports:** identify same-radius, same-period port and non-port
coastal plant catalogues with sampling effort, native/introduced status,
and general plant-community denominator. Compare within-island (i) ports,
(ii) urbanized non-port coasts, and (iii) natural non-port coasts. Estimate
both nutrient-affinity composition among candidate host plants and the same
composition among non-host plants as an observation/background control.
A strong positive trait gradient in all plants could be unrelated to
butterfly-specific host selection.

Independently layer **nitrogen deposition**, vegetation/cropland, and
wastewater loading only after plant sampling alignment. HydroWASTE is a
global set of 58,502 WWTPs with reported/estimated treatment, effluent,
and estimated watercourse outfall, not a direct measurement of local soil N:
the original publication reports 6.5 ± 3.1 km average distance between
plant site and estimated outfall. Small islands and coastal discharge
are particularly uncertain. Nearest treatment plant or raw count within
1 km must **not** be interpreted as nitrogen supplied to port plants.

For post-2011 butterfly records, retain prior occurrence, climate,
geographic connectivity, species–host availability, background and
recording effort. Do not use 2019 port records as an observed historical
shipping network or classify 2001–2024 pooled hosts as a pre-2011 cause.

**Publishing gate:** keep this in exploratory provenance until replicated
with geographically matched port–non-port botanical samples and prior
butterfly records. N indicator values alone are not sufficient to
promote sewage/nutrient effects into a GEB/JBI mechanistic claim.
