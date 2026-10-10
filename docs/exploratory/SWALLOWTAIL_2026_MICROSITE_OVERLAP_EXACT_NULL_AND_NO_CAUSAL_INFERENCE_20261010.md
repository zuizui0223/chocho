# Butterfly shared-host encounter: exact public-data nulls versus sunlight-only counterexample

**Date:** 2026-10-10  
**Status:** SOURCE-BASED MATHEMATICAL COMPARISON VERIFIED; REAL ECOLOGICAL CAUSALITY UNIDENTIFIABLE.  
**Exploration repository branch:** `exploration/realized-host-window-v01`; GEB manuscript PR #38 unchanged.

## Scientific origin and source type

[Jang et al. (2026), *Journal of Ecology and Environment* 50:15, DOI 10.5141/jee.26.017](https://doi.org/10.5141/jee.26.017) studied *Sericinus montela* and *Atrophaneura alcinous* feeding on ***Aristolochia contorta*** along four riverbanks in Korea. The article already published light-related microhabitat segregation and larval preference/performance results. The first author's 2024 SNU master's thesis also studied the same butterfly niche-partitioning system: **not an independent replication merely because it preceded the 2026 journal paper**.

We authenticated the source only as **published aggregate counts**, not original ramet/quadrat/leaf observations. 2026 paper: Methods `255 quadrats`; written site denominators sum `50+51+57+57=215`; the `S. montela only` table header gives `155` but site and month breakdown each sum `115`. No author-verified correction or raw 215/255 observation file was available; **all quadrat calculations below are conditional on the internally consistent 215 breakdown**.

## Published observed versus fixed-marginal random placement: exact calculation

**Host ramet records:** total **875**; 542 neither, 204 *S* only, 124 *A* only, 5 both. Fixed two-species margins are **S-positive=209** and **A-positive=129**. Under a deliberately naïve null that randomly overlays the two species' positive labels among all 875 identical, interchangeable observations:

- expected overlap `209 × 129 / 875 = 30.81257`, observed **5**;
- **observed/expected = 0.1623**; **83.77%** below naïve reference;
- exact one-sided hypergeometric lower tail **P(X≤5)=4.58×10⁻¹¹**; central 95% reference interval **22–40**.

**Published quadrat source totals under 215 interpretation:**

| Stratify the *same* printed quadrat table by | Observed both | Fixed-margin expected | Observed/expected | Exact one-sided reference tail | Central 95% null overlap interval |
| --- | ---: | ---: | ---: | ---: | ---: |
| All four **sites** separately (AY, CJ, JM, PT) | 11 | **48.275** | 0.228 | **8.35×10⁻²⁸** | 41–55 |
| All five **months** separately (May–September) | 11 | **44.503** | 0.247 | **2.33×10⁻²⁵** | 38–51 |

Site-specific and month-specific exact tails are obtained by **convolution** of hypergeometric distributions fixing S/A marginal counts within each site or within each month. They are two **alternative partitions of the SAME published table**, not separate studies or independent confirmation. The 2026 paper's field survey was deliberately set up around where larvae appeared; it's not a simple random sample of all possible host sites.

**Very important inferential limit:** These probabilities are mathematically exact **conditional on assumptions that are biologically unjustified** for random mixing across heterogenous observations (identical plant availability, light, host size, detection, sampling independence, fixed sampling frame). They are *not valid confirmatory p-values for ecological exclusion or interspecific competition*. Ramet observations may be nested/revisited, survey effort is selected, and site/month controls do not adjust local irradiance within quadrats. The original site×month×light×ramet data would be necessary to validate a biologically informative conditional null.

## Constructive nonidentifiability proof: no direct butterfly interaction is necessary

We explicitly constructed a **hypothetical** (NOT from the paper's source records) two-stratum population of 875 equivalent host records, respecting the paper's exact species marginal counts. Within each sunlight stratum, **the two butterfly species' presences are independent**:

| HYPOTHETICAL class | Host records | *S. montela* occurrences | *A. alcinous* occurrences | Expected both if independent **within** class |
| --- | ---: | ---: | ---: | ---: |
| Sunny | 400 | 200 | 5 | 200×5/400 = 2.50 |
| Shady | 475 | 9 | 124 | 9×124/475 = 2.35 |
| **Total** | **875** | **209** | **129** | **4.85** |

**Stronger constructive check: replicate the paper's entire four-cell contingency table exactly, not only its margins or an approximate overlap expectation.** Choose one possible integer realization in the same **hypothetical** two light strata:

| Invented sunlight class | Neither observed | S only | A only | Both observed |
| --- | ---: | ---: | ---: | ---: |
| Sunny (400) | 198 | 197 | 2 | **3** |
| Shady (475) | 344 | 7 | 122 | **2** |
| **Total** | **542** | **204** | **124** | **5** |

These are **exactly the four observed category totals published by Jang et al. (2026)**. Within each invented stratum the co-occurrence of 3 versus independence expectation 2.50 and 2 versus 2.35 is entirely compatible with the no-competition independence example. Thus **even exact agreement with all four published contingency cells cannot exclude environmental sorting as a sufficient explanation**. This example is a constructive mathematical demonstration and **is NOT a reconstruction of Korean ramet-level irradiance**; no actual light strata were retrieved.

Thus **expected co-occurrence of ≈4.85** appears while the globally pooled naïve independence expectation is **30.81** and observed Korean overlap is **5**. This toy counterexample **does not fit, verify or even estimate the actual light classes in Korea**; it rigorously demonstrates that the published pooled margins, even with a large formal random-placement deficit, **cannot identify competition versus environmental sorting**.

A light/leaf-quality axis causing both habitat choice and larval presence is biologically plausible from the 2026 published study itself (relative light ORs in opposite directions), but whether it explains all five actual co-records remains untested.

## Consequences for chocho ecological claims

1. **Publishable novelty is NOT** “butterflies that can eat the same plant have surprisingly low overlap,” which is both an already published finding and readily explainable by light heterogeneity. Nor is proving an impossible zero ecological contact while both species occupy a shared host's range.
2. The global GEB WCVP/HOSTS study measures **host-plant supply opportunity**, not same-ramet encounter, contact frequency, enemy transmission or viable adult recruitment. A GEB host map cannot provide those rates; don't silently reframe species-range coincidence as competitive interaction.
3. An independent ecological paper must test a deeper ecological result: does **microhabitat-conditioned overlap** differ from an appropriately site×date×light×leaf-area conditional null, AND if the butterflies are forced/allowed to truly meet, how does experimentally randomized competitor density alter native butterfly adult outcome under independently controlled food-access regimes?
4. Require public or author-provided row-level 2026 quadrat/ramet data to estimate a joint site×month×light test. Confirm original 215 vs 255 article denominator with corresponding author; don't invent a correction or bootstrap nonexistent ramet observations.
5. The Korean study uses *A. contorta*, while the Kyoto source study uses native *A. debilis* and introduced-to-Japan *S. montela*. Their values are not geographically or botanically interchangeable.

## Reproducibility and STOP decision

- Frozen independent mathematical question: `SWALLOWTAIL_MICROSITE_OVERLAP_NULL_PROTOCOL_20261010.json`.
- Complete Python script: `scripts/analyze_2026_swallowtail_overlap_nulls.py`; original printed-count fixture in `scripts/audit_2026_swallowtail_microhabitat_published_counts.py`.
- [GitHub Actions run 38056591090](https://github.com/zuizui0223/chocho/actions/runs/38056591090): **5 mathematical/source-limit tests**, including a constructive check that **all four** published ramet totals are exactly reproduced by hypothetical light strata; JSON receipt uploaded.
- **No new biological/causal effect claimed**. This is a valuable mathematical check revealing that published aggregate segregation is underdetermined, not proof of competition and not a second empirical replication.
- **STOP** optimizing null assumptions to claim a strong interspecific interaction on published count margins. Without matched actual light and time observations, a causal comparison is not identified. Keep GEB submission PR #38 isolated.
