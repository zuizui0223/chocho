# LGM host-specific fragmentation: v0.8 botanical negative control (independent audit v0.9)

**2026-10-22 — exploratory ecological finding, not ancient occupancy or genetic evidence**

## Current verdict

Across the original source-bound modern butterfly/moth sampling geometry, a source-fixed 0 BP- and 21 ka potential climatic larval-host patch can sometimes split. However, **we did not detect a systematic additional glacial patch-fragmentation effect specific to the insect's recorded larval host identity**, compared with *other real plant taxa* previously matched on their modern 0 BP observed climate moments.

This is a **failure to support an additional recorded-host-specific signal**, not a formal equivalence claim, and not a demonstration that larval-host resource constraints never cause historical genetic differentiation.

The earlier TTF v0.5 matched-real-other-plant climate analogue comparison was also negative for host-identity specificity. This new 2D climate/PMIP4 land-ice graph control does NOT retune or rescue the TTF result.

## Why the source evidence is informative

The source data are not simulated species selected after seeing the 21 ka response:

- **335** insects originally had alternative *real* plants selected in a TTF-only **0 BP** climate-signature matching plan. Each list excluded the insect's **recorded host IDs and host genera**. These selected plants are *unlisted*, not guaranteed physiologically unsuitable as larval hosts.
- Native WCVP region-screened 0 BP climatic niche models were available for all recorded hosts in 294/335 insects. **292/335** retained at least four preselected alternative plant groups with original source climate models; exactly **2,108** alternative groups available before opening any v0.8 LGM-patch response.
- Every actual-vs-other-plant group comparison uses the **same modern sampling geographic edges** on which all four conditions are contactable: known host 0 BP, known host 21 ka, other plant 0 BP, other plant 21 ka. PMIP4 ICE-6G_C land/ice is source-pinned (majority ≥5/9 10′ subcells), 0.5° grid and source plant SDMs are unchanged, and the nearest suitable patch must be on the same four-neighbor physically ice-free land component.
- The per-group comparison is the *difference in within-pair host-climate patch contact change*: (real-host 21 ka−0 BP share of edges on the same potential patch) − (other-plant 21 ka−0 BP share on the same edges). **Negative** values mean the real known host became more split. If no original edge is comparable under all four conditions, the donor group has **no numerical contrast** (not 0).

The source 2,108 groups are correlated technical negative controls. They are **not** 2,108 independent insect species.

## Empirical source result

| Prespecified source case | Prequalified insects | Any four-way matched edges | Real-host more split | Other plant more split | Same change | Median species contrast |
|---|---:|---:|---:|---:|---:|---:|
| Analog-only, 50 km | 292 | 203 | 29 | 43 | 131 | 0 |
| **Analog-only, 100 km (main)** | **292** | **226** | **62** | **79** | **85** | **0** |
| Analog-only, 150 km | 292 | 235 | 79 | 86 | 70 | 0 |
| Include novel climate, 50 km | 292 | 208 | 28 | 43 | 137 | 0 |
| Include novel climate, 100 km | 292 | 227 | 61 | 81 | 85 | 0 |
| Include novel climate, 150 km | 292 | 236 | 78 | 86 | 72 | 0 |

**Main source-specific evaluability:** 66/292 were NOT evaluable with even one matched real-other-plant group with a single original edge contactable in all four conditions. The main 226 evaluable insects represent just **173 distinct recorded host combinations** (181 insect genera). Sharing hostsets and other-plant donors across insects means that source modelled resource landscape outcomes are not independent across all 226 taxa.

**Frozen stability reporting condition:** With ≥4 distinct original donor plant groups, each offering ≥5 four-way matched original edges, **176 insects** qualify for a descriptive stricter readout: **56 real-host more split, 67 other-plant more split, 53 equal**, median contrast still **0**. Including extrapolated novel climate yields 177 stable insects, signs 56/68/53, median 0. **This is a prespecified denominator check, not a fresh independent confirmatory cohort.**

### Is current climate matching perfect?

No. The species-level median *signed* 0 BP same-patch baseline gap is zero, but:
- **145/226** evaluable species have exactly zero **median** modern baseline gap between real and other-plant landscape same-patch contact.
- **52/226** evaluable species have a **median absolute current same-patch gap >0.1 (10 percentage points)**.
- The other-real-plant source donor ranking matched modern observation climate **mean/SD**, *not* the modern native-WCVP-trained patch topology. These differences are not equivalent sampling controls for historic plant range or observed host use.

A median of zero can arise from a point mass at zero caused by coarse 0.5° patch classification and scarce four-way comparable edges. **It is not a confidence interval for equivalence or proof of no effect.** The difficult-to-evaluate 66 species are not automatically neutral.

## Independent audit v0.9

A separate source reaggregation read the immutable selected v0.8 records at their actual grains:
- Original frozen v0.8 0 BP-only eligibility: 292 insects, 2,108 available real-plant groups;
- All **12,648** group × scenario × radius records, **8,668** with a four-way matched original edge and **3,980 without**;
- All **1,752** species × scenario × radius rows;
- Every donor ID and 0 BP climate matching score against its original pre-LGM fixed plant identities;
- Per donor group arithmetic rechecking 0 BP and 21 ka same-component fractions, split and merge counts, modern baseline gap and exact Δ(real)−Δ(other) difference;
- Independently recomputing per-insect medians, the four-control-five-edge predeclared readout, every available/nonavailable sign count and all source summaries.

The recount agrees with v0.8. Audited compact JSON source:
`results/lgm_v08_otherplant_independent_audit_v09/summary.json`
with frozen SHA256 source identity digests. Its audit has **no refitting of 416 plant models, no change to climate or PMIP4 source masks, and no genetic response**. The test script and source conditions are in `scripts/audit_v08_otherplant_lgm_null_independent_v09.py` and `tests/test_v08_otherplant_independent_audit_v09.py`. The independent recount described here was executed by a separate read-only model-data aggregation; no claim is made that those new Python tests passed GitHub CI.

Original source v0.8 [successful run #39408984271](https://github.com/zuizui0223/chocho/actions/runs/39408984271) **PASS**, its full result data and audit receipt remain untouched.

## Ecological interpretation: what survives and what does not

1. Historical physical PMIP4 land/ice and CHELSA potential climatic suitability often **rearrange** the modelled distribution of climatically passable larval-host plant habitat. This is an ecological *model prediction*, not independent ancient data.
2. The source v0.6 exact-modern-geography pair comparison shows more hypothetical split than merge transitions within some regions but is strongly conditioned on modern anchors that happen to be contactable in both epochs.
3. When the relevant control is another real plant previously chosen using modern climate only, **the additional signature specific to the recorded host plant identity is not supported** across the source-qualified, four-way-contactable panel (this v0.8).

The robust ecological question for a genuinely new programme is whether **independently reconstructed/validated occupied LGM plant distributions**, with taxon-specific dispersal and refuge geography, correspond to independently sourced present-day within-species genetic discontinuities. This cannot be answered by repeatedly retuning the current climate-mask cutoffs, radius thresholds or posthoc taxonomic subsets. A credible next decisive dataset would use verified palaeoecological records / regionally native palaeoniches alongside a **fresh, outcome-unseen genetic panel**. If neither is available, the appropriate publishable object is a carefully labelled comparative model-based glacial climatic-resource *opportunity* and geography study, not an inferred genetic barrier.

**Inferential stop rule:** Do not relabel this null result as host-specific LGM support, do not recode unevaluable as zero, do not confuse 173 distinct source hostsets with 226 independent host graph realizations, and do not claim occupied host corridors or actual glacial insect ancestral populations.
