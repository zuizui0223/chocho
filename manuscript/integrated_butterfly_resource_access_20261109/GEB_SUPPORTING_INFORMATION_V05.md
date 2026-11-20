# Supporting Information — GEB Research Article

**Ports, the sea and two layers of island floras reshape butterfly larval-host opportunities**

**20 November 2026; journal-facing coauthor-review working version, unsubmitted.** Eight consecutively numbered appendices; figures S3.1 and S5.1. The underlying complete old v0.3 source archives remain versioned and separate; this is a source-qualified summary, not a replacement for downloadable data and code.

## Appendix S1. Source populations, taxonomic identities and geographic grain

| Analysis | Cohort | Measured response |
|---|---|---|
| Fixed known-host opportunity | 239 HOSTS butterflies × 369 WGSRPD3 botanical regions | Presence of ≥1 reported native or native+introduced larval host plant, **not butterfly occurrence** |
| Island–mainland host opportunity | 126 insular and 237 mainland botanical regions (6 mixed separate) | Fixed-source host availability and new gap filling |
| Port × whole plant flora | The **same 126** botanical islands, WPI ports × WCVP accepted species | Adjusted contemporary native/introduced plant richness and 7,875 paired island floristic similarities |
| Marine resistance | Source-matched **104** island botanical regions × **5,355** reused flora pairs | Native vs introduced plant Jaccard and undirected sea electrical resistance (not measured currents) |
| Independent later butterfly indexing | **155** butterflies, 208,482 species×cell test cases | Retrospective later indexed GBIF records, not population change |
| Independent GBIF community | **224** butterflies, 169 original primary regions, 14,196 reused region pairs | 2018–2024 **recorded** butterfly Jaccard conditional on earlier records and observation proxies |
| External island endemic | Independent 41 extant Macaronesian butterflies incl. 18 endemic | Named host records, independent wild/garden larval use, catalogue coverage |

Source HOSTS DOI **10.5519/havt50xw**, pinned data mirror 808e0b869f9ec1adf8efff87cf6a395adda103e0; Kew WCVP **v13 21 May 2024** DOI **10.34885/0yex-xv26**. Independent 224 species source GBIF download DOI **10.15468/dl.umqmsf** and annual positive/effort aggregates pinned to commit dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0. HOSTS is incompletely documented and WCVP contemporary regional plant status is not a dated arrival history. Other GBIF index and port acquisition families cannot borrow the independent 224-species DOI.

## Appendix S2. Native baseline accounting and resource concentration

| Source fixed host coverage | Island (126 botanical regions) | Mainland (237) |
|---|---:|---:|
| Butterfly×region potential cells | 30,114 | 56,643 |
| Native known-host covered | 4,829 | 21,215 |
| Added by introduced known hosts | **5,741** | **8,686** |
| Current known-host covered | **10,570** | **29,901** |
| Native baseline coverage | **16.04%** | **37.45%** |
| Fraction of initially missing resource combinations newly filled | **22.71%** | **24.52%** |
| Relative gain over native host coverage | **+118.89%** | **+40.94%** |
| Introduced-only source coverage as fraction of current covered combinations | **54.31%** | 29.05% |

For a geography class with native coverage B and initially empty fraction filled F, relative gain = F(1−B)/B. This arithmetic shows why larger proportional island gains do not demonstrate better island gap filling. The matched area×Level1 test (+10.7 percentage points introduced-attributable contemporary share) has a **different denominator** from raw gap filling, which is 22.7% island vs 24.5% mainland.

Among **670** source-accepted introduced host contributors, **38 (5.7%)** provided half of all 14,553 added butterfly×region opportunities using equal split of credits within each newly covered combination (plant Gini 0.7774, top-ten cumulative contribution 25.1%). The original 239 butterfly source native-host geographical breadth predicted absolute island additions (tied-rank Spearman ρ=+0.578896); broadest 59 species gained 2,271/5,741 combinations (39.56%), narrowest 60 gained 644/5,741 (11.22%). **82/97** one-host-family butterflies gained new global potential regions; excluding FAO crop hosts still left +46–47% global resource increase. These are mapped opportunities, not larval abundance, food quality or individual fitness.

Source-constructed butterfly-resource island-pair Jaccard increased +0.23667, compared with mainland +0.15044, but the Level1×island-class × fixed species/region gain-margin null left island-specific excess only **+0.00119**. This is not observed butterfly community homogenization.

## Appendix S3. Port-associated plant supply and sea-related geography

Matched 126-botanical-island WPI whole-flora associations, adjusting area and TDWG Level1: introduced β=**+0.2854** (HC3 +0.0622,+0.5087); native β=**+0.0954** (−0.0774,+0.2683); directly paired introduced-minus-native β=**+0.1900** (−0.0187,+0.3987). Nine overlapping leave-Level1-out fits retained positive introduced slopes +0.202 to +0.380. The introduced result is not a significant *between-origin slope contrast* or historical ship-cargo effect.

Across 7,875 pairs from the identical 126 regions, whole introduced flora median Jaccard was **0.046030**, native **0.004487** (~10.26×); this is contemporary composition of different species pools, not tenfold change through time. Separate 104-island original single-island holdout electrical resistance yielded introduced **+4.239%**, native **−0.110%** plant-Jaccard predictive MSE reduction. Strict both-island-endpoints-unseen 1,089 dyads in five source folds yielded introduced **+2.735%** (descriptive −1.649,+5.508%) versus native **−0.044%** (−0.984,+0.431%); paired difference interval crosses zero. A *different* shortest water path helps both origin classes. Alternative Simpson similarity can reverse origin ordering. The electric resistance is undirected model geometry, not observed currents.

Native-minus-introduced share of plants with EuDiS recorded seawater-dispersal morphology was **+2.9515 percentage points** among 22 biased/mostly European island regions. The independent published Tonga Table 1's broader *water dispersal* category points in the other direction (native 137/364 vs post-European introduced 137/309); **not** the same trait ontology or a matched replication.

**Figure S3.1. Original historical port contributor plot (same source data, preserved).**

![Figure S3.1. Ten influential recorded host plants at original 38 sampled botanical port buffers, not 1,596-port historical arrivals.](../../figures/integrated_island_port_butterfly_v03/Fig2_global_contributors_at_ports.svg)

## Appendix S4. Same-source HOSTS×WCVP potential resource controls

On the original 104 botanical islands, H=668 introduced WCVP plant taxa registered as hosts of the 239 butterflies; O=7,219 introduced taxa **not registered** as hosts (not confirmed true non-hosts). Under original geography E controls, H-after-O improved prediction of *WCVP/HOSTS-derived potential butterfly resource composition* by +17.609% MSE. With 199 plant-family×island-occupancy-matched label draws, observed source H-after-O +16.365% versus median matched null +7.699% (9/199 null scores at least as large, descriptive rank .050). Contemporary whole alien flora similarity improved source-potential butterfly-resource prediction +26.363% against a 225 serial-state fixed source-gain-margin null median +18.974% (excess +7.389 points). These models share WCVP/HOSTS taxonomic input between predictors and potential response and **are not independently observed butterfly evidence**.
