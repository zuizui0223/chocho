# v07b. A source-stable introduced-only region can hide documented native larval hosts

**2026-10-12 JST. Targeted post-hoc source-completeness diagnostic, NOT realized breeding evidence.** The frozen GEB PR #38 was not updated.

## A more important warning than the ~0.1 percentage-point v07 residual

Across 32 butterflies and 1,081 observed butterfly × WGSRPD3 cells, the four crossed host-source variants classify **64 adult-observed species×regions across 18 butterflies** as introduced-only larval host resource, totaling 2,638 mapped adult occurrence records. Such agreement can reflect **shared missing source coverage**, not strong independent evidence of native host absence.

*Nacaduba biocellata* contributes 9 such cells and 1,335 adult records, and *Pyrgus communis* 22 cells and 736 adult records (together 31/64 source-stable cells, 2,071/2,638 adult GBIF occurrence records). These observations are not individual larvae or breeding records.

The pinned original HOSTS resource represents Nacaduba with only one host ID, **Vachellia farnesiana 2452405**, in all four scenarios. However [South Australian Butterflies and Moths](https://sabutterflies.org.au/home/food-plants.html) explicitly lists *N. biocellata biocellata* larval feeding on *Acacia cyclops* (ID **2611847**, pinned WCVP native SOA, WAU), *A. pycnantha* (**2613563**, NSW, SOA, VIC), *A. melanoxylon* (**2612971**, NSW, QLD, SOA, TAS, VIC), and *A. saligna* (**2613720**, WAU). [Museums Victoria](https://collections.museumsvictoria.com.au/species/11366) likewise documents Australian Nacaduba larvae on wattles.

**Targeted sensitivity only:** add those four previously unrepresented, literature-documented links without changing plant ranges, other species or any original source files. Nacaduba's potential resource regions become:
- Native: **50 → 56**;
- Contemporary: **131 → 165**;
- Introduced-only additions: **81 → 109**.

Four adult-observed Australian regions NSW (272 records), QLD (295), SOA (296), WAU (300) change **introduced-only → native host present** (1,163 adult GBIF records). VIC (356) changes **no known host → native host present**. These numbers represent *mapped adult occurrence records*, not unique individuals or verified on-site larval feeding. The apparent 64 source-stable adult-observed introduced-only cells are **60** after the targeted links and their record count 2,638→1,475.

In the 32-species occurrence panel, introduced-only recovery among adult-observed cells outside the native mapped host envelope changes **64/113 → 60/108**; unexplained adult occurrence cells change **49 → 48**. Across all 239 butterflies, in the marker+BCE reconstruction native/contemporary/added geographic resource cells change **29,855/44,927/15,072 → 29,861/44,961/15,100**. This handpicked source sensitivity is NOT a proposed global correction or new manuscript value.

## Boundaries

The cited Australian larval-host compilation is species/subspecies-level, not a GPS/date matched larval cohort for each WGSRPD3 region. WCVP records are broad botanical presence. We cannot assert any specific adult butterfly fed on these hosts or that populations colonized due to introduced plants. Australian subspecies data should not automatically be transported to every global subspecies.

The separate high-count *Pyrgus communis* warrants a per-region native-host audit: [UC Davis Art Shapiro](https://butterfly.ucdavis.edu/butterfly/pyrgus/communis) documents both introduced *Malva* and native *Malvella/Sidalcea* host usage, but no given eastern US region is proven misclassified by that general information.

**Ecological decision:** Source-stable introduced-only larval-host geography should be filtered for native-host completeness before being treated as independent validation of larval colonization. Further spatial null p-values should not substitute for field-confirmed larvae, eggs, host-specific feeding, and development to adults.

Evidence: `provenance/reviewer_defenses/results/butterfly_nacaduba_host_gap_v07b_20261012.json`. Reproducible v07b ZIP with fixed parent v06 source, 64-cell ledger, independent matrix-check, 3 passing regression tests and the original **169-test** parent-suite log is available separately; no cumulative v06/v07b code has been merged into this branch.
