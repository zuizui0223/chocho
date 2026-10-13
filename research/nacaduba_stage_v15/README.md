# v15 — Intraspecific flowering provenance adds a *time dimension* to butterfly host resources

**2026-10-14. Post-hoc source-verified plant flowering analysis. Not butterfly demographic data, not a causal phenology test, and not a replacement for the frozen GEB PR #38.**

## Direct primary source

[Kato et al. (2012), `Flowering phenology and germination ability of pollens for Acacia mangium and A. auriculiformis`](https://doi.org/10.1515/sg-2012-0029), *Silvae Genetica* 61(6):228–236. Original [open-access PDF](https://reference-global.com/download/article/10.1515/sg-2012-0029.pdf) (raw SHA256 `f66086e45aa7ce5bc6e8e2f4be94d8d7e1b21e161285b7f6411f9399f5a59aa66`). Tables 4–5–7 were checked against rendered pages, *without OCR*. Primary acquisition: [GitHub Actions 38424833145](https://github.com/zuizui0223/chocho/actions/runs/38424833145).

At Iriomote, Japan, **23 clones, 39 ramets, 8 source groups** were observed every 3–4 days in three April–March flowering seasons. Authors coded flowering onset and period in **10-day bins**, not daily flower-bud quantity. Thailand-labelled source g comprises 4 clones / 6 ramets, and other seven sources 19 clones / 33 ramets.

The original paper's three **clone-level** mean percentages of any flowering `65.2 ± 9.7`, `73.9 ± 8.8` and `60.9 ± 9.4` percent (mean ± SE) have been reproduced from 23 transcribed clones. All **69 clone-years** of Table 7 flowering/not-flowering flags reconcile exactly against Table 5; overall 79 flowering ramet×year events and three completely non-flowering ramets.

## Reconstructed flowering coverage

| April–March season | g flowered / 6 | other sources flowered / 33 | Combined union, 10-day bins | g-only bins | others-only bins |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2008–09 | 6 | 19 | **19** | **9** | 2 |
| 2009–10 | 6 | 23 | **17** | **6** | 4 |
| 2010–11 | 6 | 19 | **21** | **8** | 3 |

In this planted experiment, the g seed source was first in flower in **late July, early September and early August**, respectively; the non-g union started in **late October, early November and late October**. The addition of its cultivated trees changes pooled flowering-month categories, equivalent to **90, 60, and 80 nominal days of extra 10-day coverage**. These are an arithmetic summary of observed plant flowering *categories*, not 90/60/80 additional days of independently measured nectar, bud biomass, larval survival, insect colonization or flowering in Singapore. Seed source, ramet number, genotype/planting site and other aspects were not randomized for this new contrast; do not compute an independent effect size or p-value from pseudo-replicated ramet-years.

**Original-source discrepancy that must not be concealed:** the 2009 maximum duration is stated as *8 × 10 days* in the article's p.232 prose, while Table 7 p.232 contains source-g clone 146, `early Sep., 13` and clone 149, `early Sep., 10`. Our primary analysis reproduces Table 7 verbatim and flags the discordance. Artificially capping those two durations at 8 retains the **17** all-source bins and **6** g-exclusive bins in 2009.

## Biological hypothesis and rivals

This is the same known caterpillar host plant *Acacia auriculiformis* across all eight seed sources; thus *species-level plant occupancy does not fix flowering opportunity or plant intraspecific temporal resource traits*. That new botanical quantitative result can motivate a **provenance-mediated phenological niche** hypothesis for florivorous butterflies. However, it does not show the butterfly tracks this flowering and cannot justify extending local Singapore butterfly larval feeding to other host species/regions.

- [ButterflyCircle 2010](https://butterflycircle.blogspot.com/2010/01/life-history-of-two-spotted-line-blue.html) documents *Nacaduba biocellata* larvae consuming *A. auriculiformis* flower buds/flowers for approximately **9.5–11 days** in captive rearing, not wild cohort survival.
- [PROSEA, 1997](https://prosea.prota4u.org/view.aspx?id=2982) explicitly reports the same plant can flower year-round with location-dependent pronounced peak months; [Singapore NParks](https://www.nparks.gov.sg/florafaunaweb/flora/2/6/2693) calls it free-flowering, **not** a daily flower-bud-density time series.
- [Museums Victoria](https://collections.museumsvictoria.com.au/species/11366) documents mass *Nacaduba* emergence after rains in Australian populations, suggesting an alternative rainy-season pupal/adult emergence gate separate from larval floral resources.

A genuine test needs replicate trees/provenances, daily or weekly bud density and new flowers, egg cohorts tracked to pupae/adults, rain / temperature / predators / ant attendance, and sampling effort. Distinguish floral limitation of larval recruitment from adult appearance determined by emergence/dispersal.

## Reproduce GitHub source files

```sh
cd research/nacaduba_stage_v15
python audit.py
python -m unittest discover -s . -p 'test_audit.py' -v
```

The more extensive downloadable v15 archive includes 9 local regression tests, an independent NumPy auditor, rendered 10-day flowering-season figure (PNG/SVG), full per-ramet audit/mean±SE, primary metadata and checksum manifest. The repository subset reproduces key 10-day coverage from identical, blob-SHA-verified CSV source data. **Neither output changes the frozen GEB submission.**
