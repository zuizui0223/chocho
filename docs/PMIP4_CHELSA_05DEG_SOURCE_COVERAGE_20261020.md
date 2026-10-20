# ICE-6G_C / CHELSA 0.5° source compatibility (20 October 2026)

**Verified real-data result, pre-biological screening.** The geography-compatibility numerical calculation [GitHub Actions #39366228063](https://github.com/zuizui0223/chocho/actions/runs/39366228063) completed source checks and the 259,200-cell 0.5° availability census and uploaded source output artifact **12239910484**, digest SHA256 `706d482c5e139675951a170e538ceee11249281c77e38da99f37f2b3b7524947`. The workflow's final branch commit step **failed because the branch advanced concurrently**; therefore the overall workflow status was failed although the source math and artifact production steps passed. This document cites the original immutable artifact rather than claiming a successful archive commit.

The original native 10 arcminute PMIP4 ICE-6G_C `sftlf` and `sftgif` files for ages 0 and 21 ka were source MD5-pinned independently by [Actions #39356535149](https://github.com/zuizui0223/chocho/actions/runs/39356535149), and original 0BP/21ka CHELSA-TraCE21K four bioclimatic variables were ZIP-SHA verified before raster access. ICE-6G_C PMIP4 physical 0.5° cells are formed by grouping the source 10' pixels into 3×3 native cells, lat flipped, longitude rolled so the output is north-up −180°..+180°; source terrain is `land=100 AND ice=0`. The climate grid is read at 0.5° nearest-neighbor from EPSG:4326; no climate values were imputed.

| Ice-free land threshold in 3×3 source cells | PMIP 0ka ice-free half-degree cells | PMIP 21ka ice-free half-degree cells | Newly ice-free 21ka vs 0ka | 21ka climate values available for newly ice-free cells | 21ka climate missing on continuously ice-free cells |
|---|---:|---:|---:|---:|---:|
| Any 1/9 | 63,340 | 60,913 | 9,104 | **9,104/9,104** | 311 |
| Majority >=5/9 (primary, source-frozen) | **58,765** | **57,186** | **8,356** | **8,356/8,356** | **311** |
| All 9/9 | 54,707 | 54,270 | 7,953 | **7,953/7,953** | 311 |

**Important distinction:** `newly ice-free` counts include both land emerging from the modern ocean and land losing ice. They are not exclusively paleocontinental shelf; 0.5° coarsening further obscures which subcell changed, and does not guarantee fine-scale continuous cross-coastal passability. Missing climate values on continuously ice-free regions remain explicitly missing (not assigned absence or poor suitability). The original 4-var climatology is finite in 243,746/259,200 present cells and 242,372/259,200 LGM cells; full-grid climate includes some ocean values which are then gated by physical substrate.

The source-only release consists of `summary.json`, `coarse_masks_finite_05deg.npz`, and a provenance receipt in artifact 12239910484. The exact source files SHA256 are `2c18201c522906e7872617488d7382ad5d64ddb57741237024bd7fb488c1a15d` (JSON) and `d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db` (NPZ).

**Ecological boundary:** 0.5° connected potential climatic habitat built from a modern native-WCVP-trained SDM is not proven plant occupancy in an ancient coastline, species-specific migration, larval food presence or an independently confirmed genetic barrier. The TTF historical host-specific climate negative control remains frozen negative.

The next runner replays the original 416 native plant 0BP SDMs **at all original 6,337 scoring sites** with fixed maximum absolute error 1e-5 before projecting to 0.5° potential habitat. An early real execution found an exceptional source-reproduction mismatch (WCVP plant ID 2543867, *Paullinia pinnata*, maximum difference 1.205×10⁻⁴). We are explicitly diagnosing source identity/numeric fidelity and will not relax the frozen 1e-5 constraint or silently omit the species to generate connectivity results.
