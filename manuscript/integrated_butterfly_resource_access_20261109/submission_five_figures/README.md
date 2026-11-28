# Five-main-figure manuscript synchronization

Current author-local candidate: *Port-associated introduced floras reveal a weaker oceanic isolation gradient in mapped butterfly host resources*. The previous four-main-figure Markdown and exporter are preserved as a historical, separately reproducible version.

This directory connects the revised five-main/three-SI figures to existing public repository source files, rather than publishing another copy of the internal 240-file archive.

| New figure | Former location | Evidence |
|---|---|---|
| 1 | 1 | Fixed 239-species resource baseline and denominator accounting |
| 2 | 2A–B | Port coefficients and floristic Jaccard (7,875 pairs; 7,874 defined introduced comparisons) |
| 3 | 2D–E | 126-region/239-species isolation, 123-region introduced-only share, taxonomic rank and paired slope contrast |
| 4 | 3A, B, D | 670 contributors, all 239 scatter coordinates, four gap-filling groups |
| 5 | 4A | Independent 224-species records from 84 regions |
| S3.1 | New presentation of existing sensitivity results | Five directly paired port-origin models |
| S3.2 | 2C | 104-region marine electrical resistance; adverse strict transfer retained |
| S5.1 | S5.1–S5.2 | Original regional/transfer/within- and cross-South-America scores |

Old 3C (82/97), 4B (monarch/endemic examples), 4C (155-species prediction gains), 2F (seawater morphology) and the historic port-contributor figure are reported as text or tables in the revised Word manuscript; original evidence is retained. Do not treat them as deleted scientific results.

## Rebuild

From this directory, with Python 3.13 and the pinned packages:

```bash
python -m pip install -r requirements.txt
python rebuild_plot_data.py --download
python plot_figures.py
python verify_render.py
python build_current_word.py --outdir build
```

Alternatively use `--git-repo /path/to/full/chocho/checkout` to extract the seven pinned Git blobs. All seven SHA-1 Git blob IDs and all 17 canonical numeric fingerprints are checked before rendering. PNG/PDF output filenames use the renderer stems in word_candidate_manifest.json; map figS33 to Figure S3.1. The original local environment was Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Raster bytes can vary with renderer environment; the stored embedded PNG digests are exact reference hashes, not a promise of cross-platform bit identity.

## Current editable source and export

- [Manuscript.md](Manuscript.md): 116 text paragraphs and five figure positions; updated 28 November 2026 to clarify unequal resource gains in section 3.3 and its Discussion.
- [Supporting_Information.md](Supporting_Information.md): all 155 original text paragraphs, nine tables (314 cells), and three figure positions.
- [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md): eight captions extracted from these canonical documents, not independent edited text.
- [word_candidate_manifest.json](word_candidate_manifest.json): original Word/image digests, editable-source digests, and the bounded comparison.

The initial Markdown was byte-identical to independently extracting the author-local candidate; the current main text includes a subsequent editorial clarification (see the version history in word_candidate_manifest.json). The exporter preserves paragraph and table order, text, emphasis and embedded figure bytes; locally regenerated Word files were compared against the originals (22 and 27 rendered pages). The exact page layout of new exports may differ. The original submission candidate and internal archive remain private; do not publish them or treat CI previews as a submitted manuscript. For future intentional text/figure changes, update the source version and appropriate frozen count/digest checks in a reviewed commit.

No new raw ecological model is fitted here. No data licence is assigned and no stable reviewer access or author approval is asserted. PR #189 remains a separate internal source-audit branch. The three RF numerical mismatches and historical port-source attribution remain open.
