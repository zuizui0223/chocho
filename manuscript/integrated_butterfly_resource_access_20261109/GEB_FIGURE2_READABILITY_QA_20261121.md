# GEB Figure 2 readability and document QA — 21 November 2026

## Current state

The Figure 2 readability preparation is complete. The article retains four
main figures and separate Supporting Information with eight appendices and
two figures. This is a preparation update for Draft PR #179, not a journal
submission or completion of the outstanding source-access and citation gates.

Baseline paper commit: `4439d3b0afa9707dabe0d59aa30a518ebd4ac08b`.
Immutable numerical source: `3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250`.

## Changes

1. Replaced the active Figure 2 display with a source-driven 2 × 2 panel
   layout. The original v05 SVG remains byte-identical. At the manuscript's
   6.6-inch width, the minimum label size increases from 5.43 to 8.208 pt.
   Figure numbers and scientific estimates are unchanged.
2. Corrected the caption from three to four panels, specified the actual
   log(1 + distance) axis, and made the nested native/contemporary host counts
   explicit. The accessible SVG description correctly identifies 239
   butterfly species. The revised heading describes associations without
   adding a causal claim.
3. Kept short Supporting Information table rows on one page with repeated
   column headers. In Appendix S2, the complete denominator label now stays
   with the 54.31% and 29.05% values. Removed internal PR labels from the
   journal-facing SI export; the original source SI remains unchanged.
4. Resolved the Figure 4D legend to Figure S5.1 in Appendix S5, with a mapping
   for the unchanged inset's historical Supplement S.G label. Preserved
   sentence-final punctuation when redacting a bare author-repository URL.

## Source fidelity

An independent coordinate audit inverted the plotted scales and recovered
all three A estimates and six interval endpoints, both B medians, and all
four C estimates and eight interval endpoints. Panel D copies the two
archived paths and all six anchors unchanged, using only an affine display
transform. Its five x ticks and three y ticks retain the original mapping;
anchor counts agree with the receipt within the original SVG rounding
(less than 0.003 butterfly species).

The entire numerical receipt is unchanged except for its Figure 2 display
path. The original Figure 2 SVG, raw SI and 115,972-character legacy
reference/supplement archive are byte-identical to the baseline.

Source populations remain distinct: A/B use 126 islands; C uses a separate
104-island cohort, 5,355 available pairs and 1,089 strict test pairs; D uses
126 islands and 239 butterflies, with a separate share model on 123 active
islands. The paired count-slope contrast retains nine-Level1 p = 0.0677 and
95% interval −0.0172 to +0.3907. The share-model p = 0.014 and rounded
55%/66%/72% predictions are not presented as that count contrast or as ratios
of the count curves. No ecological model was refitted.

## Verification completed

- `render_GEB_v05_figure2_readable.py --check`: pass.
- `verify_integrated_geb_v04_editorial.py`, against the immutable source
  checkout and frozen v03 manuscript: pass.
- `build_GEB_v05_blinded_main_and_SI.py`: pass; abstract 292 words, main
  body 4,936 words, 39 references, four main and two SI embedded images.
- Both DOCX files rendered with the bundled document renderer. All 22 main
  and 10 SI pages were inspected as full-page PNGs. No clipping, overlap,
  blank pages or stranded headings were found. Following the final locator
  and punctuation edits, unchanged-page hashes were checked and the six
  changed pages were inspected again.
- The final PDF text contains no author repository identifier or internal
  PR labels. The sentence before “No butterfly presence” retains its full
  stop; the Figure 4D reference resolves to the distributed SI.
- Independent scientific, code and visual review closed all findings.

Page counts above describe the delivered local PDFs; other rendering
environments may paginate differently. The standalone Figure 2 PDF uses
the same 6.6-inch physical width as the manuscript.

## Remaining preparation scope

The source-access limits remain explicit, including the unavailable
original 155-species training panel and the absent independent refit. The
negative South America results remain in the main and supplementary figures.
The next main preparation objective is source-data and citation closure:
original acquisition/contributor citations, reviewer-accessible anonymous
archiving, source and figure rights, the identifying title page and
declarations, and the remaining author-list/citation-style audit. Existing
Stage A restrictions on review requests, merge and submission remain.
