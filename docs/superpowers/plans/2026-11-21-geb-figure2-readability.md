# GEB Figure 2 Readability Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans. The user's action-first project instructions authorize reversible preparation and verification without a further design approval round.

**Goal:** Make the existing Figure 2 readable at its 6.6-inch manuscript width and deliver verified Word/PDF files.

**Architecture:** A deterministic SVG renderer reads the existing numeric receipt. It repositions the two archived model-fitted isolation curves without refitting or interpolating new ecological values. Four main figure numbers and the standalone eight-appendix SI are retained.

**Tech Stack:** Python standard library, SVG, existing CairoSVG/python-docx build, bundled document renderer.

**Spec:** The design and constraints below implement the user's continuation of the previous Figure 2 readability proposal.

## Design and alternatives

- Selected by the implementer: two rows of separately labelled A-D panels, with short chart labels and full methods/qualifications in the existing caption.
- Alternative: split Figure 2 into two numbered figures. This requires renumbering Figures 3-4 and all manuscript, SI and production references; the current readability issue can be solved within four figures.
- Alternative: retain the tall stacked layout. At 6.6 inches it reduces 12 px labels to 5.4 pt, so it does not solve the observed print-size issue.
- Preserve the original v05 SVG byte-for-byte and add the new readable SVG at a versioned path.

## Global constraints

- Source paper baseline: 4439d3b0afa9707dabe0d59aa30a518ebd4ac08b; immutable numerical evidence: 3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250.
- Preserve all numeric estimates, signed uncertainty intervals and source cohort definitions.
- Panels A/B: 126 botanical islands. C: 104 islands, 5,355 available dyads, 1,089 strictly evaluated dyads. D: 239 butterflies / 126 islands, separate introduced-only-share model on 123 active islands.
- Distinguish D's paired count-slope contrast p=0.0677 from share-model p=0.014; preserve the count curves' original log1p(distance) geometry.
- Minimum visible Figure 2 font 19 px on a 1,100 px canvas (8.2 pt at 6.6 inches).
- No model refitting, newly invented uncertainty, changes to source data, coauthor messages, merge or journal submission.

## Review focus

- A-C points and interval endpoints remain the source values on correctly scaled axes.
- D's native and contemporary counts are nested opportunities, not mutually exclusive groups or observed butterfly occurrence.
- The D share values come from a separate model and are not computed as curve ratios.
- Captions and accessible SVG description cover all four panels; manuscript/SI references resolve.
- Rendered Figure 2 labels and manuscript figures have no clipping, overlap or page overflow.

## Task 1: Rebuild Figure 2 and connect current sources

- [x] Create `scripts/render_GEB_v05_figure2_readable.py`; consume the current source receipt and archived v05 SVG, emit a new four-panel SVG.
- [x] Update the article figure path, caption path/four-panel wording, source-receipt display path and validator wording without changing numerical source fields.
- [x] Run the existing immutable-source editorial validator and compare all original numeric receipt content and curve coordinates.
- [x] Render and inspect the new figure at its delivered dimensions.

## Task 2: Verify manuscript files and preserve the result

- [x] Run the current v0.5 main/SI build and bundled DOCX page renderer.
- [x] Inspect all resulting pages; verify four main/two SI embedded images, anonymity and word counts.
- [x] Obtain independent code/scientific/visual review and fix material findings.
- [x] Record outcome and next blocking preparation item; prepare the verified change for the existing Draft paper branch. Delivery uses a normal fast-forward push only if the upstream baseline remains unchanged.

## Execution record

- Baseline sources fetched into an isolated worktree; original analysis worktrees remain untouched.
- Verification uses existing scientific source checks and direct rendered-page review. This presentation-only change does not add redundant implementation-mirroring unit tests.
- Full visual QA found and fixed one split SI table row, two export-only PR labels, a legacy supplement locator and a swallowed sentence stop. Numerical source content remains unchanged.
- Final local delivery: 22-page main and 10-page SI, 4,936 main words and 292 abstract words. Independent review closed all findings. See `manuscript/integrated_butterfly_resource_access_20261109/GEB_FIGURE2_READABILITY_QA_20261121.md` for evidence and remaining Stage A scope.
