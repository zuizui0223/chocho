# BCE/Clarke partial host-inventory sensitivity for resource novelty (2026-10-17)

**Status:** post-hoc source-sensitivity design, written after the primary HOSTS-based 239-butterfly nulls were seen, but **before** inspecting the augmented-source resource-absorption endpoint. This is not a confirmatory biological study and must not alter existing GEB submission PR #38.

## Biological question

Does the accounting result that most introduced host × butterfly × region candidate incidences coincide with a butterfly's native-host-supported resource geography depend strongly on missing larval-host records in the frozen HOSTS inventory? Import the **same, already source-validated** European BCE/Clarke accepted WCVP-ID candidate host links used in Supplementary Table S9; do not select hosts by their effect on the novel endpoint.

## Source identity

- Independently derived, post-hoc Europe-specific foodplant checklist based on Butterfly Conservation Europe / Clarke (2024), DOI 10.1002/ece3.10834, evidence classes 1–3.
- Source: successful chocho GitHub workflow run `37795847333`, artifact `11557274784`, name `bce-clarke-anonymous-s9-review-addendum-v01`. Inner zip `bce_clarke_anonymous_s9_reproducibility_v01.zip` SHA256 `168d02efee9ece4c8841120be5a8cba712c1a02cb7cff9568ce251dc922bd7ab`; exact included crosswalk `inputs/bce_accepted_species_host_edges.json` SHA256 `f46d18703fe5cc8d97699af6f084f652a704e46b0f9e8c68e32dd70decc32c53`.
- All 1,027 distinct WCVP-ID candidate host links for 81 butterflies must be added **together**. Other 158 focal butterflies unchanged. Exactly 239 baseline butterflies must be present; assert 92 European matches.
- Frozen WCVP v13 commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`; HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`; reconstruct all original and new host distributions under identical extant/non-doubtful rules.
- Source is **not globally representative**; candidate host records are not independent original feeding experiments and WCVP geographic presence is not local butterfly utilization.

## Endpoint and required identities

For each butterfly (b), let original or augmented larval-host inventory be (S_b), each host's WCVP native regions (N_h) and contemporary regions (C_h). Define `native_b` as the union of (N_h), and `contemporary_b` as the union of (C_h). Introduced host incidences are `(h,region)` for every (h \in S_b) with (region \in C_h \setminus N_h).

Count separately:
- `H`: all candidate introduced host × butterfly × region incidences;
- `A`: their incidences within the **butterfly's current source scenario's native-host** union;
- `U`: unique regions in contemporary but not in native union;
- `R = (H - A) - U`: additional introduced host incidences within novel opportunity cells.
- `H = A + R + U` (integer identity). Note the BCE augmentation **can add native hosts**, so some previously new cells are **reclassified as native**.

Reproduce the archived baseline exactly: `(H,A,R,U) = (75346,46387,14406,14553)`; 239 butterflies, 206 with expansion, 26,530 native union cells and 41,083 contemporary cells. The augmented union geography MUST replicate earlier BCE Table S9: 29,852 native cells, 44,819 contemporary, 14,967 added, 213 species with expansion. Otherwise fail and write no positive result.

Report full-239 descriptive comparison and the 81 changed-inventory species separately; classify changes in unique new regions as genuinely newly supplied vs reclassified native. Report fixed sample exposure and geographic union counts, not invented source-independent effects.

## Stopping and claim boundary

This is a **single descriptive source-robustness** check. Do **not** treat the old +271–276 spatial conditional residual as confirmed by the augmented-source descriptive proportions, because a new matched spatial null would be needed for that stronger inference. Do not tune new candidates or drop European taxa to retain the effect. A direction change, large denominator change, or reproducibility failure must be reported openly. No further hypothesis-driven searching for more impressive p-values follows from this test.
