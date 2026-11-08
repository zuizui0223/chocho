# Frozen geometry-only current-map specification — 8 November 2026

## Main question
Where would unit electrical current concentrate on the existing source-frozen ocean graph if each of the original 1,275 unordered island pairs were equally weighted?

## Fixed inputs
Use the exact original 51 island groups, 194 WPI2019 ports, 182 marine-qualified terminal ports, GSHHG 2.3.7 topology with 0.0625-degree fine water cells and 0.5-degree conductivity blocks, and unchanged frozen uniform-ocean conductivity. No source geography or conductivity is chosen using the butterfly or botanical response. The full existing 189,550 marine graph nodes and 370,996 edges must be retained.

## Algorithm
Use a grounded sparse Laplacian and island 0 as the arbitrary reference. Solve 50 independent island-minus-reference unit-current potential columns, each preserving the original distributed port weights. For any pair (i,j), the difference of columns gives the exact model pair flow. For every graph edge e, calculate mean absolute electric current over all 1,275 equally weighted unordered original pairs. Compute the absolute difference mean with a sorted potential-drop identity to avoid 1,275 independent solves. Convert edge flows to water-node through-current by assigning half each edge's absolute flow to both endpoints, then sum only for display within 0.5-degree geographic grid cells. Do not alter the biological results.

## Independent validations and deliverables
Unit tests: sorted-pair mean against enumeration, invariance to potential gauge shifts, exact 3-island series circuit, and 4-node parallel circuit with distributed ports. Full run must reproduce three fixed pair energies from the original archived 51-by-51 effective-resistance file within relative tolerance 2e-6; max normalized grounded linear residual must not exceed 1e-7. Export machine-readable node current, original verified centroids, provenance receipt, PNG, SVG and PDF.

## Inferential boundary
A high-current water corridor is a modelled connectivity feature, NOT evidence of inter-island butterfly movement, ship propagule transport, or an island stepping-stone hub. The frozen network models no transport through an island; ports are source/sink allocations and ports on the same island are not short-circuited. Do not label an island as a relay solely from coastal current. AIS and wind fields may be run separately with the same program as explicitly labelled sensitivities but the registered principal map uses the uniform geometry field.
