# GEB reviewer-defense map — manuscript v0.2

This file records post-hoc reviewer-defense work separately from the active primary-analysis surface.

| Reviewer concern | Analysis / action | Result | Manuscript consequence |
|---|---|---|---|
| Same-family host null is biased by plant commonness / recording intensity | Native-range-matched, other-Lepidoptera-use-weighted, and combined strict host nulls; direct plant-level network analysis | Simple host-identity excess disappears under the combined null: 12,853 observed added units vs null median 12,971 (p=0.584). Across 8,909 plants, network degree remains positively associated with expansion within plant-family × native-breadth strata (r=0.295, p=0.0002). | Do **not** claim butterfly-specific selection of unusually globalized hosts. Positive result: network-prominent host plants are disproportionately redistributed, while degree remains a joint proxy for ecology and recording intensity. |
| Original portfolio architecture is structurally constrained by host count | Exact host-species-count adjustment; joint host-count + butterfly-Family adjustment | Raw rho 0.492/-0.486 becomes 0.059/-0.071; joint adjustment -0.039/+0.015 | Architecture removed from headline and retained only as a structural diagnostic. |
| WGSRPD finite ceiling could hide a generalist effect | Native-breadth restrictions and rank adjustment | Host-family breadth vs proportional expansion remains near zero across thresholds; adjusted correlation ~0.012 | No hidden broad-generalist advantage. |
| Geographic bias in HOSTS/LepTraits | Dominant native-resource-region stratification | Near-zero relation in Northern America, Africa, Temperate Asia, Tropical Asia and Southern America; modest positive Europe signal | Report regional robustness and retain uneven sampling as a limitation. |
| Species are phylogenetically non-independent | Butterfly-Family random intercept; Kawahara et al. (2023) exact-species PGLS | Family mixed model beta=0.021, p=0.752. Exact tree matches 124/239; Brownian rank-PGLS beta=0.140, p=0.119; Pagel lambda=0.033, beta=-0.051, p=0.581 | Species-level phylogenetic correction does not reveal a broad-generalist advantage. |
| Occurrence validation may be dominated by *Pyrgus communis* and pooled pseudo-replication | Species-level null, species-cluster bootstrap, all leave-one-out analyses | Mean species recovery 0.583 vs region-matched null median 0.385; 0/99,999 exceedances. Removing *P. communis*: 44/93 recovered vs null median 30; 0/99,999 exceedances | Occurrence result retained, but explicitly labeled secondary reuse of the climate-stratified panel. |
| Non-native butterflies such as *Pieris brassicae* complicate host-facilitation interpretation | Leave-one-out exclusion | Removing *P. brassicae*: 62/111 recovered vs null median 46; conclusion unchanged | Mention documented non-native populations and avoid causal host-facilitation language. |
| Negative climate test is underpowered | Distance matching + bootstrap interval + MDE | Distance-matched median filtering 0.714–0.750; host-breadth effect partial rho=-0.166, bootstrap 95% interval -0.583 to 0.261; ~|rho|=0.50 needed for 80% power | Climate moved to Supplementary Figure S1 and described as secondary / imprecise. |
| Monte Carlo p-value wording | Manuscript wording audit | 0/1,999 exceedances reported as Monte Carlo p<0.001 rather than treating 0.0005 as continuous precision | Corrected in abstract, Results and Figure 1 legend. |
| Procedure-heavy prose | Manuscript compression | Main text ~4.4k words; abstract <300 words; execution terms removed from main prose | Execution history lives under provenance rather than the manuscript. |

## Current defensible story

Human redistribution of host plants has substantially expanded reconstructed butterfly resource geography, but this gain is not preferentially larger in broad family-level generalists. A simple butterfly host-identity excess is absorbed when plant native-range breadth and network-wide Lepidoptera host prominence / recording intensity are jointly controlled. Direct plant-level analysis shows that network-prominent host plants are disproportionately redistributed even within plant-family and native-breadth strata. Secondary occurrence analyses show that contemporary introduced-host geography aligns with butterfly observations beyond structural overlap expectations.

## Claim boundary

The data do not establish local larval use of every introduced host population, causal host-facilitated butterfly range expansion, or an evolutionary preference by butterflies for globally redistributed plants.
