# Chocho next ecology: causal sign reversal is the remaining distinct question

**2026-10-10** | **Evidence-backed decision note, NOT a new ecological result** | `exploration/realized-host-window-v01` | GEB PR #38 untouched

## Why the previous ideas alone are not new

Static expansion of known host ranges is partly true by construction. Exact-host shared geography is only **potential** co-exposure, not same-plant use or competition. Repeated observational butterfly range/opportunity tests, host-lag probes and the source-verified 19-year *Melitaea cinxia* host-decline model all failed independent biological or heldout prediction gates.

The previous GloBI exact-*Plantago lanceolata* source search for *Junonia coenia* (217 eat rows) and *Anartia jatrophae* (33 rows) yielded **zero** original exact-plant, larval feeding events for either species. The search was finite and non-exhaustive; this is **failed evidence support for a specific wild shared-host system**, not biological proof of absent consumption. Do not claim natural JcDV cross-species transmission from geographically overlapping distributions.

**Important literature overlap discovered 2026-10-10:**

- van Nouhuys & Kraft (2012), *Population Ecology*, https://doi.org/10.1007/s10144-011-0302-5: a field experiment already showed a **butterfly–butterfly indirect interaction mediated by a shared pupal parasitoid**. *Melitaea cinxia* had *lower*, not higher, parasitism with *M. athalia* present, attributed to host accessibility. Apparent competition and short-term apparent commensalism in butterflies are **prior art**; the local and landscape signs can oppose.
- Cuny et al. (2023), *Oecologia*, https://doi.org/10.1007/s00442-023-05465-z: sequential heterospecific herbivory and parasitoid attack on shared *Brassica* already caused **asymmetric plant-mediated indirect parasitoid effects**. In its comparison the hosts were a moth (*Mamestra brassicae*) and a butterfly (*Pieris rapae*)—not two butterflies. Field-collected host-site responses and adult-emergence effects cannot be inferred from this paper alone.
- Christensen et al. (2024), *Ecology*, https://doi.org/10.1002/ecy.4282: plant-mediated horizontal/vertical **within-species** densovirus transmission already established in *Euphydryas phaeton*.
- Muchoney et al. (2024), *Journal of Invertebrate Pathology*, PMID 39159850: same JcDV × exotic *Plantago* experimental diet across *E. phaeton* and *Anartia jatrophae* already compared infection/survival after **separate experimental inoculation**. Different species' susceptibility is prior art; it does not certify actual interspecific transfer.
- Ragonese et al. (2025), *Ecological Entomology*, https://doi.org/10.1111/een.70010: host diet × temperature × OE infection effects on monarch survival and parasite tolerance already published.

**Corrected novelty floor:** A new community ecology paper must identify whether **human introduction of host plants reverses the *causal sign* of interactions among butterfly consumers**, rather than merely add host opportunities, reproduce a single-species parasite effect, or rediscover host sharing / enemy-mediated indirect interactions.

## A sharply identifiable causal question

> **At equal available edible host biomass, do introduced versus native host plants change whether the presence of a second butterfly causes a positive, null or negative *enemy-mediated* effect on viable recruitment of a focal butterfly?**

Use biological **interaction direction** (per-capita viable adult offspring of a recipient after a heterospecific neighbor is introduced) as the response, rather than number of available butterfly × WGSRPD3 cells.

A minimal balanced experiment is a **host plant species × neighbor presence × enemy access** factorial. Quantifying a general introduced-versus-native host origin effect requires **at least two independently chosen accepted host species per origin**, each able to support the focal recipient's larval development. Holding a single tropical species against one swamp species can only identify a *species identity difference*, because origin is perfectly confounded with identity.

**Causal contrasts** (randomize plant/enclosure as the unit, stratify by host plant species and site block):

- (Y(h,n,e)) = flight-capable focal adults produced per original focal cohort on host plant (h), with a heterospecific neighbor (n\in\{0,1\}) and enemy access (e\in\{0,1\}). All assigned animals, including deaths, remain in the denominator.
- Neighbor effect on plant (h) when enemies present: (N_{h,1}=E[Y(h,1,1)]-E[Y(h,0,1)]).
- Enemy-mediated neighbor effect: (M_h=N_{h,1}-(E[Y(h,1,0)]-E[Y(h,0,0)])).
- **Main introduced-host sign-shift contrast:** mean (M_h) across introduced hosts minus mean (M_h) across native hosts, with uncertainty across **independent plant species/populations and experimental sites**, not merely larvae.
- Host-origin inference is not valid if only one botanical species represents a status. Consider origin as a descriptive grouping of chosen host species, not a randomized native/exotic property.

**Non-enemy controls:** Equalize focal density, initial plant biomass and *accessibly edible biomass through development*; define comparable plant developmental stage; measure resource depletion and plant microclimate. Sham cages are necessary because enemy exclusion changes temperature/humidity. Record plant chemistry or pathogen status where feasible but do not invent a mechanistic covariate.

**Enemy identity:** A negative effect under enemy access can arise through predation, parasitoids, pathogens, changed plant quality or resource depletion. Separately measure confirmed parasite identity and acquisition / rearing outcomes; do not label it a pathogen-mediated cross-species link without causal source verification. A population-only result without species-specific infection evidence supports conditional neighbor effects, not a viral transfer route.

### Opposing predictions

- **Resource facilitation:** neighbor access does not lower viable adult production or is positive when shared hosts provide local resource concentration without serious depletion; enemy-mediated term near zero.
- **Enemy amplification:** neighbor attracts, amplifies or carries natural enemies; recipient adult production falls despite equivalent host biomass, in enemy-open but not enemy-excluded arms.
- **Enemy dilution / accessibility buffering:** adding a more attractive or accessible alternative butterfly causes **higher** focal adult production (less attack), possibly only on particular plants.
- **Environment-dominated:** apparent sign switches vanish with matched microclimate/plant quality/enemy exposure. Static origin labels explain little or nothing after exact plant species context is controlled.

Any sign is meaningful **only** if replicate populations and biological units, source identity and field co-use are genuine.

## Source-verified published-data calibration, completed

Original dataset from author GitHub [commit 74e3e2c](https://github.com/IRagonese/MilkweedWarming2024/commit/74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2), `MWwarming_comp_May16.csv` blob SHA1 `e7a625911d8dc9b871a52d7657314a5d05472216`, passed exact source-hash validation and **3/3 tests** in [GitHub action 38016318973](https://github.com/zuizui0223/chocho/actions/runs/38016318973).

- Exactly **240 assigned monarch caterpillars**, 120 plants (2 larvae per plant), **30 independent field plots** (8 larvae and four factorial plant/parasite combinations per plot), 15 temperature plots per arm. **Plot**, not larvae, is the independent temperature replicate.
- 239 known adult-eclosion fates, **219 adult emergences**; original injury `22b` on ambient tropical control has `Surv_adult=NA`, not experimental death.
- Descriptive post-hoc host × temperature contrast, collapsing OE conditions: `(tropical-swamp survival in elevated) - (tropical-swamp survival in ambient)` = **−0.1333333** (−13.3 percentage points).
- **9,999 stratified-by-temperature plot bootstrap draws** with seed 20261010: empirical 95% interval **[−0.3000,+0.0333]**, zero within interval.
- If the single accidental injury is **conservatively counted as failed adult eclosion**, the descriptive contrast becomes **−0.1166667**, direction unchanged.
- These are **new audit calculations on an already published experiment**, not a new ecology finding. Ragonese et al. already tested host × temperature adult survival and reported a marginal overall interaction; the plot-cluster interval is a different post-hoc uncertainty specification, not an independent confirmatory replication.
- Adult eclosion is not automatically demonstrated mating, egg laying or population growth. In the original trial, survival is high, leaving limited ability to estimate subtle plant × enemy survival effects with 15 temperature plots.

## Scientific feasibility and publication stop

1. **No validated field population using a common introduced host for both experimental butterfly species** has yet passed the frozen same-site/season larval co-use gate. The specific `J. coenia × A. jatrophae` pair is **not a field-confirmed bridge**; no outcome selected to rescue it.
2. Existing 240-larva monarch data contain only **one butterfly species**, no heterospecific presence/removal factorial, and OE inoculation rather than donor-mediated interspecific infection. They cannot test (M_h).
3. A planned experiment could test (M_h), but is not completed. Obtain field/cohort feasibility and institutional approvals first; do not propose unapproved pathogen propagation or ecological introductions.
4. **GO criterion:** independently demonstrated local co-use by at least two butterfly species and enough plant-species/plot/population replication to estimate the *enemy-mediated neighbor effect* on independently observed **viable adult production**.
5. **NO-GO criterion:** if the only evidence is HOSTS/WCVP geographical overlap, published mono-specific diet effects, or one exotic host species per comparison, a generalized community-level introduction effect is not identified.
6. Main Global Ecology and Biogeography paper PR #38 stays on its frozen, honest **potential-resource** claim. Do not merge exploratory mechanistic figures or p-values.

**Current research state:** plausible, sharply falsifiable experimental question with quantitative replication insight, but **not yet a newly demonstrated ecological mechanism**.
