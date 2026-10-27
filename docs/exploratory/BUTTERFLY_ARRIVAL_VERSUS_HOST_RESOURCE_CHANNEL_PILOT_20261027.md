# Butterfly arrival vectors and relocated larval-host opportunities are separate ecological channels

**2026-10-27 — completed, source-backed posthoc pilot, NOT a verified invasion-arrival model.** [Original-source full CI #40478988726](https://github.com/zuizui0223/chocho/actions/runs/40478988726): **SUCCESS**, 6 no-network scientific controls passed. [Source-only machine-readable result](https://github.com/zuizui0223/chocho/actions/runs/40478988726/artifacts/12722375513). This draft is stacked on exploratory PR #80; frozen main/GEB paper unchanged.

## 1. Correct a false arrival presence/absence outcome before using the source

The GEB resource-validation source `butterfly_occurrence_resource_validation_v0.1.json` has exactly **115 butterfly-present species×WGSRPD3 units** OUTSIDE the butterfly's old native-host resource range: **66** units overlap a newly introduced host WCVP source region, and **49** units remain OUTSIDE the entire mapped contemporary larval host resource envelope. Original butterfly field presence is recorded for **all 115**. Therefore 66/49 is NOT a present/absent comparison and is incapable of estimating arrival, absence or colonization success. An apparent 66 successes vs 49 failures test would be conditioned on the outcome. Both original numeric sources and all species-level sums were explicitly reverified.

## 2. Reconstruct original quantitative resource geography without route assumptions

From original 11-member source-SHA-verified GEB author-only HOSTS/WCVP accepted taxon records:

- **239 original butterflies**, 2,596 original WCVP-accepted butterfly×larval-host pairs;
- **206/239** gained new resource geography after documented plant introductions;
- **14,553** extra butterfly×WGSRPD3 units, one fractional unit per species×new region split between all source accepted hosts capable of creating it;
- **670** original WCVP-accepted plant taxa contributed, top ten provide **25.10%** and top 38 provide **50.57%**.

This replay reproduces the original GEB resource geography and plant concentration *without* fitting an arrival model. Its main value is establishing which real plants must be assigned an independent transport *commodity* variable, not treating host resource credit itself as an insect-arrival predictor.

### Ten original highest resource-contribution plants: independent propagation and nursery modes

| Host plant | Percent of all original added resource units | Source-supported propagation / nursery mode (NONEXCLUSIVE) |
|---|---:|---|
| *Medicago sativa* | 4.55% | sown as forage seed |
| *Poa pratensis* | 3.29% | lawn seed **and** living sod |
| *Oxalis corniculata* | 3.01% | seed reproduction **and** nursery-container weed |
| *Cynodon dactylon* | 2.60% | seeded grass **and** live sod/sprigs |
| *Dactylis glomerata* | 2.25% | seed-established pasture grass |
| *Zea mays* | 2.03% | planted corn seed |
| *Avena sativa* | 2.03% | oat seeds |
| *Urtica urens* | 1.82% | seed-propagating annual weed |
| *Oryza sativa* | 1.80% | direct sowing **and local** transplanting of seedlings |
| *Digitaria sanguinalis* | 1.72% | seed-reproducing nursery-container weed |

Source-supported planting/propagation evidence among the top ten overlaps:

- all ten have evidence of seed planting or seed spread, covering a **25.10%** original resource-credit fraction;
- turf live sod/sprigs demonstrated for *Poa* and *Cynodon*: **5.89%** of source resource credit;
- container-nursery weed association *Oxalis* and *Digitaria*: **4.73%**;
- rice seedlings transplanted locally: **1.80%**.

**These percentages OVERLAP and do not add to a single fraction.** Each mode is a propagation / nursery-occurrence fact, **NOT the form used in historical international plant shipments or an interception of a butterfly egg/larva**. For ALL top ten, no source-linked country-year shipment volume or butterfly-stage interception was verified in this pilot: numeric arrival effects and an 'orphan resource fraction' are intentionally reported **NOT IDENTIFIABLE**, not zero. The remaining **74.90%** of resource credit across other 660 contributing host plants is not yet mode audited; do not extrapolate the ten-host mechanism to it.

Plant horticultural/biological sources include [US Forest Service on alfalfa seed](https://research.fs.usda.gov/feis/species-reviews/medsat), [Illinois extension on Kentucky bluegrass seed+sod](https://web.extension.illinois.edu/hortanswers/plantdetail.cfm?PlantID=20&PlantTypeID=12), [North Carolina Extension on *Oxalis* in nurseries](https://plants.ces.ncsu.edu/plants/oxalis-corniculata/), [Mississippi State Extension on seed and sod/sprigs for *Cynodon*](https://cto.extension.msstate.edu/bermudagrass-cynodon-dactylon), [Penn State orchardgrass seed](https://plantscience.psu.edu/outreach/plant-id/grasses/orchardgrass), [Japanese NARO oats](https://www.naro.go.jp/laboratory/carc/organic/ryokuhi-sakumotsu/enbaku.html), [RHS *Urtica* propagation](https://www.rhs.org.uk/plants/119107/urtica-urens/details), [rice direct sowing and transplant](https://rkb-odisha.in/rice-in-odisha/step-by-step-production/growth/planting/), [NC State *Digitaria* nursery weeds](https://herbarium.ncsu.edu/containerWeeds/Digitaria_sanguinalis.htm). These support categorical **possibilities**, not importer/shipment exposure.

## 3. Evidence of real arrival mechanisms is very different across butterfly examples

This difference makes the hypothesis biologically worth testing, and also prevents treating all ornamental-host invasions as independently confirmed live-plant shipments.

| Historical butterfly example | Source-supported arrival route / inference | Independently recorded postarrival resource | Interpretation |
|---|---|---|---|
| *Cacyreus marshalli*, UK 1978 | **Two larvae reported intercepted on imported South African live *Pelargonium***, historical primary study discussion; infestation eradicated in that incident, no establishment demonstrated **from that lot** | Ornamental *Pelargonium* in Europe, established in Spanish islands later | Arrival via living plant can be observed and still NOT itself imply establishment; introduced route and plant resource can be tightly coupled |
| *C. marshalli*, Balearics 1989 | Horticultural plant imports considered likely source; no original 1989 shipment lot linked to resident population by a causal trace | Ornamental *Pelargonium* confirmed as larval host in European reports | Route likely, precise founding vector remains unverified |
| *Thymelicus lineola*, Ontario by 1910 | **Eggs on hay used to pack imported goods** hypothesized in Canadian 1968 historical entomological account; initial adults found near refuse; not a confirmed egg interception | Introduced pasture/timothy and grasses, direct later field larval reports | Original host resources can arise via grass **seed**, while an insect could arrive on a different form of plant **product** such as hay rather than growing grass or seed |
| *Papilio demoleus*, Dominican Republic 2004 | Genetics traces origin to SE Asia; original genetic study explicitly says **arrival route unknown**. EPPO lists citrus consignment/early-stage transport as a **hypothesis** | Citrus and other plants support larvae in recipient countries | Genetic geographic origin must never be misread as observed commodity pathway |
| *Chilades pandava*, Guam 2005 | Route of specific Guam founder butterfly not directly verified by the cycad nursery study; Philippine nursery population likely introduced on plants is a distinct regional observation | Guam-native *Cycas micronesica* confirmed with wild-origin larvae reared to pupation in 2014 | Postarrival resource can be a **recipient-native plant**, distinct from hypothetical horticultural carrier |

Sources: [Sarto i Monteys historical account of UK *Cacyreus* interception, 1992](https://www.researchgate.net/publication/273380945_Spread_of_the_Southern_African_Lycaenid_butterfly_Cacyreus_marshalli_Butler_1898_LEP_Lycaenidae_in_the_Balearic_Archipelago_Spain_and_considerations_on_its_likely_introduction_to_continental_Europe), [2022 *Biological Invasions* repeat plant-import egg evidence](https://link.springer.com/article/10.1007/s10530-022-02895-3), [Canadian 1968 historical transport hypothesis](https://publications.gc.ca/collections/collection_2021/aac-aafc/A49-2/A49-2-1968-1.pdf), [2006 *P. demoleus* genetic provenance](https://farrell.oeb.harvard.edu/publications/provenance-old-world-swallowtail-butterflies-papilio-demoleus-lepidoptera), [EPPO speculation on citrus imports](https://gd.eppo.int/reporting/article-1282), [Guam *Chilades* indigenous host use](https://pmc.ncbi.nlm.nih.gov/articles/PMC5022404/).

The often-cited European review covered **97 non-native Lepidoptera, including MANY MOTHS**, not 97 strictly butterflies, and found almost all introduced by plants or PLANT PRODUCTS, not necessarily living potted host material. **52.5% in gardens/parks** is a postestablishment habitat association, NOT a proportion moved internationally on live plants. [Lopez-Vaamonde et al., BioRisk](https://doi.org/10.3897/biorisk.4.50); [Augustin et al., EPPO Bulletin](https://doi.org/10.1111/epp.2600). Hence these percentages cannot supply the missing live-plant shipment exposure denominator.

## Additional source-ontology limitation: wild WCVP distributions omit purely cultivated hosts

The original GEB WCVP resource envelope includes accepted host **native, extant and introduced wild** botanical distribution units, not every planted garden/crop/nursery plant record. An insect carried in a living **horticultural consignment** may find a cultivated host in a city or nursery whose population is NOT naturalized, and therefore is *missing* from a WCVP-introduced wild resource map. Conversely, an introduced wild grass in the WCVP map may supply abundant larval leaves despite an arrival via adult dispersal, hay packing or another unrelated commodity. Thus using WCVP wild records as a complete measure of botanical arrival vectors or urban/crop establishment resources would be a second source-ontology error. A verified study must preserve separately **trade of live planting material**, **planted-only (non-wild) host landscapes**, and **naturalized wild host footprints**. Neither the 25.10% seed propagation coverage nor the original 54.9% opportunity increase is a direct cross-border trade fraction.

## 4. Next identifiable macro experiment

A valid species×destination×year model must distinguish:

`potential transported insect propagules` (source-region butterfly occurrence, verified live-plant and plant-product consignment route/volume and insect-stage interception, adult natural movement, season);

`potential postarrival host resources` (the original WCVP accepted host repertoire's distribution/quantity and climate); and

`detected establishment` (repeated, preferably wild immature records, detection effort, not a single adult GBIF sighting).

The most telling comparison is neither 'pasture versus ornamental' nor the invalid 66 vs 49 cells. It is whether **same resource geography supports different establishment outcomes conditional on the arrival channel and survey exposure**, or whether an observed propagule shipment fails to generate a persistent population even in a host-rich recipient. A seed-borne host can be a large establishment opportunity yet travel without larvae; the insect could still arrive independently on hay, other live plants or by flight. The outcome is empirically open.

**Status:** original resource geography reconstructed and documented channels sampled, but no global species×country×date plant consignment/interception dataset joined. No estimated channel coefficient, orphan opportunity percentage or direct suppression/fitness result; those would currently be inventions.

## Reproducibility

- [Frozen two-channel protocol](../../provenance/reviewer_defenses/protocols/ARRIVAL_RESOURCE_CHANNEL_DECOUPLING_20261027.md)
- [Script reproducing the exact 14,553 species×region credits and correcting the observed115 source error](../../scripts/audit_butterfly_arrival_vs_host_resource_routes_20261027.py)
- [Ten exact top-ranked original host taxa with source practice URLs and UNVERIFIED foreign shipment/insect-stage fields](../../data/exploratory/arrival_resource_top10_propagation_route_pilot_20261027.csv)
- [Six synthetic checks and successful original data replay](https://github.com/zuizui0223/chocho/actions/runs/40478988726), [archive of source summary](https://github.com/zuizui0223/chocho/actions/runs/40478988726/artifacts/12722375513)

No raw author-only HOSTS accepted butterfly-host interaction rows were exported. GEB main and existing PR #75–#80 remain unchanged.
