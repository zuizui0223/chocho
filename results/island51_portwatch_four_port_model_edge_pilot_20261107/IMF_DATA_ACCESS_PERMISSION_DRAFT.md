# IMF PortWatch connection-network data access: proposed academic permission clarification

**Prepared 2026-11-07. Not sent.** This is a suggested request for the data copyright/contact team, requiring the researcher's personal review and approval before external transmission. It does not authorize any new automated bulk download and should not be interpreted as a license already granted.

**Proposed recipient:** copyright@imf.org  
**Subject:** Permission clarification — noncommercial ecological research using IMF PortWatch Maritime Connections and port registry

Dear IMF Copyright and Data Team,

I am conducting noncommercial academic macroecological research on how the historical availability of larval host plants, island geography and marine transport opportunities relate to indexed butterfly community similarity. I would like to clarify the applicable IMF PortWatch permissions before accessing the following publicly hosted ArcGIS FeatureServer data beyond small examples:

1. **PortWatch Spillover Simulator Maritime Connections**: https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/Spillover_Simulator_Maritime_Connections/FeatureServer/0 (the public server currently reports 155,521 connection records, with from/to PortWatch port IDs, positions, modeled distance/transit days and simulator-input columns).
2. **PortWatch port registry**: https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0 (the public source currently reports 2,065 port records, including port name, country/ISO3, location and PortWatch port ID).

We have already performed a strictly bounded proof-of-concept using four independently verified port identities and 12 directed connection conditions; **we have not downloaded the full global connection or registry tables**. We cite IMF PortWatch and do not treat a simulator edge as proof of an observed voyage or biological transportation.

For research reproducibility, we would like to (a) obtain a complete, source-dated snapshot of the above **publicly available** port registry and directed connection table, ideally as an official authorized bulk file or via rate-limited API queries; (b) cross-reference the registry with our original 194 World Port Index port entries, without distributing any sensitive vessel identifiers; (c) derive and publish limited **aggregated** island-pair connectivity statistics, their model-input provenance, source dates, code and reproducible comparison methodology; and (d) retain a local raw source snapshot strictly for research verification, without republishing IMF's full raw dataset.

Could you please clarify:

- Whether these two particular PortWatch layers fall under the IMF Copyright and Usage **special terms for published statistical Data** and, in particular, whether the prohibition on automated bulk downloading applies to their academic use.
- Whether written permission is required for one bounded/rate-limited full snapshot of approximately 155,521 connection records and 2,065 port-registry rows, and, if so, whether permission may be granted for the research and retention described above.
- The appropriate required credit, dataset/version citation, restrictions on redistribution of derived aggregates and source snapshots, and whether an officially prepared bulk CSV/Parquet is available.
- The construction/reference period for the static maritime connection network, its relation to AIS vessel-visit evidence, and whether the records indicate observed movements, inferred possible links, or model-specific network capacity. We will not label the source as individual observed voyages without documentation.

The planned research is observational and noncommercial. We will distinguish source connectivity from actual ship-carried plants, butterfly movement and colonization, and will explicitly state data and temporal limitations. We can provide a more detailed project description or institutional affiliation if needed.

Thank you for your guidance.

Sincerely,  
[Researcher name and academic affiliation to be supplied by the sender]

---

## Research-side constraints before any permission decision

- This file is only an **unsent draft**, not consent for external communication, not permission from IMF, and not a user credential.
- The current [source-pinned four-port test](README.md) establishes WPI-to-PortWatch identifier compatibility only for the deliberately enriched four sampled ports, not for the 194-port population.
- Even if permission is granted, first freeze a port identity and source cohort protocol; validate all candidate original WPI matches by original name, ISO3 and geodesic coordinates, and retain no-match/ambiguous statuses. Then report source island-pair coverage before any butterfly model.
- A static model edge is not an actual voyage, nor a verified 2001–2010 shipping history. The earlier global GEB larval-host resource results and island/port conditional AIS prediction findings are distinct and unchanged.
