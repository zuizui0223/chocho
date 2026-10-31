#!/usr/bin/env python3
"""Audit historical 5-km-compatible host indexing against new fixed 1-km botanical sites.

The original historical source returned specific compatible host speciesKey
within 5km of a 2019-WPI source port. This recheck asks whether such an
association remains indexed at the port's own 1km circle or matched nonport
coastal circles at the same years. Non-detection is never inferred absence.
Plant taxonomy is accepted GBIF species only, not local-alien or larval feeding.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCES={
    "focal":(
        "provenance/reviewer_defenses/results/dated_focal_host_taxonomy_20261031/focal_host_prior_to_butterfly_20261031.json",
        "94a82adf74eff89ee073fb2016613dda2aee7a26"
    ),
    "sites":(
        "provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_plant_site_facets.csv",
        "eb09a3f3ec053c00b43e8bc0678d6643ce02697b"
    ),
    "cohort":(
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_links.csv",
        "88db773c74dea9cd2d58bde0b67a394d93226869"
    )
}
def load(name):
    path,known=SOURCES[name]
    raw=(ROOT/path).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if actual!=known:raise ValueError("Source public Git blob changed "+name+" "+actual)
    if path.endswith(".json"):return json.loads(raw)
    return list(csv.DictReader(raw.decode("utf8").splitlines()))

def analyze():
    focal=load("focal")
    sites=load("sites")
    links=load("cohort")
    assert focal["butterfly_links_checked"]==3
    assert focal["island_taxon_endpoint_types"]==4
    assert len(sites)==114
    assert len(links)==5
    assert sum(x["query_status"]!="COMPLETE" for x in sites)==0
    byisland=defaultdict(list)
    for site in sites:byisland[str(site["landmass_id"])].append(site)
    transitions=[]
    grouped={}
    for rel in focal["focal_links"]:
        species=rel["butterfly"]
        host=rel["source_catalogue_compatible_host"]
        plant_key=str(rel["gbif_accepted_plant_speciesKey"])
        assert (species,host) in {(x["species"],x["host_name"]) for x in links}
        for endpoint in rel["island_record_comparisons"]:
            island=str(endpoint["island"])
            assert 2001<=endpoint["host_first_indexed_plant_year"]<=2010
            assert endpoint["host_first_indexed_plant_year"] < endpoint["butterfly_first_indexed_year"]
            assert endpoint["lead_years_of_host_index_before_butterfly_index"]==endpoint["butterfly_first_indexed_year"]-endpoint["host_first_indexed_plant_year"]
            if island not in byisland:raise ValueError("No 1km GBIF controls for focal island")
            port=[s for s in byisland[island] if s["site_type"]=="port"]
            off=[s for s in byisland[island] if s["site_type"]=="offport"]
            assert len(off)==2*len(port)
            def source_match(arr):
                return [s["site_id"] for s in arr if plant_key in s["species_keys"].split(";")]
            result={
                "butterfly":species,
                "host":host,
                "GBIF_accepted_plant_speciesKey":plant_key,
                "island":island,
                "indexed_5km_WPI_source_host_2001_2010":True,
                "original_host_first_indexed_year":endpoint["host_first_indexed_plant_year"],
                "butterfly_first_indexed_year":endpoint["butterfly_first_indexed_year"],
                "host_index_lead_years":endpoint["lead_years_of_host_index_before_butterfly_index"],
                "port_1km_sites":len(port),
                "port_1km_indexed_speciesKey_sites":len(source_match(port)),
                "port_1km_matched_site_ids":source_match(port),
                "offport_1km_sites":len(off),
                "offport_1km_indexed_speciesKey_sites":len(source_match(off)),
                "offport_1km_matched_site_ids":source_match(off),
                "one_km_nonindex_NOT_biological_absence":True,
            }
            transitions.append(result)
            grouped[(island,plant_key)]=result
    assert len(transitions)==6
    assert len(grouped)==4
    assert len(set(x["butterfly"] for x in transitions))==3
    assert len(set(x["host"] for x in transitions))==2
    summary=[{
        "island":i,"host":val["host"],"speciesKey":k,
        "5km_prior_host_source_indexed":True,
        "1km_port_site_count":val["port_1km_sites"],
        "1km_port_site_indexed_count":val["port_1km_indexed_speciesKey_sites"],
        "1km_offport_site_count":val["offport_1km_sites"],
        "1km_offport_site_indexed_count":val["offport_1km_indexed_speciesKey_sites"],
    } for (i,k),val in sorted(grouped.items(),key=lambda x:(int(x[0][0]),int(x[0][1])))]
    assert sum(q["1km_port_site_indexed_count"] for q in summary)==0
    assert sum(q["1km_offport_site_indexed_count"] for q in summary)==2
    return {
        "schema":"chocho_focal_hosts_earlier_5km_vs_1km_harbor_scale_20261031_v1",
        "status":"TEMPORALLY_PRIOR_5KM_HOST_INDEX_DOES_NOT_ESTABLISH_PORT_LOCAL_1KM_HOST",
        "source_git_blobs":{k:{"path":path,"sha":sha} for k,(path,sha) in SOURCES.items()},
        "study_window_plant_index":"2001-2010",
        "subsequent_butterfly_index_range":[2015,2021],
        "5km_source": "previously verified independent WPI-centered botanical speciesKey facet",
        "1km_source":"fixed 2001-2010 matched port vs offport 1km botanical speciesKey facets",
        "unique_host_by_island_endpoints":len(summary),
        "butterfly_species_island_endpoints_counting_shared_hosts":len(transitions),
        "named_focal_butterfly_species":3,
        "distinct_candidate_hosts":2,
        "site_index_summary":summary,
        "butterfly_island_endpoints":transitions,
        "limits":[
            "The 2001-2010 5km and 1km queries share GBIF indexing; they are not independent biological arrival dates.",
            "SpeciesKey absence from 1km GBIF index does not imply host plant absent.",
            "Old 5km buffers are centered at 2019 WPI positions and do not measure historical vessel cargo.",
            "The nonport controls are sampled from 2019 WPI-excluded coasts; not necessarily pristine or unplanted.",
            "Alnus incana offport observations do not establish native versus historically introduced status.",
            "The same Medicago host islands are used by Colias croceus and Colias hyale; the six endpoints are NOT independent replicates.",
            "This is a selected focal subset, not a random 225-butterfly host-port survey.",
            "No observed larval consumption, establishment, independent island-to-island flight, or causal stepping-stone effect.",
        ]
    }

def main():
    obj=analyze()
    (HERE/"focal_host_5km_to_1km_source_audit.json").write_text(
        json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":obj["status"],"island_host_support":obj["site_index_summary"]},
        indent=2,ensure_ascii=False))

if __name__=="__main__":main()
