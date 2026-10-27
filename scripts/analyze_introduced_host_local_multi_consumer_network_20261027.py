#!/usr/bin/env python3
"""Geographically reported shared butterfly larval hosts, original 239 taxa.

Fixed source: original HOSTS Location and WCVP accepted plant ranges.
A plant×HOSTS Location with >=2 distinct butterfly consumer reports is an
INDEXED CO-USE NODE, NOT proof of co-occurring larvae, wild breeding,
resource competition, ecological rewiring or a consumer invasion.

Artificial host-origin exchanges only where native and introduced node
labels exist within identical Location × plant-family strata; controls
original source raw HOSTS report effort, WCVP-unlisted status and source
country-area threshold without publishing original interaction pairs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    SOURCE_ZIP_SHA, EXPECTED_HOSTS_ROWS, EXPECTED_SPECIES_HOST_ROWS,
    accepted_pairs, classify_location, load_country_atoms, read_rows, sha256,
    wcvp_range,
)

SCHEMA = "chocho_239_geographically_documented_butterfly_local_shared_hosts_20261027_v01"
STATUSES = ("NATIVE_CONFIRMED", "INTRODUCED_CONFIRMED", "WCVP_UNLISTED",
            "MIXED_OR_PARTIAL", "COUNTRY_MAPPING_UNSAFE",
            "SUBREGION_MAPPING_UNSAFE", "UNMAPPED_OR_UNREVIEWED_LOCATION")
PRIMARY_CUTOFF = ".95"
CUTS = ((.95, ".95"),(.80, ".80"))
REPORTING_BINS = ("1", "2", "3-5", "6+")
MIN_PERMUTABLE_LOCALITIES = 10
MIN_PERMUTABLE_STRATA = 15
MIN_INTRODUCED_MATCHED_NODES = 20
MIN_TOTAL_MATCHED_SHARED_NODES = 10


def row_bin(n):
    if not n > 0:
        raise ValueError("No indexed accepted HOSTS source rows for observed plant")
    return "1" if n == 1 else "2" if n == 2 else "3-5" if n <= 5 else "6+"


def hypergeom(N, K, m):
    N,K,m = int(N),int(K),int(m)
    if not 0 < N or not (0<=K<=N) or not (0<=m<=N):
        raise ValueError("Invalid source node hypergeometric margins")
    start = max(0,m-(N-K))
    end=min(m,K)
    denom=math.comb(N,m)
    pmf={k:math.comb(K,k)*math.comb(N-K,m-k)/denom
         for k in range(start,end+1)}
    if not math.isclose(sum(pmf.values()),1.,rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError("Source node hypergeom normalization failed")
    return pmf


def exact_fixed_margin_source_reference(strata, target):
    """Null labels assigned to observed plant nodes; not ecological randomization."""
    pmf={0:1.}
    mean=0.
    var=0.
    source_locations=set()
    permutable_source_locations=set()
    informative=0
    participating=0
    count=Counter()
    across_labels=Counter()
    for cell in strata:
        n=cell["N"]
        m=cell["target"]
        k=cell["shared"]
        observed=cell["target_shared"]
        if not 0<=m<=n or not 0<=k<=n or not max(0,m+k-n)<=observed<=min(m,k):
            raise ValueError("Invalid observed source plant-node shared status margins")
        # Matched stratum requires BOTH native and one target origin source nodes.
        if not 0<m<n:
            continue
        participating+=1
        loc=cell["location"]
        source_locations.add(loc)
        count["source_nodes"]+=n
        count["target_nodes"]+=m
        count["native_nodes"]+=n-m
        count["shared_all"]+=k
        count["target_shared"]+=observed
        count["native_shared"]+=k-observed
        across_labels[cell["family"]]+=m
        distribution=hypergeom(n,k,m)
        mean += m*k/n
        var += m*(k/n)*(1-k/n)*((n-m)/(n-1)) if n>1 else 0.
        if len(distribution)>1:
            informative+=1
            permutable_source_locations.add(loc)
        next_ = defaultdict(float)
        for i,p in pmf.items():
            for v,q in distribution.items():
                next_[i+v]+=p*q
        pmf=dict(next_)
    if not math.isclose(sum(pmf.values()),1.,abs_tol=1e-9):
        raise ValueError("Convolved source probabilities lost mass")
    possible=sorted(pmf)
    obs=count["target_shared"]
    def pct(q):
        c=0.
        for k in possible:
            c+=pmf[k]
            if c+1e-12>=q:
                return k
        return possible[-1]
    ok=(len(permutable_source_locations)>=MIN_PERMUTABLE_LOCALITIES and
        informative>=MIN_PERMUTABLE_STRATA and
        count["target_nodes"]>=MIN_INTRODUCED_MATCHED_NODES and
        count["shared_all"]>=MIN_TOTAL_MATCHED_SHARED_NODES)
    return {
        "reference_target_status":target,
        "source_matched_location_family_strata":participating,
        "source_permutable_location_family_strata":informative,
        "source_matched_location_count":len(source_locations),
        "source_permutable_location_count":len(permutable_source_locations),
        "source_matched_plant_nodes_all":count["source_nodes"],
        "source_matched_target_plant_nodes":count["target_nodes"],
        "source_matched_native_plant_nodes":count["native_nodes"],
        "source_shared_plant_nodes_in_matched_strata":count["shared_all"],
        "source_observed_target_shared_nodes":obs,
        "source_observed_native_shared_nodes":count["native_shared"],
        "source_target_shared_node_fraction_in_matched_strata":
            obs/count["target_nodes"] if count["target_nodes"] else None,
        "source_native_shared_node_fraction_in_matched_strata":
            count["native_shared"]/count["native_nodes"] if count["native_nodes"] else None,
        "source_exact_conditional_expected_target_shared_nodes":mean,
        "source_exact_conditional_variance":var,
        "source_exact_conditional_sd":math.sqrt(max(0.,var)),
        "source_central_95pct_reference_interval":[pct(.025),pct(.975)],
        "source_conditional_left_tail_not_ecological_p":
            sum(p for k,p in pmf.items() if k<=obs),
        "source_conditional_right_tail_not_ecological_p":
            sum(p for k,p in pmf.items() if k>=obs),
        "source_conditional_permutable_support":[possible[0],possible[-1]],
        "information_gate":(
            "SOURCE_REFERENCE_WITH_EXCHANGEABLE_LABELS" if ok
            else "INSUFFICIENT_MATCHED_ORIGIN_GEOGRAPHY_INFORMATION"
        ),
        "red_line":(
            "The source plant-native/introduced label is NOT randomly assigned. "
            "A mathematical conditional tail is not an ecological p value, "
            "fitness effect, butterfly species competition, or wild sympatry."
        ),
    }


def source_nodes(hosts_path, accepted, family_by_plant, focal):
    indexed=defaultdict(lambda:{"butterflies":set(),"raw_rows":0})
    totals=Counter()
    for r in read_rows(hosts_path):
        totals["original_all"]+=1
        if not all(str(r.get(k) or "").strip() for k in
                   ("Insect Genus","Insect Species",
                    "Hostplant Genus","Hostplant Species")):
            continue
        totals["complete_binomials"]+=1
        butterfly=r["Insect Genus"].strip()+" "+r["Insect Species"].strip()
        if butterfly not in focal:
            continue
        totals["focal_raw"]+=1
        host=r["Hostplant Genus"].strip()+" "+r["Hostplant Species"].strip()
        pid=accepted.get((butterfly,host))
        if pid is None:
            totals["focal_unmatched_original_WCVP_pair"]+=1
            continue
        totals["focal_accepted_original_raw"]+=1
        loc=(r.get("Location") or "").strip()
        entry=indexed[pid,loc]
        entry["butterflies"].add(butterfly)
        entry["raw_rows"]+=1
    if (totals["original_all"],totals["complete_binomials"],
        totals["focal_raw"],totals["focal_accepted_original_raw"])!=(
        EXPECTED_HOSTS_ROWS,EXPECTED_SPECIES_HOST_ROWS,4695,3898):
        raise ValueError("Original source HOSTS row numbers have changed")
    if sum(len(v["butterflies"]) for v in indexed.values()) != 3586:
        raise ValueError("Original accepted source HOSTS plant-node edges changed")
    result=[]
    for (pid,loc),data in sorted(indexed.items()):
        fam=family_by_plant.get(pid,())
        if len(fam)!=1 or not fam[0]:
            raise ValueError("Original botanical accepted plant family ambiguous")
        genus={sp.split()[0] for sp in data["butterflies"]}
        result.append({
            "plant_id":pid,
            "original_Location":loc,
            "plant_family":fam[0],
            "raw_original_HOSTS_rows":data["raw_rows"],
            "report_intensity_bin":row_bin(data["raw_rows"]),
            "distinct_original_239_butterfly_consumers":len(data["butterflies"]),
            "documented_at_least_two_butterfly_species":len(data["butterflies"])>=2,
            "documented_at_least_two_butterfly_genera":len(genus)>=2,
        })
    return result,dict(totals)


def classify_nodes(nodes, cutoff, atoms, atom_meta, regions, native, current, *, require_original_source=True):
    source_qc=Counter()
    cache={}
    out=[]
    for node in nodes:
        key=(node["plant_id"],node["original_Location"])
        if key not in cache:
            cache[key]=classify_location(
                key[0],key[1],threshold=cutoff,
                atoms=atoms,atom_metadata=atom_meta,original_regions=regions,
                native=native,current=current)
        status=cache[key]
        if status not in STATUSES:
            raise ValueError("Unregistered WCVP original plant source status")
        source_qc[status]+=node["distinct_original_239_butterfly_consumers"]
        out.append({**node,"source_plant_region_origin":status})
    if require_original_source and cutoff==.95 and any(source_qc[k]!=v for k,v in
                           (("NATIVE_CONFIRMED",399),
                            ("INTRODUCED_CONFIRMED",179),
                            ("WCVP_UNLISTED",303))):
        raise ValueError("Original WCVP accepted HOSTS source geography changed")
    if require_original_source and sum(source_qc.values())!=3586:
        raise ValueError("Source butterfly-host×Location edges not conserved")
    return out,dict(source_qc)


def summarise_status(nodes):
    sums=defaultdict(Counter)
    locations=defaultdict(set)
    families=defaultdict(set)
    for v in nodes:
        k=v["source_plant_region_origin"]
        s=sums[k]
        s["accepted_plant_Location_nodes"]+=1
        s["documented_butterfly_Location_host_edges"]+=v["distinct_original_239_butterfly_consumers"]
        s["original_raw_HOSTS_rows"]+=v["raw_original_HOSTS_rows"]
        s["shared_plant_nodes_2plus_butterflies"]+=v["documented_at_least_two_butterfly_species"]
        s["cross_genus_shared_plant_nodes"]+=v["documented_at_least_two_butterfly_genera"]
        if v["documented_at_least_two_butterfly_species"]:
            s["butterfly_host_edges_on_shared_nodes"]+=v["distinct_original_239_butterfly_consumers"]
        locations[k].add(v["original_Location"])
        families[k].add(v["plant_family"])
    return {
        k:{
            **dict(val),
            "original_report_locations":len(locations[k]),
            "distinct_plant_families":len(families[k]),
            "source_fraction_plant_nodes_shared_by_two_butterflies":
                val["shared_plant_nodes_2plus_butterflies"]/
                val["accepted_plant_Location_nodes"],
            "source_fraction_edges_on_shared_plant_nodes":
                val["butterfly_host_edges_on_shared_nodes"]/
                val["documented_butterfly_Location_host_edges"],
        } for k,val in sorted(sums.items())
    }


def make_strata(nodes,target,*,report_bin=False):
    """Only two explicitly distinguished WCVP statuses enter each matched source null."""
    if target not in ("INTRODUCED_CONFIRMED","WCVP_UNLISTED"):
        raise ValueError("Unknown target source origin")
    bykey=defaultdict(list)
    for node in nodes:
        if node["source_plant_region_origin"] not in (
                target,"NATIVE_CONFIRMED"):
            continue
        pair=(node["original_Location"],node["plant_family"])
        key=(*pair,node["report_intensity_bin"]) if report_bin else pair
        bykey[key].append(node)
    out=[]
    for key,rows in sorted(bykey.items()):
        n=len(rows)
        target_units=sum(r["source_plant_region_origin"]==target for r in rows)
        shared=sum(r["documented_at_least_two_butterfly_species"] for r in rows)
        both=sum(r["source_plant_region_origin"]==target and
                 r["documented_at_least_two_butterfly_species"] for r in rows)
        out.append({
            "location":key[0],"family":key[1],
            "N":n, "target":target_units,"shared":shared,"target_shared":both,
        })
    return out


def source_locality_concentration(nodes,status):
    samples=defaultdict(Counter)
    for v in nodes:
        if v["source_plant_region_origin"]!=status:
            continue
        loc=v["original_Location"]
        samples[loc]["documented_host_plant_nodes"]+=1
        samples[loc]["source_shared_plant_nodes"]+=bool(
            v["documented_at_least_two_butterfly_species"])
        samples[loc]["source_butterfly_host_edges"]+=v[
            "distinct_original_239_butterfly_consumers"]
    # Source names are localities/regions already contained in HOSTS; no host link identities.
    return [
        {"source_location_name":name,**dict(q)}
        for name,q in sorted(samples.items(),
                            key=lambda p:-p[1]["source_butterfly_host_edges"])[:12]
    ]


def one_cutoff(nodes, *, cutoff, atoms, atom_meta, regions, native, current):
    rows,edgecounts=classify_nodes(nodes,cutoff,atoms,atom_meta,regions,native,current)
    summary=summarise_status(rows)
    refs={}
    for target in ("INTRODUCED_CONFIRMED","WCVP_UNLISTED"):
        x={}
        for label,reportbin in (
            ("same_Location_and_plant_family",False),
            ("same_Location_and_plant_family_and_HOSTS_raw_report_bin",True)
        ):
            strata=make_strata(rows,target,report_bin=reportbin)
            x[label]=exact_fixed_margin_source_reference(strata,target)
        refs[target]=x
    return {
        "source_original_butterfly_host_Location_edge_counts_by_status":edgecounts,
        "plant_Location_source_nodes_by_origin":summary,
        "source_origin_reference":refs,
        "top_original_indexed_introduced_host_locations":
            source_locality_concentration(rows,"INTRODUCED_CONFIRMED"),
        "top_original_indexed_unlisted_host_locations":
            source_locality_concentration(rows,"WCVP_UNLISTED"),
        "interpretation":"Observed overlap of butterfly host identities at a literature-reported HOSTS Location, not direct site/time co-use or interspecific competition",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",type=Path,required=True)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--country-atoms",type=Path,required=True)
    p.add_argument("--country-atoms-manifest",type=Path,required=True)
    p.add_argument("--country-summary",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    root=args.author_source
    original=json.loads((root/"SOURCE_MANIFEST.json").read_text())
    if original.get("inputs_count")!=11:
        raise ValueError("Not exact original author-only HOSTS/WCVP input")
    for k,v in original["sha256"].items():
        if sha256(root/k)!=v:
            raise ValueError("Original source member SHA changed")
    folder=root/"inputs"
    species_table=read_rows(folder/"anthropogenic_species_metrics.csv")
    focal={r["species"] for r in species_table}
    if len(species_table)!=239 or len(focal)!=239:
        raise ValueError("Original 239 butterfly cohort changed")
    accepted,joinQC,families=accepted_pairs(
        folder/"original/insect_host_accepted.csv",focal)
    if joinQC["source_ambiguous_accepted_id_names"]:
        raise ValueError("Original accepted WCVP source pair ambiguity")
    native,current=wcvp_range(
        folder/"original/native_extant_nondoubtful_wgsrpd3.csv",
        folder/"original/extant_nondoubtful_wgsrpd3.csv")
    atoms,atom_meta,regions,geo=load_country_atoms(
        args.country_atoms,args.country_atoms_manifest,args.country_summary)
    nodes,rawstats=source_nodes(args.hosts_resource,accepted,families,focal)
    result={}
    for cutoff,label in CUTS:
        result[label]=one_cutoff(nodes,cutoff=cutoff,atoms=atoms,
                                atom_meta=atom_meta,regions=regions,
                                native=native,current=current)
    report={
        "schema":SCHEMA,
        "status":"COMPLETE_ORIGINAL_HOSTS_SOURCE_CO_REPORTED_CONSUMERS_ONLY",
        "source":{
            "original_butterfly_species":239,
            "original_total_unique_accepted_butterfly_host_Location_edges":3586,
            "original_WCVP_author_only_source_zip_sha256":
                SOURCE_ZIP_SHA,
            "original_HOSTS_commit":
                "globalbioticinteractions/HOSTS@808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "original_HOSTS_csv_sha256":sha256(args.hosts_resource),
            "source_original_accepted_host_pair_sha256":
                sha256(folder/"original/insect_host_accepted.csv"),
            "original_geography_source":geo,
            "original_source_accepted_join_qc":joinQC,
            "source_unfiltered_raw_rows":rawstats,
            "accepted_plant_Location_nodes_before_WCVP_origin":len(nodes),
        },
        "original_WCVP_country_dominance":result,
        "no_independent_field_competition_claim":True,
        "limits":[
            "HOSTS Location is a publication geography, often a whole country, and co-reporting is not co-occurrence of butterflies in the same habitat, time or plant individual.",
            "A shared plant node is not resource competition; prior Japanese butterfly research supports resource filtering rather than displacement as one plausible mechanism.",
            "WCVP host origin is a regional source distribution, not butterfly indigenous/alien population status or dated host introduction.",
            "Only 239 preselected butterflies; report/country distribution is highly source- and literature-biased, with broad biogeographic labels unclassifiable.",
            "Fixed-label hypergeom exchanges are artificial source controls: the distribution of native and introduced plants is not exchangeable in ecology.",
            "Conditioning on original source HOSTS raw record count may adjust part of ascertainment but cannot recover unobserved local host use.",
            "Original source butterfly×accepted plant×Location raw edges and copyrighted source records were not exported.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    primary=result[".95"]
    summary=primary["plant_Location_source_nodes_by_origin"]
    print(json.dumps({
        "status":report["status"],
        "original239":239,
        "source_accepted_butterfly_host_location_edges":3586,
        "source_total_accepted_plant_location_nodes":len(nodes),
        "native_origin":summary.get("NATIVE_CONFIRMED"),
        "introduced_origin":summary.get("INTRODUCED_CONFIRMED"),
        "unlisted_origin":summary.get("WCVP_UNLISTED"),
        "matched_origin":primary["source_origin_reference"]["INTRODUCED_CONFIRMED"],
        "unlisted_control":primary["source_origin_reference"]["WCVP_UNLISTED"],
        "area80_introduced_origin":result[".80"]["plant_Location_source_nodes_by_origin"].get("INTRODUCED_CONFIRMED"),
    },ensure_ascii=False,sort_keys=True))


if __name__ == "__main__":
    main()
