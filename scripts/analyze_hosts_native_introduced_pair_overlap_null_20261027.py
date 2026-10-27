#!/usr/bin/env python3
"""Frozen conditional native/introduced HOSTS documentation overlap audit.

Question: does the very low number of same accepted butterfly-host pairs
reported in both native and introduced WCVP plant areas exceed or fall below
a simple *source-reporting* exchangeability reference? This is NOT a
biological randomized null for genuine host adoption/niche expansion.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    EXPECTED_HOSTS_ROWS, EXPECTED_SPECIES_HOST_ROWS,
    SOURCE_ZIP_SHA, accepted_pairs, classify_location,
    load_country_atoms, read_rows, sha256, wcvp_range, observations,
)

SCHEMA = "chocho_239_hosts_location_pair_overlap_conditional_null_20261027_v01"
STATUS = ("NATIVE_CONFIRMED", "INTRODUCED_CONFIRMED", "WCVP_UNLISTED")
Q_SPECS = ((.95, "0.95"), (.80, "0.8"))
MIN_INFORMATIVE_BUTTERFLIES = 20
MIN_EXPECTED_OVERLAP = 5.0


def exact_hypergeom(universe, marked, draws):
    """P[overlap=k] for uniformly chosen native-marked source plant IDs."""
    u, m, d = int(universe), int(marked), int(draws)
    if not 0 <= m <= u or not 0 <= d <= u or u <= 0:
        raise ValueError("Impossible original accepted-host source margins")
    lo, hi = max(0, d - u + m), min(m, d)
    den = math.comb(u, d)
    dist = {
        k: math.comb(m,k)*math.comb(u-m,d-k)/den
        for k in range(lo,hi+1)
    }
    if abs(sum(dist.values())-1)>1e-10:
        raise ValueError("Source-host hypergeometric null not normalized")
    return dist


def exact_overlap_reference(parts, observed):
    """Exact DP convolution of source-matched species/family hypergeom margins."""
    pmf={0:1.0}
    mean=0.
    var=0.
    n_informative=0
    species_informative=set()
    for row in parts:
        u, native, target=int(row["U"]),int(row["N"]),int(row["T"])
        overlap=int(row["O"])
        if not max(0,native+target-u) <= overlap <= min(native,target):
            raise ValueError("Source overlap outside possible fixed margin")
        distribution=exact_hypergeom(u,target,native)
        mean += native*target/u
        var += (
            native * (target/u) * (1-target/u) * ((u-native)/(u-1))
            if u>1 else 0.
        )
        if len(distribution)>1:
            n_informative+=1
            species_informative.add(row["species"])
        updated=defaultdict(float)
        for existing,p in pmf.items():
            for k,q in distribution.items():
                updated[existing+k] += p*q
        pmf=dict(updated)
    total_mass=sum(pmf.values())
    if not math.isclose(total_mass,1.0,rel_tol=1e-8,abs_tol=1e-8):
        raise ValueError("Exact source-reference mixture failed probability identity")
    support=sorted(pmf)
    def quantile(q):
        cdf=0.
        for k in support:
            cdf+=pmf[k]
            if cdf+1e-12>=q:
                return int(k)
        return support[-1]
    return {
        "observed_native_target_shared_pairs": int(observed),
        "expected_shared_pairs_source_reference":float(mean),
        "observed_minus_expected_shared_pairs":float(observed-mean),
        "null_variance":float(var),
        "null_sd":float(math.sqrt(max(var,0))),
        "fixed_margin_possible_support":[support[0] if support else 0,
                                         support[-1] if support else 0],
        "source_mathematical_left_tail_not_biological_p":
            float(sum(p for k,p in pmf.items() if k<=observed)),
        "source_mathematical_right_tail_not_biological_p":
            float(sum(p for k,p in pmf.items() if k>=observed)),
        "source_reference_95pct_interval":[quantile(.025), quantile(.975)],
        "source_strata":len(parts),
        "permutable_strata":n_informative,
        "individually_permutable_butterfly_species":len(species_informative),
        "source_recording_information_gate":(
            "SOURCE_REFERENCE_INFORMATIVE" if len(species_informative)>=
            MIN_INFORMATIVE_BUTTERFLIES and mean>=MIN_EXPECTED_OVERLAP
            else "SOURCE_REFERENCE_TOO_SPARSE_FOR_INTERPRETATION"
        ),
        "conditional_model_is_artificial_not_biological_randomization":True,
    }


def family_strata(pairs, families, native, target):
    result=[]
    groups=defaultdict(set)
    for pid in pairs:
        val=families.get(pid,["<UNKNOWN>"])
        if len(val)!=1:
            raise ValueError("Original accepted WCVP host plant family is ambiguous")
        groups[val[0]].add(pid)
    for family,units in sorted(groups.items()):
        N=len(native & units)
        T=len(target & units)
        O=len(native & target & units)
        result.append({"U":len(units),"N":N,"T":T,"O":O,
                       "family":family})
    return result


def compare_documentation(sets_by_species, plant_families, target_type):
    if target_type not in ("INTRODUCED_CONFIRMED","WCVP_UNLISTED"):
        raise ValueError("Unfrozen introduced/native documentation comparison")
    coarse=[]
    conditional=[]
    totals=Counter()
    species_w_both=0
    for species in sorted(sets_by_species):
        groups=sets_by_species[species]
        U=groups["all"]
        N=groups["NATIVE_CONFIRMED"]
        T=groups[target_type]
        if not N<=U or not T<=U:
            raise ValueError("Documented accepted plant pair absent from source universe")
        overlap=len(N&T)
        coarse.append({"species":species,"U":len(U),"N":len(N),
                       "T":len(T),"O":overlap})
        for r in family_strata(U,plant_families,N,T):
            conditional.append({"species":species,**r})
        totals["accepted_host_universe_species_pair_sum"]+=len(U)
        totals["source_native_documented_species_host_pairs"]+=len(N)
        totals["source_target_documented_species_host_pairs"]+=len(T)
        totals["source_native_and_target_documented_pairs"]+=overlap
        totals["butterfly_species_both_native_and_target_source_docs"]+=int(
            bool(N) and bool(T))
        totals["target_positive_butterflies"]+=int(bool(T))
        totals["native_positive_butterflies"]+=int(bool(N))
    out={
        "source_margins":dict(totals),
        "within_butterfly_uniform_host_list_null":
            exact_overlap_reference(coarse,totals["source_native_and_target_documented_pairs"]),
        "within_butterfly_and_host_family_null":
            exact_overlap_reference(conditional,totals["source_native_and_target_documented_pairs"]),
        "ecological_interpretation_prohibited": (
            "No native HOSTS report does not establish a newly adopted host. "
            "The source host list is not a random ecological pool; even exact "
            "hypergeometric tails are artificial fixed-margin source references."
        )
    }
    return out


def source_host_sets(hosts_path, accepted, focal, *, native, current, atoms,
                     atom_meta, regions, threshold):
    all_hosts=defaultdict(set)
    for (species, raw),pid in accepted.items():
        all_hosts[species].add(pid)
    triples=set()
    raw=Counter()
    for r in read_rows(hosts_path):
        raw["total"]+=1
        fields=("Insect Genus","Insect Species","Hostplant Genus","Hostplant Species")
        if not all(str(r.get(k) or "").strip() for k in fields):
            continue
        raw["full_binomials"]+=1
        sp=r["Insect Genus"].strip()+" "+r["Insect Species"].strip()
        if sp not in focal:
            continue
        raw["focal"]+=1
        host=r["Hostplant Genus"].strip()+" "+r["Hostplant Species"].strip()
        pid=accepted.get((sp,host))
        if pid is None:
            raw["unmatched"]+=1
            continue
        raw["accepted"]+=1
        label=(r.get("Location") or "").strip()
        triples.add((sp,pid,label))
    if (raw["total"],raw["full_binomials"],raw["focal"],raw["accepted"],
        len(triples))!=(EXPECTED_HOSTS_ROWS,EXPECTED_SPECIES_HOST_ROWS,
                        4695,3898,3586):
        raise ValueError("Original 239 source accounting drifted")
    cache={}
    out={sp:{"all":set(all_hosts[sp]),**{s:set() for s in STATUS}}
         for sp in focal}
    hist=Counter()
    for sp,pid,label in sorted(triples):
        k=(pid,label)
        if k not in cache:
            cache[k]=classify_location(
                pid,label,threshold=threshold,atoms=atoms,
                atom_metadata=atom_meta,original_regions=regions,
                native=native,current=current)
        source_status=cache[k]
        hist[source_status]+=1
        if source_status in STATUS:
            out[sp][source_status].add(pid)
    return out, dict(raw), dict(hist)


def run(args):
    source=args.author_source
    m=json.loads((source/"SOURCE_MANIFEST.json").read_text())
    if m.get("inputs_count")!=11:
        raise ValueError("Not exact 11 member original author source manifest")
    for name, expected in m["sha256"].items():
        if sha256(source/name)!=expected:
            raise ValueError("Original author-only HOSTS/WCVP member hash changed")
    root=source/"inputs"
    roster=read_rows(root/"anthropogenic_species_metrics.csv")
    focal={r["species"] for r in roster}
    if len(roster)!=239 or len(focal)!=239:
        raise ValueError("Wrong original 239 cohort")
    accepted,meta,families=accepted_pairs(root/"original/insect_host_accepted.csv",focal)
    if meta["source_ambiguous_accepted_id_names"]:
        raise ValueError("Source accepted host pair ambiguous")
    native, current=wcvp_range(root/"original/native_extant_nondoubtful_wgsrpd3.csv",
                               root/"original/extant_nondoubtful_wgsrpd3.csv")
    atoms,atom_meta,regions,geo=load_country_atoms(
        args.country_atoms,args.country_atoms_manifest,args.country_summary)
    outcomes={}
    for q,label in ((.95,"0.95"),(.80,"0.8")):
        data,flow,status=source_host_sets(args.hosts_resource,accepted,focal,
            native=native,current=current,atoms=atoms,atom_meta=atom_meta,
            regions=regions,threshold=q)
        if label=="0.95" and (status.get("INTRODUCED_CONFIRMED"),
                              status.get("NATIVE_CONFIRMED"),
                              status.get("WCVP_UNLISTED"))!=(179,399,303):
            raise ValueError("Original previous WCVP report geography source margin drift")
        introduced=compare_documentation(data,families,"INTRODUCED_CONFIRMED")
        unlisted=compare_documentation(data,families,"WCVP_UNLISTED")
        if label=="0.95" and (
            introduced["source_margins"]["source_target_documented_species_host_pairs"],
            introduced["source_margins"]["source_native_and_target_documented_pairs"],
            introduced["source_margins"]["target_positive_butterflies"])!=(156,11,62):
            raise ValueError("Original introduced documented pair source outcome drifted")
        outcomes[label]={
            "source_unique_original_species_plant_location_triples":sum(status.values()),
            "source_origin_location_counts":status,
            "introduced_plant_area_source_reports":introduced,
            "WCVP_unlisted_plant_area_report_control":unlisted,
        }
    result={
        "schema":SCHEMA,
        "status":"SOURCE_HOST_LIST_DOC_OVERLAP_CONDITIONAL_REFERENCE_ONLY",
        "source":{
            "focal_original_butterflies":239,
            "HOSTS_resource_sha256":sha256(args.hosts_resource),
            "original_author_source_zip_sha256":SOURCE_ZIP_SHA,
            "original_host_pair_accepted_sha256":
                sha256(root/"original/insect_host_accepted.csv"),
            "original_country_geography":geo,
            "original_WCVP_source_sha256":{
                "native":sha256(root/"original/native_extant_nondoubtful_wgsrpd3.csv"),
                "current":sha256(root/"original/extant_nondoubtful_wgsrpd3.csv"),
            },
            "original_plant_family_resolution":meta,
        },
        "geographic_area_thresholds":outcomes,
        "not_identified":[
            "The butterfly's native versus invaded range and arrival time.",
            "Native-area HOSTS record coverage and independent evidence of native non-use.",
            "Local oviposition, realized larval host acceptance, wild recruitment and fitness.",
            "A randomized geography/host-plant source mechanism satisfying exchangeability.",
            "Crop identity, habitat and climate-linked true niche breadth of the butterfly.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":result["status"],
        "original239":239,
        "original_source95_triples":outcomes["0.95"]["source_unique_original_species_plant_location_triples"],
        "primary95":{
            key:outcomes["0.95"]["introduced_plant_area_source_reports"][key]
            for key in ("source_margins","within_butterfly_uniform_host_list_null",
                        "within_butterfly_and_host_family_null")},
        "unlisted95":{
            key:outcomes["0.95"]["WCVP_unlisted_plant_area_report_control"][key]
            for key in ("source_margins","within_butterfly_uniform_host_list_null",
                        "within_butterfly_and_host_family_null")},
        "intro80":outcomes["0.8"]["introduced_plant_area_source_reports"][
            "source_margins"],
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",type=Path,required=True)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--country-atoms",type=Path,required=True)
    p.add_argument("--country-atoms-manifest",type=Path,required=True)
    p.add_argument("--country-summary",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
