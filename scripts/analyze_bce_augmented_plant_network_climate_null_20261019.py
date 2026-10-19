#!/usr/bin/env python3
"""BCE-augmented butterfly climate analogue analysis with plant-network nulls.

Data: the exact author-only WCVP/HOSTS/CHELSA archive previously verified
for the original 239 butterflies, supplemented by exactly 1,027 accepted-ID
host links from the pinned partial European BCE/Clarke supplement.

This analysis models potential regional larval host opportunity ONLY, NOT
observed butterfly distribution, host acceptance, larval fitness or breeding.
The 52 otherwise missing host-family names are pinned via a GBIF taxonomy
backbone sidecar. All other plant families use source HOSTS plant IDs/genus.
No raw WCVP/HOSTS geographic table is written to the result artifact.

Both MCMC nulls preserve:
 - original butterfly-to-host plant identity (including the BCE additions);
 - native plant and butterfly resource geography, fixed throughout draws;
 - the plant's introduced number of regions by plant family / continent;
 - region-specific introduced plant-family degree;
The second mode also fixes each plant's coarse BIO1×BIO12 climatic region
strata within continents (2×2 quantile bins).
Because butterfly-added resource unions change under plant swaps, compare
*aligned fraction* O/A with each draw's own A, never raw O.

Experimental status: retrospective exploratory structural null; not an
independent randomized intervention and not a confirmatory p-value.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from analyze_butterfly_climate_analog_macro_20261019 import (
    qbounds, read_csv, source_model,
)
from analyze_climate_analog_fixedmargin_null_20261019 import ess_approx

ROOT=Path(__file__).resolve().parents[1]
TAXON=ROOT/"data/exploratory/bce_additional_plant_family_gbif_crosswalk_20261019.json"
SEEDS=(20261019,20261020)
MODES=("plant_family_continent","plant_family_continent_BIO1_BIO12_2x2")
QS=(.05,.10)
EXPECTED_BASE=(26530,41083,14553)
EXPECTED_AUGMENTED=(29852,44819,14967)
EXPECTED_CLIMATE={(0.05):(9949,14961),(0.10):(7370,14961)}
EXPECTED_MATCHED_ORIGINAL_ELIGIBLE={
    (0.05):(9903,14903),
    (0.10):(7331,14903),
}
ALL_RESPONSE_KEYS=("0.05","0.1","matched_0.05","matched_0.1")
EXPECTED_EXTRA_HOST_EDGES=1027
EXPECTED_EXTRA_ACCEPTED_PLANTS=679
EXPECTED_INTRODUCED_PLANT_REGION_EDGES=37869
EXPECTED_OUTSIDE_CLIMATE={"Coenonympha oedippus":["ANT","HMD"]}
DEFAULT_DRAWS=199
DEFAULT_BURN=10
DEFAULT_SPACING=4

def audited_extra_host_families(source_root,extra_ids):
    """Never guess a plant family from id order; fail on ambiguity."""
    sys.path.insert(0,str(source_root/"code"))
    from analyze_bce_resource_novelty_source_sensitivity import (
        rows as source_rows,
    )
    all_hosts=source_rows(source_root/"inputs/original/insect_host_accepted.csv")
    by_id=defaultdict(set)
    by_genus=defaultdict(set)
    for row in all_hosts:
        taxon_id=row["accepted_plant_name_id"].strip()
        family=row["family"].strip()
        genus=row["accepted_name"].strip().split()[0]
        if family and taxon_id and genus:
            by_id[taxon_id].add(family)
            by_genus[genus].add(family)
    candidates=json.loads(
        (source_root/"inputs/bce/accepted_host_edges.json").read_text(encoding="utf-8"))
    accepted_name={}
    for species_rec in candidates["per_butterfly"]:
        for item in species_rec["details"]:
            if item["match_status"]!="ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY":
                continue
            host=item["accepted_plant_name_id"]
            name=item["accepted_name"]
            if host in accepted_name and accepted_name[host]!=name:
                raise ValueError("Ambiguous WCVP accepted botanical taxon name")
            accepted_name[host]=name
    assert set(accepted_name)==set(extra_ids) and len(extra_ids)==679
    sidecar=json.loads(TAXON.read_text(encoding="utf-8"))
    assert sidecar["schema"]=="chocho_bce_partial_unseen_genera_GBIF_plant_family_crosswalk_v0.1"
    assert len(sidecar["crosswalk"])==52
    tax={
        item["wcvp_accepted_plant_name_id"]:item
        for item in sidecar["crosswalk"]
    }
    assert len(tax)==52
    families={}
    evidence=Counter()
    exceptions=[]
    for h in sorted(extra_ids):
        direct=by_id[h]
        genus=by_genus[accepted_name[h].split()[0]]
        if len(direct)==1:
            family=next(iter(direct))
            evidence["existing_accepted_plant_id"]+=1
        elif direct:
            raise ValueError("Conflicting source accepted plant family "+h)
        elif len(genus)==1:
            family=next(iter(genus))
            evidence["unique_original_HOSTS_genus"]+=1
        elif genus:
            raise ValueError("Conflicting source botanical genus family "+h)
        else:
            rec=tax.get(h)
            if not rec or rec["wcvp_accepted_name"]!=accepted_name[h]:
                raise ValueError("No GBIF-checked taxonomic family for "+h)
            family=rec["verified_family"]
            exact=(rec["match_type"]=="EXACT" and rec["confidence"]>=90 and
                   rec["gbif_matched_canonical"]==accepted_name[h])
            single_fuzzy=(
                accepted_name[h]=="Galactites tomentosus"
                and family=="Asteraceae" and rec["match_type"]=="FUZZY"
                and rec["gbif_matched_canonical"]=="Galactites tomentosa"
                and "Galactites&rank=GENUS" in rec["manual_genus_corroboration"]
            )
            if not (exact or single_fuzzy):
                raise ValueError("Unverified crosswalk family, no rescue "+h)
            if single_fuzzy:exceptions.append(h)
            evidence["GBIF_backbone_checked"]+=1
        families[h]=family
    assert dict(evidence)=={
        "existing_accepted_plant_id":323,
        "unique_original_HOSTS_genus":304,
        "GBIF_backbone_checked":52,
    },dict(evidence)
    assert exceptions==["3043225"]
    assert set(tax)=={
        h for h in extra_ids if not by_id[h] and not by_genus[accepted_name[h].split()[0]]
    }
    return families,dict(evidence),{
        "pinned_crosswalk_file":str(TAXON.relative_to(ROOT)),
        "GBIF_taxa_resolved":len(tax),
        "cross_database_botanical_spelling_exception":exceptions,
        "attribution":"GBIF taxonomy family only, NOT GBIF range or host interactions"
    }


def source_augmented_model(source_root):
    sys.path.insert(0,str(source_root/"code"))
    from analyze_bce_resource_novelty_source_sensitivity import (
        load_candidate_host_links,host_geography,combine_host_ranges,
    )
    base,regions,climate,orig,original_provenance=source_model(source_root)
    assert len(regions)==365 and set(regions)==set(climate)
    data=source_root/"inputs"
    bce,source=load_candidate_host_links(data/"bce/accepted_host_edges.json")
    assert len(bce)==81 and sum(len(v) for v in bce.values())==1027
    additional_ids=set().union(*bce.values())
    assert len(additional_ids)==679
    mapped,source_types,gbif=audited_extra_host_families(source_root,additional_ids)
    original_focal_hosts=set().union(*[set(h) for h in base["species_hosts"]])
    links=[
        sorted(set(hosts)|bce.get(s,set()))
        for s,hosts in zip(base["species"],base["species_hosts"])
    ]
    assert len(links)==239 and len(set().union(*map(set,links)))==2228
    # Derive WCVP accepted-ID plant native/current ranges, not new observations.
    n0=host_geography(data/"original/native_extant_nondoubtful_wgsrpd3.csv")
    c0=host_geography(data/"original/extant_nondoubtful_wgsrpd3.csv")
    n_extra=host_geography(data/"bce/additional_native_wgsrpd3.csv")
    c_extra=host_geography(data/"bce/additional_contemporary_wgsrpd3.csv")
    native_geo=combine_host_ranges(n0,n_extra,additional_ids)
    contemporary_geo=combine_host_ranges(c0,c_extra,additional_ids)
    index={code:i for i,code in enumerate(regions)}
    all_hosts=set().union(*map(set,links))
    level1_lookup=json.loads(
        (data/"level1_lookup.json").read_text(encoding="utf-8"))["units"]
    level1={i:level1_lookup[code] for code,i in index.items()}
    full_family={**base["host_families"],**mapped}
    assert set(all_hosts)<=set(full_family)
    native_sets={}
    introduced_sets={}
    full_native_total=full_contemporary_total=full_added_total=0
    out_of_climate={}
    baseline_climate=[]
    added_climate=[]
    added_climate_all=0
    for host in all_hosts:
        native=native_geo.get(host,set())
        current=contemporary_geo.get(host,set())
        if not native<=current:
            raise ValueError("WCVP host native not subset of current "+host)
        native_sets[host]={index[u] for u in native if u in index}
        introduced_sets[host]={index[u] for u in current-native if u in index}
    for sp,hosts in zip(base["species"],links):
        native=set().union(*(native_geo.get(h,set()) for h in hosts))
        current=set().union(*(contemporary_geo.get(h,set()) for h in hosts))
        added=current-native
        if not native<=current:raise ValueError("source butterfly graph drift")
        full_native_total+=len(native)
        full_contemporary_total+=len(current)
        full_added_total+=len(added)
        outside=sorted(added-set(index))
        if outside:out_of_climate[sp]=outside
        nm=0
        introduced=0
        for h in hosts:
            for k in native_sets[h]:nm |= 1<<k
            for k in introduced_sets[h]:introduced|=1<<k
        new=introduced & ~nm
        if nm.bit_count()!=len(native & set(index)):
            raise ValueError("native source vs 365 universe mismatch "+sp)
        if new.bit_count()!=len(added & set(index)):
            raise ValueError("added source vs 365 universe mismatch "+sp)
        baseline_climate.append(nm)
        added_climate.append(new)
        added_climate_all+=new.bit_count()
    assert (full_native_total,full_contemporary_total,full_added_total)==EXPECTED_AUGMENTED
    assert added_climate_all==14965
    assert out_of_climate==EXPECTED_OUTSIDE_CLIMATE
    original_eligible={
        row["species"] for row in orig if row["native"]>=10
    }
    augmented_eligible={
        sp for sp,mask in zip(base["species"],baseline_climate)
        if mask.bit_count()>=10
    }
    assert len(original_eligible)==234 and len(augmented_eligible)==235
    assert augmented_eligible-original_eligible=={"Satyrium pruni"}
    assert not (original_eligible-augmented_eligible)
    assert sum(len(v) for v in introduced_sets.values())==37869
    # Validate original source-model geography as independent anti-drift.
    assert (sum(x["native"] for x in orig),
            sum(x["contemporary"] for x in orig),
            sum(x["added"] for x in orig))==EXPECTED_BASE
    output={
        "species":base["species"],
        "species_hosts":links,
        "native_species":baseline_climate,
        "native_sets":native_sets,
        "introduced_sets":introduced_sets,
        "host_families":full_family,
        "level1":level1,
    }
    return output,regions,climate,{
        "frozen_source":original_provenance,
        "BCE_1027_crosswalk_SHA256":source["bce_accepted_edges_sha256"],
        "extra_host_family_method_counts":source_types,
        "external_taxonomy":gbif,
        "whole_source_native_contemporary_added":list(EXPECTED_AUGMENTED),
        "source_added_regions_missing_climate":out_of_climate,
        "climate_covered_added_all_239":added_climate_all,
        "introduced_plant_region_links":sum(map(len,introduced_sets.values())),
        "climate_eligible_original_HOSTS_butterflies":len(original_eligible),
        "newly_eligible_under_BCE":sorted(augmented_eligible-original_eligible)
    }


def analog_mask(model,regions,climate):
    good=[]
    result={q:{} for q in QS}
    for i,mask in enumerate(model["native_species"]):
        native=[j for j in range(len(regions)) if mask>>j & 1]
        if len(native)<10:continue
        good.append(i)
        for q in QS:
            a,b=qbounds([climate[regions[j]][0] for j in native],q)
            c,d=qbounds([climate[regions[j]][2] for j in native],q)
            bmask=0
            for j,code in enumerate(regions):
                bi=climate[code]
                if a<=bi[0]<=b and c<=bi[2]<=d:bmask|=1<<j
            result[q][i]=bmask
    assert len(good)==235
    return result,good


def metrics(model,plants,analogs):
    # Rebuild newly supported resource union from exact fixed host identities,
    # NEVER multiply raw introduced plant incidences by nominal butterfly count.
    plantbits={
        h:sum(1<<i for i in cells) for h,cells in plants.items()
    }
    out={}
    newly={}
    for i,hosts in enumerate(model["species_hosts"]):
        native=model["native_species"][i]
        union=0
        for h in hosts:union |= plantbits[h]
        newly[i]=union & ~native
    matched_exclude=model["species"].index("Satyrium pruni")
    for q,eligible in analogs.items():
        A=sum(newly[i].bit_count() for i in eligible)
        O=sum((newly[i]&analog).bit_count() for i,analog in eligible.items())
        if A<=0 or O>A:raise ValueError("invalid added resource denominator")
        out[str(q)]={"aligned":O,"evaluable_added":A,"aligned_fraction":O/A}
        assert matched_exclude in eligible
        matched_A=A-newly[matched_exclude].bit_count()
        matched_O=O-(newly[matched_exclude]&eligible[matched_exclude]).bit_count()
        assert matched_A>0 and 0<=matched_O<=matched_A
        out["matched_"+str(q)]={
            "aligned":matched_O,"evaluable_added":matched_A,
            "aligned_fraction":matched_O/matched_A,
        }
    return out,sum(x.bit_count() for x in newly.values())


def marginal_signature(model,plants,classes):
    by_plant=Counter();by_region_family=Counter()
    for h,regions in plants.items():
        assert not (regions & model["native_sets"].get(h,set())),h
        fam=model["host_families"][h]
        for idx in regions:
            key=(model["level1"][idx],
                 classes.get(idx) if classes is not None else None)
            by_plant[h,key]+=1
            by_region_family[idx,fam]+=1
    return by_plant,by_region_family,sum(map(len,plants.values()))


def null_chain(model,regions,climate,analogs,source_plants,mode,seed,*,draws,burn,spacing):
    from analyze_butterfly_shared_host_structural_null import (
        attempt_range_swaps,make_groups
    )
    from analyze_butterfly_climate_matched_plant_range_null import (
        climate_strata,grouped_introduced_edges
    )
    plants={h:set(cells) for h,cells in source_plants.items()}
    E=sum(map(len,plants.values()))
    assert E==EXPECTED_INTRODUCED_PLANT_REGION_EDGES
    if mode=="plant_family_continent":
        groups,cumulative,total=make_groups(
            model["native_sets"],plants,model["host_families"],model["level1"])
        assert total==E
        classes=None
        structure={"total_introduced_plant_region_edges":E,
            "climate_covered_fraction":1.0,
            "viable_edge_fraction":(cumulative[-1] if cumulative else 0)/E}
    else:
        classes,_=climate_strata(
            climate,regions,model["level1"],"temperature_precipitation_2x2")
        groups,cumulative,structure=grouped_introduced_edges(
            plants,model["host_families"],model["level1"],classes)
        # Source API returns viable_fraction_of_total, while this script
        # uses one common key across both null modes. Preserve the original.
        structure={
            **structure,
            "viable_edge_fraction":structure["viable_fraction_of_total"],
        }
    assert structure["total_introduced_plant_region_edges"]==E
    assert structure["viable_edge_fraction"]>=0.75
    assert structure["climate_covered_fraction"]>=0.95
    source_sig=marginal_signature(model,source_plants,classes)
    assert source_sig[2]==E
    rng=random.Random(seed)
    accepted_burn=attempt_range_swaps(
        groups,cumulative,model["native_sets"],plants,burn*E,rng)
    moved=sum(len(plants[h]^source_plants[h]) for h in plants)/(2*E)
    assert marginal_signature(model,plants,classes)==source_sig
    history={key:[] for key in ALL_RESPONSE_KEYS}
    sampled_accept=[]
    for draw in range(draws):
        n=attempt_range_swaps(
            groups,cumulative,model["native_sets"],plants,spacing*E,rng)
        sampled_accept.append(n)
        values,total=metrics(model,plants,analogs)
        for key,datum in values.items():
            history[key].append({
                **datum,"total_all_butterflies_new_regions":total,
            })
        if draw in (0,draws//2,draws-1):
            if marginal_signature(model,plants,classes)!=source_sig:
                raise ValueError("plant-family or plant-stratum margin drift")
    summaries={}
    for k in ALL_RESPONSE_KEYS:
        values=[z["aligned_fraction"] for z in history[k]]
        ds=[z["evaluable_added"] for z in history[k]]
        summaries[k]={
            "ess_diagnostic":ess_approx(values),
            "null_fraction_median":float(np.median(values)),
            "null_fraction_ci_95":[float(np.quantile(values,p)) for p in (.025,.975)],
            "null_added_regions_median":float(np.median(ds)),
            "null_added_regions_range":[min(ds),max(ds)]
        }
    return {
        "seed":seed,
        "structure":structure,
        "moved_after_burn":moved,
        "accepted_burn":accepted_burn,
        "accepted_per_sample":float(np.mean(sampled_accept)),
        "samples":history,
        "summaries":summaries,
        "source_invariants":"plant by continent(/coarse-climate) and recipient plant family; checked at burn, beginning, middle, last"
    }


def qualify(observed,chains,mode):
    result={};gates={}
    for k in ALL_RESPONSE_KEYS:
        stats=[c["summaries"][k] for c in chains]
        a,b=[z["null_fraction_median"] for z in stats]
        ess=[z["ess_diagnostic"]["ess_approx"] for z in stats]
        mid_agree=abs(a-b)<=2*max(z["ess_diagnostic"]["std"] for z in stats)
        ok=(all(e>=30 for e in ess) and mid_agree
            and all(x["moved_after_burn"]>=0.10 for x in chains)
            and all(x["accepted_burn"]>=100 and x["accepted_per_sample"]>=100
                    for x in chains))
        gates[k]={"qualified":ok,"ess_each_ge_30":all(e>=30 for e in ess),
                  "medians_agree":mid_agree,
                  "moved_ge_10pct":all(x["moved_after_burn"]>=.10 for x in chains),
                  "sample_accepted_ge_100":all(x["accepted_per_sample"]>=100 for x in chains)}
        frac=observed[k]["aligned_fraction"]
        result[k]={
            "observed":observed[k],
            "null_fraction_median_by_chain":[a,b],
            "null_fraction_95ci_by_chain":[x["null_fraction_ci_95"] for x in stats],
            "observed_minus_null_median_pp_by_chain":[100*(frac-x) for x in (a,b)],
            "null_added_regions_median_by_chain":[z["null_added_regions_median"] for z in stats],
            "null_added_regions_range_by_chain":[z["null_added_regions_range"] for z in stats],
            "ess_by_chain":ess,
            "right_tail_add_one_serial_not_exact_randomization_p_by_chain":[
                (1+sum(z["aligned_fraction"]>=frac for z in chain["samples"][k]))/
                (len(chain["samples"][k])+1)
                for chain in chains
            ] if ok else None,
            "numerically_qualified":ok,
        }
    return {
        "mode":mode,
        "status":"QUALIFIED_EXPLORATORY_SOURCE_CONDITIONAL_MCMC" if
        all(v["qualified"] for v in gates.values()) else
        "NOT_EVALUABLE_MCMC_QUALITY_GATE",
        "quality":gates,
        "result":result,
        "chains":chains,
    }


def run(source_root,*,draws=DEFAULT_DRAWS,burn=DEFAULT_BURN,spacing=DEFAULT_SPACING,modes=MODES):
    model,regions,climate,provenance=source_augmented_model(source_root)
    masks,included=analog_mask(model,regions,climate)
    observed,all_added=metrics(model,model["introduced_sets"],masks)
    assert all_added==14965 and len(included)==235
    for q,expect in EXPECTED_CLIMATE.items():
        assert (
            observed[str(q)]["aligned"],observed[str(q)]["evaluable_added"]
        )==expect,(q,observed[str(q)])
    for q,expect in EXPECTED_MATCHED_ORIGINAL_ELIGIBLE.items():
        key="matched_"+str(q)
        assert (observed[key]["aligned"],observed[key]["evaluable_added"])==expect,            (key,observed[key])
    modes_out={}
    for mode in modes:
        chains=[
            null_chain(model,regions,climate,masks,model["introduced_sets"],
                       mode,seed,draws=draws,burn=burn,spacing=spacing)
            for seed in SEEDS
        ]
        modes_out[mode]=qualify(observed,chains,mode)
    return {
        "schema":"chocho_239_butterfly_bce_augmented_host_network_climate_null_v0.1",
        "status":"QUALIFIED_EXPLORATORY_ALL_TESTED_MODES" if
        all(v["status"]=="QUALIFIED_EXPLORATORY_SOURCE_CONDITIONAL_MCMC"
            for v in modes_out.values()) else "PARTIAL_OR_NOT_EVALUABLE",
        "source_scope":{
            "total_butterflies":239,
            "climate_eligible_butterflies":235,
            "original_eligibility_matched_butterflies":234,
            "accepted_extra_host_links":1027,
            "unique_extra_accepted_plant_ids":679,
            "introduced_plant_region_edges":37869,
            "all_239_new_resource_cells":14967,
            "climate_supported_new_resource_cells_all":14965,
            "climate_eligible_new_resource_cells":14961,
            "climate_missing_added_cells":2,
            "missing_region_codes_by_species":EXPECTED_OUTSIDE_CLIMATE
        },
        "provenance":provenance,
        "observed":observed,
        "draws_per_chain":draws,
        "burn_attempts_per_plant_edge":burn,
        "spacing_attempts_per_plant_edge":spacing,
        "modes":modes_out,
        "boundaries":[
            "Europe-biased partial host catalogue; not a representative global source replication.",
            "Families resolved using frozen HOSTS accepted ID/unique genus and 52 pinned GBIF backbone taxonomies; 1 GBIF/WCVP botanical spelling mismatch audited.",
            "Two subantarctic WCVP introduced-added cells outside CHELSA 365-region universe are retained in source geography, excluded only from climatic response denominator.",
            "No host-plant native abundance, exact climatic physiological niche of butterflies, butterfly occupancy or recruitment, host quality or time series measured.",
            "Plant-null changes unique added butterfly×region count; report fraction O/A with each draw's own denominator.",
            "Serial MCMC add-one right-tail fractions are not exact permutation p-values or evidence of butterfly population causation.",
            "Null conditions coarse plant climate 2×2; finer plant spatial/climate autocorrelation still unmodelled.",
            "Source WCVP/HOSTS accepted IDs and geographic sidecars are author-only; no raw provider redistribution authorised.",
            "Both a 235-butterfly augmented-eligibility response and a fixed-original-eligibility 234-butterfly response are evaluated from identical null draws, excluding only newly eligible Satyrium pruni from the matched response.",
            "This is exploratory, post hoc; original frozen GEB submission manuscript PR38 unchanged."
        ]
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--draws",type=int,default=DEFAULT_DRAWS)
    ap.add_argument("--burn-factor",type=int,default=DEFAULT_BURN)
    ap.add_argument("--spacing-factor",type=int,default=DEFAULT_SPACING)
    ap.add_argument("--mode",choices=MODES,nargs="*",default=list(MODES))
    args=ap.parse_args()
    res=run(args.source_root,draws=args.draws,burn=args.burn_factor,
            spacing=args.spacing_factor,modes=args.mode)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(res,indent=2,sort_keys=True)+"\n",
                                encoding="utf-8")
    print(json.dumps({
        "status":res["status"],"source_scope":res["source_scope"],
        "family_method_counts":res["provenance"]["extra_host_family_method_counts"],
        "observed":res["observed"],
        "modes":{
            mode:{
                "status":x["status"],
                "quality":x["quality"],
                "result":x["result"],
                "moved":[z["moved_after_burn"] for z in x["chains"]],
                "accepted_per_sample":[z["accepted_per_sample"] for z in x["chains"]]
            } for mode,x in res["modes"].items()
        }
    },sort_keys=True))
