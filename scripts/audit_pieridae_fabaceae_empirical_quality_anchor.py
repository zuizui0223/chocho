#!/usr/bin/env python3
"""Audit independent Phoebis-Senna rearing outcomes and HOSTS-link sensitivity.

Published observational data (not novel field data):
Koptur et al. 2024 Insects 15:123, DOI 10.3390/insects15020123, Table 1.
Reported "Counted" = adult emerged OR parasitoid emerged.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
from statsmodels.stats.contingency_tables import StratifiedTable, Table2x2

DOI="10.3390/insects15020123"
# Records: locality, Senna host, local origin, Found, Counted, Parasitized
TABLE1=[
 ("Pinecrest","Senna chapmanii","native",16,13,1),
 ("South Miami","Senna chapmanii","native",67,55,6),
 ("Westchester","Senna chapmanii","native",26,21,0),
 ("Pinecrest","Senna ligustrina","native",21,17,2),
 ("South Miami","Senna ligustrina","native",1,1,0),
 ("Westchester","Senna ligustrina","native",174,128,7),
 ("Pinecrest","Senna polyphylla","introduced",5,5,2),
 ("South Miami","Senna polyphylla","introduced",45,32,5),
 ("Westchester","Senna polyphylla","introduced",8,7,1),
 ("Pinecrest","Senna surattensis","introduced",17,16,3),
 ("South Miami","Senna surattensis","introduced",83,65,14),
 ("Westchester","Senna surattensis","introduced",231,158,24),
]
SITES=["Pinecrest","South Miami","Westchester"]
SUPPORT_SPECIES=["Phoebis philea","Phoebis sennae"]
EXTERNAL_HOSTS=["Senna ligustrina","Senna surattensis"]

def read_rows(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))

def union_ranges(ids,distribution):
    result=set()
    for key in ids:result|=distribution.get(key,set())
    return result

def summarize_rearing():
    by_site=defaultdict(lambda:defaultdict(lambda:Counter()))
    for site,plant,origin,found,counted,parasitized in TABLE1:
        assert 0<=parasitized<=counted<=found
        c=by_site[site][origin]
        c.update({"found":found,"counted":counted,"parasitized":parasitized,"adult_confirmed":counted-parasitized})
    expected={
        "Pinecrest":{"native":(37,30,3),"introduced":(22,21,5)},
        "South Miami":{"native":(68,56,6),"introduced":(128,97,19)},
        "Westchester":{"native":(200,149,7),"introduced":(239,165,25)}
    }
    for site in SITES:
        for origin in ("native","introduced"):
            c=by_site[site][origin]
            if (c["found"],c["counted"],c["parasitized"])!=expected[site][origin]:
                raise RuntimeError(f"Table 1 count drift for {site} {origin}")
    totals={origin:Counter() for origin in ("native","introduced")}
    for site in SITES:
        for origin in totals:totals[origin].update(by_site[site][origin])
    if (totals["native"]["counted"],totals["native"]["parasitized"],totals["introduced"]["counted"],totals["introduced"]["parasitized"])!=(235,17,283,49):
        raise RuntimeError("Table 1 summed fate count drift")
    mat=[
        [
            [by_site[site]["introduced"]["parasitized"],by_site[site]["introduced"]["adult_confirmed"]],
            [by_site[site]["native"]["parasitized"],by_site[site]["native"]["adult_confirmed"]],
        ]
        for site in SITES
    ]
    strat=StratifiedTable(mat)
    unadjusted=Table2x2([
        [totals["introduced"]["parasitized"],totals["introduced"]["adult_confirmed"]],
        [totals["native"]["parasitized"],totals["native"]["adult_confirmed"]],
    ])
    results={
        "source":{"doi":DOI,"title":"Pierid Butterflies, Legume Hostplants, and Parasitoids in Urban Areas of Southern Florida","table":"Table 1; three urban Miami localities","year":2024},
        "conditional_denominator":"Only individuals counted with adult emergence or parasitoid emergence; 176 found individuals had unaccounted fates and are excluded from parasite probability",
        "by_site":[{
            "site":site,
            "native":dict(by_site[site]["native"]),
            "introduced":dict(by_site[site]["introduced"]),
            "conditional_parasitism_or":float(Table2x2(mat[i]).oddsratio)
        } for i,site in enumerate(SITES)],
        "totals":{k:dict(v) for k,v in totals.items()},
        "summary":{
            "native_parasitoid_fraction_among_counted":17/235,
            "introduced_parasitoid_fraction_among_counted":49/283,
            "raw_odds_ratio":float(unadjusted.oddsratio),
            "raw_or_ci95":[float(x) for x in unadjusted.oddsratio_confint()],
            "site_stratified_mantel_haenszel_or":float(strat.oddsratio_pooled),
            "site_stratified_mh_or_ci95":[float(x) for x in strat.oddsratio_pooled_confint()],
            "site_stratified_mh_p":float(strat.test_null_odds().pvalue),
            "between_site_odds_heterogeneity_p":float(strat.test_equal_odds().pvalue),
            "unaccounted_outcomes_native":305-235,
            "unaccounted_outcomes_introduced":389-283
        },
        "source_limitations":[
            "This is a reanalysis of already published counts, not an independent ecological discovery.",
            "Plant taxa are not randomized to origin, so Senna species, setting, season and origin are confounded.",
            "Weekly observations and larvae on the same shrubs induce dependence; individual-level contingency p-values may be anti-conservative.",
            "Outcome-conditioned Counted excludes lost and unresolved juveniles; ratios are NOT whole-cohort parasitism or survival rates.",
            "Parasitoids were identified across Phoebis hosts; species-specific parasitoid risks cannot be inferred from Table 1.",
            "The abstract misnames the orange-barred Phoebis; the text/figures refer to Phoebis philea. Do not use species-specific numeric estimates without original data."
        ]
    }
    return results

def test_external_link_augmentation(interaction,native,contemporary,panel):
    by_sp=defaultdict(set); fam_by_id={}; name_to_ids=defaultdict(set)
    for r in interaction:
        sp=str(r["insect_species"]).strip();pid=str(r["accepted_plant_name_id"]).strip()
        by_sp[sp].add(pid)
        fam_by_id[pid]=str(r["family"]).strip()
        name_to_ids[str(r["accepted_name"]).strip()].add(pid)
    for name in EXTERNAL_HOSTS:
        if len(name_to_ids[name])!=1:raise RuntimeError(f"ambiguous or missing plant name: {name}")
    ids={name:next(iter(name_to_ids[name])) for name in EXTERNAL_HOSTS}
    if ids!={"Senna ligustrina":"2490490","Senna surattensis":"2897969"}:
        raise RuntimeError(f"WCVP frozen ID mismatch: {ids}")
    nat=defaultdict(set);cur=defaultdict(set)
    for r in native:nat[r["accepted_plant_name_id"]].add(r["area_code_l3"])
    for r in contemporary:cur[r["accepted_plant_name_id"]].add(r["area_code_l3"])
    if "FLA" not in nat[ids["Senna ligustrina"]]:raise RuntimeError("native Florida Senna ligustrina status drift")
    if "FLA" in nat[ids["Senna surattensis"]] or "FLA" not in cur[ids["Senna surattensis"]]:
        raise RuntimeError("introduced Florida Senna surattensis status drift")
    def summary(sp,hosts):
        n=union_ranges(hosts,nat);c=union_ranges(hosts,cur)
        return n,c,c-n
    rows=[]
    for sp in SUPPORT_SPECIES:
        base=by_sp[sp]
        new=base|set(ids.values())
        assert len(new-base)==2, f"expected both independently published Senna edges to be missing for {sp}"
        oldn,oldc,olda=summary(sp,base)
        newn,newc,newa=summary(sp,new)
        rows.append({
            "species":sp,
            "added_external_exact_host_links":2,
            "original_native_known_resource_regions":len(oldn),
            "augmented_native_known_resource_regions":len(newn),
            "original_contemporary_known_resource_regions":len(oldc),
            "augmented_contemporary_known_resource_regions":len(newc),
            "original_introduced_only_known_resource_regions":len(olda),
            "augmented_introduced_only_known_resource_regions":len(newa),
            "native_regions_gained":";".join(sorted(newn-oldn)),
            "contemporary_regions_gained":";".join(sorted(newc-oldc)),
            "introduced_only_regions_removed":";".join(sorted(olda-newa)),
            "introduced_only_regions_gained":";".join(sorted(newa-olda))
        })
    expected={"Phoebis philea":(83,208,125,85,208,123),"Phoebis sennae":(104,226,122,106,227,121)}
    for r in rows:
        got=(r["original_native_known_resource_regions"],r["original_contemporary_known_resource_regions"],
             r["original_introduced_only_known_resource_regions"],r["augmented_native_known_resource_regions"],
             r["augmented_contemporary_known_resource_regions"],r["augmented_introduced_only_known_resource_regions"])
        if got!=expected[r["species"]]:raise RuntimeError(f"external link regional sensitivity drift {r['species']}: {got}")
    panel_count={};sum_original=sum_new=0
    for sp in panel:
        oldn,oldc,olda=summary(sp,by_sp[sp])
        hosts=by_sp[sp]|set(ids.values()) if sp in SUPPORT_SPECIES else by_sp[sp]
        newn,newc,newa=summary(sp,hosts)
        if sp in SUPPORT_SPECIES:panel_count[sp]={"before":len(olda),"after":len(newa)}
        sum_original+=len(olda)
        sum_new+=len(newa)
    result={
        "status":"SENSITIVITY_TO_TWO_EXTERNAL_PUBLISHED_SENNA_LINKS_PER_BUTTERFLY",
        "panel_species":len(panel),
        "source_study_doi":DOI,
        "frozen_WCVP_ids":ids,
        "published_locality":"South Florida; larvae reared on the two exotics Senna polyphylla and Senna surattensis, plus native Senna chapmanii and Senna ligustrina",
        "added_links":[{"butterfly":sp,"host":name,"accepted_plant_id":ids[name]} for sp in SUPPORT_SPECIES for name in EXTERNAL_HOSTS],
        "source_locality_note":"Local rearing establishes use at the study site; applying a link to all WCVP regions is the ORIGINAL model's unverified range-wide assumption.",
        "species_results":rows,
        "frozen_24_species_introduced_only_units_before":sum_original,
        "frozen_24_species_introduced_only_units_after":sum_new,
        "two_species_total_change":sum_new-sum_original,
        "limitations":[
            "Host Senna polyphylla was not in the frozen accepted HOSTS/WCVP sidecar; it is explicitly omitted rather than assigned a guessed ID.",
            "Missing HOSTS links can understate known local trophic options even when regional resource envelope sizes barely change due to other plants.",
            "Newly inferred native range opportunity in AND and MDV is a model extrapolation, not a documented butterfly occurrence or successful larvae there."
        ]
    }
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ("insect_host_csv","native_distribution_csv","contemporary_distribution_csv","panel_protocol_json","output_json","output_species_csv"):
        p.add_argument("--"+name.replace("_","-"),required=True,type=Path)
    a=p.parse_args()
    panel=json.loads(a.panel_protocol_json.read_text(encoding="utf-8"))["species"]
    if len(panel)!=24:raise RuntimeError("expected 24 butterfly panel")
    parasitism=summarize_rearing()
    link_validation=test_external_link_augmentation(
        read_rows(a.insect_host_csv),read_rows(a.native_distribution_csv),read_rows(a.contemporary_distribution_csv),panel
    )
    results={
        "schema":"chocho_pieridae_fabaceae_empirical_quality_anchor_v0.1",
        "status":"EXPLORATORY_PUBLISHED_DATA_REANALYSIS_AND_HOST_LINK_COVERAGE_DIAGNOSTIC",
        "rearing":parasitism,"host_link_sensitivity":link_validation,
        "inference":"The observed spatial network is potential resource availability. This independent literature example documents both successful rearing on some introduced hosts and greater parasitoid emergence odds on introduced hosts. Both can coexist; neither implies a global causal effect or generalizable butterfly fitness outcome.",
        "do_not_promote_to_main_GEB":True
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    with a.output_species_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(link_validation["species_results"][0]))
        w.writeheader();w.writerows(link_validation["species_results"])
    print(json.dumps({"rearing":parasitism["summary"],"link_sensitivity":link_validation["species_results"],
                      "24species_delta":link_validation["two_species_total_change"]},indent=2))

if __name__=="__main__":main()
