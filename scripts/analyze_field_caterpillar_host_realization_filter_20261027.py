#!/usr/bin/env python3
"""External field caterpillar detection vs fixed global butterfly HOSTS links.

External: Whitener & Hurlbert Ecosphere 2026 DOI 10.1002/ecs2.70615,
fixed public hurlbertlab/caterpillars-on-plants commit 2c5690ac...,
their field host plant Table S2, 100 >=40 survey woody plant species.

Local: original author's accepted HOSTS/WCVP 239-butterfly taxonomy.
This is not a study of 239-butterfly larval field occurrence: field outcome
contains larvae from moths AND butterflies. Zero exact HOSTS matches is
documentary uncertainty, NOT plant inability to support butterflies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs, read_rows, sha256,
)

SCHEMA = "chocho_original239_host_opportunity_external_field_larval_filter_20261027_v01"
EXTERNAL_COMMIT="2c5690ac7c6dba7eca506da886c49781c37617bd"
EXTERNAL_BLOB="bc0d66a2fa2e44c4c1d959fa4700f4256cbb238e"
EXTERNAL_PAPER="10.1002/ecs2.70615"
ORIGINAL_SOURCE_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
ORIGINAL_HOSTS_COMMIT="808e0b869f9ec1adf8efff87cf6a395adda103e0"
SEED=20261027
BOOTSTRAPS=499
RIDGE=1e-3


def git_blob_sha(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def normalize_binomial(text):
    pieces=" ".join(str(text or "").strip().split()).split(" ")
    if len(pieces)!=2:
        return None
    genus,species=pieces
    if not (genus[0:1].isupper() and genus[1:].isalpha()
            and genus.isalpha() and species.isalpha()
            and species.islower() and len(species)>=2):
        return None
    if species in ("sp","spp","species","group"):
        return None
    return " ".join(pieces)


def original_survey_plant_panel(path):
    raw=Path(path).read_bytes()
    if git_blob_sha(raw)!=EXTERNAL_BLOB:
        raise ValueError("External exact-sourced original published Table_S2 changed")
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    expected={
        "Plant species", "Family", "USDA Native Status","Origin",
        "n_Surveys","n_Branches","n_Sites","Pct with caterpillars"
    }
    if len(rows)!=100 or not expected<=set(rows[0]):
        raise ValueError("Not original 100 survey-plant source published field subset")
    out=[]
    for row in rows:
        name=normalize_binomial(row["Plant species"])
        if name is None:
            raise ValueError("Source 100-plant species binomial unverified")
        if any(name==x["field_plant_binomial"] for x in out):
            raise ValueError("Duplicate original field plant species")
        origin=row["Origin"].strip().lower()
        if origin not in ("native","alien",""):
            raise ValueError("Unrecognized original USDA field plant origin")
        n=int(row["n_Surveys"])
        if n<40:
            raise ValueError("Original Table S2 40-survey source inclusion changed")
        frac=float(row["Pct with caterpillars"])/100.
        if not 0<=frac<=1:
            raise ValueError("Impossible field larval occurrence estimate")
        larvae=round(n*frac)
        if abs(100*larvae/n-100*frac)>.0051:
            raise ValueError("Rounded published source fraction cannot reconstruct binomial response")
        out.append({
            "field_plant_binomial":name,
            "field_plant_family":row["Family"].strip(),
            "field_plant_genus":name.split()[0],
            "original_USDA_origin":origin or "source_origin_unspecified",
            "field_original_n_surveys":n,
            "field_original_n_branches":int(row["n_Branches"]),
            "field_original_n_sites":int(row["n_Sites"]),
            "field_original_caterpillar_pct":float(row["Pct with caterpillars"]),
            "field_reconstructed_positive_survey_events":larvae,
        })
    original_status=Counter(x["original_USDA_origin"] for x in out)
    if original_status!={"native":85,"alien":14,"source_origin_unspecified":1}:
        raise ValueError("Original public source plant native/alien/unknown composition changed")
    return out,{
        "public_original_table_sha256":hashlib.sha256(raw).hexdigest(),
        "public_original_table_git_blob":git_blob_sha(raw),
        "original_100_plants_by_USDA_origin":dict(original_status),
        "original100_surveys":sum(r["field_original_n_surveys"] for r in out),
    }


def accepted_global_links(path, butterfly_roster):
    accepted,meta,_=accepted_pairs(path,butterfly_roster)
    by_name=defaultdict(set)
    for (sp, raw), pid in accepted.items():
        name=normalize_binomial(raw)
        if name is not None:
            by_name[name].add(pid)
    ambiguity={name for name,ids in by_name.items() if len(ids)>1}
    by_accepted_id=defaultdict(set)
    for (sp,raw),pid in accepted.items():
        by_accepted_id[pid].add(sp)
    total_species_accepted_pairs=sum(len(s) for s in by_accepted_id.values())
    if total_species_accepted_pairs!=2596:
        raise ValueError("Original 2596 accepted 239-species host links changed")
    return by_name,by_accepted_id,ambiguity,{
        **meta,
        "original_239_accepted_species_plant_pairs":total_species_accepted_pairs,
        "original_named_WCVP_accepted_host_ids":len(by_accepted_id),
        "original_unambiguous_source_host_binomial_names":len(by_name)-len(ambiguity),
        "original_ambiguous_source_host_binomial_names":len(ambiguity),
    }


def broad_catalogue_host_links(raw_hosts_path, field_names):
    """All Lepidoptera (moths included), exact source names. Not accepted 239 links."""
    out=defaultdict(set)
    count=Counter()
    with Path(raw_hosts_path).open(encoding="utf-8-sig",newline="") as f:
        it=csv.DictReader(f)
        required={"Insect Genus","Insect Species","Hostplant Genus","Hostplant Species"}
        if not required<=set(it.fieldnames or []):
            raise ValueError("Missing original HOSTS species column names")
        for r in it:
            count["raw_rows"]+=1
            if any(not str(r[k] or "").strip() for k in required):
                continue
            count["species_level_rows"]+=1
            host=normalize_binomial(r["Hostplant Genus"]+" "+r["Hostplant Species"])
            if host not in field_names:
                continue
            insect=normalize_binomial(r["Insect Genus"]+" "+r["Insect Species"])
            if insect is None:
                continue
            out[host].add(insect)
    if (count["raw_rows"],count["species_level_rows"])!=(140485,86826):
        raise ValueError("Original HOSTS full Lepidoptera source row count changed")
    return out,dict(count)


def join(survey,by_name,by_pid,ambiguous,other):
    new=[]
    qc=Counter()
    for row in survey:
        name=row["field_plant_binomial"]
        matches=by_name.get(name,set())
        status="EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH" if (
            len(matches)==1 and name not in ambiguous) else (
                "AMBIGUOUS_SOURCE_ACCEPTED_PLANT_MATCH" if matches
                else "NO_EXACT_239_SOURCE_ACCEPTED_PLANT_MATCH"
            )
        count=len(by_pid[next(iter(matches))]) if (
            status=="EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH") else None
        qc[status]+=1
        new.append({
            **row,"original239_accepted_HOSTS_match_status":status,
            "documented_original_239_butterfly_host_species_count":count,
            "source_broader_original_HOSTS_lepidopteran_consumers":len(other.get(name,set())),
        })
    return new,dict(qc)


def rank_average(vals):
    indexed=sorted(range(len(vals)),key=lambda k:vals[k])
    out=[0.0]*len(vals)
    i=0
    while i<len(indexed):
        j=i+1
        while j<len(indexed) and vals[indexed[j]]==vals[indexed[i]]:
            j+=1
        rank=(i+1+j)/2.0
        for k in indexed[i:j]:
            out[k]=rank
        i=j
    return np.array(out,float)


def spearman(a,b):
    if len(a)<3 or len(a)!=len(b):
        return None
    x=rank_average(a);y=rank_average(b)
    if np.std(x)<1e-12 or np.std(y)<1e-12:
        return None
    return float(np.corrcoef(x,y)[0,1])


def expit(z):
    z=np.asarray(z,float)
    return 1/(1+np.exp(-np.clip(z,-36.,36.)))


def solve_logit(X,successes,trials,*,ridge=RIDGE):
    X=np.asarray(X,float)
    y=np.asarray(successes,float)
    n=np.asarray(trials,float)
    if X.ndim!=2 or X.shape[0]!=len(y) or len(y)!=len(n) or not len(y):
        raise ValueError("Invalid binomial source observation group design")
    if X.shape[0]<X.shape[1]+2 or np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Field plant model rank insufficient")
    if np.any(n<=0) or np.any(y<0) or np.any(y>n):
        raise ValueError("Invalid original rounded source binomial counts")
    if np.sum(y)<3 or np.sum(n-y)<3:
        raise ValueError("Insufficient independent field detection outcomes")
    pen=np.diag([0]+[ridge]*(X.shape[1]-1))
    b=np.zeros(X.shape[1])
    global_p=np.clip(y.sum()/n.sum(),1e-7,1-1e-7)
    b[0]=math.log(global_p/(1-global_p))
    for _ in range(150):
        p=expit(X@b)
        w=n*p*(1-p)
        H=(X.T*w)@X+pen
        grad=X.T@(y-n*p)-pen@b
        step=np.linalg.solve(H,grad)
        def ll(beta):
            prob=expit(X@beta)
            return float(np.sum(y*np.log(np.maximum(prob,1e-15)) +
                         (n-y)*np.log(np.maximum(1-prob,1e-15)))
                         -.5*beta@pen@beta)
        cur=ll(b)
        scale=1.
        for _line in range(25):
            nxt=b+scale*step
            if ll(nxt)>=cur-1e-8:
                break
            scale*=.5
        else:
            raise ValueError("Grouped field logistic Newton step not improving")
        b=nxt
        if np.linalg.norm(scale*step)<1e-8:
            return b
    raise ValueError("Grouped binomial plant model failed convergence")


def model_records(rows,predictor,*,include_zero=False):
    available=[r for r in rows if r["original_USDA_origin"] in ("native","alien")
               and (include_zero or r["original239_accepted_HOSTS_match_status"]==
                    "EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH")]
    if predictor=="original239":
        vals=[float(r["documented_original_239_butterfly_host_species_count"]
                    or 0.) for r in available]
    elif predictor=="all_lepidoptera":
        vals=[float(r["source_broader_original_HOSTS_lepidopteran_consumers"]) for r in available]
    else:
        raise ValueError("Original linked source predictor not among predeclared sources")
    return available,np.log1p(np.array(vals,dtype=float))


def conditional_model(rows, predictor, *,include_zero=False, boot=BOOTSTRAPS):
    selected,vals=model_records(rows,predictor,include_zero=include_zero)
    n=len(selected)
    n_native=sum(r["original_USDA_origin"]=="native" for r in selected)
    n_alien=sum(r["original_USDA_origin"]=="alien" for r in selected)
    genera={r["field_plant_genus"] for r in selected}
    gate=(n>=15 and n_native>=3 and n_alien>=3 and len(genera)>=10
          and np.std(vals,ddof=1)>1e-12)
    info={
        "population":"all100_zero_as_documented_source_only" if include_zero
                      else "only_exact_original_239_source_linked_plants",
        "predictor":predictor,
        "plants":n,"native_plants":n_native,"alien_plants":n_alien,
        "plant_genera":len(genera),
        "source_documented_zero_is_not_ecological_absence":True,
        "descriptive_spearman_source_host_count_vs_caterpillar_pct":
            spearman([float(x) for x in vals],
                     [r["field_original_caterpillar_pct"] for r in selected]),
        "original_site_reported_positive_source_events":sum(
            r["field_reconstructed_positive_survey_events"] for r in selected),
        "original_survey_events":sum(r["field_original_n_surveys"] for r in selected),
        "prespecified_information_gate":(
            "SOURCE_LINKED_FIELD_MODEL_INFORMATION_OK" if gate else
            "INSUFFICIENT_EXACT_MATCH_NATIVE_ALIEN_SUPPORT"),
    }
    if not gate:
        return info
    mean=float(vals.mean());sd=float(vals.std(ddof=1))
    xx=(vals-mean)/sd
    native=np.array([r["original_USDA_origin"]=="native" for r in selected],float)
    successes=np.array([r["field_reconstructed_positive_survey_events"] for r in selected])
    counts=np.array([r["field_original_n_surveys"] for r in selected])
    X=np.column_stack([np.ones(n),native,xx])
    beta=solve_logit(X,successes,counts)
    beta_native=solve_logit(X[:,:2],successes,counts)
    info.update({
        "log_host_consumer_predictor_mean_across_plant_species":mean,
        "log_host_consumer_predictor_SD_across_plant_species":sd,
        "source_native_origin_only_logistic_beta":float(beta_native[1]),
        "source_host_count_per_SD_adjusted_log_odds_beta":float(beta[2]),
        "source_host_count_per_SD_adjusted_odds_ratio":float(math.exp(beta[2])),
        "source_native_vs_alien_logistic_beta_adjusted_for_host_count":float(beta[1]),
        "source_native_vs_alien_adjusted_odds_ratio":float(math.exp(beta[1])),
        "source_native_vs_host_count_correlation":float(np.corrcoef(native,xx)[0,1]),
    })
    # Out-of-genus validation, no raw survey pseudorep inferential p values.
    holdouts=Counter()
    for genus in sorted(genera):
        test=np.array([r["field_plant_genus"]==genus for r in selected])
        train=~test
        if len(set(native[train]))<2 or len(set(vals[train]))<3 or \
           np.linalg.matrix_rank(X[train])<3:
            holdouts["skipped_genera_due_source_rank_or_support"]+=1
            continue
        try:
            f=solve_logit(X[train],successes[train],counts[train])
            z=solve_logit(X[train,:2],successes[train],counts[train])
        except (ValueError,np.linalg.LinAlgError):
            holdouts["skipped_genera_due_fit_failure"]+=1
            continue
        pf=expit(X[test]@f)
        pz=expit(X[test,:2]@z)
        y=successes[test];q=counts[test]
        sf=float(np.sum(y*np.log(np.maximum(pf,1e-15))+
                        (q-y)*np.log(np.maximum(1-pf,1e-15))))
        sz=float(np.sum(y*np.log(np.maximum(pz,1e-15))+
                        (q-y)*np.log(np.maximum(1-pz,1e-15))))
        holdouts["test_log_score_full"]+=sf
        holdouts["test_log_score_origin_only"]+=sz
        holdouts["test_surveys"]+=int(q.sum())
        holdouts["validated_heldout_genera"]+=1
    info["leave_one_plant_genus_out"]={
        **dict(holdouts),
        "incremental_log_predictive_gain_per_1000_surveys":(
            1000*(holdouts["test_log_score_full"]-
                  holdouts["test_log_score_origin_only"])/holdouts["test_surveys"]
            if holdouts["test_surveys"] else None),
        "gain_positive_means_source_host_count_improves_out_of_genus_prediction":True,
    }
    # Source genus block bootstrap is descriptive, NOT 68k independent events.
    rng=np.random.default_rng(SEED+("all_lepidoptera"==predictor)+int(include_zero)*10)
    by_genus=defaultdict(list)
    for i,r in enumerate(selected):
        by_genus[r["field_plant_genus"]].append(i)
    names=sorted(by_genus)
    effect=[];failed=0
    for _ in range(boot):
        groupdraw=rng.integers(0,len(names),len(names))
        idx=np.array([i for k in groupdraw for i in by_genus[names[int(k)]]],
                     dtype=int)
        if len(set(native[idx]))<2 or np.linalg.matrix_rank(X[idx])<3:
            failed+=1
            continue
        try:
            b=solve_logit(X[idx],successes[idx],counts[idx])
            effect.append(float(b[2]))
        except (ValueError,np.linalg.LinAlgError):
            failed+=1
    info["plant_genus_bootstrap_499"]={
        "successful_draws":len(effect),"failed_draws":failed,
        "per_SD_log_odds_beta_percentile95":[float(x) for x in np.quantile(
            effect,[.025,.975])] if len(effect)>=.9*boot else None,
        "per_SD_OR_percentile95":[float(math.exp(x)) for x in np.quantile(
            effect,[.025,.975])] if len(effect)>=.9*boot else None,
        "generally_qualified":len(effect)>=.9*boot,
    }
    return info


def run(args):
    m=json.loads((args.author_source/"SOURCE_MANIFEST.json").read_text())
    if m.get("inputs_count")!=11:
        raise ValueError("Wrong original author-only source archive")
    for name,digest in m["sha256"].items():
        if sha256(args.author_source/name)!=digest:
            raise ValueError("Author-only original WCVP/HOSTS source member changed")
    focal_table=read_rows(args.author_source/"inputs/anthropogenic_species_metrics.csv")
    focal={r["species"] for r in focal_table}
    if len(focal_table)!=239 or len(focal)!=239:
        raise ValueError("Original 239 source butterfly identity not recovered")
    survey,paper=original_survey_plant_panel(args.external_plant_table)
    by_name,by_id,ambiguous,accepted=accepted_global_links(
        args.author_source/"inputs/original/insect_host_accepted.csv",focal)
    catalogue,sourceflow=broad_catalogue_host_links(
        args.original_hosts_resource,{r["field_plant_binomial"] for r in survey})
    rows,qc=join(survey,by_name,by_id,ambiguous,catalogue)
    breakdown={}
    for origin in ("native","alien","source_origin_unspecified"):
        rs=[r for r in rows if r["original_USDA_origin"]==origin]
        matched=[r for r in rs if r["original239_accepted_HOSTS_match_status"]==
                 "EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH"]
        totals=sum(r["field_original_n_surveys"] for r in rs)
        pos=sum(r["field_reconstructed_positive_survey_events"] for r in rs)
        breakdown[origin]={
            "field_survey_plant_species":len(rs),
            "field_survey_events":totals,
            "field_positive_larval_survey_events_rounded_from_published_pct":pos,
            "field_caterpillar_detection_pct_preselected_100plant_subset":
                100*pos/totals if totals else None,
            "exact_original239_host_linked_field_plants":len(matched),
            "source_239_host_linked_plants_frac":len(matched)/len(rs) if rs else None,
            "broader_HOSTS_exact_name_linked_field_plants":sum(
                r["source_broader_original_HOSTS_lepidopteran_consumers"]>0 for r in rs),
        }
    # The subset source-native vs alien replication is a PREVIOUSLY PUBLISHED
    # source result. Its inclusion is a cross-data QA, not an independent claim.
    result={
        "schema":SCHEMA,
        "status":"EXTERNAL_FIELD_LARVAE_VS_ORIGINAL_HOSTS_CATALOGUE_LINKAGE_AUDIT",
        "source":{
            "original_author_only_host_source_inner_SHA256":ORIGINAL_SOURCE_SHA,
            "original_HOSTS_commit":ORIGINAL_HOSTS_COMMIT,
            "original_HOSTS_original_source_SHA256":sha256(args.original_hosts_resource),
            "external_field_study_DOI":EXTERNAL_PAPER,
            "external_public_repo":"hurlbertlab/caterpillars-on-plants",
            "external_public_commit":EXTERNAL_COMMIT,
            "external_public_Table_S2":paper,
            "accepted_original_239":accepted,
            "source_original_HOSTS_complete_counts":sourceflow,
        },
        "original_100_field_plants_by_origin":breakdown,
        "source_link_quality":qc,
        "external_field_plant_rows":[
            {k:r[k] for k in (
                "field_plant_binomial","field_plant_family","field_plant_genus",
                "original_USDA_origin","field_original_n_surveys",
                "field_original_n_sites","field_original_caterpillar_pct",
                "field_reconstructed_positive_survey_events",
                "original239_accepted_HOSTS_match_status",
                "documented_original_239_butterfly_host_species_count",
                "source_broader_original_HOSTS_lepidopteran_consumers")}
            for r in rows
        ],
        "statistical_source_model":{
            "original_239_exact_host_linked_plant_test":
                conditional_model(rows,"original239"),
            "original_239_zero_documented_source_sensitivity":
                conditional_model(rows,"original239",include_zero=True),
            "all_original_HOSTS_lepidopteran_broad_coverage_benchmark":
                conditional_model(rows,"all_lepidoptera",include_zero=True),
        },
        "biological_red_lines":[
            "External field larval outcome contains Lepidoptera moths and butterflies; original 239 butterfly source host edges are not the field species identities.",
            "Original 100 plants are field-survey-preselected >40 events, not all species or all 68392 field surveys; full-sample published 8.1% vs 2.2% is a separate denominator.",
            "The native/alien larval detection comparison is already published by Whitener & Hurlbert 2026 and must not be claimed novel here.",
            "Zero exact original HOSTS plant name matches is not proof that no butterfly could or does use the plant; accepted plant synonyms and source work effort matter.",
            "Original field study USDA plant status is not the globally averaged WCVP country source plant native/introduced classification.",
            "Branch surveys have site/year/plant dependence, and binomial event scores are not independent biological experiments.",
            "Published whole-study field caterpillar occurrence is not individual larval survival, adult recruitment, invasive butterfly establishment, or interspecific competition.",
            "Field plant sample does not supply contemporary overlap of multiple identified butterfly larvae on exact same plant individual.",
            "Only source derived 100 public plant names, counts and aggregate original HOSTS links are exported; no author-only full accepted host matrix."
        ],
        "computed_utc":datetime.now(timezone.utc).isoformat(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print(json.dumps({
        "status":result["status"],"breakdown":breakdown,"join":qc,
        "primary":result["statistical_source_model"]["original_239_exact_host_linked_plant_test"],
        "source_documented_zero":result["statistical_source_model"]["original_239_zero_documented_source_sensitivity"],
        "all_lepidoptera":result["statistical_source_model"]["all_original_HOSTS_lepidopteran_broad_coverage_benchmark"],
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--original-hosts-resource",type=Path,required=True)
    p.add_argument("--external-plant-table",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
