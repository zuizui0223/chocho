#!/usr/bin/env python3
"""Stream original pinned HOSTS resource.csv, resolve only the frozen 239
butterflies' host binomials against WCVP v13 accepted species IDs.

Rationale: original generic R sidecar resolves ALL insects in HOSTS before
taking 239 butterflies, an expensive and unnecessary full-catalogue operation.
This is a source-compatible restricted parser, with the same unambiguous
accepted-ID rules and no island-outcome-based host selection.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import defaultdict,Counter
from pathlib import Path
import pandas as pd
import pyreadr

def ident(v):
    t=str(v).strip()
    if re.fullmatch(r"[0-9]+(?:\.0+)?",t):return t.split(".")[0]
    return t

def stream_hosts(path,allowed):
    edges=set();hits=Counter();rows=0
    with Path(path).open(newline="",encoding="utf-8") as f:
        it=csv.DictReader(f)
        req={"Insect Genus","Insect Species","Hostplant Genus","Hostplant Species"}
        if not req.issubset(it.fieldnames or []):
            raise ValueError("Pinned HOSTS source schema missing insect/host binomial columns")
        for r in it:
            rows+=1
            ins=" ".join(x.strip() for x in (r["Insect Genus"] or "",r["Insect Species"] or "")).strip()
            if ins not in allowed:continue
            genus=(r["Hostplant Genus"] or "").strip()
            epithet=(r["Hostplant Species"] or "").strip()
            if not genus or not epithet:continue
            plant=genus+" "+epithet
            edges.add((ins,plant))
            hits[ins]+=1
    return edges,{"HOSTS_resource_rows_read":rows,
                 "original_239_panel_butterflies_with_host_names":len(hits),
                 "unique_butterfly_host_binomial_edges":len(edges)}

def accepted_mapping(path,requested):
    v=pyreadr.read_r(str(path))["wcvp_names"]
    need={"plant_name_id","taxon_rank","taxon_status","taxon_name","accepted_plant_name_id","family"}
    if not need.issubset(v.columns):raise ValueError("WCVP accepted species schema changed")
    # Keep exact printed host binomials; alias/synonym must resolve to ONE
    # species-level accepted WCVP ID, not genus-level or fuzzy.
    exact=v[(v.taxon_name.isin(requested))&
            (v.taxon_rank.astype(str).str.strip().str.casefold()=="species")].copy()
    status=exact.taxon_status.astype(str).str.casefold()
    exact["plant_id_accepted"]=exact.accepted_plant_name_id.map(ident)
    exact.loc[(status=="accepted") & exact.plant_id_accepted.isin(("","nan","None","NA")),
              "plant_id_accepted"]=exact.loc[(status=="accepted") &
                  exact.plant_id_accepted.isin(("","nan","None","NA")),"plant_name_id"].map(ident)
    exact=exact[~exact.plant_id_accepted.isin(("","nan","None","NA"))]
    choices=defaultdict(set)
    for r in exact.itertuples(index=False):
        choices[str(r.taxon_name)].add(str(r.plant_id_accepted))
    all_names=v.assign(plant_name_id=v.plant_name_id.map(ident))
    table=all_names.drop_duplicates("plant_name_id").set_index("plant_name_id")
    resolved={};ambig=0
    for name,ids in choices.items():
        if len(ids)!=1:
            ambig+=1;continue
        idx=next(iter(ids))
        if idx not in table.index:continue
        row=table.loc[idx]
        if isinstance(row,pd.DataFrame):raise ValueError("WCVP accepted ID has duplicate taxonomy")
        resolved[name]={"plant_id":idx,"accepted_name":str(row.taxon_name),
                        "family":str(row.family)}
    return resolved,{"WCVP_accepted_name_rows":len(v),
                     "exact_matched_host_binomials":len(resolved),
                     "ambiguous_host_binomials":ambig}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("hosts","wcvp-names","panel","out","audit"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    panel=pd.read_csv(a.panel)
    if len(panel)!=239 or panel.species.nunique()!=239:
        raise ValueError("Frozen scientific 239-butterfly panel differs")
    edges,src=stream_hosts(a.hosts,set(panel.species))
    if len(edges)<400:
        raise ValueError("HOSTS panel sampled too few source host relations")
    resolutions,tax=accepted_mapping(a.wcvp_names,{x[1] for x in edges})
    rows=[]
    for ins,host in sorted(edges):
        if host not in resolutions:continue
        row=resolutions[host]
        rows.append({"insect_species":ins,"input_host_name":host,
                     "accepted_plant_name_id":row["plant_id"],
                     "accepted_name":row["accepted_name"],"family":row["family"]})
    side=pd.DataFrame(rows).drop_duplicates(["insect_species","accepted_plant_name_id"])
    if side.insect_species.nunique()<210 or side.accepted_plant_name_id.nunique()<200:
        raise ValueError("Restricted HOSTS×WCVP source matching unexpectedly low")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    side.to_csv(a.out,index=False)
    audit={"schema":"chocho_239_butterfly_pinned_HOSTS_WCVP_full_diet_exact_source_v1",
           "source_HOSTS_pinned":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
           "source_WCVP_pinned":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
           "WCVP_names_sha256":hashlib.sha256(a.wcvp_names.read_bytes()).hexdigest(),
           "HOSTS_source_sha256":hashlib.sha256(a.hosts.read_bytes()).hexdigest(),
           "butterfly_panel_sha256":hashlib.sha256(a.panel.read_bytes()).hexdigest(),
           **src,**tax,
           "resolved_source_butterfly_host_edges":len(side),
           "butterflies_with_accepted_hosts":int(side.insect_species.nunique()),
           "distinct_accepted_plant_ids":int(side.accepted_plant_name_id.nunique()),
           "host_diet_defined_without_source_island_gains":True,
           "unmatched_host_aliases_are_unavailable_not_zero_traits":True}
    a.audit.parent.mkdir(parents=True,exist_ok=True)
    a.audit.write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("WCVP_HOSTS_FULL239_STREAMED_SIDE_CAR",json.dumps(audit,sort_keys=True),flush=True)

if __name__=="__main__":main()
