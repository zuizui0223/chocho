#!/usr/bin/env Rscript
# Reconstruct ORIGINAL-GEB-style resource envelopes for a source-selected panel
# of independently classified island-endemic butterflies. Retains original
# frozen plant taxonomy/WGSRPD3 grain, but bypasses LepTraits and 339/239
# selection. NOT an estimator for all 18 island endemics.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript analyze_independent_endemic_host_only_WCVP_v13_20261118.R WCVP_REPO SOURCE_EDGES.csv OUTPUT")
rrepo<-args[[1L]]
panel_file<-args[[2L]]
out<-args[[3L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
env<-new.env(parent=emptyenv())
load(file.path(rrepo,"data","wcvp_names.rda"),envir=env)
load(file.path(rrepo,"data","wcvp_distributions.rda"),envir=env)
n<-get("wcvp_names",envir=env)
d<-get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","taxon_name","taxon_status","taxon_rank") %in% names(n)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(d)))
edge<-read.csv(panel_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c(""))
stopifnot(nrow(edge)==13L,length(unique(edge$butterfly))==8L,
          length(unique(edge$plant_species))==11L)
stopifnot(all(edge$source_region_TDWG3 %in% c("CNY","AZO","MDR")))
stopifnot(all(edge$in_original_GEB239==0L))
tier_allow<-c("larval_reported","oviposition_only","captivity_only","nonlocal_database")
stopifnot(all(edge$evidence_tier %in% tier_allow))
n$plant_name_id<-as.character(n$plant_name_id)
n$taxon_name<-as.character(n$taxon_name)
n$taxon_rank<-tolower(trimws(as.character(n$taxon_rank)))
n$taxon_status<-tolower(trimws(as.character(n$taxon_status)))
d$plant_name_id<-as.character(d$plant_name_id)
accepted<-n[!is.na(n$taxon_name) & n$taxon_rank=="species" &
              n$taxon_status=="accepted",c("plant_name_id","taxon_name"),drop=FALSE]
valid_d<-d[!is.na(d$area_code_l3) &
  as.character(d$extinct)=="0" &
  as.character(d$location_doubtful)=="0" &
  as.character(d$introduced) %in% c("0","1"),
  c("plant_name_id","area_code_l3","introduced"),drop=FALSE]
valid_d$area_code_l3<-as.character(valid_d$area_code_l3)
valid_d$introduced<-as.character(valid_d$introduced)
valid_d<-unique(valid_d)
taxon_result<-lapply(sort(unique(edge$plant_species)),function(plant){
   matches<-accepted[accepted$taxon_name==plant,,drop=FALSE]
   id<-if(nrow(matches)==1L)matches$plant_name_id[[1L]] else ""
   dv<-valid_d[valid_d$plant_name_id==id & nzchar(id),,drop=FALSE]
   native<-sort(unique(dv$area_code_l3[dv$introduced=="0"]))
   introduced<-sort(unique(dv$area_code_l3[dv$introduced=="1"]))
   conflict<-intersect(native,introduced)
   status<-if(nrow(matches)==0L)"no_accepted_name" else if(nrow(matches)>1L) "multiple_accepted_names" else
           if(length(conflict))"conflicting_geographic_origin" else if(!nrow(dv))"no_valid_distribution" else "resolved"
   list(plant=plant,accepted_id=id,matching_accepted_species=nrow(matches),
        native=native,introduced=introduced,conflict=conflict,
        status=status)
})
names(taxon_result)<-vapply(taxon_result,function(z)z$plant,character(1))
resolved<-lapply(seq_len(nrow(edge)),function(i){
  row<-edge[i,,drop=FALSE]
  t<-taxon_result[[row$plant_species]]
  src<-row$source_region_TDWG3
  cbind(row,
    accepted_WCVP_ID=t$accepted_id,
    accepted_plant_status=t$status,
    source_region_native=as.integer(t$status=="resolved" && src %in% t$native),
    source_region_introduced=as.integer(t$status=="resolved" && src %in% t$introduced),
    source_region_uncertain=as.integer(t$status!="resolved" || (!src %in% c(t$native,t$introduced))),
    global_native_regions=paste(t$native,collapse=";"),
    global_introduced_regions=paste(t$introduced,collapse=";"),
    stringsAsFactors=FALSE)
})
resolved<-do.call(rbind,resolved)
if(nrow(resolved)!=13L)stop("Source edges lost before taxonomic resolution")
if(any(resolved$accepted_plant_status!="resolved")) {
  cat("WARN_UNRESOLVED_SOURCE_PLANTS=",paste(unique(resolved$plant_species[resolved$accepted_plant_status!="resolved"]),collapse=";"),"\n",sep="")
}
if(any(resolved$source_region_uncertain==1L))
  cat("WARN_SOURCE_REGION_PLANT_STATUS_UNRESOLVED=",paste(unique(resolved$plant_species[resolved$source_region_uncertain==1L]),collapse=";"),"\n",sep="")
write.csv(resolved,file.path(out,"source_13_butterfly_plant_species_WCVPv13_origin.csv"),row.names=FALSE,na="")
# Modes are nested evidence source classes and explicitly show how
# reliance on HOSTS-only or oviposition records affects results.
modes<-list(
  "wild_and_documented_larval"=c("larval_reported"),
  "plus_published_oviposition_and_captivity"=c("larval_reported","oviposition_only","captivity_only"),
  "all_named_sources_including_NHM_nonlocal"=tier_allow,
  "NHM_original_nonlocal_only"=c("nonlocal_database"))
outcomes<-list()
idb<-sort(unique(edge$butterfly))
for(mode in names(modes)) {
  for(butterfly in idb) {
    z<-resolved[resolved$butterfly==butterfly & resolved$evidence_tier %in% modes[[mode]],,drop=FALSE]
    source_region<-unique(edge$source_region_TDWG3[edge$butterfly==butterfly])
    if(length(source_region)!=1L)stop("Butterfly assigned more than one archipelago source")
    host_ids<-sort(unique(z$accepted_WCVP_ID[z$accepted_plant_status=="resolved"]))
    native<-character(0); contemporary<-character(0)
    if(length(host_ids)){
      for(id in host_ids){
        pp<-valid_d[valid_d$plant_name_id==id,,drop=FALSE]
        nat<-sort(unique(pp$area_code_l3[pp$introduced=="0"]))
        intro<-sort(unique(pp$area_code_l3[pp$introduced=="1"]))
        conflict<-intersect(nat,intro)
        if(length(conflict))stop("WCVP geographic source conflict was not quarantined")
        native<-union(native,nat)
        contemporary<-union(contemporary,c(nat,intro))
      }
    }
    native<-sort(native);contemporary<-sort(contemporary)
    added<-sort(setdiff(contemporary,native))
    if(length(setdiff(native,contemporary)))stop("Original native ranges not subset contemporary")
    no_data<-nrow(z)==0L
    source_native<-sum(z$source_region_native)
    source_intro<-sum(z$source_region_introduced)
    unresolved<-sum(z$source_region_uncertain)
    # If no source rows, counts are MISSING, not zero. Never infer
    # 'no host' from no source in an evidence-tier subset.
    outcomes[[length(outcomes)+1L]]<-data.frame(
      source_butterfly=butterfly,
      source_region_TDWG3=source_region,
      mode=mode,
      source_host_species_n=nrow(z),
      unresolved_host_species_or_local_status_n=unresolved,
      source_resolved_hosts_native_n=source_native,
      source_resolved_hosts_introduced_n=source_intro,
      source_region_native_host_species_present=as.integer(source_native>0L),
      source_region_any_host_species_present=as.integer((source_native+source_intro)>0L),
      source_region_added_host_species_when_native_already_present=as.integer(source_native>0L && source_intro>0L),
      source_region_newly_host_present_after_intro=as.integer(source_native==0L && source_intro>0L && unresolved==0L && !no_data),
      source_region_newly_host_present_identifiable=as.integer(!no_data && unresolved==0L),
      potential_global_native_regions=if(no_data)NA_integer_ else length(native),
      potential_global_contemporary_regions=if(no_data)NA_integer_ else length(contemporary),
      potential_global_added_regions=if(no_data)NA_integer_ else length(added),
      potential_global_contemporary_over_native_ratio=if(no_data||!length(native)) NA_real_ else length(contemporary)/length(native),
      host_species_without_source_record=as.integer(no_data),
      geographic_unit="WCVP WGSRPD3 botanical region, NOT individual island or butterfly occurrence",
      stringsAsFactors=FALSE)
  }
}
output<-do.call(rbind,outcomes)
stopifnot(nrow(output)==length(modes)*length(idb))
# A source-recorded native+introduced host can augment resources but
# not create the first host-present CNY region for P. cheiranthi.
p<-output[output$source_butterfly=="Pieris cheiranthi" &
          output$mode=="all_named_sources_including_NHM_nonlocal",,drop=FALSE]
if(nrow(p)!=1L || p$source_resolved_hosts_native_n<1L || p$source_resolved_hosts_introduced_n<1L ||
   p$source_region_newly_host_present_after_intro!=0L)
  stop("Focal Pieris known portfolio expansion was misread as a new host-present island-region")
# The two source-only NHM Azores butterflies share Festuca jubata
# and the original zero introduced-range experiment.
for(b in c("Hipparchia azorina","Hipparchia miguelensis")) {
  rr<-output[output$source_butterfly==b &
    output$mode=="all_named_sources_including_NHM_nonlocal",,drop=FALSE]
  if(nrow(rr)!=1L || rr$potential_global_added_regions!=0L ||
     rr$potential_global_native_regions!=2L ||
     rr$potential_global_contemporary_regions!=2L)
    stop("Original fixed Festuca source calculation changed for "+b)
}
write.csv(output,file.path(out,"eight_independent_endemics_host_only_resource_by_evidence.csv"),row.names=FALSE,na="")
cat("INDEPENDENT_EIGHT_END_EMICS_WCVP_SOURCE_ONLY_RESOURCE_COMPLETED\n")
for(i in seq_len(nrow(output))){
  r<-output[i,,drop=FALSE]
  if(r$mode=="all_named_sources_including_NHM_nonlocal")
    cat(sprintf("END_EMIC_HOST_RESOURCE=%s,source_region=%s,native_hosts=%d,introduced_hosts=%d,global_native=%s,global_contemporary=%s,global_added=%s,ratio=%s,uncertain=%d\n",
       r$source_butterfly,r$source_region_TDWG3,r$source_resolved_hosts_native_n,r$source_resolved_hosts_introduced_n,
       r$potential_global_native_regions,r$potential_global_contemporary_regions,
       r$potential_global_added_regions,r$potential_global_contemporary_over_native_ratio,
       r$unresolved_host_species_or_local_status_n))
}
cat("SCIENTIFIC_BOUNDARY: source-selected host-only potential geography, not representative of 18 island endemics; local co-occurrence and butterfly colonization unmeasured\n")
