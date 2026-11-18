#!/usr/bin/env Rscript
# Literature-contested Festuca jubata versus Festuca petraea azorina-group
# larval host names, original WCVP v13, native vs introduced distribution.
# This is a hypothetical source-declared HOST alternative sensitivity, not
# a new observed butterfly-host pair and not diet breadth imputation.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript <original_WCVP_REPO> <out_dir>")
wcvp_repo<-args[[1L]]
out<-args[[2L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
env<-new.env(parent=emptyenv())
load(file.path(wcvp_repo,"data","wcvp_names.rda"),envir=env)
load(file.path(wcvp_repo,"data","wcvp_distributions.rda"),envir=env)
nm<-get("wcvp_names",envir=env)
ds<-get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","taxon_name","taxon_rank","taxon_status")%in%names(nm)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")%in%names(ds)))
source_candidates<-c("Festuca jubata","Festuca petraea")
nm$plant_name_id<-as.character(nm$plant_name_id)
nm$taxon_name<-as.character(nm$taxon_name)
nm$taxon_rank<-tolower(trimws(as.character(nm$taxon_rank)))
nm$taxon_status<-tolower(trimws(as.character(nm$taxon_status)))
ds$plant_name_id<-as.character(ds$plant_name_id)
candidates<-lapply(source_candidates,function(plant){
  matched<-nm[nm$taxon_name==plant&nm$taxon_rank=="species"&
              nm$taxon_status=="accepted",,drop=FALSE]
  if(nrow(matched)!=1L)stop(paste("Original WCVP v13 unique accepted plant does not exist:",plant))
  id<-matched$plant_name_id[[1L]]
  raw<-ds[ds$plant_name_id==id &
           as.character(ds$extinct)=="0" &
           as.character(ds$location_doubtful)=="0" &
           as.character(ds$introduced)%in%c("0","1") &
           !is.na(ds$area_code_l3)&nzchar(as.character(ds$area_code_l3)),
           c("area_code_l3","introduced"),drop=FALSE]
  native<-sort(unique(as.character(raw$area_code_l3)[as.character(raw$introduced)=="0"]))
  introduced<-sort(unique(as.character(raw$area_code_l3)[as.character(raw$introduced)=="1"]))
  if(!"AZO"%in%native)stop(paste("Source regional native AZO absent:",plant))
  list(plant=plant,id=id,native=native,introduced=introduced)
})
names(candidates)<-source_candidates
rows<-lapply(candidates,function(item){
  data.frame(plant=item$plant,original_WCVP_accepted_id=item$id,
             original_native_area_codes=paste(item$native,collapse=";"),
             original_introduced_area_codes=paste(item$introduced,collapse=";"),
             source_native_AZO=1L,native_unit_n=length(item$native),
             introduced_region_n=length(setdiff(item$introduced,item$native)),
             stringsAsFactors=FALSE)
})
rows<-do.call(rbind,rows)
scenario_specs<-list(
    only_original_NHM_F_jubata=c("Festuca jubata"),
    only_historical_Fuchs_F_petraea=c("Festuca petraea"),
    both_reported_Festuca_candidates=c("Festuca jubata","Festuca petraea")
)
scenario<-lapply(names(scenario_specs),function(sc){
  names<-scenario_specs[[sc]]
  native<-sort(unique(unlist(lapply(names,function(n)candidates[[n]]$native),use.names=FALSE)))
  introduced<-sort(unique(unlist(lapply(names,function(n)candidates[[n]]$introduced),use.names=FALSE)))
  contemporary<-sort(unique(c(native,introduced)))
  added<-setdiff(contemporary,native)
  data.frame(
    scenario=sc,
    supported_host_source=if(sc=="only_original_NHM_F_jubata")"original_NHM_HOSTS" else
      if(sc=="only_historical_Fuchs_F_petraea")"historical_Fuchs_1993_review_alternative" else
      "unverified_union_of_literature_alternatives",
    assumed_host_species=paste(names,collapse=";"),
    native_WGSRPD3_units=length(native),
    contemporary_WGSRPD3_units=length(contemporary),
    introduced_added_WGSRPD3_units=length(added),
    ratio=if(length(native))length(contemporary)/length(native) else NA_real_,
    accepted_WCVP_AZO_native=as.integer("AZO"%in%native),
    field_evidence_that_each_butterfly_uses_all_assumed_hosts=0L,
    stringsAsFactors=FALSE)
})
scenario<-do.call(rbind,scenario)
write.csv(rows,file.path(out,"Festuca_two_literature_host_candidates_WCVPv13.csv"),row.names=FALSE,na="")
write.csv(scenario,file.path(out,"Festuca_source_host_sensitivity_three_scenarios_WCVPv13.csv"),row.names=FALSE,na="")
cat("FESTUCA_AZORES_2018_HOST_CONCEPT_SENSITIVITY_SOURCE_COMPLETED\n")
for(i in seq_len(nrow(rows))){
  r<-rows[i,,drop=FALSE]
  cat(sprintf("FESTUCA_CANDIDATE=%s,native=%s,introduced=%s\n",
    r$plant,r$original_native_area_codes,r$original_introduced_area_codes))
}
for(i in seq_len(nrow(scenario))){
  r<-scenario[i,,drop=FALSE]
  cat(sprintf("FESTUCA_SCENARIO=%s,native_n=%d,contemporary_n=%d,added_n=%d\n",
    r$scenario,r$native_WGSRPD3_units,r$contemporary_WGSRPD3_units,
    r$introduced_added_WGSRPD3_units))
}
