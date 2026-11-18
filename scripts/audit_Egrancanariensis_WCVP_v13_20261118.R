#!/usr/bin/env Rscript
# Verify new endemic butterfly larval host exact species in original WCVP v13.
# The Canary region CNY is island-group-level, not the site of actual feeding.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript audit_Egrancanariensis_WCVP_v13.R WCVP_REPO SOURCE_CSV OUTPUT_DIR")
wcvp_repo<-args[[1L]]
source_file<-args[[2L]]
output<-args[[3L]]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
env <- new.env(parent=emptyenv())
load(file.path(wcvp_repo,"data","wcvp_names.rda"),envir=env)
load(file.path(wcvp_repo,"data","wcvp_distributions.rda"),envir=env)
n<-get("wcvp_names",envir=env)
d<-get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","taxon_name","taxon_rank","taxon_status") %in% names(n)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(d)))
input<-read.csv(source_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings="")
stopifnot(nrow(input)==4L,all(c("butterfly","host_plant","plant_name_resolution","evidence_class") %in% names(input)))
to_verify<-sort(unique(input$host_plant[input$plant_name_resolution=="species" &
                  input$butterfly=="Euchloe grancanariensis"]))
stopifnot(identical(to_verify,c("Descurainia preauxiana","Hirschfeldia incana")))
n$plant_name_id<-as.character(n$plant_name_id)
n$taxon_name<-as.character(n$taxon_name)
n$taxon_status<-tolower(trimws(as.character(n$taxon_status)))
n$taxon_rank<-tolower(trimws(as.character(n$taxon_rank)))
d$plant_name_id<-as.character(d$plant_name_id)
result<-lapply(to_verify,function(plant){
    m<-n[!is.na(n$taxon_name)&n$taxon_name==plant&
        n$taxon_rank=="species"&n$taxon_status=="accepted",,drop=FALSE]
    key<-if(nrow(m)==1L)m$plant_name_id[[1L]] else NA_character_
    dat<-if(!is.na(key))d[d$plant_name_id==key&
         as.character(d$area_code_l3)=="CNY"&
         as.character(d$extinct)=="0"&
         as.character(d$location_doubtful)=="0"&
         as.character(d$introduced)%in%c("0","1"),,drop=FALSE]else d[FALSE,,drop=FALSE]
    origin<-sort(unique(as.character(dat$introduced)))
    status<-if(nrow(m)!=1L)"not_unique_accepted_species" else
       if(identical(origin,"0"))"native_CNY" else
       if(identical(origin,"1"))"introduced_CNY" else
       if(length(origin)==0L)"no_CNY_distribution" else "CNY_origin_conflict"
    data.frame(
      butterfly="Euchloe grancanariensis", plant=plant,
      source_region="Canary Islands CNY", original_WCVP_accepted_ID=if(is.na(key))"" else key,
      exact_accepted_species_n=nrow(m), source_v13_origin=status,
      direct_wild_larva_on_species=as.integer(plant=="Hirschfeldia incana"),
      direct_oviposition_on_species=as.integer(plant=="Descurainia preauxiana"),
      confirmed_larval_development_to_adult=0L,
      stringsAsFactors=FALSE)
})
result<-do.call(rbind,result)
write.csv(result,file.path(output,"Egrancanariensis_two_hosts_original_WCVP_v13_CNY.csv"),row.names=FALSE)
cat("EUCHLOE_GRANCANARIENSIS_SOURCE_WCVP_SPECIES_CROSSWALK_COMPLETE\n")
for(i in seq_len(nrow(result))){
 r<-result[i,,drop=FALSE]
 cat(sprintf("EG_SOURCE_HOST=%s,original_WCVP_CNY=%s,accepted_matches=%d\n",
    r$plant,r$source_v13_origin,r$exact_accepted_species_n))
}
if(!all(result$source_v13_origin=="native_CNY"))stop(
 "Source frozen WCVP does not label both source-named plant species as native CNY; require manual source audit")
