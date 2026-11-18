#!/usr/bin/env Rscript
# Check original GEB WCVP v13 species-level Canary native/introduced status
# of three named wild/captive larval host plants from Acosta-Fernandez 2010.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript audit_Euchloe_field_lab_WCVP_v13_20261118.R WCVP_REPO SOURCE_CSV OUT_DIR")
source_repo<-args[[1L]]
datafile<-args[[2L]]
outdir<-args[[3L]]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
env<-new.env(parent=emptyenv())
load(file.path(source_repo,"data","wcvp_names.rda"),envir=env)
load(file.path(source_repo,"data","wcvp_distributions.rda"),envir=env)
names<-get("wcvp_names",envir=env)
dist<-get("wcvp_distributions",envir=env)
source<-read.csv(datafile,stringsAsFactors=FALSE,check.names=FALSE,na.strings="")
stopifnot(nrow(source)==4L)
plants<-sort(unique(source$offered_or_field_larval_host))
stopifnot(identical(plants,c("Hirschfeldia incana","Sisymbrium erysimoides","Sisymbrium irio")))
names$taxon_rank<-tolower(trimws(as.character(names$taxon_rank)))
names$taxon_status<-tolower(trimws(as.character(names$taxon_status)))
names$plant_name_id<-as.character(names$plant_name_id)
dist$plant_name_id<-as.character(dist$plant_name_id)
records<-lapply(plants,function(nm){
  m<-names[
    as.character(names$taxon_name)==nm & names$taxon_rank=="species" &
    names$taxon_status=="accepted",
    c("plant_name_id","taxon_name"),drop=FALSE]
  id<-if(nrow(m)==1L)m$plant_name_id[[1L]] else NA_character_
  d<-if(!is.na(id))dist[
    dist$plant_name_id==id &
    as.character(dist$area_code_l3)=="CNY" &
    as.character(dist$extinct)=="0" &
    as.character(dist$location_doubtful)=="0" &
    as.character(dist$introduced) %in% c("0","1"),
    ,drop=FALSE] else dist[FALSE,,drop=FALSE]
  states<-sort(unique(as.character(d$introduced)))
  origin<-if(nrow(m)!=1L)"unresolved_accepted_taxon" else
          if(identical(states,"0"))"native_CNY" else
          if(identical(states,"1"))"introduced_CNY" else
          if(length(states)==0L)"no_source_CNY_distribution" else
          "conflicting_CNY_status"
  data.frame(plant_species=nm,accepted_WCVP_plant_id=if(is.na(id))"" else id,
             n_accepted_species_matches=nrow(m),WCVP_v13_CNY_origin=origin,
             stringsAsFactors=FALSE)
})
rows<-do.call(rbind,records)
write.csv(rows,file.path(outdir,"Euchloe_three_host_plants_WCVPv13_CNY.csv"),
          row.names=FALSE,na="")
for(i in seq_len(nrow(rows)))
  cat(sprintf("WCVP2015_EUCHLOE_PLANT=%s,CNY_ORIGIN=%s\n",rows$plant_species[i],rows$WCVP_v13_CNY_origin[i]))
stopifnot(nrow(rows)==3L,all(rows$n_accepted_species_matches==1L))
if(!all(rows$WCVP_v13_CNY_origin=="native_CNY"))
   stop("Original pinned WCVP disagrees on native/introduced classification, do not promote paper claims")
cat("PINNED_EUCHLOE_TWO_ISLAND_HOST_SPECIES_SOURCE_CNY_VALIDATED\n")
