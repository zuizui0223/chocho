#!/usr/bin/env Rscript
# Original WCVP v13 single-plant host-only resource geography for two independently
# classified Azores endemic Hipparchia. Critical: original 339 excludes BOTH taxa.
# This is a taxonomic feasibility/sensitivity probe, NOT a GEB239 extension.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript <original_wcvp_repo> <pinned_hosts_resource_csv> <output_dir>")
wcvp_repo <- args[[1L]]
hosts_file <- args[[2L]]
out <- args[[3L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
env <- new.env(parent=emptyenv())
load(file.path(wcvp_repo,"data","wcvp_names.rda"),envir=env)
load(file.path(wcvp_repo,"data","wcvp_distributions.rda"),envir=env)
nm <- get("wcvp_names",envir=env)
di <- get("wcvp_distributions",envir=env)
stopifnot(all(c("taxon_name","taxon_rank","taxon_status","plant_name_id") %in% names(nm)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(di)))
h <- read.csv(hosts_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings="")
stopifnot(nrow(h)==140485L)
stopifnot(all(c("Insect Genus","Insect Species","Insect Family","Hostplant Genus","Hostplant Species")%in%names(h)))
i <- h[!is.na(h[["Insect Genus"]]) & h[["Insect Genus"]]=="Hipparchia" &
  !is.na(h[["Insect Species"]]) & h[["Insect Species"]] %in% c("azorina","miguelensis"),,drop=FALSE]
i <- i[!is.na(i[["Insect Family"]]) & i[["Insect Family"]]=="Nymphalidae",,drop=FALSE]
i$butterfly <- paste(i[["Insect Genus"]],i[["Insect Species"]])
i$plant <- paste(trimws(i[["Hostplant Genus"]]),trimws(i[["Hostplant Species"]]))
for (b in c("Hipparchia azorina","Hipparchia miguelensis")) {
  source_rows <- i[i$butterfly==b,,drop=FALSE]
  stopifnot(nrow(source_rows)==1L)
  if(!identical(source_rows$plant,"Festuca jubata"))stop(paste("Original HOSTS named species changed",b))
}
# NHM original raw association is NOT locality-specific observed larval feeding.
nm$plant_name_id <- as.character(nm$plant_name_id)
nm$taxon_name <- as.character(nm$taxon_name)
nm$taxon_rank <- tolower(trimws(as.character(nm$taxon_rank)))
nm$taxon_status <- tolower(trimws(as.character(nm$taxon_status)))
plants <- nm[nm$taxon_name=="Festuca jubata" &
  nm$taxon_rank=="species" & nm$taxon_status=="accepted",,drop=FALSE]
if(nrow(plants)!=1L)stop("Festuca jubata cannot be uniquely matched to accepted WCVP v13 species")
accepted_id <- plants$plant_name_id[[1L]]
di$plant_name_id <- as.character(di$plant_name_id)
d <- di[di$plant_name_id==accepted_id,,drop=FALSE]
# The distinction is native-only versus contemporary (=native + introduced).
d <- d[!is.na(d$area_code_l3) & nzchar(as.character(d$area_code_l3)) &
  as.character(d$extinct)=="0" & as.character(d$location_doubtful)=="0" &
  as.character(d$introduced)%in%c("0","1"),,drop=FALSE]
d <- unique(data.frame(
  original_WCVP_accepted_plant_name_id=d$plant_name_id,
  accepted_plant_species="Festuca jubata",
  area_code_l3=as.character(d$area_code_l3),
  introduced=as.character(d$introduced),
  stringsAsFactors=FALSE))
if(nrow(d)==0L)stop("No verified native or introduced WCVP distribution for Festuca jubata")
origin_by_area <- split(d$introduced,d$area_code_l3)
region <- do.call(rbind,lapply(sort(names(origin_by_area)),function(area){
  status<-sort(unique(origin_by_area[[area]]))
  data.frame(
    host_plant="Festuca jubata",
    accepted_WCVP_plant_id=accepted_id,
    area_code_l3=area,
    has_native=as.integer("0"%in%status),
    has_introduced=as.integer("1"%in%status),
    status=if(identical(status,"0"))"native" else if(identical(status,"1"))"introduced" else "conflicting_same_region",
    stringsAsFactors=FALSE)
}))
stopifnot("AZO"%in%region$area_code_l3)
azores<-region[region$area_code_l3=="AZO",,drop=FALSE]
if(nrow(azores)!=1L || azores$has_native[[1L]]!=1L)stop("Azores origin is not reliably native; hold")
native<-sort(region$area_code_l3[region$has_native==1L])
contemporary<-sort(region$area_code_l3[region$has_native==1L|region$has_introduced==1L])
added<-setdiff(contemporary,native)
if(length(contemporary)<length(native))stop("Contemporary source cannot remove native regions")
result<-do.call(rbind,lapply(sort(unique(i$butterfly)),function(b){
  data.frame(
    independently_classified_Azores_endemic_butterfly=b,
    original_NHM_HOSTS_plant="Festuca jubata",
    original_source_HOSTS_host_species_records=nrow(i[i$butterfly==b,,drop=FALSE]),
    source_WCVP_accepted_plant_id=accepted_id,
    native_WGSRPD3_units=length(native),
    contemporary_WGSRPD3_units=length(contemporary),
    introduced_added_WGSRPD3_units=length(added),
    contemporary_over_native_ratio=if(length(native))length(contemporary)/length(native) else NA_real_,
    native_Azores_plant_source=as.integer("AZO"%in%native),
    original_GEB339_source_panel_included=0L,
    original_GEB239_final_sample_included=0L,
    environmental_trait_selection_bypassed=1L,
    host_quality_or_actual_larval_use_measured=0L,
    independent_evolutionary_replication_for_host_plant=0L,
    stringsAsFactors=FALSE)
}))
if(nrow(result)!=2L)stop("Not exactly the two original independent Azores endemics")
write.csv(region,file.path(out,"Festuca_jubata_originalWCVPv13_native_introduced_regions.csv"),row.names=FALSE,na="")
write.csv(result,file.path(out,"Azores_two_endemic_host_only_WCVPv13_resource_envelope.csv"),row.names=FALSE,na="")
cat("AZORES_HIPPARCHIA_PINNED_NHM_WCVP_HOST_ONLY_RESOURCE_VALIDATED\n")
cat(sprintf("PLANT_WCVP_ACCEPTED_ID=%s\n",accepted_id))
cat(sprintf("PLANT_NATIVE_REGIONS=%s\n",paste(native,collapse=";")))
cat(sprintf("PLANT_INTRODUCED_ONLY_REGIONS=%s\n",paste(added,collapse=";")))
cat(sprintf("PLANT_NATIVE_UNITS=%d,CONTEMPORARY_UNITS=%d,ADDED_UNITS=%d\n",length(native),length(contemporary),length(added)))
cat(sprintf("SOURCE_TWO_ENDEMICS_SHARE_IDENTICAL_ONE_HOST_PLANT=%s\n",
  if(identical(sort(unique(i$plant)),"Festuca jubata"))"TRUE" else "FALSE"))
cat("IMPORTANT: hypothetical known-host potential geography only; no field larvae or causal dispersal inferred\n")
