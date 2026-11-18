#!/usr/bin/env Rscript
# Check source-named larval host PLANT SPECIES against the original WCVP v13
# native/introduced species origin for the Canary TDWG3 region CNY.
#
# A regional taxon does NOT automatically imply presence on each named island.
# Source-local photo and independently verified EXOS all-islands T. majus
# presence are the only admissible local alien-host confirmations here.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript audit_focal_island_plant_WCVP_v13_20261118.R WCVP_REPO EVIDENCE_CSV OUT_DIR")
wcvp_repo <- args[[1L]]
events_csv <- args[[2L]]
output <- args[[3L]]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
env <- new.env(parent=emptyenv())
load(file.path(wcvp_repo,"data","wcvp_names.rda"),envir=env)
load(file.path(wcvp_repo,"data","wcvp_distributions.rda"),envir=env)
stopifnot(exists("wcvp_names",envir=env),exists("wcvp_distributions",envir=env))
nm <- get("wcvp_names",envir=env)
ds <- get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","taxon_name","taxon_rank","taxon_status") %in% names(nm)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(ds)))
events <- read.csv(events_csv,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c(""))
stopifnot(nrow(events)==13L,length(unique(events$event_id))==13L)
nm$plant_name_id <- as.character(nm$plant_name_id)
nm$taxon_name <- as.character(nm$taxon_name)
nm$taxon_rank <- tolower(trimws(as.character(nm$taxon_rank)))
nm$taxon_status <- tolower(trimws(as.character(nm$taxon_status)))
accepted <- nm[
  nm$taxon_rank=="species" & nm$taxon_status=="accepted" &
  !is.na(nm$taxon_name) & nzchar(nm$taxon_name),
  c("plant_name_id","taxon_name"),drop=FALSE]
d <- ds[!is.na(ds$plant_name_id)&as.character(ds$area_code_l3)=="CNY" &
   as.character(ds$extinct)=="0"&as.character(ds$location_doubtful)=="0" &
   as.character(ds$introduced)%in%c("0","1"),
   c("plant_name_id","introduced"),drop=FALSE]
d$plant_name_id <- as.character(d$plant_name_id)
d$introduced <- as.character(d$introduced)
d <- unique(d)
valid_types=c("native","introduced")
audit <- lapply(seq_len(nrow(events)),function(i){
  row <- events[i,,drop=FALSE]
  candidate <- row$reported_host_binomial
  human_origin <- row$reported_plant_origin
  stopifnot(human_origin %in% valid_types)
  trusted <- row$name_resolution_grade %in% c("high","medium")
  matches <- accepted[accepted$taxon_name==candidate,,drop=FALSE]
  accepted_id <- if(nrow(matches)==1L)matches$plant_name_id[[1]] else NA_character_
  dist <- if(!is.na(accepted_id))d[d$plant_name_id==accepted_id,,drop=FALSE]else d[FALSE,,drop=FALSE]
  states <- sort(unique(dist$introduced))
  source_status <- if(nrow(matches)==0L)"no_exact_accepted_species" else
     if(nrow(matches)>1L)"multiple_accepted_species" else
     if(length(states)==0L)"accepted_species_no_CNY_distribution" else
     if(identical(states,"0"))"native_CNY" else
     if(identical(states,"1"))"introduced_CNY" else "conflicting_native_and_introduced_CNY"
  conclusion <- if(!trusted)"quarantined_plant_taxonomic_concept" else
     if(source_status=="native_CNY")"native_CNY_exact_species" else
     if(source_status=="introduced_CNY")"introduced_CNY_exact_species" else
     "species_identity_or_origin_not_confirmed"
  data.frame(
     event_id=row$event_id,
     butterfly_name=row$butterfly_name,
     reported_host_binomial=candidate,
     named_island=if(is.na(row$named_island))"" else row$named_island,
     direct_larva_use=as.integer(row$strict_direct_named_island_larva_use),
     historic_extinct=as.integer(row$historically_extinct_butterfly),
     published_origin=human_origin,
     source_taxonomic_grade=row$name_resolution_grade,
     original_WCVP_accepted_plant_name_id=if(is.na(accepted_id))"" else accepted_id,
     accepted_matching_n=nrow(matches),
     wcvp_CNY_status=source_status,
     admissible_species_origin=conclusion,
     trusted_island_specific_direct_introduced_use=as.integer(
       trusted && source_status=="introduced_CNY" &&
       human_origin=="introduced" &&
       row$strict_direct_named_island_larva_use==1L &&
       row$historically_extinct_butterfly==0L &&
       !is.na(row$named_island) && nzchar(row$named_island)
     ),
     stringsAsFactors=FALSE)
})
audit <- do.call(rbind,audit)
audit <- audit[order(audit$event_id),,drop=FALSE]
write.csv(audit,file.path(output,"focal_island_named_host_plant_species_WCVP_CNY_v13.csv"),
          row.names=FALSE,na="")
positive <- audit[audit$trusted_island_specific_direct_introduced_use==1L,,drop=FALSE]
cat("ORIGINAL_WCVP_NAMED_PLANT_SPECIES_CNY_CROSSWALK_COMPLETED\n")
cat(sprintf("EVENTS=%d\n",nrow(audit)))
cat(sprintf("DIRECT_ISLAND_ALIEN_LARVAE_WCVP_SPECIES_CONFIRMED=%d\n",nrow(positive)))
for(i in seq_len(nrow(audit))){
 a <- audit[i,,drop=FALSE]
 cat(sprintf("PLANT_EVENT=%s,butterfly=%s,plant=%s,island=%s,wcvp=%s,admissible=%s\n",
   a$event_id,a$butterfly_name,a$reported_host_binomial,a$named_island,
   a$wcvp_CNY_status,a$admissible_species_origin))
}
# A species-level proof of any recorded exotic caterpillar host on the target
# islands MUST originate from exact plant taxonomy + regional origin and a
# source field larva. Do not count genus-level congeners as host species.
if(nrow(positive)>0L){
  stopifnot(all(positive$reported_host_binomial=="Tropaeolum majus"))
  stopifnot(all(positive$butterfly_name=="Pieris cheiranthi"))
  stopifnot(identical(sort(unique(positive$named_island)),c("La Palma","Tenerife")))
}
