#!/usr/bin/env Rscript
# WCVP v13 source-frozen Canary whole-flora host-genus plausibility audit.
# A genus having introduced taxa on CNY DOES NOT identify any actual host
# species used by a butterfly. This is a taxonomic-grain non-identification test.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: audit_canary_host_genus_WCVP_v13_20261118.R <wcvp-repo> <uniform-host-csv> <output-dir>")
wcvp_repo <- args[[1L]]
butterfly_file <- args[[2L]]
output <- args[[3L]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo,"data","wcvp_names.rda"),envir=e)
load(file.path(wcvp_repo,"data","wcvp_distributions.rda"),envir=e)
stopifnot(exists("wcvp_names",envir=e),exists("wcvp_distributions",envir=e))
names_df <- get("wcvp_names",envir=e)
dist_df <- get("wcvp_distributions",envir=e)
stopifnot(all(c("plant_name_id","taxon_rank","taxon_status","taxon_name") %in% names(names_df)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(dist_df)))
butterflies <- read.csv(butterfly_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings="")
stopifnot(nrow(butterflies)==33L)
host_by_butterfly <- lapply(seq_len(nrow(butterflies)),function(i){
    names <- butterflies$larval_host_genera_uniform[i]
    if (is.na(names) || !nzchar(names)) return(NULL)
    taxa <- trimws(strsplit(names,";",fixed=TRUE)[[1L]])
    stopifnot(all(nzchar(taxa)),length(unique(taxa))==length(taxa))
    data.frame(
      butterfly=rep(butterflies$independent_checklist_species[i],length(taxa)),
      butterfly_origin=rep(butterflies$canary_origin[i],length(taxa)),
      genus=taxa,stringsAsFactors=FALSE)
})
host_by_butterfly <- do.call(rbind,host_by_butterfly)
host_by_butterfly <- unique(host_by_butterfly)
stopifnot(nrow(host_by_butterfly)==86L)
genus_set <- sort(unique(host_by_butterfly$genus))
stopifnot(length(genus_set)==52L)
names_df$plant_name_id <- as.character(names_df$plant_name_id)
names_df$taxon_name <- as.character(names_df$taxon_name)
names_df$taxon_rank <- tolower(as.character(names_df$taxon_rank))
names_df$taxon_status <- tolower(as.character(names_df$taxon_status))
# Include only accepted SPECIES; do not inherit a genus-wide status from a
# synonym, accepted genus-level taxa, infraspecies, doubtful/extinct occurrences.
accepted <- names_df[
  names_df$taxon_rank=="species" &
  names_df$taxon_status=="accepted" &
  !is.na(names_df$plant_name_id) &
  !is.na(names_df$taxon_name),
  c("plant_name_id","taxon_name"),drop=FALSE]
accepted$genus <- sub("[[:space:]].*$","",trimws(accepted$taxon_name))
accepted <- accepted[accepted$genus %in% genus_set,,drop=FALSE]
stopifnot(nrow(accepted)>0L)
region <- dist_df[
    as.character(dist_df$area_code_l3)=="CNY" &
    as.character(dist_df$extinct)=="0" &
    as.character(dist_df$location_doubtful)=="0" &
    as.character(dist_df$introduced) %in% c("0","1"),
    c("plant_name_id","introduced"),drop=FALSE]
region$plant_name_id <- as.character(region$plant_name_id)
region$introduced <- as.character(region$introduced)
region <- unique(region)
m <- merge(accepted, region, by="plant_name_id")
m <- unique(m[,c("plant_name_id","taxon_name","genus","introduced"),drop=FALSE])
stopifnot(nrow(m)>50L)
# If the same accepted species is assigned both native and introduced flags
# on the same CNY level-3 area, classify its status as conflict, not evidence.
species_status <- aggregate(
   introduced~plant_name_id+taxon_name+genus,
   data=m,FUN=function(v) paste(sort(unique(v)),collapse=";"))
species_status$status <- ifelse(species_status$introduced=="0","native",
  ifelse(species_status$introduced=="1","introduced","conflicting"))
out_genus <- lapply(genus_set,function(genus){
  ss <- species_status[species_status$genus==genus,,drop=FALSE]
  native <- sum(ss$status=="native")
  introduced <- sum(ss$status=="introduced")
  conflict <- sum(ss$status=="conflicting")
  data.frame(
    host_genus=genus,
    accepted_wcvp_native_species_on_CNY=native,
    accepted_wcvp_introduced_species_on_CNY=introduced,
    accepted_wcvp_conflicting_origin_species_on_CNY=conflict,
    accepted_wcvp_total_species_on_CNY=nrow(ss),
    region_status=if(native>0L && introduced>0L) "both_native_and_introduced" else
                  if(native>0L) "native_only" else
                  if(introduced>0L) "introduced_only" else
                  if(conflict>0L) "conflicted_only" else "no_species_recorded",
    stringsAsFactors=FALSE)
})
out_genus <- do.call(rbind,out_genus)
stopifnot(nrow(out_genus)==52L, length(unique(out_genus$host_genus))==52L)
out_pairs <- merge(host_by_butterfly,out_genus,by.x="genus",by.y="host_genus",
                   all.x=TRUE,sort=FALSE)
stopifnot(nrow(out_pairs)==86L, !any(is.na(out_pairs$region_status)))
out_pairs <- out_pairs[order(out_pairs$butterfly,out_pairs$genus),,drop=FALSE]
summary_rows <- lapply(split(out_pairs,out_pairs$butterfly),function(p){
  data.frame(
    butterfly=p$butterfly[1L],
    butterfly_origin=p$butterfly_origin[1L],
    known_larval_host_genera=nrow(p),
    host_genera_with_any_native_CNY_species=sum(p$accepted_wcvp_native_species_on_CNY>0L),
    host_genera_with_any_introduced_CNY_species=sum(p$accepted_wcvp_introduced_species_on_CNY>0L),
    host_genera_mixed_native_introduced_CNY=sum(p$region_status=="both_native_and_introduced"),
    host_genera_with_no_extant_nondoubtful_accepted_CNY_species=sum(p$region_status=="no_species_recorded"),
    stringsAsFactors=FALSE)
})
out_butterfly <- do.call(rbind,summary_rows)
stopifnot(nrow(out_butterfly)==24L)
write.csv(out_genus,file.path(output,"canary_CNY_WCVPv13_host_genus_native_introduced_inventory.csv"),row.names=FALSE,na="")
write.csv(out_pairs,file.path(output,"canary_CNY_WCVPv13_butterfly_genus_context_pairs.csv"),row.names=FALSE,na="")
write.csv(out_butterfly,file.path(output,"canary_CNY_WCVPv13_butterfly_genus_context_summary.csv"),row.names=FALSE,na="")
origin_profiles <- lapply(split(out_butterfly,out_butterfly$butterfly_origin),function(sub){
  data.frame(
    butterfly_origin=sub$butterfly_origin[[1L]],
    larval_genus_documented_butterflies=nrow(sub),
    has_any_CNY_accepted_native_congener=sum(sub$host_genera_with_any_native_CNY_species>0L),
    has_any_CNY_accepted_introduced_congener=sum(sub$host_genera_with_any_introduced_CNY_species>0L),
    has_any_CNY_mixed_native_introduced_host_genus=sum(sub$host_genera_mixed_native_introduced_CNY>0L),
    all_larval_host_genera_no_CNY_species=sum(
      sub$host_genera_with_no_extant_nondoubtful_accepted_CNY_species==sub$known_larval_host_genera
    ),
    stringsAsFactors=FALSE)
})
origin_profiles <- do.call(rbind,origin_profiles)
write.csv(origin_profiles,file.path(output,"canary_CNY_WCVPv13_origin_profile_nonbiological.csv"),row.names=FALSE,na="")
for (i in seq_len(nrow(origin_profiles))) {
  x<-origin_profiles[i,,drop=FALSE]
  cat(sprintf("WCVP_ORIGIN_PROFILE=%s,documented_butterflies=%d,native_genus_congener=%d,introduced_genus_congener=%d,mixed_genus=%d\\n",
  x$butterfly_origin, x$larval_genus_documented_butterflies, x$has_any_CNY_accepted_native_congener,
  x$has_any_CNY_accepted_introduced_congener,x$has_any_CNY_mixed_native_introduced_host_genus))
}
cat("PINNED_WCVP_CANARY_GENUS_CONTEXT_PASSED\n")
cat(sprintf("ORIGINAL_CNY_WCVP_ACCEPTED_SPECIES=%d\n",nrow(species_status)))
cat(sprintf("HOST_GENERA_WITH_CNY_NATIVE=%d\n",sum(out_genus$accepted_wcvp_native_species_on_CNY>0L)))
cat(sprintf("HOST_GENERA_WITH_CNY_INTRODUCED=%d\n",sum(out_genus$accepted_wcvp_introduced_species_on_CNY>0L)))
cat(sprintf("HOST_GENERA_WITH_BOTH_CNY_ORIGINS=%d\n",sum(out_genus$region_status=="both_native_and_introduced")))
cat(sprintf("HOST_GENERA_WITHOUT_ANY_CNY_ACCEPTED_SPECIES=%d\n",sum(out_genus$region_status=="no_species_recorded")))
cat(sprintf("BUTTERFLIES_WITH_ANY_HOST_GENUS_CONTAINING_CNY_INTRODUCED_SPECIES=%d\n",sum(out_butterfly$host_genera_with_any_introduced_CNY_species>0L)))
cat("INTERPRETATION: genus co-listing is NOT evidence that any introduced congener is a real larval host\n")
