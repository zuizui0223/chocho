#!/usr/bin/env Rscript
# Rebuild ONLY additional BCE/Clarke accepted-plant host WCVP ranges.
# Uses frozen WCVP v13 source; no live BCE network requests and no R packages.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L)
  stop("usage: Rscript scripts/build_bce_additional_wcvp_ranges.R <WCVP_repo> <IDs_csv> <output_dir>")
root <- args[[1]]
ids_file <- args[[2]]
out <- args[[3]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
requested <- read.csv(ids_file,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot("accepted_plant_name_id" %in% names(requested))
ids <- unique(as.character(requested$accepted_plant_name_id))
ids <- ids[!is.na(ids)&nzchar(ids)]
if(length(ids)<100L) stop("Unexpectedly small BCE accepted host ID list")
env <- new.env(parent=emptyenv())
load(file.path(root,"data","wcvp_distributions.rda"),envir=env)
stopifnot(exists("wcvp_distributions",envir=env))
d <- get("wcvp_distributions",envir=env)
required <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required %in% names(d)))
keep <- as.character(d$plant_name_id) %in% ids &
  !is.na(d$extinct) & as.integer(d$extinct)==0L &
  !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L
rows <- d[keep, required, drop=FALSE]
rows$plant_name_id <- as.character(rows$plant_name_id)
rows$area_code_l3 <- as.character(rows$area_code_l3)
rows <- rows[!is.na(rows$area_code_l3)&nzchar(rows$area_code_l3),,drop=FALSE]
contemporary <- unique(rows[c("plant_name_id","area_code_l3")])
native <- unique(rows[!is.na(rows$introduced)&as.integer(rows$introduced)==0L,
  c("plant_name_id","area_code_l3"),drop=FALSE])
names(contemporary)[1] <- "accepted_plant_name_id"
names(native)[1] <- "accepted_plant_name_id"
contemporary <- contemporary[order(contemporary$accepted_plant_name_id,contemporary$area_code_l3),]
native <- native[order(native$accepted_plant_name_id,native$area_code_l3),]
write.csv(contemporary,file.path(out,"additional_contemporary_wgsrpd3.csv"),row.names=FALSE,na="")
write.csv(native,file.path(out,"additional_native_wgsrpd3.csv"),row.names=FALSE,na="")
cat(sprintf(
  "requestedIDs=%s WCVPacceptedIDsWithPresentRegion=%s nativeRows=%s presentRows=%s\n",
  length(ids),length(unique(contemporary$accepted_plant_name_id)),nrow(native),nrow(contemporary)
))
