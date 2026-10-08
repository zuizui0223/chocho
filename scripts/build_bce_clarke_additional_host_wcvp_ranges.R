#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L)
  stop("usage: build_bce_clarke_additional_host_wcvp_ranges.R <pinned_WCVP_repo> <candidate_ids_csv> <out_dir>")
root <- args[[1]]
file <- args[[2]]
out <- args[[3]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input <- read.csv(file,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot("accepted_plant_name_id" %in% names(input))
ids <- unique(as.character(input$accepted_plant_name_id))
ids <- ids[!is.na(ids)&nzchar(ids)]
if (!length(ids)) stop("No accepted plant IDs in frozen source output")
env <- new.env(parent=emptyenv())
load(file.path(root,"data","wcvp_distributions.rda"),envir=env)
stopifnot(exists("wcvp_distributions",envir=env))
d <- get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(d)))
available <- d[
  as.character(d$plant_name_id) %in% ids &
  !is.na(d$extinct)&as.integer(d$extinct)==0L &
  !is.na(d$location_doubtful)&as.integer(d$location_doubtful)==0L,
  c("plant_name_id","area_code_l3","introduced"),drop=FALSE
]
available$plant_name_id <- as.character(available$plant_name_id)
available$area_code_l3 <- as.character(available$area_code_l3)
available <- available[!is.na(available$area_code_l3)&nzchar(available$area_code_l3),,drop=FALSE]
current <- unique(available[c("plant_name_id","area_code_l3")])
native <- unique(available[
   !is.na(available$introduced)&as.integer(available$introduced)==0L,
   c("plant_name_id","area_code_l3"),drop=FALSE
])
names(current)[[1]] <- "accepted_plant_name_id"
names(native)[[1]] <- "accepted_plant_name_id"
write.csv(current,file.path(out,"additional_contemporary_wgsrpd3.csv"),row.names=FALSE,na="")
write.csv(native,file.path(out,"additional_native_wgsrpd3.csv"),row.names=FALSE,na="")
cat(sprintf("requestedAcceptedIDs=%s withCurrentRegion=%s nativePairs=%s contemporaryPairs=%s\n",
   length(ids),length(unique(current$accepted_plant_name_id)),nrow(native),nrow(current)))
