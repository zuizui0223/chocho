#!/usr/bin/env Rscript
# Extract source-locked, strictly matching accepted plant x WGSRPD3 native /
# introduced flags for the original HOSTS literature row analysis.
# These WCVP region labels are not historical local host introduction dates.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("usage: Rscript extract_HOSTS_location_origin_20261116.R <wcvp_repo> <host_name_resolution.csv> <strict_locality_crosswalk.csv> <out.csv>")
}
wcvp <- args[[1]]
resolution <- read.csv(args[[2]], check.names=FALSE, stringsAsFactors=FALSE)
localities <- read.csv(args[[3]], check.names=FALSE, stringsAsFactors=FALSE)
need <- c("match_status","accepted_plant_name_id")
if (!all(need %in% names(resolution))) stop("Missing host-name-resolution fields")
if (!all(c("WGSRPD3_codes","comparison_group","island_origin_class") %in% names(localities))) {
  stop("Missing curated geography crosswalk fields")
}

wanted_ids <- sort(unique(as.character(resolution$accepted_plant_name_id[resolution$match_status=="matched"])))
wanted_ids <- wanted_ids[!is.na(wanted_ids) & nzchar(wanted_ids)]
wanted_codes <- sort(unique(unlist(strsplit(as.character(localities$WGSRPD3_codes), "|", fixed=TRUE))))
wanted_codes <- wanted_codes[!is.na(wanted_codes) & nzchar(wanted_codes)]

e <- new.env(parent=emptyenv())
load(file.path(wcvp,"data","wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_distributions",envir=e))
d <- get("wcvp_distributions",envir=e)
required <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
if (!all(required %in% names(d))) stop("Unexpected WCVP v13 distribution columns")
keep <- !is.na(d$plant_name_id) & as.character(d$plant_name_id) %in% wanted_ids &
  !is.na(d$area_code_l3) & as.character(d$area_code_l3) %in% wanted_codes &
  !is.na(d$extinct) & as.integer(d$extinct)==0L &
  !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L &
  !is.na(d$introduced) & as.integer(d$introduced) %in% c(0L,1L)
o <- unique(d[keep,c("plant_name_id","area_code_l3","introduced")])
colnames(o) <- c("accepted_plant_name_id","area_code_l3","introduced")
o$accepted_plant_name_id <- as.character(o$accepted_plant_name_id)
o$area_code_l3 <- as.character(o$area_code_l3)
o$introduced <- as.integer(o$introduced)
o <- o[order(o$accepted_plant_name_id,o$area_code_l3,o$introduced),]
dir.create(dirname(args[[4]]),showWarnings=FALSE,recursive=TRUE)
write.csv(o,args[[4]],row.names=FALSE,na="")
cat(sprintf("WCVP_HOST_REGION_ORIGIN ids=%d regions=%d rows=%d native=%d introduced=%d\n",
 length(unique(o$accepted_plant_name_id)), length(unique(o$area_code_l3)),nrow(o),
 sum(o$introduced==0),sum(o$introduced==1)))
