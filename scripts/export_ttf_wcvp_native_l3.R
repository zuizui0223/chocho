#!/usr/bin/env Rscript
# Export original WCVP native WGSRPD-3 ranges + geometry, without any LGM score.
suppressPackageStartupMessages(library(sf))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: export_ttf_native_wcvp.R <rWCVPdata_checkout> <ttf_host_manifest.csv> <outdir>")
}
repo <- args[[1]]; manifest <- args[[2]]; out <- args[[3]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)

e <- new.env(parent=emptyenv())
load(file.path(repo, "data", "wcvp_distributions.rda"), envir=e)
load(file.path(repo, "data", "wgsrpd3.rda"), envir=e)
stopifnot(exists("wcvp_distributions",e), exists("wgsrpd3",e))
d <- get("wcvp_distributions",e)
g <- get("wgsrpd3",e)
stopifnot(inherits(g,"sf"), nrow(g)>=360, nrow(g)<=400)
m <- read.csv(manifest, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(m)==1037L, all(c("wcvp_id_as_frozen","accepted_host_name") %in% names(m)))
m$accepted_plant_name_id <- sprintf("%.0f",as.numeric(m$wcvp_id_as_frozen))
taxa <- unique(m[c("accepted_plant_name_id","accepted_host_name")])
stopifnot(nrow(taxa)==518L, !anyDuplicated(taxa$accepted_plant_name_id))
required <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required %in% names(d)))
lookup <- d[as.character(d$plant_name_id) %in% taxa$accepted_plant_name_id,required,drop=FALSE]
stopifnot(nrow(lookup)>0)
native <- lookup[!is.na(lookup$introduced) & !is.na(lookup$extinct) &
 !is.na(lookup$location_doubtful) &
 as.integer(lookup$introduced)==0L & as.integer(lookup$extinct)==0L &
 as.integer(lookup$location_doubtful)==0L,,drop=FALSE]
native <- unique(data.frame(
 accepted_plant_name_id=as.character(native$plant_name_id),
 area_code_l3=trimws(as.character(native$area_code_l3)),
 stringsAsFactors=FALSE))
native <- native[nzchar(native$area_code_l3),,drop=FALSE]
codes <- c("LEVEL3_COD","LEVEL3_CODE","level3_cod","level3_code","area_code_l3")
code <- codes[codes %in% names(g)][1]
if (is.na(code)) stop("no L3 code on original regional geometry")
names(g)[names(g)==code] <- "area_code_l3"
g$area_code_l3 <- trimws(as.character(g$area_code_l3))
if (anyDuplicated(g$area_code_l3))stop("duplicate WGSRPD3 geometry unit")
if (is.na(st_crs(g)))stop("unknown source CRS: no implicit coordinate assignment")
g <- st_transform(g,4326)
stopifnot(all(native$area_code_l3 %in% g$area_code_l3))
listed <- sort(unique(native$accepted_plant_name_id))
stopifnot(length(listed)==517L)
missing <- taxa[!taxa$accepted_plant_name_id %in% listed,,drop=FALSE]
stopifnot(nrow(missing)==1L,missing$accepted_host_name[[1]]=="Brachypodium pinnatum")
native <- merge(native,taxa,by="accepted_plant_name_id",all.x=TRUE,sort=FALSE)
native <- native[order(native$accepted_plant_name_id,native$area_code_l3),]
write.csv(native,file.path(out,"native_host_wgsrpd3.csv"),row.names=FALSE,na="")
write.csv(missing,file.path(out,"unresolved_native_host.csv"),row.names=FALSE,na="")
geo <- g["area_code_l3"]
geo <- geo[order(geo$area_code_l3),]
st_write(geo,file.path(out,"wgsrpd3_regions.geojson"),driver="GeoJSON",delete_dsn=TRUE,quiet=TRUE)
summary <- list(schema="chocho_wcvp_native_l3_v01",
 source="rWCVPdata pinned 65bed76bae9d644ccb6ad200c05f9f5071d89e05; WCVP snapshot version 13 2024-05-21",
 original_manifest_rows=1037L, original_hosts=518L,
 hosts_with_native_records=517L,missing_host_names=as.character(missing$accepted_host_name),
 native_L3_links=nrow(native),polygon_regions=nrow(g),
 rule="introduced==0 AND extinct==0 AND location_doubtful==0; current native L3 is NOT LGM occupied habitat",
 genetic_response_used=FALSE,paleoclimate_response_used=FALSE)
writeLines(jsonlite::toJSON(summary,pretty=TRUE,auto_unbox=TRUE),
           file.path(out,"native_source_summary.json"))
cat(jsonlite::toJSON(summary,auto_unbox=TRUE),"\n")
