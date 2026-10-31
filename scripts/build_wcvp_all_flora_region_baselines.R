# Rebuild the comparison background for host-resource composition.
# Use the SAME fixed WCVP v13 source as the existing GEB resource analysis.
# Output is source-derived all accepted vascular-plant species × WGSRPD3 region.
# Do not treat "contemporary" as a date-stamped historical flora.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("usage: Rscript scripts/build_wcvp_all_flora_region_baselines.R <wcvp_repo> <output_dir>")
}
wcvp_repo <- args[[1]]
out <- args[[2]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_names", envir=e),
          exists("wcvp_distributions", envir=e))
names_df <- get("wcvp_names", envir=e)
dist_df <- get("wcvp_distributions", envir=e)
req_names <- c("plant_name_id", "taxon_rank", "taxon_status")
req_dist <- c("plant_name_id", "area_code_l3", "introduced",
              "extinct", "location_doubtful")
stopifnot(all(req_names %in% names(names_df)))
stopifnot(all(req_dist %in% names(dist_df)))
# Filtering directly to accepted species prevents synonym double-counting.
valid <- tolower(trimws(as.character(names_df$taxon_rank))) == "species" &
         tolower(trimws(as.character(names_df$taxon_status))) == "accepted"
accepted <- unique(as.character(names_df$plant_name_id[!is.na(valid) & valid]))
accepted <- accepted[!is.na(accepted) & nzchar(accepted)]
keep <- !is.na(dist_df$plant_name_id) &
        as.character(dist_df$plant_name_id) %in% accepted &
        !is.na(dist_df$extinct) & as.integer(dist_df$extinct) == 0L &
        !is.na(dist_df$location_doubtful) &
        as.integer(dist_df$location_doubtful) == 0L &
        !is.na(dist_df$area_code_l3) &
        nzchar(trimws(as.character(dist_df$area_code_l3)))
d <- dist_df[keep, req_dist, drop=FALSE]
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- as.character(d$area_code_l3)
d$introduced <- as.integer(d$introduced)
d <- d[!is.na(d$introduced),,drop=FALSE]
# Native is a historical counterfactual based on present-day native flags,
# NOT reconstructed plant locations at any particular calendar year.
all <- unique(d[c("plant_name_id", "area_code_l3")])
native <- unique(d[d$introduced == 0L,c("plant_name_id","area_code_l3"),drop=FALSE])
names(all) <- c("accepted_plant_name_id", "area_code_l3")
names(native) <- c("accepted_plant_name_id", "area_code_l3")
all <- all[order(all$accepted_plant_name_id, all$area_code_l3),,drop=FALSE]
native <- native[order(native$accepted_plant_name_id, native$area_code_l3),,drop=FALSE]
stopifnot(!anyDuplicated(all), !anyDuplicated(native))
stopifnot(nrow(native) <= nrow(all), all(native$accepted_plant_name_id %in% all$accepted_plant_name_id))
write.csv(native, file.path(out,"allflora_native_accepted_wgsrpd3.csv"),
          row.names=FALSE, na="")
write.csv(all, file.path(out,"allflora_contemporary_accepted_wgsrpd3.csv"),
          row.names=FALSE, na="")
receipt <- data.frame(
  metric=c("accepted_species_global","native_plant_region_rows",
           "contemporary_plant_region_rows","native_regions","contemporary_regions",
           "native_distinct_species","contemporary_distinct_species"),
  value=c(length(accepted),nrow(native),nrow(all),
          length(unique(native$area_code_l3)),length(unique(all$area_code_l3)),
          length(unique(native$accepted_plant_name_id)),
          length(unique(all$accepted_plant_name_id))),
  stringsAsFactors=FALSE
)
write.csv(receipt, file.path(out,"allflora_source_counts.csv"),row.names=FALSE,na="")
print(receipt,row.names=FALSE)
