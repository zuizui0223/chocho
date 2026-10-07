args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("usage: build_wcvp_introduced_plant_pool.R <wcvp_repo> <output_csv>")
}
wcvp_repo <- args[[1]]
output_csv <- args[[2]]

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
names_df <- get("wcvp_names", envir=e)
dist_df <- get("wcvp_distributions", envir=e)

required_names <- c("plant_name_id","taxon_rank","taxon_status","family","taxon_name")
required_dist <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required_names %in% names(names_df)))
stopifnot(all(required_dist %in% names(dist_df)))

accepted_species <- names_df[
  tolower(trimws(as.character(names_df$taxon_rank))) == "species" &
  tolower(trimws(as.character(names_df$taxon_status))) == "accepted",
  required_names,
  drop=FALSE
]
accepted_species$plant_name_id <- as.character(accepted_species$plant_name_id)
accepted_species <- accepted_species[!duplicated(accepted_species$plant_name_id),,drop=FALSE]
accepted_ids <- accepted_species$plant_name_id

introduced_flag <- suppressWarnings(as.integer(as.character(dist_df$introduced)))
extinct_flag <- suppressWarnings(as.integer(as.character(dist_df$extinct)))
doubt_flag <- suppressWarnings(as.integer(as.character(dist_df$location_doubtful)))

keep <-
  as.character(dist_df$plant_name_id) %in% accepted_ids &
  introduced_flag == 1L &
  extinct_flag == 0L &
  doubt_flag == 0L &
  !is.na(dist_df$area_code_l3) &
  nzchar(trimws(as.character(dist_df$area_code_l3)))

d <- dist_df[keep, c("plant_name_id","area_code_l3"), drop=FALSE]
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- trimws(as.character(d$area_code_l3))
d <- unique(d)

breadth <- aggregate(area_code_l3 ~ plant_name_id, d, function(x) length(unique(x)))
names(breadth)[2] <- "introduced_wgsrpd3_count"
meta <- accepted_species[,c("plant_name_id","taxon_name","family"),drop=FALSE]
names(meta) <- c("plant_name_id","accepted_name","family")
out <- merge(d, breadth, by="plant_name_id", all.x=TRUE, sort=FALSE)
out <- merge(out, meta, by="plant_name_id", all.x=TRUE, sort=FALSE)
out <- out[,c("plant_name_id","accepted_name","family","area_code_l3","introduced_wgsrpd3_count")]
out <- out[order(out$area_code_l3,out$accepted_name,out$plant_name_id),]

dir.create(dirname(output_csv),recursive=TRUE,showWarnings=FALSE)
write.csv(out, output_csv, row.names=FALSE, na="")
cat(sprintf(
  "introduced_rows=%d introduced_species=%d regions=%d\n",
  nrow(out), length(unique(out$plant_name_id)), length(unique(out$area_code_l3))
))
