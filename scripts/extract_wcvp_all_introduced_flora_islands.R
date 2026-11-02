args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) {
  stop("usage: extract_wcvp_all_introduced_flora_islands.R <pinned_WCVP_repo> <Q1_geographic_regions.csv> <output_dir>")
}
source_repo <- args[[1]]
region_csv <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

e <- new.env(parent=emptyenv())
load(file.path(source_repo,"data","wcvp_names.rda"),envir=e)
load(file.path(source_repo,"data","wcvp_distributions.rda"),envir=e)
stopifnot(exists("wcvp_names",envir=e),exists("wcvp_distributions",envir=e))
nms <- get("wcvp_names",envir=e)
dist <- get("wcvp_distributions",envir=e)
required_nms <- c("plant_name_id","taxon_rank","taxon_status")
required_dist <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required_nms %in% names(nms)), all(required_dist %in% names(dist)))

geo <- read.csv(region_csv, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("area_code_l3","class","native","contemporary") %in% names(geo)))
isles <- sort(as.character(geo$area_code_l3[geo$class=="island"]))
stopifnot(length(isles)==126L, "CNY" %in% isles, !("RHO" %in% isles))

species <- tolower(trimws(as.character(nms$taxon_rank)))=="species"
accepted <- tolower(trimws(as.character(nms$taxon_status)))=="accepted"
valid_ids <- unique(as.character(nms$plant_name_id[species & accepted]))
valid_ids <- valid_ids[!is.na(valid_ids) & nzchar(valid_ids)]
if(length(valid_ids)<200000L) {
  stop(sprintf("Too few accepted WCVP species: %d",length(valid_ids)))
}
d <- dist[as.character(dist$plant_name_id) %in% valid_ids &
          as.character(dist$area_code_l3) %in% isles &
          as.integer(dist$extinct)==0L &
          as.integer(dist$location_doubtful)==0L,
          required_dist,drop=FALSE]
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- as.character(d$area_code_l3)
d$introduced <- as.integer(d$introduced)
d <- unique(d)
if(nrow(d)<10000L)stop("Expected many contemporary island WCVP species×area distributions")
if(anyNA(d$introduced))stop("Unknown introduction statuses in full flora; must preserve rather than silently classify")
native <- unique(d[d$introduced==0L,c("area_code_l3","plant_name_id"),drop=FALSE])
alien <- unique(d[d$introduced==1L,c("area_code_l3","plant_name_id"),drop=FALSE])
all <- unique(d[,c("area_code_l3","plant_name_id"),drop=FALSE])
joined <- merge(native,alien,by=c("area_code_l3","plant_name_id"))
if(nrow(joined)>0L)stop(sprintf("Conflicting WCVP native and alien same taxon×region labels: %d",nrow(joined)))
count_by_area <- function(dat) table(factor(as.character(dat$area_code_l3),levels=isles))
out <- data.frame(
  area_code_l3=isles,
  native_flora_species_count=as.integer(count_by_area(native)),
  introduced_flora_species_count=as.integer(count_by_area(alien)),
  contemporary_flora_species_count=as.integer(count_by_area(all)),
  stringsAsFactors=FALSE
)
stopifnot(nrow(out)==126L)
stopifnot(all(out$native_flora_species_count+out$introduced_flora_species_count==
              out$contemporary_flora_species_count))
write.csv(out,file.path(outdir,"wcvp_all_flora_native_introduced_by_island.csv"),
          row.names=FALSE,na="")
audit <- data.frame(
  metric=c("WCVP_accepted_species_global","strict_island_regions",
           "extant_nondoubtful_accepted_species_island_records",
           "native_species_island_records","introduced_species_island_records",
           "regions_with_alien_species","regions_with_zero_alien_species",
           "conflicting_native_alien_statuses"),
  value=c(length(valid_ids),length(isles),nrow(all),nrow(native),nrow(alien),
          sum(out$introduced_flora_species_count>0L),
          sum(out$introduced_flora_species_count==0L),nrow(joined)),
  stringsAsFactors=FALSE
)
write.csv(audit,file.path(outdir,"WCVP_all_flora_island_source_audit.csv"),
          row.names=FALSE,na="")
cat(sprintf("WCVP全植物相 accepted=%d island=%d native_cell=%d alien_cell=%d\n",
            length(valid_ids),nrow(out),nrow(native),nrow(alien)))
