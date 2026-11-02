args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("usage: source_wcvp_repo matched_host_pairs_csv output_taxonomy_csv")
source_repo <- args[[1]]
pairs <- read.csv(args[[2]],stringsAsFactors=FALSE,check.names=FALSE,
                  colClasses=c("plant_id"="character","island_code"="character"))
stopifnot(all(c("plant_id","island_code") %in% colnames(pairs)))
ids <- sort(unique(as.character(pairs$plant_id)))
stopifnot(length(ids)==44L, nrow(pairs)==100L)
e <- new.env(parent=emptyenv())
load(file.path(source_repo,"data","wcvp_names.rda"),envir=e)
stopifnot(exists("wcvp_names",envir=e))
w <- get("wcvp_names",envir=e)
req <- c("plant_name_id","taxon_name","taxon_rank","taxon_status")
stopifnot(all(req %in% names(w)))
w$plant_name_id <- as.character(w$plant_name_id)
v <- w[w$plant_name_id %in% ids,req,drop=FALSE]
if (nrow(v)!=length(ids) || anyDuplicated(v$plant_name_id)) {
  stop(sprintf("WCVP accepted ID mapping incomplete: %d rows vs %d IDs",nrow(v),length(ids)))
}
v <- v[match(ids,v$plant_name_id),,drop=FALSE]
stopifnot(identical(as.character(v$plant_name_id),ids))
eligible <- tolower(as.character(v$taxon_rank))=="species" &
            tolower(as.character(v$taxon_status))=="accepted"
out <- data.frame(
  plant_id=ids,
  wcvp_accepted_name=as.character(v$taxon_name),
  wcvp_taxon_rank=as.character(v$taxon_rank),
  wcvp_taxon_status=as.character(v$taxon_status),
  gbif_species_query_eligible=ifelse(eligible,"True","False"),
  stringsAsFactors=FALSE
)
if (any(!nzchar(out$wcvp_accepted_name)) || anyNA(out$wcvp_accepted_name))
  stop("WCVP accepted plant name missing")
dir.create(dirname(args[[3]]),showWarnings=FALSE,recursive=TRUE)
write.csv(out,args[[3]],row.names=FALSE,na="")
cat(sprintf("WCVP 44 source plants mapped; species-level accepted=%d, other=%d\n",
            sum(eligible),sum(!eligible)))
