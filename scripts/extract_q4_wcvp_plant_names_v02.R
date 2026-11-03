args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop("usage: extract_q4_wcvp_plant_names.R <pinned_WCVP_repo> <pairs_csv> <out_csv>")
source_repo <- args[[1]]
pairs_path <- args[[2]]
out <- args[[3]]
p <- read.csv(pairs_path, stringsAsFactors=FALSE, check.names=FALSE,
              colClasses=c("plant_id"="character"))
stopifnot(all(c("plant_id","island_code") %in% names(p)))
ids <- sort(unique(as.character(p$plant_id)))
stopifnot(length(ids)==564L)
env <- new.env(parent=emptyenv())
load(file.path(source_repo,"data","wcvp_names.rda"),envir=env)
stopifnot(exists("wcvp_names",envir=env))
w <- get("wcvp_names",envir=env)
req <- c("plant_name_id","taxon_rank","taxon_status","taxon_name","accepted_plant_name_id")
stopifnot(all(req %in% names(w)))
w$plant_name_id <- as.character(w$plant_name_id)
take <- w[w$plant_name_id %in% ids, req, drop=FALSE]
if (nrow(take)!=length(ids) || anyDuplicated(take$plant_name_id)) {
 stop(sprintf("WCVP accepted plant ID mapping not unique/complete: got %d for %d",nrow(take),length(ids)))
}
take <- take[match(ids,take$plant_name_id),,drop=FALSE]
if (!identical(as.character(take$plant_name_id),ids)) stop("WCVP matched ID order drift")
if (any(is.na(take$taxon_name) | !nzchar(trimws(as.character(take$taxon_name))))) stop("Missing plant names")
rank <- tolower(trimws(as.character(take$taxon_rank)))
status <- tolower(trimws(as.character(take$taxon_status)))
eligible <- rank=="species" & status=="accepted"
outdf <- data.frame(
  plant_id=ids, wcvp_accepted_name=as.character(take$taxon_name),
  wcvp_taxon_rank=rank, wcvp_taxon_status=status,
  gbif_species_query_eligible=ifelse(eligible,"True","False"),
  stringsAsFactors=FALSE
)
dir.create(dirname(out),showWarnings=FALSE,recursive=TRUE)
write.csv(outdf,out,row.names=FALSE,na="")
cat(sprintf("Q4 v02 WCVP taxonomy crosswalk PASS: total=%d species-eligible=%d infraspecific_or_ambiguous=%d\n",nrow(outdf),sum(eligible),sum(!eligible)))
print(table(rank,status,useNA="ifany"))
