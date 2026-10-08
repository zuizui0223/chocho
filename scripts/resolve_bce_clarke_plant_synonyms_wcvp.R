#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L)
  stop("usage: resolve_bce_clarke_plant_synonyms_wcvp.R <pinned_WCVP_repo> <candidate_names_csv> <output_csv>")
repo <- args[[1]]
candidates <- read.csv(args[[2]], stringsAsFactors=FALSE, check.names=FALSE)
stopifnot("plant_binomial" %in% names(candidates))
queries <- sort(unique(trimws(as.character(candidates$plant_binomial))))
if (length(queries)<5L || any(!grepl("^[A-Z][A-Za-z-]+ [a-z][a-z-]+$", queries)))
  stop("Non-species-name inputs encountered")
env <- new.env(parent=emptyenv())
load(file.path(repo,"data","wcvp_names.rda"),envir=env)
if (!exists("wcvp_names",envir=env))stop("Missing pinned WCVP names")
n <- get("wcvp_names",envir=env)
required <- c("taxon_name","taxon_rank","taxon_status","plant_name_id","accepted_plant_name_id")
stopifnot(all(required %in% names(n)))
eligible <- n[
  as.character(n$taxon_name) %in% queries &
  tolower(trimws(as.character(n$taxon_rank)))=="species",
  ,drop=FALSE
]
eligible$accepted_plant_name_id <- as.character(eligible$accepted_plant_name_id)
eligible$plant_name_id <- as.character(eligible$plant_name_id)
accepted <- tolower(trimws(as.character(eligible$taxon_status)))=="accepted"
blank <- is.na(eligible$accepted_plant_name_id) |
  !nzchar(eligible$accepted_plant_name_id)
eligible$accepted_plant_name_id[accepted & blank] <-
  eligible$plant_name_id[accepted & blank]
byname <- split(seq_len(nrow(eligible)),as.character(eligible$taxon_name))
name_to_acc <- setNames(as.character(n$taxon_name),as.character(n$plant_name_id))
out <- lapply(queries,function(q) {
  idx <- byname[[q]]
  if (is.null(idx)) idx <- integer(0)
  source <- eligible[idx,,drop=FALSE]
  valid <- !is.na(source$accepted_plant_name_id) &
    nzchar(source$accepted_plant_name_id)
  chosen <- source[valid,,drop=FALSE]
  primary <- chosen[
    tolower(trimws(as.character(chosen$taxon_status)))=="accepted",,
    drop=FALSE]
  ids <- unique(primary$accepted_plant_name_id)
  priority <- "UNIQUE_EXACT_ACCEPTED"
  if (length(ids)==0) {
    ids <- unique(chosen$accepted_plant_name_id)
    priority <- "UNIQUE_SYNONYM_TO_ACCEPTED"
  }
  status <- if (length(ids)==0) "UNMATCHED_EXACT_SPECIES_NAME" else
    if (length(ids)>1) "AMBIGUOUS_ACCEPTED_ID" else priority
  id <- if (length(ids)==1) ids[[1]] else ""
  nm <- if (nzchar(id) && id %in% names(name_to_acc))
    name_to_acc[[id]] else ""
  data.frame(plant_binomial=q,match_status=status,
     accepted_plant_name_id=id,accepted_name=nm,
     rows_in_wcvp=nrow(source),competing_accepted_ids=length(ids),
     stringsAsFactors=FALSE)
})
result <- do.call(rbind,out)
write.csv(result,args[[3]],row.names=FALSE,na="")
cat("WCVP matching", nrow(result),"\n")
print(table(result$match_status))
