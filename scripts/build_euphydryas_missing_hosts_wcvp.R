#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: build_euphydryas_missing_hosts_wcvp.R <pinned_wcvp_repo> <output_dir>")
wcvp <- args[[1]]
output <- args[[2]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(file.path(wcvp, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp, "data", "wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_names", envir=e), exists("wcvp_distributions", envir=e))
n <- get("wcvp_names", envir=e)
d <- get("wcvp_distributions", envir=e)
stopifnot(all(c("taxon_name","taxon_rank","taxon_status","plant_name_id",
                "accepted_plant_name_id") %in% names(n)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct",
                "location_doubtful") %in% names(d)))
focal <- c("Castilleja hispida", "Castilleja levisecta")
taxa <- data.frame(
  host_name=focal, match_status=rep("",length(focal)),
  accepted_plant_name_id=rep("",length(focal)),
  accepted_name=rep("",length(focal)),
  exact_accepted_species_candidates=integer(length(focal)),
  stringsAsFactors=FALSE
)
for (i in seq_along(focal)) {
  hit <- n[
    as.character(n$taxon_name) == focal[[i]] &
      tolower(trimws(as.character(n$taxon_rank))) == "species" &
      tolower(trimws(as.character(n$taxon_status))) == "accepted",
    , drop=FALSE
  ]
  taxa$exact_accepted_species_candidates[[i]] <- nrow(hit)
  ids <- unique(as.character(hit$accepted_plant_name_id))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  # For accepted records with no separately populated accepted ID, use plant_name_id.
  if (!length(ids) && nrow(hit)) {
    ids <- unique(as.character(hit$plant_name_id))
    ids <- ids[!is.na(ids) & nzchar(ids)]
  }
  if (length(ids) != 1L) {
    taxa$match_status[[i]] <- if (nrow(hit)==0L) "NOT_EXACT_ACCEPTED" else "AMBIGUOUS"
    next
  }
  taxa$match_status[[i]] <- "UNIQUE_EXACT_ACCEPTED"
  taxa$accepted_plant_name_id[[i]] <- ids[[1]]
  taxa$accepted_name[[i]] <- focal[[i]]
}
if (anyDuplicated(taxa$accepted_plant_name_id[taxa$match_status=="UNIQUE_EXACT_ACCEPTED"])) {
  stop("Two distinct focal accepted taxa resolve to the same accepted ID")
}
write.csv(taxa, file.path(output, "focal_taxon_matches.csv"), row.names=FALSE, na="")
accepted <- taxa[taxa$match_status=="UNIQUE_EXACT_ACCEPTED",,drop=FALSE]
source <- d[
  as.character(d$plant_name_id) %in% accepted$accepted_plant_name_id &
    !is.na(d$extinct) & as.integer(d$extinct)==0L &
    !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L,
  c("plant_name_id","area_code_l3","introduced"),drop=FALSE
]
source$plant_name_id <- as.character(source$plant_name_id)
source$area_code_l3 <- as.character(source$area_code_l3)
source <- source[!is.na(source$area_code_l3)&nzchar(source$area_code_l3),,drop=FALSE]
contemporary <- unique(source[c("plant_name_id","area_code_l3")])
native <- unique(source[
  !is.na(source$introduced) & as.integer(source$introduced)==0L,
  c("plant_name_id","area_code_l3"),drop=FALSE
])
names(native)[[1]] <- "accepted_plant_name_id"
names(contemporary)[[1]] <- "accepted_plant_name_id"
write.csv(native, file.path(output, "castilleja_native.csv"), row.names=FALSE, na="")
write.csv(contemporary, file.path(output, "castilleja_contemporary.csv"), row.names=FALSE, na="")
cat("Exact accepted matches and distribution records\n")
print(taxa, row.names=FALSE)
for (id in accepted$accepted_plant_name_id) {
  cat(sprintf("%s native=%d contemporary=%d\n", id,
      sum(native$accepted_plant_name_id==id),
      sum(contemporary$accepted_plant_name_id==id)))
}
