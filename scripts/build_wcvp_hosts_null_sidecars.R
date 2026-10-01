args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("usage: build_wcvp_hosts_null_sidecars.R <wcvp_repo> <hosts_repo> <native_output_dir> <contemporary_output_dir>")
}
wcvp_repo <- args[[1]]
hosts_repo <- args[[2]]
native_output_dir <- args[[3]]
contemporary_output_dir <- args[[4]]
dir.create(native_output_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(contemporary_output_dir, recursive=TRUE, showWarnings=FALSE)

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_names", envir=e), exists("wcvp_distributions", envir=e))
names_df <- get("wcvp_names", envir=e)
dist_df <- get("wcvp_distributions", envir=e)

hosts <- read.csv(
  file.path(hosts_repo, "resource.csv"),
  stringsAsFactors=FALSE,
  check.names=FALSE,
  na.strings=c("")
)
required_hosts <- c("Hostplant Genus","Hostplant Species")
stopifnot(all(required_hosts %in% names(hosts)))

genus <- trimws(ifelse(is.na(hosts[["Hostplant Genus"]]), "", hosts[["Hostplant Genus"]]))
species <- trimws(ifelse(is.na(hosts[["Hostplant Species"]]), "", hosts[["Hostplant Species"]]))
keep <- nzchar(genus) & nzchar(species)
host_names <- sort(unique(paste(genus[keep], species[keep])))

required_names <- c(
  "plant_name_id","taxon_rank","taxon_status","family","taxon_name",
  "accepted_plant_name_id"
)
stopifnot(all(required_names %in% names(names_df)))

rank_species <- tolower(trimws(as.character(names_df$taxon_rank))) == "species"
exact <- names_df[
  rank_species & as.character(names_df$taxon_name) %in% host_names,
  required_names,
  drop=FALSE
]
exact$plant_name_id <- as.character(exact$plant_name_id)
exact$accepted_plant_name_id <- as.character(exact$accepted_plant_name_id)

status_lower <- tolower(trimws(as.character(exact$taxon_status)))
missing_acc <- is.na(exact$accepted_plant_name_id) | !nzchar(exact$accepted_plant_name_id)
exact$accepted_plant_name_id[missing_acc & status_lower=="accepted"] <-
  exact$plant_name_id[missing_acc & status_lower=="accepted"]

id_to_name <- setNames(
  as.character(names_df$taxon_name),
  as.character(names_df$plant_name_id)
)
id_to_family <- setNames(
  as.character(names_df$family),
  as.character(names_df$plant_name_id)
)

split_ids <- split(exact$accepted_plant_name_id, as.character(exact$taxon_name))
resolved <- lapply(host_names, function(nm) {
  ids <- unique(split_ids[[nm]])
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids)==0) {
    return(data.frame(
      input_host_name=nm, match_status="unmatched",
      accepted_plant_name_id="", accepted_name="", family="",
      stringsAsFactors=FALSE
    ))
  }
  if (length(ids)>1) {
    return(data.frame(
      input_host_name=nm, match_status="ambiguous",
      accepted_plant_name_id=paste(sort(ids), collapse=";"),
      accepted_name="", family="", stringsAsFactors=FALSE
    ))
  }
  id <- ids[[1]]
  data.frame(
    input_host_name=nm, match_status="matched",
    accepted_plant_name_id=id,
    accepted_name=ifelse(is.na(id_to_name[[id]]),"",id_to_name[[id]]),
    family=ifelse(is.na(id_to_family[[id]]),"",id_to_family[[id]]),
    stringsAsFactors=FALSE
  )
})
resolved <- do.call(rbind, resolved)

insect_genus <- trimws(ifelse(is.na(hosts[["Insect Genus"]]), "", hosts[["Insect Genus"]]))
insect_species_epithet <- trimws(ifelse(is.na(hosts[["Insect Species"]]), "", hosts[["Insect Species"]]))
insect_keep <- keep & nzchar(insect_genus) & nzchar(insect_species_epithet)
interaction <- data.frame(
  insect_species=paste(insect_genus[insect_keep], insect_species_epithet[insect_keep]),
  input_host_name=paste(genus[insect_keep], species[insect_keep]),
  stringsAsFactors=FALSE
)
interaction <- unique(interaction)
interaction <- merge(
  interaction,
  resolved[resolved$match_status=="matched",
           c("input_host_name","accepted_plant_name_id","accepted_name","family"),
           drop=FALSE],
  by="input_host_name",
  all=FALSE,
  sort=FALSE
)
interaction <- unique(interaction)
interaction <- interaction[
  order(interaction$insect_species, interaction$accepted_plant_name_id),
  ,
  drop=FALSE
]
write.csv(
  interaction,
  file.path(native_output_dir, "insect_host_accepted.csv"),
  row.names=FALSE, na=""
)

matched_ids <- unique(resolved$accepted_plant_name_id[resolved$match_status=="matched"])
required_dist <- c(
  "plant_name_id","area_code_l3","introduced","extinct","location_doubtful"
)
stopifnot(all(required_dist %in% names(dist_df)))

base_keep <-
  as.character(dist_df$plant_name_id) %in% matched_ids &
  as.integer(dist_df$extinct)==0L &
  as.integer(dist_df$location_doubtful)==0L

contemporary <- dist_df[base_keep, required_dist, drop=FALSE]
contemporary$plant_name_id <- as.character(contemporary$plant_name_id)
contemporary$area_code_l3 <- as.character(contemporary$area_code_l3)
contemporary <- unique(contemporary[c("plant_name_id","area_code_l3")])
contemporary <- contemporary[nzchar(contemporary$area_code_l3),,drop=FALSE]
names(contemporary)[1] <- "accepted_plant_name_id"
contemporary$accepted_name <- unname(id_to_name[contemporary$accepted_plant_name_id])
contemporary <- contemporary[order(contemporary$accepted_plant_name_id, contemporary$area_code_l3),]
write.csv(
  contemporary,
  file.path(contemporary_output_dir, "extant_nondoubtful_wgsrpd3.csv"),
  row.names=FALSE, na=""
)

native <- dist_df[
  base_keep & as.integer(dist_df$introduced)==0L,
  required_dist,
  drop=FALSE
]
native$plant_name_id <- as.character(native$plant_name_id)
native$area_code_l3 <- as.character(native$area_code_l3)
native <- unique(native[c("plant_name_id","area_code_l3")])
native <- native[nzchar(native$area_code_l3),,drop=FALSE]
names(native)[1] <- "accepted_plant_name_id"
native$accepted_name <- unname(id_to_name[native$accepted_plant_name_id])
native <- native[order(native$accepted_plant_name_id, native$area_code_l3),]
write.csv(
  native,
  file.path(native_output_dir, "native_extant_nondoubtful_wgsrpd3.csv"),
  row.names=FALSE, na=""
)

counts <- data.frame(
  key=c(
    "resolved_insect_host_pairs",
    "native_host_x_wgsrpd3_rows",
    "contemporary_host_x_wgsrpd3_rows",
    "contemporary_minus_native_rows"
  ),
  value=c(
    nrow(interaction),
    nrow(native),
    nrow(contemporary),
    nrow(contemporary)-nrow(native)
  )
)
write.table(
  counts,
  file.path(native_output_dir, "null_sidecar_counts.tsv"),
  row.names=FALSE, quote=FALSE, sep="\t"
)
cat(sprintf(
  "pairs=%d native_rows=%d contemporary_rows=%d\n",
  nrow(interaction), nrow(native), nrow(contemporary)
))
