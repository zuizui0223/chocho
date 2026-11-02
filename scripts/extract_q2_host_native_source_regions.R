args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("usage: extract_q2_host_native_source_regions.R <pinned_wcvp_repo> <q2_introduced_edges_csv> <output_dir>")
}
wcvp_repo <- args[[1]]
edges_csv <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

edges <- read.csv(edges_csv, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("plant_id", "species", "island_code") %in% names(edges)))
ids <- sort(unique(as.character(edges$plant_id)))
stopifnot(length(ids) == 560L)

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_distributions", envir=e))
d <- get("wcvp_distributions", envir=e)
required <- c("plant_name_id", "area_code_l3", "introduced", "extinct", "location_doubtful")
stopifnot(all(required %in% names(d)))

subset <- d[
  as.character(d$plant_name_id) %in% ids &
  as.integer(d$introduced)==0L &
  as.integer(d$extinct)==0L &
  as.integer(d$location_doubtful)==0L,
  c("plant_name_id","area_code_l3"), drop=FALSE
]
subset$plant_name_id <- as.character(subset$plant_name_id)
subset$area_code_l3 <- as.character(subset$area_code_l3)
subset <- unique(subset)
subset <- subset[!is.na(subset$area_code_l3) & nzchar(subset$area_code_l3),,drop=FALSE]
names(subset)[1] <- "plant_id"
subset <- subset[order(subset$plant_id,subset$area_code_l3),]
write.csv(subset,file.path(output_dir,"q2_introduced_host_native_regions.csv"),
          row.names=FALSE,na="")

counts <- data.frame(
  key=c("distinct_intro_host_ids", "native_host_region_rows", "hosts_with_native_data"),
  value=c(length(ids),nrow(subset),length(unique(subset$plant_id)))
)
write.table(counts,file.path(output_dir,"source_counts.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
cat(sprintf("requested_plants=%d native_rows=%d with_native=%d\n",
            length(ids),nrow(subset),length(unique(subset$plant_id))))
