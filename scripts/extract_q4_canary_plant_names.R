args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("usage: extract_q4_canary_plant_names.R <pinned_wcvp> <CNY_three_pairs.csv> <out.csv>")
src <- args[1]; pairs <- read.csv(args[2],stringsAsFactors=FALSE,check.names=FALSE,colClasses=c(plant_id="character"))
ids <- sort(unique(as.character(pairs$plant_id)))
if (length(ids)!=3L || nrow(pairs)!=3L || any(pairs$island_code!="CNY")) stop("Not fixed 3 CNY plant×island pairs")
e <- new.env(parent=emptyenv())
load(file.path(src,"data","wcvp_names.rda"),envir=e)
if (!exists("wcvp_names",envir=e)) stop("No pinned WCVP names")
w <- get("wcvp_names",envir=e)
stopifnot(all(c("plant_name_id","taxon_name","taxon_rank","taxon_status") %in% names(w)))
w$plant_name_id <- as.character(w$plant_name_id)
x <- w[w$plant_name_id %in% ids,c("plant_name_id","taxon_name","taxon_rank","taxon_status"),drop=FALSE]
if (nrow(x)!=3L || anyDuplicated(x$plant_name_id)) stop("Missing/duplicated accepted plant IDs")
x <- x[match(ids,x$plant_name_id),,drop=FALSE]
rank <- tolower(trimws(as.character(x$taxon_rank)))
status <- tolower(trimws(as.character(x$taxon_status)))
out <- data.frame(plant_id=ids,wcvp_accepted_name=as.character(x$taxon_name),
                  wcvp_taxon_rank=rank,wcvp_taxon_status=status,
                  gbif_species_query_eligible=ifelse(rank=="species" & status=="accepted","True","False"))
dir.create(dirname(args[3]),showWarnings=FALSE,recursive=TRUE)
write.csv(out,args[3],row.names=FALSE,na="")
cat("CNY WCVP 3 PLANT NAME CANDIDATES",paste(out$wcvp_accepted_name,collapse=" | "),"\n")
