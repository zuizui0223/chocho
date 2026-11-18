#!/usr/bin/env Rscript
# Extract accepted, extant, nondoubtful WCVP vascular flora in the same
# location-coded WGSRPD3 regions used for the NHM HOSTS reports.
# Family-specific denominators allow taxon-composition-matched counterfactuals.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=6L)stop("Usage: Rscript extract_WCVP_hostreport_family_flora_20261116.R <WCVPdir> <strict_localities.csv> <FAO_exact.txt> <FAO_wildcard_genera.txt> <FAO_all_genera.txt> <outdir>")
wp <- args[[1]]; geopath <- args[[2]]; outdir <- args[[6]]
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
loc <- read.csv(geopath, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("HOSTS_location","WGSRPD3_codes") %in% names(loc)),nrow(loc)>=90L)
regions <- sort(unique(unlist(strsplit(loc$WGSRPD3_codes,"|",fixed=TRUE))))
regions <- regions[nzchar(regions)]
exact <- unique(trimws(readLines(args[[3]],warn=FALSE,encoding="UTF-8")))
wild <- unique(trimws(readLines(args[[4]],warn=FALSE,encoding="UTF-8")))
gen_all <- unique(trimws(readLines(args[[5]],warn=FALSE,encoding="UTF-8")))

e <- new.env(parent=emptyenv())
load(file.path(wp,"data","wcvp_names.rda"),envir=e)
load(file.path(wp,"data","wcvp_distributions.rda"),envir=e)
nms <- get("wcvp_names",envir=e)
d <- get("wcvp_distributions",envir=e)
needn<-c("plant_name_id","taxon_rank","taxon_status","taxon_name","family")
needd<-c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(needn %in% names(nms)),all(needd %in% names(d)))
valid <- !is.na(nms$taxon_rank) & !is.na(nms$taxon_status) &
  tolower(trimws(as.character(nms$taxon_rank)))=="species" &
  tolower(trimws(as.character(nms$taxon_status)))=="accepted"
accepted <- nms[valid,needn,drop=FALSE]
accepted$plant_name_id<-as.character(accepted$plant_name_id)
accepted<-accepted[!duplicated(accepted$plant_name_id),]
if(nrow(accepted)<200000L)stop("Far too few accepted WCVP species")
keep <- as.character(d$area_code_l3) %in% regions &
  as.character(d$plant_name_id) %in% accepted$plant_name_id &
  !is.na(d$extinct) & as.integer(d$extinct)==0L &
  !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L &
  !is.na(d$introduced) & as.integer(d$introduced) %in% c(0L,1L)
q <- unique(d[keep,c("plant_name_id","area_code_l3","introduced"),drop=FALSE])
q$plant_name_id<-as.character(q$plant_name_id)
q$area_code_l3<-as.character(q$area_code_l3)
q$introduced<-as.integer(q$introduced)
if(nrow(q)<30000L)stop("Too few accepted WCVP species-by-region records")
# Retain all genuine accepted flora for whole-flora denominators, but
# exclude ambiguous native+introduced labels for the same plant/region.
key<-paste(q$plant_name_id,q$area_code_l3,sep="|")
mixed <- unique(key[duplicated(key) | duplicated(key,fromLast=TRUE)])
invalid_n <- length(mixed)
if(invalid_n>0L)q<-q[!(key %in% mixed),,drop=FALSE]
index<-match(q$plant_name_id,accepted$plant_name_id)
if(anyNA(index))stop("Unmatched accepted WCVP taxon after filtering")
q$taxon_name<-as.character(accepted$taxon_name[index])
q$family<-trimws(as.character(accepted$family[index]))
q$family[is.na(q$family)|!nzchar(q$family)]<-"__UNKNOWN_FAMILY__"
# Only the first two words, as done in the existing chocho FAO parser.
q$binomial<-sub("^([A-Z][A-Za-z-]+ [a-z][a-z-]+).*","\\1",trimws(q$taxon_name))
q$genus<-sub(" .*","",trimws(q$taxon_name))
q$crop_exact<-q$binomial %in% exact
q$crop_wild<-q$genus %in% wild
q$crop_genus<-q$genus %in% gen_all
if(sum(q$crop_exact)==0L)stop("FAO exact crop species did not match WCVP")
modes<-c("all_hosts","exclude_FAO_exact","exclude_FAO_exact_plus_genus_wildcards","exclude_FAO_all_crop_genera")
out<-list()
for(m in modes){
 keep<-switch(m,
  all_hosts=rep(TRUE,nrow(q)),
  exclude_FAO_exact= !q$crop_exact,
  exclude_FAO_exact_plus_genus_wildcards=!(q$crop_exact|q$crop_wild),
  exclude_FAO_all_crop_genera= !q$crop_genus)
 v<-q[keep,c("area_code_l3","family","introduced"),drop=FALSE]
 overall<-v
 overall$family<-"__ALL_FAMILIES__"
 v<-rbind(v,overall)
 v$cropmode<-m
 agg<-aggregate(
   rep.int(1L,nrow(v)),
   by=list(area_code_l3=v$area_code_l3,family=v$family,cropmode=v$cropmode,introduced=v$introduced),
   FUN=sum)
 names(agg)[5]<-"n"
 out[[m]]<-agg
}
all<-do.call(rbind,out)
all<-all[order(all$area_code_l3,all$cropmode,all$family,all$introduced),]
write.csv(all,file.path(outdir,"wcvp_flora_region_family_origin_FAO.csv"),row.names=FALSE,na="")
audit<-list(
  WCVP_original_accepted_species=nrow(accepted),
  source_region_count=length(regions),
  filtered_source_species_region_origin_rows=nrow(q),
  contradictory_native_introduced_species_region_keys_excluded=invalid_n,
  missing_family_rows=sum(q$family=="__UNKNOWN_FAMILY__"),
  FAO_taxa_source=list(exact_binomials=length(exact),explicit_wildcard_genera=length(wild),all_fao_genera=length(gen_all)),
  family_rows=nrow(all),
  status="CONTEMPORARY_WCVP_FLORA_SOURCE_NOT_HOST_AVAILABILITY",
  caveat="WCVP vascular flora is not a census of plants physically accessible to butterfly larvae; family-matched counterfactuals only adjust taxonomic opportunity and contemporary source geography"
)
writeLines(jsonlite_not_needed<-paste(capture.output(dput(audit)),collapse="\n"),
 file.path(outdir,"floral_source_dput.txt"))
cat("WCVP_FAMILY_FLORA_SUMMARY", "accepted",nrow(accepted),"regions",length(regions),
 "accepted_species_regions",nrow(q),"contradictory",invalid_n,"count_rows",nrow(all),
 "missing_family_rows",audit$missing_family_rows,"\n")
