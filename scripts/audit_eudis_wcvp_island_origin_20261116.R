#!/usr/bin/env Rscript
# SOURCE COVERAGE ONLY: crosswalk EuDiS v5 morphological thalassochory
# with accepted WCVP v13 species and region-specific native/introduced flora.
# This does NOT measure seed movement by ocean currents or ship traffic.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=5L) stop("Usage: Rscript <script> WCVPdir EuDiS.csv Q1_regions.csv original_WCVP98.csv output_dir")
dir.create(a[[5]], recursive=TRUE, showWarnings=FALSE)

normalize_binomial <- function(s) {
  s <- trimws(as.character(s))
  m <- regexec("^([A-Z][A-Za-z-]+)[[:space:]]+([a-z][a-z-]+)(?:[[:space:]].*)?$",s,perl=TRUE)
  x <- regmatches(s,m)
  vapply(x,function(z) if(length(z)>2L) tolower(paste(z[2],z[3])) else NA_character_,character(1))
}
eu <- read.csv(a[[2]], stringsAsFactors=FALSE,check.names=FALSE,fileEncoding="UTF-8-BOM")
names(eu) <- trimws(names(eu))
required <- c("Species_POWO","Thalassochorous","Further_clarification_needed")
if(!all(required %in% names(eu))) stop("Required EuDiS published fields absent")
if(nrow(eu)<9000L || nrow(eu)>15000L) stop("Unexpected EuDiS source size")
eu$key <- normalize_binomial(eu$Species_POWO)
eu$trait <- trimws(as.character(eu$Thalassochorous))
eu$uncertain_tha <- grepl("THA",as.character(eu$Further_clarification_needed),fixed=TRUE)
eu$trait[!eu$trait %in% c("0","1") | eu$uncertain_tha] <- NA_character_
eu <- eu[!is.na(eu$key),c("key","trait"),drop=FALSE]
trait_groups <- split(eu$trait,eu$key)
trait_lookup <- vapply(trait_groups,function(z) {
  u <- unique(z)
  if(length(u)==1L && !is.na(u[[1]])) as.integer(u[[1]]) else NA_integer_
},integer(1))
if(length(trait_lookup)<7000L || sum(trait_lookup==1L,na.rm=TRUE)<100L)
  stop("EuDiS trait identity or positive-thalassochory coverage unexpectedly low")

geo <- read.csv(a[[3]],stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c("area_code_l3","level1","class") %in% names(geo))) stop("Wrong Q1 geography source")
geo <- geo[!duplicated(geo$area_code_l3),]
if(sum(geo$class=="island")!=126L) stop("Original 126-island geography changed")
regions <- sort(as.character(geo$area_code_l3))
envir <- new.env(parent=emptyenv())
load(file.path(a[[1]],"data","wcvp_names.rda"),envir=envir)
load(file.path(a[[1]],"data","wcvp_distributions.rda"),envir=envir)
nm <- get("wcvp_names",envir=envir)
ds <- get("wcvp_distributions",envir=envir)
names_need <- c("plant_name_id","taxon_name","taxon_rank","taxon_status")
dist_need <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
if(!all(names_need %in% names(nm)) || !all(dist_need %in% names(ds)))
  stop("Pinned WCVP structure has changed")
ok <- !is.na(nm$taxon_rank) & !is.na(nm$taxon_status) &
  tolower(trimws(as.character(nm$taxon_rank)))=="species" &
  tolower(trimws(as.character(nm$taxon_status)))=="accepted"
names_ok <- nm[ok,c("plant_name_id","taxon_name"),drop=FALSE]
names_ok$plant_name_id <- as.character(names_ok$plant_name_id)
names_ok <- names_ok[!duplicated(names_ok$plant_name_id),]
if(nrow(names_ok)<200000L) stop("WCVP accepted-species population too small")
names_ok$bin <- normalize_binomial(names_ok$taxon_name)
name_freq <- table(names_ok$bin, useNA="no")
ambiguous_bins <- names(name_freq)[name_freq>1L]

# Count every accepted contemporary species by origin (not only EuDiS species)
valid <- as.character(ds$plant_name_id) %in% names_ok$plant_name_id &
  as.character(ds$area_code_l3) %in% regions &
  !is.na(ds$extinct) & as.integer(ds$extinct)==0L &
  !is.na(ds$location_doubtful) & as.integer(ds$location_doubtful)==0L &
  !is.na(ds$introduced) & as.integer(ds$introduced) %in% 0:1
q <- unique(ds[valid,c("plant_name_id","area_code_l3","introduced"),drop=FALSE])
q$plant_name_id <- as.character(q$plant_name_id)
q$area_code_l3 <- as.character(q$area_code_l3)
q$introduced <- as.integer(q$introduced)
if(nrow(q)<30000L) stop("Contemporary WCVP species-by-region source unexpectedly small")
idarea <- paste(q$plant_name_id,q$area_code_l3,sep="|")
if(anyDuplicated(idarea)) stop("Conflicting WCVP native/introduced status for one species-region")
idindex <- match(q$plant_name_id,names_ok$plant_name_id)
if(anyNA(idindex)) stop("WCVP accepted ID missing after inclusion filter")
q$bin <- names_ok$bin[idindex]
q$thalassochorous <- unname(trait_lookup[match(q$bin,names(trait_lookup))])
q$thalassochorous[q$bin %in% ambiguous_bins] <- NA_integer_

# Original port-flora WCVP counts are an independent fail-closed crosscheck.
def_table <- function(d) table(factor(d$area_code_l3,levels=regions),
                               factor(d$introduced,levels=0:1))
all_t <- def_table(q)
matched <- q[!is.na(q$thalassochorous),,drop=FALSE]
match_t <- def_table(matched)
positive_t <- def_table(matched[matched$thalassochorous==1L,,drop=FALSE])
port_orig <- read.csv(a[[4]],stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c("code","native","introduced") %in% names(port_orig))) stop("WCVP port flora control missing")
shared <- intersect(as.character(port_orig$code),regions)
if(length(shared)<70L) stop("Too few original island-port WCVP control regions")
port_index <- match(shared,as.character(port_orig$code))
if(any(all_t[shared,"0"]!=port_orig$native[port_index]) ||
   any(all_t[shared,"1"]!=port_orig$introduced[port_index]))
  stop("WCVP original island-port native/introduced totals do not reproduce")

rows <- do.call(rbind,lapply(0:1,function(k) {
  data.frame(area_code_l3=regions,introduced=k,
    WCVP_accepted_species=as.integer(all_t[,as.character(k)]),
    EuDiS_matched_species=as.integer(match_t[,as.character(k)]),
    thalasso_trait_1=as.integer(positive_t[,as.character(k)]),
    stringsAsFactors=FALSE)
}))
rows$thalasso_trait_0 <- rows$EuDiS_matched_species-rows$thalasso_trait_1
rows$trait_coverage <- ifelse(rows$WCVP_accepted_species>0,
          rows$EuDiS_matched_species/rows$WCVP_accepted_species,NA_real_)
rows$thalasso_trait_fraction_among_EuDiS_matched <-
    ifelse(rows$EuDiS_matched_species>0,
           rows$thalasso_trait_1/rows$EuDiS_matched_species,NA_real_)
rows$class <- geo$class[match(rows$area_code_l3,geo$area_code_l3)]
rows$level1 <- geo$level1[match(rows$area_code_l3,geo$area_code_l3)]
rows$eligible <- rows$EuDiS_matched_species>=20L & rows$trait_coverage>=0.05
write.csv(rows,file.path(a[[5]],"wcvp_region_origin_eudis_thalasso_coverage.csv"),
          row.names=FALSE,na="")

summary_rows <- do.call(rbind,lapply(c("island","mainland"),function(g) {
  do.call(rbind,lapply(0:1,function(k) {
    z <- rows[rows$class==g & rows$introduced==k,]
    e <- z[z$eligible,,drop=FALSE]
    data.frame(group=g,origin=if(k==0L) "native" else "introduced",
      source_regions=nrow(z), eligible_regions=nrow(e),
      total_WCVP_species_region=sum(z$WCVP_accepted_species),
      total_EuDiS_matched_species_region=sum(z$EuDiS_matched_species),
      median_EuDiS_coverage_eligible=if(nrow(e)>0) median(e$trait_coverage) else NA_real_,
      median_thalasso_trait_share_eligible=if(nrow(e)>0)
        median(e$thalasso_trait_fraction_among_EuDiS_matched) else NA_real_,
      stringsAsFactors=FALSE)
  }))
}))
write.csv(summary_rows,file.path(a[[5]],"thalasso_origin_cohort_coverage_summary.csv"),
          row.names=FALSE,na="")
write.csv(data.frame(binomial=names(trait_lookup),thalassochorous=unname(trait_lookup)),
          file.path(a[[5]],"eudis_pinned_binomial_trait_crosswalk.csv"),
          row.names=FALSE,na="")
cat("EUDIS_WCVP_SOURCE_GATE",
  "accepted_WCVP",nrow(names_ok),
  "pinned_EuDiS_rows",nrow(eu),
  "unique_EuDiS_binomials",length(trait_lookup),
  "thalassochorous_traits",sum(trait_lookup==1L,na.rm=TRUE),
  "original_port_flora_matched_regions",length(shared),
  "island_native_eligible",summary_rows$eligible_regions[summary_rows$group=="island" & summary_rows$origin=="native"],
  "island_introduced_eligible",summary_rows$eligible_regions[summary_rows$group=="island" & summary_rows$origin=="introduced"],
  "mainland_native_eligible",summary_rows$eligible_regions[summary_rows$group=="mainland" & summary_rows$origin=="native"],
  "mainland_introduced_eligible",summary_rows$eligible_regions[summary_rows$group=="mainland" & summary_rows$origin=="introduced"],"\n")
cat("SOURCE_BOUNDARY EuDiS morphology is not observed drift; European-only source creates nonrandom global trait missingness; no sea currents or vessel OD evaluated.\n")
