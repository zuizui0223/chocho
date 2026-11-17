#!/usr/bin/env Rscript
# Reconstruct ORIGINAL accepted WCVP v13 126-island native and introduced
# incidence matrices, then calculate cross-island whole-flora similarities.
# No trait selection, no HOSTS network, no ecological absence inference.
# WCVP species identity is matched by accepted plant_name_id.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)
 stop("usage: Rscript <pinned_wcvp_repo> <Q1_region_metrics.csv> <original_126_WPI_WCVP_counts.csv> <out_dir>")
wcvp<-args[[1]];q1path<-args[[2]];verify_path<-args[[3]];out<-args[[4]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(Matrix))

geo<-read.csv(q1path,stringsAsFactors=FALSE)
original<-read.csv(verify_path,stringsAsFactors=FALSE)
isles<-sort(as.character(geo$area_code_l3[geo$class=="island"]))
stopifnot(length(isles)==126L,any(isles=="CNY"),!any(isles=="RHO"))
stopifnot(nrow(original)==126L,identical(sort(as.character(original$code)),isles))

env<-new.env(parent=emptyenv())
load(file.path(wcvp,"data","wcvp_names.rda"),envir=env)
load(file.path(wcvp,"data","wcvp_distributions.rda"),envir=env)
nm<-get("wcvp_names",envir=env)
d0<-get("wcvp_distributions",envir=env)
stopifnot(all(c("plant_name_id","taxon_rank","taxon_status") %in% names(nm)),
 all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful") %in% names(d0)))
accepted<-tolower(trimws(as.character(nm$taxon_rank)))=="species" &
           tolower(trimws(as.character(nm$taxon_status)))=="accepted"
ids<-unique(as.character(nm$plant_name_id[accepted & !is.na(accepted)]))
ids<-ids[!is.na(ids) & nzchar(ids)]
stopifnot(length(ids)==360064L)
in_region<-as.character(d0$area_code_l3) %in% isles
q<-d0[in_region & as.character(d0$plant_name_id) %in% ids &
    !is.na(d0$extinct) & as.integer(d0$extinct)==0L &
    !is.na(d0$location_doubtful) & as.integer(d0$location_doubtful)==0L,
    c("plant_name_id","area_code_l3","introduced"),drop=FALSE]
q$plant_name_id<-as.character(q$plant_name_id)
q$area_code_l3<-as.character(q$area_code_l3)
q$introduced<-as.integer(q$introduced)
if(anyNA(q$introduced) || any(!q$introduced %in% c(0L,1L)))
 stop("Original WCVP source has unclassified regional native/introduction status")
q<-unique(q)
if(nrow(q)<100000L)stop("Implausibly few original WCVP accepted island species×region records")
if(anyDuplicated(paste(q$plant_name_id,q$area_code_l3,sep="|")))
 stop("Original WCVP species×island has conflicting origin statuses; abort")

# Sparse 126×360064 matrices, exact species identities, not EuDiS matches.
i<-match(q$area_code_l3,isles)
j<-match(q$plant_name_id,ids)
stopifnot(!anyNA(i),!anyNA(j))
mk<-function(mask) {
 Matrix::sparseMatrix(i=i[mask],j=j[mask],x=rep.int(1L,sum(mask)),
  dims=c(length(isles),length(ids)),dimnames=list(isles,NULL))
}
N<-mk(q$introduced==0L)
A<-mk(q$introduced==1L)
native_counts<-as.integer(Matrix::rowSums(N))
alien_counts<-as.integer(Matrix::rowSums(A))
both_counts<-native_counts+alien_counts
ix<-match(isles,original$code)
stopifnot(!anyNA(ix),all(native_counts==original$native[ix]),
 all(alien_counts==original$introduced[ix]),
 all(both_counts==original$total[ix]))
# Reproduce source original independently verified 126-island WCVP totals.
stopifnot(sum(native_counts)==184997L,sum(alien_counts)==39024L)

NN<-as.matrix(Matrix::tcrossprod(N))
AA<-as.matrix(Matrix::tcrossprod(A))
NA<-as.matrix(N %*% Matrix::t(A))
AN<-as.matrix(A %*% Matrix::t(N))
stopifnot(all(diag(NN)==native_counts),all(diag(AA)==alien_counts))
stopifnot(all(NN==t(NN)),all(AA==t(AA)),all(NA==t(AN)))
pair_indices<-which(upper.tri(NN),arr.ind=TRUE)
stopifnot(nrow(pair_indices)==choose(126L,2L))
ii<-pair_indices[,1];jj<-pair_indices[,2]
sharedN<-as.integer(NN[pair_indices])
sharedA<-as.integer(AA[pair_indices])
unionN<-native_counts[ii]+native_counts[jj]-sharedN
unionA<-alien_counts[ii]+alien_counts[jj]-sharedA
J<-function(shared,den) ifelse(den>0,shared/den,NA_real_)
Sor<-function(shared,n1,n2) ifelse((n1+n2)>0,2*shared/(n1+n2),NA_real_)
overlap<-function(shared,n1,n2) ifelse(pmin(n1,n2)>0,shared/pmin(n1,n2),NA_real_)
p<-data.frame(
 source_region=isles[ii],destination_region=isles[jj],
 tdwg_level1_source=original$level1[ix[ii]],
 tdwg_level1_destination=original$level1[ix[jj]],
 source_wpi_2023_ports=as.integer(original$wpi_port_count[ix[ii]]),
 dest_wpi_2023_ports=as.integer(original$wpi_port_count[ix[jj]]),
 source_native_WCVP_accepted=native_counts[ii],dest_native_WCVP_accepted=native_counts[jj],
 source_alien_WCVP_accepted=alien_counts[ii],dest_alien_WCVP_accepted=alien_counts[jj],
 native_shared_accepted_taxa=sharedN,native_union_accepted_taxa=unionN,
 native_Jaccard=J(sharedN,unionN),
 native_Sorensen=Sor(sharedN,native_counts[ii],native_counts[jj]),
 native_Simpson_overlap=overlap(sharedN,native_counts[ii],native_counts[jj]),
 alien_shared_accepted_taxa=sharedA,alien_union_accepted_taxa=unionA,
 alien_Jaccard=J(sharedA,unionA),
 alien_Sorensen=Sor(sharedA,alien_counts[ii],alien_counts[jj]),
 alien_Simpson_overlap=overlap(sharedA,alien_counts[ii],alien_counts[jj]),
 # Species taxon native in SOURCE and introduced in DESTINATION. This is
 # an origin contrast, NOT observed dated natural/man-made dispersal.
 native_in_source_and_alien_in_dest=as.integer(NA[pair_indices]),
 alien_in_source_and_native_in_dest=as.integer(AN[pair_indices]),
 same_TDWG_Level1=as.integer(as.character(original$level1[ix[ii]])==
                              as.character(original$level1[ix[jj]])),
 stringsAsFactors=FALSE)
stopifnot(nrow(p)==7875L,
 all(p$native_shared_accepted_taxa<=pmin(p$source_native_WCVP_accepted,p$dest_native_WCVP_accepted)),
 all(p$alien_shared_accepted_taxa<=pmin(p$source_alien_WCVP_accepted,p$dest_alien_WCVP_accepted)))
write.csv(p,file.path(out,"WCVP_126island_7875_native_alien_flora_similarity.csv"),
          row.names=FALSE,na="")
region<-data.frame(
 code=isles,level1=original$level1[ix],
 area_km2=original$area_km2[ix],
 WPI_2023_ports=original$wpi_port_count[ix],
 native_WCVP_accepted=native_counts,
 alien_WCVP_accepted=alien_counts,
 total_WCVP_accepted=both_counts,
 stringsAsFactors=FALSE)
write.csv(region,file.path(out,"original126_WCVP_native_alien_region_counts.csv"),row.names=FALSE)
receipt<-data.frame(
 metric=c("source_WCVP_accepted_species_global","source_WCVP_island_regions",
  "source_WCVP_native_species_region_records",
  "source_WCVP_alien_species_region_records",
  "source_WCVP_all_species_region_records",
  "source_distinct_undirected_island_pairs",
  "nonempty_native_similarity_pair_count",
  "nonempty_alien_similarity_pair_count",
  "distinct_source_level1_groups",
  "native_and_alien_status_conflict_count"),
 value=c(length(ids),length(isles),sum(native_counts),sum(alien_counts),
  sum(both_counts),nrow(p),sum(is.finite(p$native_Jaccard)),
  sum(is.finite(p$alien_Jaccard)),
  length(unique(region$level1)),0L),stringsAsFactors=FALSE)
write.csv(receipt,file.path(out,"WCVP_126island_similarity_source_receipt.csv"),row.names=FALSE)
cat("WCVP126_PAIRWISE_ACCEPTED_FLORA_SOURCE_VERIFIED",
 paste(paste(receipt$metric,receipt$value,sep="="),collapse=" "),"\n")
cat("WCVP126_PAIRWISE_SIMILARITY_DESCRIPTION",
 sprintf("native Jaccard median=%.6f alien Jaccard median=%.6f native positive %d alien positive %d",
 median(p$native_Jaccard,na.rm=TRUE),
 median(p$alien_Jaccard,na.rm=TRUE),
 sum(p$native_Jaccard>0,na.rm=TRUE),
 sum(p$alien_Jaccard>0,na.rm=TRUE)),"\n")
cat("CAUSAL_STATUS only contemporaneous accepted plant occurrence co-reporting, not origin-destination arrival or shipping\n")
