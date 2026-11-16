#!/usr/bin/env Rscript
# Accepted species origin × recipient WCVP status on original AZO,CNY,CVI
# and original EuDiS thalassochory. Records contemporary plant geography
# only: does not measure seed arrival, salinity tolerance or butterfly hosts.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop("Usage: Rscript script <WCVP-repo> <EuDiS.csv> <original-island-flora.csv> <outputdir>")
wcvp <- args[[1]]; eudis_path<-args[[2]]; flora_path<-args[[3]]; out<-args[[4]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
codes<-c("AZO","CNY","CVI")
edges<-data.frame(src=c("AZO","CNY","AZO","CNY","CVI","CVI"),
                  dst=c("CNY","CVI","CVI","AZO","CNY","AZO"))
normal_binom <- function(s){
  s<-trimws(as.character(s))
  m<-regexec("^([A-Z][A-Za-z-]+)[[:space:]]+([a-z][a-z-]+)(?:[[:space:]].*)?$",s,perl=TRUE)
  x<-regmatches(s,m)
  vapply(x,function(z)if(length(z)>2L)tolower(paste(z[2],z[3])) else NA_character_,character(1))
}
e<-read.csv(eudis_path,sep=";",check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
names(e)<-gsub("[[:space:]]+","_",trimws(names(e)))
if(nrow(e)!=9874L || !all(c("Species_POWO","Thalassochorous","Further_clarification_needed") %in% names(e)))
   stop("Published EuDiS v5 source has changed")
e$key<-normal_binom(e$Species_POWO)
e$trait<-trimws(as.character(e$Thalassochorous))
e$trait[!e$trait %in% c("0","1") |
  grepl("THA",as.character(e$Further_clarification_needed),fixed=TRUE)]<-NA_character_
e<-e[!is.na(e$key),c("key","trait")]
groups<-split(e$trait,e$key)
tl<-vapply(groups,function(x){
  u<-unique(x); if(length(u)!=1L || is.na(u[1])) NA_integer_ else as.integer(u[1])
},integer(1))
if(sum(tl==1L,na.rm=TRUE)!=190L)stop("EuDiS high-confidence thalassochory positive count changed")

source_env<-new.env(parent=emptyenv())
load(file.path(wcvp,"data","wcvp_names.rda"),envir=source_env)
load(file.path(wcvp,"data","wcvp_distributions.rda"),envir=source_env)
n<-get("wcvp_names",envir=source_env)
d<-get("wcvp_distributions",envir=source_env)
needed_names<-c("plant_name_id","taxon_rank","taxon_status","taxon_name")
needed_dist<-c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
if(!all(needed_names %in% names(n)) || !all(needed_dist %in% names(d)))stop("Pinned WCVP source schema changed")
valid_names <- !is.na(n$taxon_rank) & !is.na(n$taxon_status) &
 tolower(trimws(as.character(n$taxon_rank)))=="species" &
 tolower(trimws(as.character(n$taxon_status)))=="accepted"
nam<-unique(n[valid_names,c("plant_name_id","taxon_name")])
nam$plant_name_id<-as.character(nam$plant_name_id)
nam<-nam[!duplicated(nam$plant_name_id),]
if(nrow(nam)!=360064L)stop("Pinned WCVP accepted count changed")
nam$key<-normal_binom(nam$taxon_name)
freq<-table(nam$key,useNA="no")
ambiguous<-names(freq)[freq>1L]
dkeep<-as.character(d$area_code_l3) %in% codes &
 as.character(d$plant_name_id) %in% nam$plant_name_id &
 !is.na(d$extinct) & as.integer(d$extinct)==0L &
 !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L &
 !is.na(d$introduced) & as.integer(d$introduced) %in% c(0L,1L)
q<-unique(d[dkeep,c("plant_name_id","area_code_l3","introduced")])
q$plant_name_id<-as.character(q$plant_name_id)
q$area_code_l3<-as.character(q$area_code_l3)
q$introduced<-as.integer(q$introduced)
key<-paste(q$plant_name_id,q$area_code_l3,sep="|")
if(anyDuplicated(key))stop("Multiple WCVP origin assignments for same accepted species×region")
index<-match(q$plant_name_id,nam$plant_name_id)
if(anyNA(index))stop("WCVP region accepted ID lacks name row")
q$accepted_name<-nam$taxon_name[index]
q$bin<-nam$key[index]
q$thalassochorous<-unname(tl[match(q$bin,names(tl))])
q$thalassochorous[q$bin %in% ambiguous]<-NA_integer_

# Cross-verify all accepted WCVP flora within original source regions.
flora<-read.csv(flora_path,stringsAsFactors=FALSE)
ff<-flora[flora$code %in% codes,]
if(nrow(ff)!=3L)stop("Original source control must include exactly all three regions")
for(c in codes){
  z<-q[q$area_code_l3==c,]
  f<-ff[ff$code==c,]
  if(sum(z$introduced==0L)!=f$native || sum(z$introduced==1L)!=f$introduced)
     stop(paste("Original WCVP-native/alien flora source mismatch",c))
}
# New source-side cohort counts are required to match the prior original
# 22 EuDiS × WCVP joint region source audit.
expected_n_positive<-c(AZO=24L,CNY=51L,CVI=10L)
got<-sapply(codes,function(c)sum(q$area_code_l3==c & q$introduced==0L &
                                    q$thalassochorous==1L,na.rm=TRUE))
if(!identical(as.integer(got),as.integer(expected_n_positive[codes])))
   stop(paste("EuDiS source-native thalassochory numerator drift",paste(got,collapse=",")))

out_records<-list(); species_rows<-list();nrec<-0L
for(k in seq_len(nrow(edges))){
  src<-edges$src[k]; dst<-edges$dst[k]
  source<-q[q$area_code_l3==src & q$introduced==0L &
             !is.na(q$thalassochorous),]
  target<-q[q$area_code_l3==dst,c("plant_name_id","introduced")]
  target_stat<-target$introduced[match(source$plant_name_id,target$plant_name_id)]
  source$dst_status<-ifelse(is.na(target_stat),"WCVP_NOT_LISTED",
       ifelse(target_stat==0L,"WCVP_NATIVE","WCVP_INTRODUCED"))
  for(trait in 0:1){
    grp<-source[source$thalassochorous==trait,]
    cts<-table(factor(grp$dst_status,
       levels=c("WCVP_NATIVE","WCVP_INTRODUCED","WCVP_NOT_LISTED")))
    N<-nrow(grp)
    if(N==0L)stop("No source-native species with known EuDiS trait group")
    nrec<-nrec+1L
    out_records[[nrec]]<-data.frame(
      source=src,destination=dst,source_origin="WCVP_NATIVE",
      eudis_thalassochory=trait,source_native_trait_species=N,
      destination_WCVP_native=as.integer(cts["WCVP_NATIVE"]),
      destination_WCVP_introduced=as.integer(cts["WCVP_INTRODUCED"]),
      destination_WCVP_not_listed=as.integer(cts["WCVP_NOT_LISTED"]),
      destination_WCVP_total_listed=as.integer(cts["WCVP_NATIVE"]+cts["WCVP_INTRODUCED"]),
      destination_listed_share=(cts["WCVP_NATIVE"]+cts["WCVP_INTRODUCED"])/N,
      stringsAsFactors=FALSE)
  }
  # Only the outcome-neutral original-source trait1 species IDs and
  # destination flags, to allow independent verification of overlap.
  a<-source[source$thalassochorous==1L,
            c("plant_name_id","accepted_name","thalassochorous","dst_status")]
  a$source<-src; a$destination<-dst
  species_rows[[k]]<-a
}
tab<-do.call(rbind,out_records)
if(nrow(tab)!=12L || any(tab$destination_WCVP_native+tab$destination_WCVP_introduced+
 tab$destination_WCVP_not_listed!=tab$source_native_trait_species))
 stop("Source/destination aggregation not conserved")
write.csv(tab,file.path(out,"source_native_to_other_region_WCVP_by_EuDiS.csv"),
          row.names=FALSE,na="")
write.csv(do.call(rbind,species_rows),
   file.path(out,"thalasso1_source_native_plant_species_and_destination_flags.csv"),
   row.names=FALSE,na="")
# Comparator diff within each ordered source/destination pair, not a
# year, dispersal probability, or causal regression.
comparisons<-do.call(rbind,lapply(seq(1,nrow(tab),by=2),function(i){
  z<-tab[i:(i+1),]
  if(!identical(z$eudis_thalassochory,c(0L,1L)))
    stop("Trait row ordering changed")
  data.frame(source=z$source[1],destination=z$destination[1],
    trait1_source_n=z$source_native_trait_species[2],
    trait0_source_n=z$source_native_trait_species[1],
    trait1_dest_listed_share=z$destination_listed_share[2],
    trait0_dest_listed_share=z$destination_listed_share[1],
    trait1_minus_trait0_destination_listing_pp=
      100*(z$destination_listed_share[2]-z$destination_listed_share[1]),
    stringsAsFactors=FALSE)
}))
write.csv(comparisons,file.path(out,"source_native_trait1_vs_trait0_destination_overlap.csv"),
          row.names=FALSE,na="")
cat("ATLANTIC_WCVP_SOURCE_DEST_ORIGIN_CROSSWALK",paste(capture.output(
  print(tab,row.names=FALSE)),collapse=" "),"\n")
cat("ATLANTIC_EUDIS_PLANT_DEST_RECORD_DIAGNOSTIC",
  paste(capture.output(print(comparisons,row.names=FALSE)),collapse=" "),"\n")
cat("CLAIM_BOUNDARY Contemporary WCVP region cooccurrence, NOT historical seed transport, source-specific dispersal or butterfly colonization.\n")
