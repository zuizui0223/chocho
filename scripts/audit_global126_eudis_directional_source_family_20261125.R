#!/usr/bin/env Rscript
# Source-matched WCVP x EuDiS, preserving original name/origin filters.
# Output only aggregated directional trait contrasts, never source species rows.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop("Usage: WCVP-repo EuDiS.csv island-flora.csv original-coverage.csv outdir")
wcvp<-args[1]; eudis_path<-args[2]; flora_path<-args[3]; coverage_path<-args[4]; out<-args[5]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
flora<-read.csv(flora_path,stringsAsFactors=FALSE)
codes<-sort(as.character(flora$code))
stopifnot(length(codes)==126L,!anyDuplicated(codes))
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
nam<-unique(n[valid_names,c("plant_name_id","taxon_name","family")])
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
q$family<-trimws(as.character(nam$family[index]))
q$family[is.na(q$family) | !nzchar(q$family)]<-"__UNKNOWN_FAMILY__"
q$bin<-nam$key[index]
q$thalassochorous<-unname(tl[match(q$bin,names(tl))])
q$thalassochorous[q$bin %in% ambiguous]<-NA_integer_


coverage<-read.csv(coverage_path,stringsAsFactors=FALSE)
stopifnot(nrow(coverage)==738L)
for(c in codes){
 z<-q[q$area_code_l3==c,]
 f<-flora[flora$code==c,]
 stopifnot(nrow(f)==1L,sum(z$introduced==0L)==f$native,sum(z$introduced==1L)==f$introduced)
 for(o in 0:1){
  old<-coverage[coverage$area_code_l3==c & coverage$introduced==o,]
  zz<-z[z$introduced==o,]
  stopifnot(nrow(old)==1L,nrow(zz)==old$WCVP_accepted_species,
   sum(!is.na(zz$thalassochorous))==old$EuDiS_matched_species,
   sum(zz$thalassochorous==1L,na.rm=TRUE)==old$thalasso_trait_1,
   sum(zz$thalassochorous==0L,na.rm=TRUE)==old$thalasso_trait_0)
 }
}
# Trait known accepted source-native taxa across all original botanical regions.
# Equal-source comparison via fixed family weighting: for each source and mixed
# family, weight by the number of positive-trait source species. Compare the
# same plant pool for native and introduced destination status.
known<-q[!is.na(q$thalassochorous),]
ids<-sort(unique(known$plant_name_id))
status<-matrix(-1L,nrow=length(ids),ncol=length(codes),dimnames=list(ids,codes))
status[cbind(match(known$plant_name_id,ids),match(known$area_code_l3,codes))]<-known$introduced
res<-list(); gates<-list()
for(src in codes){
 a<-known[known$area_code_l3==src & known$introduced==0L,]
 stopifnot(!anyDuplicated(a$plant_name_id))
 a$family[a$family=="__UNKNOWN_FAMILY__"]<-NA
 tt<-table(a$family,factor(a$thalassochorous,levels=0:1))
 fam<-rownames(tt)[tt[,1]>0 & tt[,2]>0]
 m<-a[!is.na(a$family) & a$family %in% fam,]
 pos<-sum(m$thalassochorous==1L); neg<-sum(m$thalassochorous==0L)
 eligible<-pos>=10L & neg>=10L & length(fam)>=3L
 old<-coverage[coverage$area_code_l3==src,]
 # Source frame membership is assigned separately from original paired CSV in Python.
 gates[[src]]<-data.frame(source=src,n_known=nrow(a),n_pos=sum(a$thalassochorous==1L),n_neg=sum(a$thalassochorous==0L),mixed_families=length(fam),matched_pos=pos,matched_neg=neg,eligible=eligible)
 if(!eligible) next
 # signed family-standardized weights sum +1 among trait1 and -1 among trait0.
 w<-numeric(nrow(m))
 for(f in fam){
  p<-which(m$family==f & m$thalassochorous==1L)
  n<-which(m$family==f & m$thalassochorous==0L)
  w[p]<-1/pos
  w[n]<- -length(p)/(pos*length(n))
 }
 stopifnot(abs(sum(w))<1e-9,abs(sum(w[w>0])-1)<1e-9)
 S<-status[match(m$plant_name_id,ids),,drop=FALSE]
 native<-as.numeric(crossprod(w,S==0L))
 alien<-as.numeric(crossprod(w,S==1L))
 stopifnot(all(abs(native)<=1+1e-10),all(abs(alien)<=1+1e-10))
 res[[src]]<-data.frame(source=src,destination=codes,matched_pos=pos,matched_neg=neg,mixed_families=length(fam),native_excess_pp=100*native,alien_excess_pp=100*alien,native_minus_alien_excess_pp=100*(native-alien))
}
g<-do.call(rbind,gates); stopifnot(nrow(g)==126L)
write.csv(g,file.path(out,"trait_source_coverage.csv"),row.names=FALSE)
z<-do.call(rbind,res); z<-z[z$source!=z$destination,]
stopifnot(!anyDuplicated(paste(z$source,z$destination)),nrow(z)==sum(g$eligible)*125L)
write.csv(z,file.path(out,"family_matched_directional_trait_contrasts.csv"),row.names=FALSE)
cat("WCVP_EUDIS_ALL126_SOURCE_GATES_PASS",nrow(g),"family_supported_sources",sum(g$eligible),"directed_trait_rows",nrow(z),"\n")
