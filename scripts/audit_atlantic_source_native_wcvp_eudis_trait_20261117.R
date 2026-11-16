#!/usr/bin/env Rscript
# Source-only WCVP accepted plants and EuDiS thalassochory at fixed
# AZO/CNY/CVI current corridors. A species shared by two regions does
# NOT prove route direction or natural drift. Never code unlisted species
# as a biological absence; report unreported_in_destination.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("usage: Rscript source_WCVP_repo EuDiS_v5.csv outdir")
wp<-args[[1]];eup<-args[[2]];outdir<-args[[3]]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
routes<-matrix(c("AZO","CNY","CNY","CVI","AZO","CVI"),
               ncol=2,byrow=TRUE,dimnames=list(NULL,c("source","destination")))
regions<-sort(unique(as.vector(routes)))
stopifnot(identical(regions,c("AZO","CNY","CVI")))

binomial<-function(x) {
 x<-trimws(as.character(x))
 m<-regexec("^([A-Z][A-Za-z-]+)[[:space:]]+([a-z][a-z-]+)(?:[[:space:]].*)?$",x,perl=TRUE)
 z<-regmatches(x,m)
 vapply(z,function(a) if(length(a)>2L)tolower(paste(a[2],a[3])) else NA_character_,character(1))
}
eu<-read.csv(eup,stringsAsFactors=FALSE,check.names=FALSE,sep=";",fileEncoding="UTF-8-BOM")
names(eu)<-gsub("[[:space:]]+","_",trimws(names(eu)))
stopifnot(nrow(eu)==9874L,all(c("Species_POWO","Thalassochorous","Further_clarification_needed") %in% names(eu)))
eu$key<-binomial(eu$Species_POWO)
eu$trait<-trimws(as.character(eu$Thalassochorous))
eu$trait[!eu$trait %in% c("0","1") |
         grepl("THA",as.character(eu$Further_clarification_needed),fixed=TRUE)]<-NA_character_
eu<-eu[!is.na(eu$key),c("key","trait"),drop=FALSE]
trait_by_bin<-vapply(split(eu$trait,eu$key),function(x) {
 u<-unique(x); if(length(u)!=1L || is.na(u[1])) NA_integer_ else as.integer(u[1])
},integer(1))
stopifnot(sum(trait_by_bin==1L,na.rm=TRUE)==190L)

e<-new.env(parent=emptyenv())
load(file.path(wp,"data","wcvp_names.rda"),envir=e)
load(file.path(wp,"data","wcvp_distributions.rda"),envir=e)
nm<-get("wcvp_names",envir=e);d<-get("wcvp_distributions",envir=e)
required_name<-c("plant_name_id","taxon_name","family","taxon_rank","taxon_status")
required_dist<-c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required_name %in% names(nm)),all(required_dist %in% names(d)))
ok<- !is.na(nm$taxon_rank) & !is.na(nm$taxon_status) &
 tolower(trimws(as.character(nm$taxon_rank)))=="species" &
 tolower(trimws(as.character(nm$taxon_status)))=="accepted"
nm<-nm[ok,c("plant_name_id","taxon_name","family"),drop=FALSE]
nm$plant_name_id<-as.character(nm$plant_name_id)
nm<-nm[!duplicated(nm$plant_name_id),]
stopifnot(nrow(nm)==360064L)
nm$key<-binomial(nm$taxon_name)
freq<-table(nm$key,useNA="no")
ambig<-names(freq)[freq>1L]
valid<-as.character(d$area_code_l3) %in% regions &
  as.character(d$plant_name_id) %in% nm$plant_name_id &
  !is.na(d$extinct) & as.integer(d$extinct)==0L &
  !is.na(d$location_doubtful) & as.integer(d$location_doubtful)==0L &
  !is.na(d$introduced) & as.integer(d$introduced) %in% c(0L,1L)
q<-unique(d[valid,c("plant_name_id","area_code_l3","introduced"),drop=FALSE])
q$plant_name_id<-as.character(q$plant_name_id)
q$area_code_l3<-as.character(q$area_code_l3)
q$introduced<-as.integer(q$introduced)
if(nrow(q)<2000L)stop("Implausibly small accepted flora in three Atlantic regions")
if(anyDuplicated(paste(q$plant_name_id,q$area_code_l3,sep="|")))
 stop("Contradictory native/introduced WCVP origin for same accepted taxon×region")
pos<-match(q$plant_name_id,nm$plant_name_id)
if(anyNA(pos))stop("WCVP name metadata missing")
q$key<-nm$key[pos];q$name<-as.character(nm$taxon_name[pos]);q$family<-as.character(nm$family[pos])
q$thalassochorous<-unname(trait_by_bin[match(q$key,names(trait_by_bin))])
q$thalassochorous[q$key %in% ambig]<-NA_integer_
observed_pos<-tapply(q$thalassochorous==1L & q$introduced==0L,
                      q$area_code_l3,sum,na.rm=TRUE)
expected_pos<-c(AZO=24L,CNY=51L,CVI=10L)
if(!all(as.integer(observed_pos[names(expected_pos)])==unname(expected_pos)))
 stop(paste("Pinned EuDiS native-thalassochory origin totals changed:",paste(names(observed_pos),observed_pos,collapse=" ")))

source_den<-q[q$introduced==0L,]
source_totals<-aggregate(rep(1L,nrow(source_den)),
                        by=list(source=source_den$area_code_l3),FUN=sum)
names(source_totals)[2]<-"all_native_accepted_WCVP_taxa"
pair_rows<-list()
species_rows<-list()
for(ii in seq_len(nrow(routes))) {
 src<-routes[ii,1];dst<-routes[ii,2]
 original<-q[q$area_code_l3==src & q$introduced==0L,]
 dest<-q[q$area_code_l3==dst,c("plant_name_id","introduced"),drop=FALSE]
 # Keep only explicit source-borne native plant names with EuDiS trait
 source<-original[!is.na(original$thalassochorous),,drop=FALSE]
 if(nrow(source)<20L)stop("Too few source-native EuDiS matched species")
 di<-match(source$plant_name_id,dest$plant_name_id)
 target_flag<-dest$introduced[di]
 target<-ifelse(is.na(di),"unreported_in_destination",
               ifelse(target_flag==0L,"destination_native","destination_introduced"))
 id<-paste(src,dst,sep="->")
 species_rows[[id]]<-data.frame(
   source=src,destination=dst,plant_name_id=source$plant_name_id,
   accepted_plant_name=source$name,family=source$family,
   source_origin="native",destination_origin=target,
   EuDiS_thalassochorous=source$thalassochorous,
   stringsAsFactors=FALSE)
 for(tr in 0:1) {
   w<-source$thalassochorous==tr
   if(sum(w)==0L)stop("Source-native EuDiS thalassochory stratum empty")
   t<-target[w]
   pair_rows[[paste(id,tr)]]<-data.frame(
     source=src,destination=dst,trait=tr,
     source_native_EuDiS_species=sum(w),
     destination_native_reported=sum(t=="destination_native"),
     destination_introduced_reported=sum(t=="destination_introduced"),
     destination_unreported_in_WCVP=sum(t=="unreported_in_destination"),
     destination_native_fraction=sum(t=="destination_native")/sum(w),
     destination_any_report_fraction=sum(t!="unreported_in_destination")/sum(w),
     source_native_accepted_total=
         source_totals$all_native_accepted_WCVP_taxa[match(src,source_totals$source)],
     stringsAsFactors=FALSE)
 }
}
tab<-do.call(rbind,pair_rows)
if(nrow(tab)!=6L)stop("All source pairs × trait rows not recovered")
stopifnot(all(tab$destination_native_reported+tab$destination_introduced_reported+
    tab$destination_unreported_in_WCVP==tab$source_native_EuDiS_species))
write.csv(tab,file.path(outdir,"wcvp_eudis_atlantic_source_native_destination_origin_counts.csv"),
          row.names=FALSE,na="")
species<-do.call(rbind,species_rows)
write.csv(species,file.path(outdir,"source_native_eudis_named_species_three_routes.csv"),
          row.names=FALSE,na="")
cat("ATLANTIC_EUDIS_WCVP_SOURCE_DESTINATION_REPORT",jsonlite_not_needed<-paste(
  apply(tab,1,function(z)paste(z[c("source","destination","trait","source_native_EuDiS_species",
             "destination_native_reported","destination_introduced_reported",
             "destination_unreported_in_WCVP")],collapse=":")),collapse=" "),"\n")
cat("INTERPRETATION species geographically shared but no dated colonization or route transport, no causal ship versus current claim\n")
