#!/usr/bin/env Rscript
# A strictly source-audited eight-butterfly host-only reconstruction in Azores.
# Original NHM HOSTS associations (unlocalized), original fixed WCVP v13
# accepted species and native vs contemporary WGSRPD3 distribution.
# Both endemic Hipparchia share Festuca jubata; all taxa treated equally
# but original GEB239 6 natives had different fixed source-panel/trait filters.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript <original_WCVP_repo> <original_NHM_hosts_resource.csv> <output_dir>")
wrepo<-args[[1L]]
source_hosts<-args[[2L]]
out<-args[[3L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
species<-c("Colias crocea","Danaus plexippus","Hipparchia azorina",
           "Hipparchia miguelensis","Lampides boeticus","Pieris brassicae",
           "Vanessa atalanta","Vanessa cardui")
origin<-setNames(c("native","native","endemic","endemic","native","native","native","native"),species)
family<-setNames(c("Pieridae","Nymphalidae","Nymphalidae","Nymphalidae",
                    "Lycaenidae","Pieridae","Nymphalidae","Nymphalidae"),species)
w<-new.env(parent=emptyenv())
load(file.path(wrepo,"data","wcvp_names.rda"),envir=w)
load(file.path(wrepo,"data","wcvp_distributions.rda"),envir=w)
n<-get("wcvp_names",envir=w)
d<-get("wcvp_distributions",envir=w)
stopifnot(all(c("plant_name_id","taxon_name","taxon_rank","taxon_status",
                "accepted_plant_name_id","family")%in%names(n)))
stopifnot(all(c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")%in%names(d)))
h<-read.csv(source_hosts,check.names=FALSE,stringsAsFactors=FALSE,na.strings="")
if(nrow(h)!=140485L)stop("Original NHM mirror drift")
stopifnot(all(c("Insect Genus","Insect Species","Insect Family",
                "Hostplant Genus","Hostplant Species")%in%names(h)))
h$butterfly<-paste(trimws(h[["Insect Genus"]]),trimws(h[["Insect Species"]]))
h<-h[h$butterfly%in%species & !is.na(h[["Insect Family"]]),,drop=FALSE]
h<-h[as.character(h[["Insect Family"]])==unname(family[h$butterfly]),,drop=FALSE]
h$plant<-paste(trimws(as.character(h[["Hostplant Genus"]])),
               trimws(as.character(h[["Hostplant Species"]])))
h$plant<-trimws(h$plant)
valid_species<-function(x)grepl("^[A-Z][A-Za-z-]+ [a-z][a-z-]+$",x) &
                !grepl(" (sp|spp|indet|unknown|cf|aff)$",x)
h<-h[valid_species(h$plant),,drop=FALSE]
if(length(setdiff(species,unique(h$butterfly)))>0L)
  stop("At least one source butterfly lacks original named host species")
h<-unique(h[c("butterfly","plant")])
for(sp in c("Hipparchia azorina","Hipparchia miguelensis")){
  if(!identical(sort(h$plant[h$butterfly==sp]),"Festuca jubata"))
    stop("A known host for endemic source changed; original source is one Festuca jubata per taxon")
}
n$plant_name_id<-as.character(n$plant_name_id)
n$accepted_plant_name_id<-as.character(n$accepted_plant_name_id)
n$taxon_rank<-tolower(trimws(as.character(n$taxon_rank)))
n$taxon_status<-tolower(trimws(as.character(n$taxon_status)))
n$taxon_name<-as.character(n$taxon_name)
accepted<-n[n$taxon_rank=="species"&n$taxon_status=="accepted",,drop=FALSE]
id_to_name<-setNames(accepted$taxon_name,accepted$plant_name_id)
# Align accepted WCVP species with the original GEB source script:
# exact source plant species binomial -> unique WCVP accepted_plant_name_id
matching<-n[n$taxon_rank=="species"&n$taxon_name%in%h$plant,
            c("plant_name_id","accepted_plant_name_id","taxon_status","taxon_name","family"),drop=FALSE]
matching$accepted_plant_name_id[matching$taxon_status=="accepted" &
  (is.na(matching$accepted_plant_name_id)|matching$accepted_plant_name_id=="")] <-
  matching$plant_name_id[matching$taxon_status=="accepted" &
  (is.na(matching$accepted_plant_name_id)|matching$accepted_plant_name_id=="")]
resolved<-lapply(sort(unique(h$plant)),function(plant){
    q<-matching[matching$taxon_name==plant,,drop=FALSE]
    ids<-sort(unique(q$accepted_plant_name_id[!is.na(q$accepted_plant_name_id)&
                                      nzchar(q$accepted_plant_name_id)]))
    accepted_id<-if(length(ids)==1L)ids[[1L]] else ""
    valid<-nzchar(accepted_id) && accepted_id%in%names(id_to_name)
    data.frame(input_plant=plant,source_candidate_ids=length(ids),
      status=if(length(ids)==0L)"unmatched" else if(length(ids)>1L)"ambiguous" else
          if(valid)"accepted_unique" else "accepted_id_unverified",
      accepted_plant_name_id=if(valid)accepted_id else "",
      accepted_plant_name=if(valid)unname(id_to_name[[accepted_id]]) else "",
      stringsAsFactors=FALSE)
})
resolved<-do.call(rbind,resolved)
d$plant_name_id<-as.character(d$plant_name_id)
distribution<-d[as.character(d$extinct)=="0" &
                 as.character(d$location_doubtful)=="0" &
                 as.character(d$introduced)%in%c("0","1") &
                 !is.na(d$area_code_l3) & nzchar(as.character(d$area_code_l3)),
                 c("plant_name_id","area_code_l3","introduced"),drop=FALSE]
distribution<-unique(distribution)
distribution$plant_name_id<-as.character(distribution$plant_name_id)
distribution$area_code_l3<-as.character(distribution$area_code_l3)
distribution$introduced<-as.character(distribution$introduced)
resolved_ok<-resolved[resolved$status=="accepted_unique",,drop=FALSE]
source_accepted<-merge(h,resolved_ok,by.x="plant",by.y="input_plant",all.x=TRUE,sort=FALSE)
if(nrow(source_accepted)!=nrow(h))stop("Lost original named host interactions in name join")
source_accepted$status[is.na(source_accepted$status)]<-"unmatched"
source_accepted$accepted_plant_name_id[is.na(source_accepted$accepted_plant_name_id)]<-""
out_rows<-lapply(species,function(sp){
  edges<-source_accepted[source_accepted$butterfly==sp,,drop=FALSE]
  ids<-unique(edges$accepted_plant_name_id[nzchar(edges$accepted_plant_name_id)])
  plant_regions<-distribution[distribution$plant_name_id%in%ids,,drop=FALSE]
  native<-sort(unique(plant_regions$area_code_l3[plant_regions$introduced=="0"]))
  contemporary<-sort(unique(plant_regions$area_code_l3[plant_regions$introduced%in%c("0","1")]))
  gained<-setdiff(contemporary,native)
  data.frame(
    independently_listed_Azores_butterfly=sp,
    independently_classified_Azores_origin=unname(origin[sp]),
    original_NHM_raw_named_host_species=length(unique(edges$plant)),
    accepted_unique_WCVP_host_species=length(ids),
    unresolved_original_host_species=length(unique(edges$plant[!nzchar(edges$accepted_plant_name_id)])),
    primary_native_regions_n=length(native),
    contemporary_host_regions_n=length(contemporary),
    introduced_added_regions_n=length(gained),
    native_coverage_in_AZO=as.integer("AZO"%in%native),
    contemporary_over_native_ratio=if(length(native))length(contemporary)/length(native) else NA_real_,
    source_contains_real_larval_use_locality_evidence=0L,
    source_GEB339_included=as.integer(origin[sp]=="native"),
    source_GEB239_included=as.integer(origin[sp]=="native"),
    original_preselected_population_not_recomputed=1L,
    stringsAsFactors=FALSE)
})
out_rows<-do.call(rbind,out_rows)
if(nrow(out_rows)!=8L)stop("Source independent eight not complete")
if(!all(out_rows$accepted_unique_WCVP_host_species<=out_rows$original_NHM_raw_named_host_species))
  stop("Accepted host count exceeds NHM raw records")
if(!all(out_rows$contemporary_host_regions_n>=out_rows$primary_native_regions_n))
  stop("Modern potential host geography cannot exclude original native ranges")
native<-out_rows[out_rows$independently_classified_Azores_origin=="native",,drop=FALSE]
endemic<-out_rows[out_rows$independently_classified_Azores_origin=="endemic",,drop=FALSE]
if(nrow(native)!=6L || nrow(endemic)!=2L)stop("Original official source origin strata changed")
if(!identical(endemic$independently_listed_Azores_butterfly,c("Hipparchia azorina","Hipparchia miguelensis")))
  stop("Endemic source identities changed")
if(any(endemic$accepted_unique_WCVP_host_species!=1L))stop("Original endemic one host WCVP accession failed")
write.csv(out_rows,file.path(out,"Azores_independent8_host_only_butterfly_WGSRPD3.csv"),row.names=FALSE,na="")
write.csv(resolved,file.path(out,"Azores_independent8_original_HOSTS_WCVP_name_resolution.csv"),row.names=FALSE,na="")
write.csv(source_accepted,file.path(out,"Azores_independent8_source_HOSTS_WCVP_plant_species_pairs.csv"),row.names=FALSE,na="")
cat("ORIGINAL_SOURCE_AZORES_EIGHT_HOST_ONLY_RECONSTRUCTION_COMPLETE\n")
for(sp in species){
  x<-out_rows[out_rows$independently_listed_Azores_butterfly==sp,,drop=FALSE]
  cat(sprintf("AZORES_SOURCE_BUTTERFLY=%s,origin=%s,original_hosts=%d,resolved_hosts=%d,unresolved=%d,native_units=%d,contemporary=%d,added=%d,ratio=%s\n",
    sp,x$independently_classified_Azores_origin,x$original_NHM_raw_named_host_species,
    x$accepted_unique_WCVP_host_species,x$unresolved_original_host_species,
    x$primary_native_regions_n,x$contemporary_host_regions_n,x$introduced_added_regions_n,
    as.character(x$contemporary_over_native_ratio)))
}
cat(sprintf("AZORES_EIGHT_TOTAL_SOURCE_HOSTS=%d\n",nrow(h)))
cat(sprintf("AZORES_SOURCE_ONLY_NATIVE_ADDED_SUM=%d,ENDEMIC_ADDED_SUM=%d\n",
    sum(native$introduced_added_regions_n),sum(endemic$introduced_added_regions_n)))
cat("INTERPRETATION: original exact NHM species-level host-only geographic opportunity; cannot infer actual local larval use or selection mechanism\n")
