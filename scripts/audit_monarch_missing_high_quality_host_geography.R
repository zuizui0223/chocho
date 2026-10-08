args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5) stop("usage: <WCVP repo> <HOSTS interaction CSV> <native WCVP sidecar> <current sidecar> <output directory>")
root<-args[[1]]; hosts_csv<-args[[2]]; native_csv<-args[[3]]; current_csv<-args[[4]]; outdir<-args[[5]]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
e<-new.env(parent=emptyenv())
load(file.path(root,"data","wcvp_names.rda"),envir=e)
load(file.path(root,"data","wcvp_distributions.rda"),envir=e)
nam<-get("wcvp_names",envir=e); dist<-get("wcvp_distributions",envir=e)
wanted<-"Asclepias curassavica"
sel<-nam[
 as.character(nam$taxon_name)==wanted &
 tolower(trimws(as.character(nam$taxon_rank)))=="species" &
 tolower(trimws(as.character(nam$taxon_status)))=="accepted",,drop=FALSE]
if(nrow(sel)!=1) stop("independent high performance Asclepias taxon unresolved")
new_id<-as.character(sel$plant_name_id[[1]])
dd<-dist[as.character(dist$plant_name_id)==new_id &
 suppressWarnings(as.integer(as.character(dist$extinct)))==0L &
 suppressWarnings(as.integer(as.character(dist$location_doubtful)))==0L,,drop=FALSE]
dd$area_code_l3<-trimws(as.character(dd$area_code_l3))
dd$introduced<-suppressWarnings(as.integer(as.character(dd$introduced)))
dd<-unique(dd[!is.na(dd$area_code_l3)&nzchar(dd$area_code_l3)&dd$introduced%in%c(0L,1L),c("area_code_l3","introduced"),drop=FALSE])
if(nrow(dd)==0) stop("new host has no pinned WCVP regions")
for(code in unique(dd$area_code_l3)){
 if(length(unique(dd$introduced[dd$area_code_l3==code]))!=1)stop(paste("plant origin ambiguous",code))
}
add_native<-sort(unique(dd$area_code_l3[dd$introduced==0L]))
add_current<-sort(unique(dd$area_code_l3))
links<-read.csv(hosts_csv,stringsAsFactors=FALSE,colClasses="character")
links<-links[links$insect_species=="Danaus plexippus",,drop=FALSE]
if(nrow(links)!=42) stop("monarch HOSTS pair count drift")
links<-unique(links[,c("accepted_plant_name_id","accepted_name"),drop=FALSE])
if(nrow(links)!=40)stop("monarch distinct accepted host count drift")
if(any(links$accepted_name==wanted))stop("proven high host already in frozen HOSTS")
disputed<-links$accepted_plant_name_id[links$accepted_name=="Gossypium arboreum"]
if(length(disputed)!=1)stop("Gossypium arboreum ambiguous or absent")
if(new_id%in%links$accepted_plant_name_id)stop("high host id already recorded under different name")
na<-read.csv(native_csv,stringsAsFactors=FALSE,colClasses="character")
co<-read.csv(current_csv,stringsAsFactors=FALSE,colClasses="character")
orig_native<-split(na$area_code_l3,na$accepted_plant_name_id)
orig_current<-split(co$area_code_l3,co$accepted_plant_name_id)
orig_native[[new_id]]<-add_native
orig_current[[new_id]]<-add_current
footprint<-function(ids,source){
 values<-unlist(source[ids],use.names=FALSE)
 sort(unique(values[!is.na(values)&nzchar(values)]))
}
calc<-function(ids){
 n<-footprint(ids,orig_native)
 c<-footprint(ids,orig_current)
 if(!all(n%in%c))stop("native opportunity not subset current")
 list(native=n,current=c,added=setdiff(c,n))
}
original<-links$accepted_plant_name_id
scenarios<-list(
 baseline=original,
 high_quality_omission_corrected=c(original,new_id),
 disputed_cotton_excluded=setdiff(original,disputed),
 both_changes=c(setdiff(original,disputed),new_id)
)
sets<-lapply(scenarios,calc)
base<-sets$baseline
if(length(base$native)!=181||length(base$current)!=263||length(base$added)!=82)stop("frozen 181/263/82 drift")
if(length(sets$disputed_cotton_excluded$added)!=79)stop("cotton exclusion regression drift")
table<-do.call(rbind,lapply(names(sets),function(nm){
 v<-sets[[nm]]
 data.frame(
 scenario=nm,
 native_regions=length(v$native),
 current_regions=length(v$current),
 introduced_only_regions=length(v$added),
 net_native_change=length(v$native)-length(base$native),
 net_current_change=length(v$current)-length(base$current),
 net_introduced_only_change=length(v$added)-length(base$added),
 originally_added_regions_reclassified_or_lost=length(setdiff(base$added,v$added)),
 newly_added_region_opportunities=length(setdiff(v$added,base$added)),
 Puerto_Rico_PUE_native=as.integer("PUE"%in%v$native),
 Puerto_Rico_PUE_current=as.integer("PUE"%in%v$current),
 stringsAsFactors=FALSE)
}))
region_codes<-sort(unique(unlist(lapply(sets,function(v)v$current))))
detail<-data.frame(area_code_l3=region_codes,stringsAsFactors=FALSE)
for(nm in names(sets)){
 x<-sets[[nm]]
 detail[[paste0(nm,"_native")]]<-as.integer(region_codes%in%x$native)
 detail[[paste0(nm,"_current")]]<-as.integer(region_codes%in%x$current)
 detail[[paste0(nm,"_introduced_only")]]<-as.integer(region_codes%in%x$added)
}
plant<-data.frame(accepted_plant_name_id=new_id,accepted_name=wanted,
   native_regions=length(add_native),
   contemporary_regions=length(add_current),
   introduced_regions=length(setdiff(add_current,add_native)),
   PUE_native=as.integer("PUE"%in%add_native),
   PUE_current=as.integer("PUE"%in%add_current),
   stringsAsFactors=FALSE)
write.csv(table,file.path(outdir,"monarch_four_quality_link_scenarios.csv"),row.names=FALSE,na="")
write.csv(detail,file.path(outdir,"monarch_quality_link_region_reclassifications.csv"),row.names=FALSE,na="")
write.csv(plant,file.path(outdir,"independently_high_performance_Asclepias_curassavica_WCVP.csv"),row.names=FALSE,na="")
write.csv(dd,file.path(outdir,"curassavica_native_and_introduced_WCVP_regions.csv"),row.names=FALSE,na="")
cat("SOURCE: PLOS ONE 2022 Table Greenstein et al DOI 10.1371/journal.pone.0269701 H3 high-quality A.curassavica\n")
cat("DISCLAIMER: Gossypium arboreum is disputed, not independently proven unsuitable; removal is counterfactual sensitivity.\n")
cat("DISCLAIMER: all counts concern regional potential resource geography, not observed larval success or butterfly occupancy.\n")
print(plant)
print(table)
