args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("usage: map_monarch_127_quality_wcvp.R <wcvp-root> <input-source-csv> <output-csv>")
root <- args[[1]]; input <- args[[2]]; output <- args[[3]]
e <- new.env(parent=emptyenv())
load(file.path(root,"data","wcvp_names.rda"),envir=e)
load(file.path(root,"data","wcvp_distributions.rda"),envir=e)
n <- get("wcvp_names",envir=e);d <- get("wcvp_distributions",envir=e)
require_names <- c("plant_name_id","taxon_rank","taxon_status","family","taxon_name","accepted_plant_name_id")
require_dist <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
if(!all(require_names%in%names(n)) || !all(require_dist%in%names(d)))stop("WCVP schema changed")
source <- read.csv(input,stringsAsFactors=FALSE)
if(nrow(source)!=127 || length(unique(source$original_species_name))!=127)stop("source table drift")
n$plant_name_id <- as.character(n$plant_name_id)
n$accepted_plant_name_id <- as.character(n$accepted_plant_name_id)
n$taxon_name <- as.character(n$taxon_name)
is_acc <- tolower(as.character(n$taxon_status))=="accepted"
no_id <- is.na(n$accepted_plant_name_id) | !nzchar(n$accepted_plant_name_id)
n$accepted_plant_name_id[no_id & is_acc] <- n$plant_name_id[no_id & is_acc]
by_id <- n[!duplicated(n$plant_name_id),c("plant_name_id","taxon_name","family")]
canon_name <- setNames(by_id$taxon_name,by_id$plant_name_id)
canon_fam <- setNames(as.character(by_id$family),by_id$plant_name_id)
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- as.character(d$area_code_l3)
d$introduced <- suppressWarnings(as.integer(as.character(d$introduced)))
d$extinct <- suppressWarnings(as.integer(as.character(d$extinct)))
d$location_doubtful <- suppressWarnings(as.integer(as.character(d$location_doubtful)))
d <- d[!is.na(d$extinct)&d$extinct==0L & !is.na(d$location_doubtful)&d$location_doubtful==0L &
       !is.na(d$area_code_l3)&nzchar(d$area_code_l3),,drop=FALSE]
results <- vector("list",nrow(source))
for(i in seq_len(nrow(source))){
  name <- source$original_species_name[i]
  hits <- n[tolower(trimws(as.character(n$taxon_rank)))=="species" & n$taxon_name==name,,drop=FALSE]
  # Exact accepted-name priority. Duplicate spelling rows can denote alternate
  # synonyms/homonyms; an actual unique accepted species record takes precedence.
  accepted_hits <- hits[
    tolower(trimws(as.character(hits$taxon_status)))=="accepted" &
    hits$plant_name_id==hits$accepted_plant_name_id,,drop=FALSE]
  accepted_ids <- unique(accepted_hits$plant_name_id)
  accepted_ids <- accepted_ids[!is.na(accepted_ids)&nzchar(accepted_ids)]
  ids <- unique(hits$accepted_plant_name_id)
  ids <- ids[!is.na(ids)&nzchar(ids)]
  if(length(accepted_ids)==1){
    id <- accepted_ids[[1]]
    status <- "matched"
    matched_as <- "accepted_name_prioritized"
  } else if(length(ids)==1){
    id <- ids[[1]]
    status <- "matched"
    matched_as <- "unambiguous_synonym"
  } else {
    id <- ""
    status <- if(length(ids)==0) "unmatched" else "ambiguous"
    matched_as <- ""
  }
  canonical <- if(nzchar(id) && id%in%names(canon_name))canon_name[[id]] else ""
  fam <- if(nzchar(id) && id%in%names(canon_fam))canon_fam[[id]] else ""
  dd <- if(nzchar(id)) d[d$plant_name_id==id,,drop=FALSE] else d[FALSE,,drop=FALSE]
  nat <- sort(unique(dd$area_code_l3[!is.na(dd$introduced)&dd$introduced==0L]))
  ali <- sort(unique(dd$area_code_l3[!is.na(dd$introduced)&dd$introduced==1L]))
  unk <- sort(unique(dd$area_code_l3[is.na(dd$introduced)]))
  results[[i]] <- data.frame(
    original_species_name=name,quality_subclass=source$quality_subclass[i],
    quality_bin=source$quality_bin[i],original_source_class_text=source$original_source_class_text[i],
    match_status=status,matched_as=matched_as,
    accepted_id=id,accepted_name=canonical,family=fam,
    native_regions=length(nat),introduced_regions=length(ali),
    unknown_status_regions=length(unk),current_regions=length(unique(c(nat,ali,unk))),
    native_codes=paste(nat,collapse=";"),introduced_codes=paste(ali,collapse=";"),
    stringsAsFactors=FALSE)
}
ans <- do.call(rbind,results)
if(nrow(ans)!=127)stop("output panel drift")
control <- ans[ans$original_species_name=="Asclepias curassavica",,drop=FALSE]
if(nrow(control)!=1 || control$match_status!="matched" ||
   control$accepted_id!="500848" ||
   control$native_regions!=41 || control$introduced_regions!=97){
  stop("Asclepias curassavica independent quality and WCVP positive-control mismatch")
}
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
write.csv(ans,output,row.names=FALSE,na="")
print(table(ans$match_status,ans$quality_bin))
print(ans[ans$original_species_name %in% c("Asclepias curassavica","Araujia sericifera"),
          c("original_species_name","quality_subclass","match_status","native_regions","introduced_regions")])
