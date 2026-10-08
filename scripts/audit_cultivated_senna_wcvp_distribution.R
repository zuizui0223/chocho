args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("usage: audit_cultivated_senna_wcvp_distribution.R <frozen_WCVP_repo> <summary_csv> <region_csv>")
root <- args[[1]]
summary_csv <- args[[2]]
region_csv <- args[[3]]
e <- new.env(parent=emptyenv())
load(file.path(root,"data","wcvp_names.rda"),envir=e)
load(file.path(root,"data","wcvp_distributions.rda"),envir=e)
n <- get("wcvp_names",envir=e)
d <- get("wcvp_distributions",envir=e)
targets <- c("Senna polyphylla","Senna surattensis","Senna ligustrina","Senna alata")
present_names <- data.frame(
  taxon_name=as.character(n$taxon_name),
  plant_name_id=as.character(n$plant_name_id),
  taxon_rank=tolower(as.character(n$taxon_rank)),
  taxon_status=tolower(as.character(n$taxon_status)),
  stringsAsFactors=FALSE
)
selected <- present_names[
  present_names$taxon_name %in% targets &
  present_names$taxon_rank=="species" &
  present_names$taxon_status=="accepted",
  ,
  drop=FALSE
]
if(any(duplicated(selected$taxon_name))) stop("ambiguous accepted Senna name")
if(!all(targets %in% selected$taxon_name)) stop(paste("unresolved accepted names:",paste(setdiff(targets,selected$taxon_name),collapse=", ")))
result <- data.frame(
  plant_name=character(),accepted_wcvp_id=character(),
  total_mapped_regions=integer(),native_regions=integer(),introduced_regions=integer(),
  FLA_present=logical(),FLA_introduced=logical(),FLA_native=logical(),
  stringsAsFactors=FALSE
)
regions <- data.frame(plant_name=character(),accepted_wcvp_id=character(),area_code_l3=character(),
                      introduced=integer(),stringsAsFactors=FALSE)
for (nm in targets){
  id <- selected$plant_name_id[selected$taxon_name==nm][1]
  idx <- as.character(d$plant_name_id)==id &
    suppressWarnings(as.integer(as.character(d$extinct)))==0L &
    suppressWarnings(as.integer(as.character(d$location_doubtful)))==0L
  dd <- d[idx,,drop=FALSE]
  if(nrow(dd)>0) {
    dd$area_code_l3 <- trimws(as.character(dd$area_code_l3))
    dd$introduced <- suppressWarnings(as.integer(as.character(dd$introduced)))
    dd <- unique(dd[!is.na(dd$area_code_l3)&nzchar(dd$area_code_l3),c("area_code_l3","introduced")])
    regions <- rbind(regions,data.frame(plant_name=nm,accepted_wcvp_id=id,area_code_l3=dd$area_code_l3,
                                      introduced=dd$introduced,stringsAsFactors=FALSE))
  }
  fla <- dd$introduced[dd$area_code_l3=="FLA"]
  result <- rbind(result,data.frame(
    plant_name=nm,accepted_wcvp_id=id,
    total_mapped_regions=nrow(dd),native_regions=sum(dd$introduced==0L,na.rm=TRUE),
    introduced_regions=sum(dd$introduced==1L,na.rm=TRUE),
    FLA_present=length(fla)>0,FLA_introduced=any(fla==1L),FLA_native=any(fla==0L),
    stringsAsFactors=FALSE
  ))
}
dir.create(dirname(summary_csv),recursive=TRUE,showWarnings=FALSE)
write.csv(result,summary_csv,row.names=FALSE,na="")
write.csv(regions,region_csv,row.names=FALSE,na="")
print(result)
