args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: analyze_butterfly_pgls_match_bias.R <anthropogenic_csv> <output_json>")
}

suppressPackageStartupMessages({
  library(jsonlite)
  library(megatrees)
})

normalize_species <- function(x) {
  trimws(gsub("_", " ", gsub("\\*+$", "", as.character(x))))
}

anth <- read.csv(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("species", "host_family_count", "log_resource_expansion", "introduced_added_units")
if (!all(required %in% names(anth))) {
  stop("anthropogenic table missing required columns")
}

data("tree_butterfly", package = "megatrees")
tree_species <- unique(normalize_species(tree_butterfly$tip.label))
anth$tree_match <- anth$species %in% tree_species

matched <- anth[anth$tree_match, , drop = FALSE]
unmatched <- anth[!anth$tree_match, , drop = FALSE]
if (nrow(matched) != 124) stop(paste("unexpected exact match count:", nrow(matched)))

wilcox_host <- wilcox.test(
  matched$host_family_count,
  unmatched$host_family_count,
  exact = FALSE
)
wilcox_expansion <- wilcox.test(
  matched$log_resource_expansion,
  unmatched$log_resource_expansion,
  exact = FALSE
)

expanded_table <- matrix(
  c(
    sum(matched$introduced_added_units > 0),
    sum(matched$introduced_added_units == 0),
    sum(unmatched$introduced_added_units > 0),
    sum(unmatched$introduced_added_units == 0)
  ),
  nrow = 2,
  byrow = TRUE
)
fisher_expanded <- fisher.test(expanded_table)

safe_spearman <- function(x, y) {
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  unname(cor(x, y, method = "spearman"))
}

payload <- list(
  schema = "chocho_butterfly_kawahara_match_bias_v0.1",
  status = "POSTHOC_PGLS_COVERAGE_DIAGNOSTIC",
  panel_species = nrow(anth),
  exact_tree_matches = nrow(matched),
  unmatched_species = nrow(unmatched),
  matched = list(
    median_host_family_count = median(matched$host_family_count),
    mean_host_family_count = mean(matched$host_family_count),
    median_log_resource_expansion = median(matched$log_resource_expansion),
    mean_log_resource_expansion = mean(matched$log_resource_expansion),
    expanded_fraction = mean(matched$introduced_added_units > 0),
    rho_host_family_vs_log_expansion = safe_spearman(
      matched$host_family_count,
      matched$log_resource_expansion
    )
  ),
  unmatched = list(
    median_host_family_count = median(unmatched$host_family_count),
    mean_host_family_count = mean(unmatched$host_family_count),
    median_log_resource_expansion = median(unmatched$log_resource_expansion),
    mean_log_resource_expansion = mean(unmatched$log_resource_expansion),
    expanded_fraction = mean(unmatched$introduced_added_units > 0),
    rho_host_family_vs_log_expansion = safe_spearman(
      unmatched$host_family_count,
      unmatched$log_resource_expansion
    )
  ),
  tests = list(
    wilcoxon_host_family_count_p = unname(wilcox_host$p.value),
    wilcoxon_log_resource_expansion_p = unname(wilcox_expansion$p.value),
    fisher_expanded_fraction_p = unname(fisher_expanded$p.value)
  ),
  interpretation_rule = paste(
    "PGLS coverage is not obviously selective with respect to the focal quantities",
    "if exact-match and unmatched species have similar host-family breadth and",
    "resource-expansion distributions."
  ),
  claim_boundary = paste(
    "This checks observed PGLS coverage bias in the two focal variables; it does",
    "not make the 124-species phylogenetic subset equivalent to the full panel."
  )
)

dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), args[[2]])
cat(toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
