args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: analyze_butterfly_pgls_kawahara.R <leptraits_csv> <anthropogenic_csv> <output_json>")
}
leptraits_path <- args[[1]]
anth_path <- args[[2]]
output_path <- args[[3]]

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
  library(jsonlite)
  library(megatrees)
})

normalize_species <- function(x) {
  trimws(gsub("_", " ", gsub("\\*+$", "", as.character(x))))
}

fit_pgls <- function(tree, dat, label) {
  tip_species <- normalize_species(tree$tip.label)
  if (anyDuplicated(tip_species)) {
    stop(paste("duplicated normalized tips in", label))
  }
  rownames(dat) <- dat$species
  keep <- tip_species[tip_species %in% dat$species]
  if (length(keep) < 30) {
    stop(paste("too few matched species in", label, length(keep)))
  }

  tip_lookup <- setNames(tree$tip.label, tip_species)
  keep_tip <- unname(tip_lookup[keep])
  sub_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tip))
  ordered_species <- normalize_species(sub_tree$tip.label)
  d <- dat[ordered_species, , drop = FALSE]
  d$tree_tip <- sub_tree$tip.label

  rank_x <- rank(d$host_family_count, ties.method = "average")
  rank_y <- rank(d$log_resource_expansion, ties.method = "average")
  d$x_rank_z <- as.numeric(scale(rank_x))
  d$y_rank_z <- as.numeric(scale(rank_y))
  d$x_log <- log1p(d$host_family_count)

  # Some dated-tree branches are zero or effectively zero. Replace only
  # non-positive lengths by a negligible positive value to keep the published
  # topology while avoiding singular covariance matrices in GLS.
  positive_lengths <- sub_tree$edge.length[sub_tree$edge.length > 0]
  if (!length(positive_lengths)) stop("tree has no positive branch lengths")
  epsilon <- min(positive_lengths) * 1e-6
  sub_tree$edge.length[sub_tree$edge.length <= 0] <- epsilon

  fit_brownian_rank <- gls(
    y_rank_z ~ x_rank_z,
    data = d,
    correlation = corBrownian(
      value = 1,
      phy = sub_tree,
      form = ~tree_tip
    ),
    method = "ML",
    control = glsControl(opt = "optim", maxIter = 1000)
  )
  fit_brownian_raw <- gls(
    log_resource_expansion ~ x_log,
    data = d,
    correlation = corBrownian(
      value = 1,
      phy = sub_tree,
      form = ~tree_tip
    ),
    method = "ML",
    control = glsControl(opt = "optim", maxIter = 1000)
  )

  extract <- function(fit, term) {
    tt <- summary(fit)$tTable
    ci <- intervals(fit, which = "coef")$coef
    list(
      coefficient = unname(tt[term, "Value"]),
      standard_error = unname(tt[term, "Std.Error"]),
      t_value = unname(tt[term, "t-value"]),
      p_value = unname(tt[term, "p-value"]),
      ci95_lower = unname(ci[term, "lower"]),
      ci95_upper = unname(ci[term, "upper"])
    )
  }

  pagel_result <- tryCatch({
    fit_pagel <- gls(
      y_rank_z ~ x_rank_z,
      data = d,
      correlation = corPagel(
        value = 0.5,
        phy = sub_tree,
        form = ~tree_tip,
        fixed = FALSE
      ),
      method = "ML",
      control = glsControl(opt = "optim", maxIter = 1000)
    )
    list(
      fit = extract(fit_pagel, "x_rank_z"),
      lambda = unname(coef(fit_pagel$modelStruct$corStruct, unconstrained = FALSE)),
      success = TRUE
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })

  list(
    label = label,
    species = nrow(d),
    brownian_rank_pgls = extract(fit_brownian_rank, "x_rank_z"),
    brownian_raw_pgls = extract(fit_brownian_raw, "x_log"),
    pagel_rank_pgls = pagel_result,
    branch_length_epsilon = epsilon
  )
}

lep <- read.csv(leptraits_path, stringsAsFactors = FALSE, check.names = FALSE)
anth <- read.csv(anth_path, stringsAsFactors = FALSE, check.names = FALSE)
required_lep <- c("Family", "Genus", "Species")
if (!all(required_lep %in% names(lep))) stop("LepTraits taxonomy columns missing")

tax <- unique(lep[, required_lep])
tax <- tax[nzchar(tax$Species), , drop = FALSE]
names(tax) <- c("butterfly_family", "genus", "species")
dat <- merge(
  anth[, c("species", "host_family_count", "log_resource_expansion")],
  tax,
  by = "species",
  all.x = TRUE,
  sort = FALSE
)
if (any(!nzchar(dat$butterfly_family)) || any(!nzchar(dat$genus))) {
  stop("missing LepTraits genus/family")
}

data("tree_butterfly", package = "megatrees")
base_tree <- tree_butterfly
base_species <- normalize_species(base_tree$tip.label)
direct_species <- intersect(dat$species, base_species)
direct_tree <- drop.tip(
  base_tree,
  base_tree$tip.label[!base_species %in% direct_species]
)
direct_result <- fit_pgls(direct_tree, dat, "Kawahara_2023_exact_species")

payload <- list(
  schema = "chocho_butterfly_kawahara_pgls_v0.1",
  status = "POSTHOC_SPECIES_LEVEL_PHYLOGENETIC_SENSITIVITY",
  source = list(
    reference = "Kawahara et al. 2023 Nature Ecology & Evolution",
    doi = "10.1038/s41559-023-02041-9",
    tree_species = 2244,
    tree_source = "megatrees::tree_butterfly"
  ),
  resource_panel_species = nrow(dat),
  exact_tree_match_species = length(direct_species),
  exact_species_result = direct_result,
  coverage_fraction = length(direct_species) / nrow(dat),
  interpretation_rule = paste(
    "The near-zero host-breadth association is phylogenetically robust if the",
    "exact-species PGLS on taxa directly present in Kawahara et al. (2023)",
    "remains small and statistically unsupported."
  ),
  claim_boundary = paste(
    "This analysis deliberately excludes panel species absent from the published",
    "Kawahara et al. (2023) tree rather than imputing their phylogenetic placement."
  )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), output_path)
cat(toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
