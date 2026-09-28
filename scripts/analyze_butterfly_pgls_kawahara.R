args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: analyze_butterfly_pgls_kawahara.R <leptraits_csv> <anthropogenic_csv> <output_json>")
}
leptraits_path <- args[[1]]
anth_path <- args[[2]]
output_path <- args[[3]]

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(megatrees)
})

normalize_species <- function(x) {
  trimws(gsub("_", " ", gsub("\\*+$", "", as.character(x))))
}

gls_from_covariance <- function(y, x, V, lambda = 1.0) {
  n <- length(y)
  X <- cbind(Intercept = 1, predictor = x)
  p <- ncol(X)

  # Pagel lambda rescales off-diagonal shared-history covariance while
  # preserving each tip's marginal variance.
  D <- diag(diag(V))
  C <- V - D
  Vl <- D + lambda * C

  # Numerical nugget is tiny relative to tree depth and only stabilizes
  # matrix factorization for near-identical covariance rows.
  nugget <- max(diag(Vl)) * 1e-10
  Vl <- Vl + diag(nugget, n)

  cholV <- chol(Vl)
  VinvX <- backsolve(cholV, forwardsolve(t(cholV), X))
  Vinvy <- backsolve(cholV, forwardsolve(t(cholV), y))
  XtVinvX <- crossprod(X, VinvX)
  beta <- solve(XtVinvX, crossprod(X, Vinvy))

  resid <- as.numeric(y - X %*% beta)
  Vinvr <- backsolve(cholV, forwardsolve(t(cholV), resid))
  rss <- sum(resid * Vinvr)
  sigma2_reml <- rss / (n - p)
  covariance_beta <- sigma2_reml * solve(XtVinvX)
  se <- sqrt(diag(covariance_beta))
  t_value <- beta / se
  p_value <- 2 * pt(abs(t_value), df = n - p, lower.tail = FALSE)
  crit <- qt(0.975, df = n - p)

  logdet <- 2 * sum(log(diag(cholV)))
  sigma2_ml <- rss / n
  loglik_ml <- -0.5 * (
    n * (log(2 * pi) + 1 + log(sigma2_ml)) + logdet
  )

  list(
    coefficient = unname(beta[2, 1]),
    standard_error = unname(se[2]),
    t_value = unname(t_value[2, 1]),
    p_value = unname(p_value[2, 1]),
    ci95_lower = unname(beta[2, 1] - crit * se[2]),
    ci95_upper = unname(beta[2, 1] + crit * se[2]),
    sigma2_reml = unname(sigma2_reml),
    loglik_ml = unname(loglik_ml),
    lambda = lambda,
    nugget = nugget
  )
}

optimize_lambda <- function(y, x, V) {
  objective <- function(lambda) {
    out <- tryCatch(
      gls_from_covariance(y, x, V, lambda = lambda),
      error = function(e) NULL
    )
    if (is.null(out) || !is.finite(out$loglik_ml)) return(1e30)
    -out$loglik_ml
  }
  opt <- optimize(objective, interval = c(0, 1), tol = 1e-7)
  fit <- gls_from_covariance(y, x, V, lambda = opt$minimum)
  fit$optimization_objective <- opt$objective
  fit
}

fit_pgls <- function(tree, dat, label) {
  tip_species <- normalize_species(tree$tip.label)
  if (anyDuplicated(tip_species)) {
    stop(paste("duplicated normalized tips in", label))
  }

  matched <- intersect(dat$species, tip_species)
  if (length(matched) < 30) {
    stop(paste("too few matched species in", label, length(matched)))
  }

  keep_original <- tree$tip.label[tip_species %in% matched]
  sub_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_original))

  # Normalize tip labels after pruning so tree and data use identical names.
  sub_tree$tip.label <- normalize_species(sub_tree$tip.label)
  if (anyDuplicated(sub_tree$tip.label)) stop("normalized pruned tips duplicated")

  if (is.null(sub_tree$edge.length)) stop("tree has no branch lengths")
  if (any(!is.finite(sub_tree$edge.length))) stop("non-finite tree branch lengths")
  positive <- sub_tree$edge.length[sub_tree$edge.length > 0]
  if (!length(positive)) stop("tree has no positive branch lengths")
  epsilon <- min(positive) * 1e-6
  sub_tree$edge.length[sub_tree$edge.length <= 0] <- epsilon

  rownames(dat) <- dat$species
  ordered_species <- sub_tree$tip.label
  d <- dat[ordered_species, , drop = FALSE]
  if (any(is.na(d$species))) stop("tree/data ordering failed")

  V <- vcv.phylo(sub_tree, corr = FALSE)
  V <- V[ordered_species, ordered_species, drop = FALSE]
  if (any(!is.finite(V))) stop("non-finite phylogenetic covariance")

  rank_x <- rank(d$host_family_count, ties.method = "average")
  rank_y <- rank(d$log_resource_expansion, ties.method = "average")
  x_rank_z <- as.numeric(scale(rank_x))
  y_rank_z <- as.numeric(scale(rank_y))
  x_log <- log1p(d$host_family_count)
  y_raw <- d$log_resource_expansion

  brownian_rank <- gls_from_covariance(y_rank_z, x_rank_z, V, lambda = 1)
  brownian_raw <- gls_from_covariance(y_raw, x_log, V, lambda = 1)
  pagel_rank <- optimize_lambda(y_rank_z, x_rank_z, V)

  list(
    label = label,
    species = nrow(d),
    matched_species = ordered_species,
    brownian_rank_pgls = brownian_rank,
    brownian_raw_pgls = brownian_raw,
    pagel_rank_pgls = pagel_rank,
    branch_length_epsilon = epsilon,
    tree_height_range = range(node.depth.edgelength(sub_tree)[seq_len(Ntip(sub_tree))])
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
  schema = "chocho_butterfly_kawahara_pgls_v0.2",
  status = "POSTHOC_SPECIES_LEVEL_PHYLOGENETIC_SENSITIVITY",
  source = list(
    reference = "Kawahara et al. 2023 Nature Ecology & Evolution",
    doi = "10.1038/s41559-023-02041-9",
    published_tree_species = 2244,
    tree_source = "megatrees::tree_butterfly"
  ),
  resource_panel_species = nrow(dat),
  exact_tree_match_species = length(direct_species),
  coverage_fraction = length(direct_species) / nrow(dat),
  exact_species_result = direct_result,
  interpretation_rule = paste(
    "The near-zero host-breadth association is phylogenetically robust if both",
    "Brownian and estimated-Pagel-lambda species-level GLS estimates remain small",
    "and statistically unsupported on exact species matches."
  ),
  claim_boundary = paste(
    "Panel species absent from the published tree are excluded rather than",
    "phylogenetically imputed."
  )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(
  toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15),
  output_path
)
cat(toJSON(payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
