#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(phylolm)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
tree_file <- file.path(base_dir, "data/raw/bacteria_globdb_r226.tree")
data_file <- file.path(base_dir, "data/processed/glo_input_binary.rds")
out_dir <- file.path(base_dir, "results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading the full phylogeny and binary repertoire data...")
tree <- read.tree(tree_file)
if (is.null(tree$edge.length)) stop("Tree branch lengths are required.")
nonpositive_edges <- sum(tree$edge.length <= 0)
edge_epsilon <- min(tree$edge.length[tree$edge.length > 0]) / 2
if (nonpositive_edges > 0) {
  tree$edge.length[tree$edge.length <= 0] <- edge_epsilon
  message(sprintf(
    "Regularised %d non-positive branch lengths to %.3g.",
    nonpositive_edges,
    edge_epsilon
  ))
}

dat <- as.data.table(readRDS(data_file))[, .(id, GH_richness, ABC_richness)]
setnames(dat, "id", "ID")
dat[, `:=`(
  GH_richness = as.numeric(GH_richness),
  ABC_richness = as.numeric(ABC_richness)
)]
dat <- dat[
  !is.na(ID) & ID != "" &
    is.finite(GH_richness) & GH_richness > 0 &
    is.finite(ABC_richness) & ABC_richness > 0
]
dat[, log10_ratio := log10(ABC_richness / GH_richness)]
dat <- unique(dat[is.finite(log10_ratio), .(ID, log10_ratio)], by = "ID")

shared <- intersect(tree$tip.label, dat$ID)
tree_use <- keep.tip(tree, shared)
trait_data <- as.data.frame(dat[match(tree_use$tip.label, ID)])
rownames(trait_data) <- trait_data$ID

stopifnot(
  nrow(trait_data) == Ntip(tree_use),
  identical(rownames(trait_data), tree_use$tip.label),
  all(is.finite(trait_data$log10_ratio))
)

message(sprintf("Fitting Pagel's lambda to the full tree (%s tips)...",
                format(Ntip(tree_use), big.mark = ",")))
fit_lambda <- phylolm(
  log10_ratio ~ 1,
  data = trait_data,
  phy = tree_use,
  model = "lambda",
  full.matrix = FALSE
)

# Lambda = 0 is a star phylogeny and is equivalent to independent Gaussian
# residuals for this intercept-only model. Compare maximum likelihoods using
# the same response and fixed-effect structure.
fit_independent <- lm(log10_ratio ~ 1, data = trait_data)
ll_lambda <- as.numeric(logLik(fit_lambda))
ll_zero <- as.numeric(logLik(fit_independent))
lr_statistic <- 2 * (ll_lambda - ll_zero)
lr_p <- pchisq(lr_statistic, df = 1, lower.tail = FALSE)

result <- data.table(
  n = Ntip(tree_use),
  trait = "log10(uptake repertoire / depolymerisation repertoire)",
  pagel_lambda = unname(fit_lambda$optpar),
  logLik_lambda = ll_lambda,
  logLik_lambda0 = ll_zero,
  likelihood_ratio = lr_statistic,
  df = 1L,
  P_value = lr_p,
  mean_log10_ratio = mean(trait_data$log10_ratio),
  sd_log10_ratio = sd(trait_data$log10_ratio)
  ,nonpositive_edges_regularised = nonpositive_edges
  ,edge_length_epsilon = edge_epsilon
)

fwrite(result, file.path(out_dir, "continuous_ratio_full_tree_pagel_lambda.csv"))
saveRDS(fit_lambda, file.path(out_dir, "continuous_ratio_full_tree_pagel_lambda_fit.rds"))

print(result)
message("Saved full-tree result to: ", out_dir)
