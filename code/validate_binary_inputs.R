#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

glo <- readRDS(file.path(base_dir, "data/processed/glo_input_binary.rds"))
assert_true(nrow(glo) == 249311, "Unexpected number of GlobDB MAGs.")
assert_true(all(glo$GH_richness > 0 & glo$ABC_richness > 0), "GlobDB contains a zero repertoire.")
assert_true(
  isTRUE(all.equal(glo$ABC_GH_ratio, glo$ABC_richness / glo$GH_richness, tolerance = 1e-12)),
  "GlobDB ratios do not reproduce from binary repertoire totals."
)

img <- readRDS(file.path(base_dir, "data/processed/img_input.rds"))
gh_cols <- names(img)[match("GH1", names(img)):match("GH98", names(img))]
uptake_cols <- names(img)[match("agaEIIA", names(img)):match("ynfM", names(img))]
img_features <- unlist(img[c(gh_cols, uptake_cols)], use.names = FALSE)
assert_true(all(is.na(img_features) | img_features %in% c(0, 1)), "IMG feature matrix is not binary.")
assert_true(
  all(rowSums(img[gh_cols], na.rm = TRUE) == img$GH_total_abundance),
  "IMG depolymerisation totals do not equal binary feature sums."
)
assert_true(
  all(rowSums(img[uptake_cols], na.rm = TRUE) == img$ABC_total_abundance),
  "IMG uptake totals do not equal binary feature sums."
)
assert_true(
  isTRUE(all.equal(img$ABC_GH_ratio, img$ABC_total_abundance / img$GH_total_abundance, tolerance = 1e-12)),
  "IMG ratios do not reproduce from binary repertoire totals."
)

loma <- read.csv(file.path(base_dir, "data/processed/loma_input.csv"), check.names = FALSE)
assert_true(all(loma$GH_total_abundance > 0 & loma$ABC_total_abundance > 0), "Loma contains a zero repertoire.")
assert_true(
  isTRUE(all.equal(loma$ABC_GH_ratio, loma$ABC_total_abundance / loma$GH_total_abundance, tolerance = 1e-12)),
  "Loma ratios do not reproduce from binary repertoire totals."
)

message("Binary input validation passed for GlobDB, IMG and Loma analysis-ready inputs.")
