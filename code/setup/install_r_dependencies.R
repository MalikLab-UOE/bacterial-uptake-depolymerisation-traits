#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "ape", "cowplot", "data.table", "dplyr", "effectsize", "emmeans",
  "ggnewscale", "ggplot2", "ggtext", "mgcv", "multcompView",
  "openxlsx", "patchwork", "phylolm", "purrr", "ragg", "readr",
  "readxl", "rstatix", "sandwich", "scales", "stringr", "tibble",
  "tidyr", "tidyverse"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioconductor_packages <- c("ggtree", "ggtreeExtra")
missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (!requireNamespace("plotbiomes", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("valentinitnelav/plotbiomes")
}

message("R dependencies are installed.")
