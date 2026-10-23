#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
root <- normalizePath(dirname(script_file), mustWork = TRUE)

scripts <- c(
  "Figure_2.R",
  "Supplementary_Figure_S2.R",
  "Supplementary_Figure_S3.R"
)

for (script in scripts) {
  message("Running ", script)
  status <- system2("Rscript", file.path(root, "code", script))
  if (status != 0) stop("Failed: ", script, call. = FALSE)
}

message("Supplementary Figures S1-S3 completed. Outputs are in results/figures/.")
