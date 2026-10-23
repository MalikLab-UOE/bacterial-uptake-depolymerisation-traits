#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
root <- normalizePath(dirname(script_file), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
include_tree <- "--include-tree" %in% args

scripts <- c(
  "Figure_1.R",
  "Figure_2.R",
  "Figure_4.R",
  "Figure_5_panels_a_c.R",
  "Figure_5_family_panel.R",
  "Figure_5_panel_b.R",
  "Figure_5_panel_c_GH_genomes.R",
  "Figure_5_panel_c_GH_transcripts.R"
)

if (include_tree) {
  scripts <- append(
    scripts,
    c("Figure_3_tree_visualisation.R", "Figure_3_Pagel_lambda.R"),
    after = 2
  )
}

for (script in scripts) {
  message("Running ", script)
  status <- system2("Rscript", file.path(root, "code", script))
  if (status != 0) stop("Failed: ", script, call. = FALSE)
}

message("Requested figure analyses completed. Outputs are in results/.")
if (!include_tree) {
  message("Figure 3 was skipped by default. Use --include-tree on a high-memory system.")
}
