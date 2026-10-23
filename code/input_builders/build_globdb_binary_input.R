#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(tools)
})

# ---------------------------
# helpers
# ---------------------------
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) paste0(label, " not found: ") else "File not found: "
    stop(msg, path, call. = FALSE)
  }
}

make_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

to_numeric_cols <- function(df, cols) {
  if (length(cols) == 0) return(df)
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    )
}

# binary richness: count gene families with value > 0
add_binary_richness <- function(df, id_col, total_col) {
  gene_cols <- setdiff(names(df), id_col)
  
  if (length(gene_cols) == 0) {
    df[[total_col]] <- 0
    return(df)
  }
  
  df <- to_numeric_cols(df, gene_cols)
  
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(gene_cols),
        ~ dplyr::if_else(is.na(.x), 0, dplyr::if_else(.x > 0, 1, 0))
      )
    ) %>%
    dplyr::mutate(
      "{total_col}" := rowSums(
        dplyr::across(dplyr::all_of(gene_cols)),
        na.rm = TRUE
      )
    )
}

# ---------------------------
# inputs
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

gh_merged_csv <- ifelse(length(args) >= 1, args[1], "Input_Data/GLO_CAZyme_GH_merged.csv")
gh_rules_xlsx <- ifelse(length(args) >= 2, args[2], "Input_Data/microtrait_GH.xlsx")

trans_matrix_csv <- ifelse(length(args) >= 3, args[3], "Input_Data/GLO_transporters_counts_matrix.csv")
abc_rules_xlsx   <- ifelse(length(args) >= 4, args[4], "Input_Data/microtrait_transporters.xlsx")
abc_sheet        <- ifelse(length(args) >= 5, as.integer(args[5]), 2)

tax_csv   <- ifelse(length(args) >= 6, args[6], "Input_Data/GLO_globdb_r226_taxonomy_1.csv")
stats_csv <- ifelse(length(args) >= 7, args[7], "Input_Data/GLO_globdb_r226_tax_plus_stats.csv")

min_comp <- ifelse(length(args) >= 8, as.numeric(args[8]), 50)
max_cont <- ifelse(length(args) >= 9, as.numeric(args[9]), 10)

out_rds <- ifelse(length(args) >= 10, args[10], "data/processed/glo_input_binary.rds")
out_csv <- ifelse(length(args) >= 11, args[11], "data/processed/glo_input_binary.csv")

message("=== Build GLO input: binary richness version ===")
message("GH merged csv:      ", gh_merged_csv)
message("GH rules xlsx:      ", gh_rules_xlsx)
message("Transporter matrix: ", trans_matrix_csv)
message("ABC rules xlsx:     ", abc_rules_xlsx, " (sheet=", abc_sheet, ")")
message("Taxonomy csv:       ", tax_csv)
message("Stats csv:          ", stats_csv)
message("QC: completeness >= ", min_comp, ", contamination <= ", max_cont)
message("Output RDS:         ", out_rds)
message("Output CSV:         ", out_csv)

stop_if_missing(gh_merged_csv, "CAZyme_GH_merged.csv")
stop_if_missing(gh_rules_xlsx, "microtrait_GH.xlsx")
stop_if_missing(trans_matrix_csv, "transporters_counts_matrix.csv")
stop_if_missing(abc_rules_xlsx, "microtrait_transporters.xlsx")
stop_if_missing(tax_csv, "taxonomy csv")
stop_if_missing(stats_csv, "stats csv")

# ---------------------------
# 1) GH binary richness
# ---------------------------
hmm_data  <- readr::read_csv(gh_merged_csv, show_col_types = FALSE)
gh_traits <- readxl::read_excel(gh_rules_xlsx)

if (!("id" %in% names(hmm_data))) stop("Column 'id' missing in GH merged CSV.", call. = FALSE)
if (!("GH" %in% names(gh_traits))) stop("Column 'GH' missing in GH rules xlsx.", call. = FALSE)

hmm_data <- hmm_data %>%
  dplyr::mutate(id = as.character(id))

gh_genes <- unique(stats::na.omit(gh_traits$GH))
gh_cols  <- intersect(names(hmm_data), gh_genes)

if (length(gh_cols) == 0) {
  stop("No GH columns matched between GH list and GH merged CSV.", call. = FALSE)
}

gh_only <- hmm_data %>%
  dplyr::select(id, dplyr::all_of(gh_cols)) %>%
  add_binary_richness(id_col = "id", total_col = "GH_richness") %>%
  dplyr::select(id, GH_richness)

message("Matched GH cols: ", length(gh_cols))
message("GH rows: ", nrow(gh_only))

# ---------------------------
# 2) ABC binary richness
# ---------------------------
micro_trans <- readr::read_csv(trans_matrix_csv, show_col_types = FALSE)
abc_traits  <- readxl::read_excel(abc_rules_xlsx, sheet = abc_sheet)

if (!("Sample" %in% names(micro_trans))) {
  stop("Column 'Sample' missing in transporter matrix.", call. = FALSE)
}

micro_trans <- micro_trans %>%
  dplyr::mutate(Sample = as.character(Sample))

abc_colname <- if ("KO" %in% names(abc_traits)) {
  "KO"
} else if ("ABC" %in% names(abc_traits)) {
  "ABC"
} else {
  NA_character_
}

if (is.na(abc_colname)) {
  stop("ABC rules sheet must contain a column named 'KO' or 'ABC'.", call. = FALSE)
}

abc_genes <- unique(stats::na.omit(abc_traits[[abc_colname]]))
abc_cols  <- intersect(names(micro_trans), abc_genes)

if (length(abc_cols) == 0) {
  stop("No ABC/KO columns matched between transporter matrix and ABC list.", call. = FALSE)
}

abc_only <- micro_trans %>%
  dplyr::select(Sample, dplyr::all_of(abc_cols)) %>%
  add_binary_richness(id_col = "Sample", total_col = "ABC_richness") %>%
  dplyr::select(Sample, ABC_richness) %>%
  dplyr::mutate(
    Sample = gsub("_KOfam_accession_count$", "", Sample)
  )

message("Matched ABC/KO cols: ", length(abc_cols))
message("ABC rows: ", nrow(abc_only))

# ---------------------------
# 3) Merge GH + ABC richness
# ---------------------------
merged_binary <- dplyr::full_join(
  gh_only,
  abc_only,
  by = c("id" = "Sample")
) %>%
  dplyr::mutate(
    GH_richness  = dplyr::coalesce(GH_richness, 0),
    ABC_richness = dplyr::coalesce(ABC_richness, 0)
  )

message("Merged binary rows: ", nrow(merged_binary))

# ---------------------------
# 4) Join taxonomy + stats
# ---------------------------
tax_meta  <- readr::read_csv(tax_csv, show_col_types = FALSE)
stat_meta <- readr::read_csv(stats_csv, show_col_types = FALSE)

if (!("MAG" %in% names(tax_meta))) {
  stop("Column 'MAG' missing in taxonomy file.", call. = FALSE)
}

if (!("ID" %in% names(stat_meta))) {
  stop("Column 'ID' missing in stats file.", call. = FALSE)
}

tax_meta <- tax_meta %>%
  dplyr::mutate(MAG = as.character(MAG))

stat_meta <- stat_meta %>%
  dplyr::mutate(ID = as.character(ID))

merged_table <- dplyr::full_join(
  merged_binary,
  tax_meta,
  by = c("id" = "MAG")
)

merged_table <- dplyr::full_join(
  merged_table,
  stat_meta,
  by = c("id" = "ID")
)

# ---------------------------
# 5) QC filter + bacteria + richness ratio
# ---------------------------
need_qc <- c("checkm2_completeness", "checkm2_contamination")
miss_qc <- setdiff(need_qc, names(merged_table))

if (length(miss_qc) > 0) {
  stop("QC columns missing: ", paste(miss_qc, collapse = ", "), call. = FALSE)
}

if (!("Kingdom" %in% names(merged_table))) {
  stop("Column 'Kingdom' missing (need d__Bacteria filter).", call. = FALSE)
}

merged_table <- merged_table %>%
  dplyr::filter(
    checkm2_completeness >= min_comp,
    checkm2_contamination <= max_cont,
    Kingdom == "d__Bacteria",
    GH_richness > 0,
    ABC_richness > 0
  ) %>%
  dplyr::mutate(
    ABC_GH_ratio = ABC_richness / GH_richness,
    ABC_GH_ratio = dplyr::if_else(
      is.finite(ABC_GH_ratio),
      ABC_GH_ratio,
      NA_real_
    )
  ) %>%
  dplyr::filter(is.finite(ABC_GH_ratio))

message("Final rows after QC + bacteria + nonzero richness: ", nrow(merged_table))

# ---------------------------
# save
# ---------------------------
make_dir(dirname(out_rds))
make_dir(dirname(out_csv))

saveRDS(merged_table, out_rds)
readr::write_csv(merged_table, out_csv)

message("Done ✅")