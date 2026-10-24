#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(mgcv)
  library(effectsize)
  library(patchwork)
  library(cowplot)
})

set.seed(20261015)

ggplot2::theme_set(
  ggplot2::theme_classic(base_size = 9.5) +
    ggplot2::theme(
      plot.title = element_text(
        face = "bold", size = 10.5, hjust = 0.5,
        margin = margin(b = 1)
      )
    )
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
in_rds <- file.path(base_dir, "data/processed/glo_input_binary.rds")
fig_dir <- file.path(base_dir, "results/figures")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(fig_dir, "Supplementary_Figure_S3.png")
out_pdf <- file.path(fig_dir, "Supplementary_Figure_S3.pdf")
out_stats <- file.path(tab_dir, "Genome_size_GAM_by_phylum_stats.csv")
out_extremes <- file.path(tab_dir, "Genome_size_top_bottom10_by_phylum.csv")

phylum_levels <- c(
  "Pseudomonadota",
  "Bacillota",
  "Bacteroidota",
  "Actinomycetota",
  "Acidobacteriota",
  "Chloroflexota",
  "Verrucomicrobiota"
)

dataset_levels <- c("All MAGs", phylum_levels)

strategy_levels <- c("Bottom 10%", "Top 10%")
delta_threshold <- 0.147

format_p <- function(p, prefix = "P") {
  if (p < 0.001) {
    paste(prefix, "< 0.001")
  } else {
    paste0(prefix, " = ", formatC(p, digits = 3, format = "f"))
  }
}

delta_magnitude <- function(delta) {
  case_when(
    abs(delta) < 0.147 ~ "negligible",
    abs(delta) < 0.330 ~ "small",
    abs(delta) < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

dat <- readRDS(in_rds) %>%
  transmute(
    Phylum = sub("^p__", "", as.character(Phylum)),
    genome_size_Mb = suppressWarnings(as.numeric(genome_size)) / 1e6,
    ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  filter(
    is.finite(genome_size_Mb),
    genome_size_Mb > 0,
    genome_size_Mb <= 10,
    is.finite(ratio),
    ratio > 0
  ) %>%
  mutate(
    log10_ratio = log10(ratio),
    Phylum = as.character(Phylum)
  ) %>%
  filter(is.finite(log10_ratio))

if (nrow(dat) == 0) stop("No valid phylum-level observations available.")

analyse_dataset <- function(dataset_name) {
  this_data <- if (dataset_name == "All MAGs") {
    dat
  } else {
    dat %>% filter(Phylum == dataset_name)
  }

  gam_fit <- gam(
    genome_size_Mb ~ s(log10_ratio, k = 10),
    data = this_data,
    method = "REML"
  )
  gam_summary <- summary(gam_fit)

  prediction_data <- tibble(
    log10_ratio = seq(
      quantile(this_data$log10_ratio, 0.005, na.rm = TRUE),
      quantile(this_data$log10_ratio, 0.995, na.rm = TRUE),
      length.out = 300
    )
  )
  prediction <- predict(gam_fit, newdata = prediction_data, se.fit = TRUE)
  prediction_data <- prediction_data %>%
    mutate(
      fit = as.numeric(prediction$fit),
      lower = fit - 1.96 * as.numeric(prediction$se.fit),
      upper = fit + 1.96 * as.numeric(prediction$se.fit)
    )

  gam_stats <- tibble(
    Dataset = dataset_name,
    n = nrow(this_data),
    edf = gam_summary$s.table[1, "edf"],
    adjusted_R2 = gam_summary$r.sq,
    deviance_explained = gam_summary$dev.expl,
    p_value = gam_summary$s.table[1, "p-value"]
  )

  gam_label <- paste0(
    "n = ", format(nrow(this_data), big.mark = ","),
    "\nedf = ", sprintf("%.2f", gam_stats$edf),
    "\nDeviance explained = ", sprintf("%.2f%%", 100 * gam_stats$deviance_explained),
    "\n", format_p(gam_stats$p_value)
  )

  point_alpha <- case_when(
    nrow(this_data) > 100000 ~ 0.12,
    nrow(this_data) > 20000 ~ 0.17,
    TRUE ~ 0.26
  )
  point_size <- case_when(
    nrow(this_data) > 100000 ~ 0.62,
    nrow(this_data) > 20000 ~ 0.68,
    TRUE ~ 0.76
  )

  p_gam <- ggplot(this_data, aes(log10_ratio, genome_size_Mb)) +
    geom_point(
      aes(color = log10_ratio),
      size = point_size,
      alpha = point_alpha,
      stroke = 0
    ) +
    geom_ribbon(
      data = prediction_data,
      aes(x = log10_ratio, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "grey65",
      alpha = 0.30
    ) +
    geom_line(
      data = prediction_data,
      aes(x = log10_ratio, y = fit),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.72
    ) +
    annotate(
      "text",
      x = -1.39,
      y = 9.88,
      label = gam_label,
      hjust = 0,
      vjust = 1,
      size = 2.55,
      lineheight = 1.0
    ) +
    scale_x_continuous(
      breaks = c(-1, -0.5, 0, 0.5, 1, 1.5),
      limits = c(-1.45, 1.58)
    ) +
    scale_y_continuous(
      breaks = c(0, 2.5, 5, 7.5, 10),
      limits = c(0, 10.2)
    ) +
    scale_color_gradient2(
      low = "#287A62",
      mid = "#F3F1EC",
      high = "#D9483B",
      midpoint = 0,
      limits = c(-2, 2),
      breaks = c(-2, -1, 0, 1, 2),
      oob = scales::squish,
      name = expression("Ratio ("*log[10]*")")
    ) +
    labs(
      title = NULL,
      x = expression("Uptake-to-depolymerisation ratio ("*log[10]*")"),
      y = "Genome size (Mb)"
    ) +
    theme_classic(base_size = 9.5) +
    theme(
      axis.title = element_text(face = "bold", size = 9.5),
      axis.text = element_text(color = "black", size = 8),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.4),
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 7.5),
      plot.margin = margin(4, 5, 4, 4)
    )

  q10 <- quantile(this_data$ratio, 0.10, na.rm = TRUE)
  q90 <- quantile(this_data$ratio, 0.90, na.rm = TRUE)
  extreme_df <- this_data %>%
    mutate(
      strategy = case_when(
        ratio <= q10 ~ "Bottom 10%",
        ratio >= q90 ~ "Top 10%",
        TRUE ~ NA_character_
      ),
      strategy = factor(strategy, levels = strategy_levels)
    ) %>%
    filter(!is.na(strategy))

  wilcox_result <- wilcox.test(
    genome_size_Mb ~ strategy,
    data = extreme_df,
    exact = FALSE
  )
  delta_result <- effectsize::cliffs_delta(
    genome_size_Mb ~ strategy,
    data = extreme_df
  )
  delta_value <- as.numeric(delta_result$r_rank_biserial)
  magnitude <- delta_magnitude(delta_value)
  practical_difference <-
    wilcox_result$p.value < 0.05 && abs(delta_value) >= delta_threshold
  strategy_letters <- if (practical_difference) c("a", "b") else c("a", "a")

  extreme_summary <- extreme_df %>%
    group_by(strategy) %>%
    summarise(
      n = n(),
      median_log10_ratio = median(log10_ratio),
      median_genome_size_Mb = median(genome_size_Mb),
      q25_genome_size_Mb = quantile(genome_size_Mb, 0.25),
      q75_genome_size_Mb = quantile(genome_size_Mb, 0.75),
      .groups = "drop"
    ) %>%
    mutate(
      Dataset = dataset_name,
      q10_ratio = q10,
      q90_ratio = q90,
      wilcoxon_p_value = wilcox_result$p.value,
      cliff_delta = delta_value,
      delta_ci_low = as.numeric(delta_result$CI_low),
      delta_ci_high = as.numeric(delta_result$CI_high),
      delta_magnitude = magnitude,
      delta_threshold = delta_threshold,
      practical_difference = practical_difference,
      group_letter = strategy_letters
    )

  extreme_df <- extreme_df %>%
    left_join(
      extreme_summary %>% select(strategy, median_log10_ratio),
      by = "strategy"
    )

  count_df <- extreme_summary %>%
    mutate(
      label = paste0(
        "n = ", format(n, big.mark = ","),
        "\nmedian = ", sprintf("%.2f", median_log10_ratio)
      ),
      y = -0.62
    )
  letter_df <- extreme_summary %>%
    transmute(
      strategy,
      label = group_letter,
      y = 8.62,
      hjust = if_else(strategy == "Bottom 10%", 1.7, -0.7)
    )
  test_label <- paste0(
    format_p(wilcox_result$p.value, "Wilcoxon P"),
    "\nCliff's delta = ", sprintf("%.2f", delta_value),
    " (", magnitude, ")"
  )

  p_extreme <- ggplot(extreme_df, aes(strategy, genome_size_Mb)) +
    geom_violin(
      aes(fill = median_log10_ratio),
      width = 0.86,
      trim = TRUE,
      alpha = 0.48,
      linewidth = 0.58,
      color = "grey25"
    ) +
    geom_jitter(
      aes(color = log10_ratio),
      width = 0.13,
      shape = 16,
      size = 0.70,
      alpha = 0.28
    ) +
    geom_boxplot(
      aes(fill = median_log10_ratio),
      width = 0.16,
      outlier.shape = NA,
      color = "grey20",
      linewidth = 0.52,
      alpha = 0.82
    ) +
    geom_text(
      data = count_df,
      aes(strategy, y, label = label),
      inherit.aes = FALSE,
      size = 2.55,
      color = "grey30",
      lineheight = 0.9
    ) +
    geom_text(
      data = letter_df,
      aes(strategy, y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = 3.6,
      color = "black",
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = 9.72,
      label = test_label,
      size = 2.55,
      fontface = "bold",
      lineheight = 1.0
    ) +
    scale_fill_gradient2(
      low = "#287A62",
      mid = "#F3F1EC",
      high = "#D9483B",
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    scale_color_gradient2(
      low = "#287A62",
      mid = "#F3F1EC",
      high = "#D9483B",
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish,
      guide = "none"
    ) +
    guides(fill = "none") +
    scale_x_discrete(labels = c(
      "Bottom 10%" = "Bottom 10%\n(depolymerisation-\noptimised)",
      "Top 10%" = "Top 10%\n(uptake-optimised)"
    )) +
    scale_y_continuous(
      breaks = c(0, 2.5, 5, 7.5, 10),
      limits = c(-1.05, 10.2)
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 9.5) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 8),
      axis.text.x = element_text(lineheight = 0.9),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.4),
      plot.margin = margin(4, 4, 4, 5)
    )

  list(
    p_gam = p_gam,
    p_extreme = p_extreme,
    gam_stats = gam_stats,
    extreme_stats = extreme_summary
  )
}

results <- lapply(dataset_levels, analyse_dataset)

write_csv(bind_rows(lapply(results, `[[`, "gam_stats")), out_stats)
write_csv(bind_rows(lapply(results, `[[`, "extreme_stats")), out_extremes)

shared_legend <- cowplot::get_legend(
  results[[1]]$p_gam + theme(legend.position = "right")
)

pair_plots <- lapply(seq_along(results), function(i) {
  pair <- (
    (results[[i]]$p_gam + theme(legend.position = "none")) |
      results[[i]]$p_extreme
    ) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(title = dataset_levels[[i]])

  cowplot::ggdraw(pair) +
    cowplot::draw_label(
      letters[[i]],
      x = 0.004,
      y = 0.996,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 14
    )
})

figure_body <- wrap_plots(
  pair_plots,
  ncol = 2
)

combined <- (
  figure_body | wrap_elements(full = shared_legend)
) +
  plot_layout(widths = c(1, 0.055))

ggsave(out_png, combined, width = 16.5, height = 18.0, dpi = 320, bg = "white")
ggsave(out_pdf, combined, width = 16.5, height = 18.0, bg = "white")

message("Saved: ", out_png)
message("Saved: ", out_pdf)
message("Saved: ", out_stats)
message("Saved: ", out_extremes)
