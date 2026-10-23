#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(effectsize)
  library(multcompView)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
panel_a_input <- file.path(base_dir, "data/processed/loma_input.csv")
panel_c_input <- file.path(base_dir, "data/raw/LOMA_merged_mags_50_10.csv")
fig_dir <- file.path(base_dir, "results/figures")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

delta_threshold <- 0.147
alpha_threshold <- 0.05
time_levels <- c("T1", "T2", "T3", "T4")
time_labels <- c(T1 = "2 months", T2 = "7 months", T3 = "14 months", T4 = "17 months")
time_cols <- c(T1 = "#D52B34", T2 = "#F39A72", T3 = "#8CC5DE", T4 = "#176AA5")
tax_cols <- c(Bacteria = "#F4A3A3", Fungi = "#59DAD8")

delta_magnitude <- function(delta) {
  abs_delta <- abs(delta)
  dplyr::case_when(
    abs_delta < 0.147 ~ "negligible",
    abs_delta < 0.330 ~ "small",
    abs_delta < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

pairwise_wilcox_delta <- function(data, response, group, level_order) {
  pairs <- combn(level_order, 2, simplify = FALSE)
  out <- lapply(pairs, function(pair) {
    group_1 <- pair[[1]]
    group_2 <- pair[[2]]
    pair_data <- data %>%
      filter(.data[[group]] %in% c(group_1, group_2)) %>%
      transmute(
        value = suppressWarnings(as.numeric(.data[[response]])),
        comparison_group = factor(.data[[group]], levels = c(group_1, group_2))
      ) %>%
      filter(is.finite(value), !is.na(comparison_group))

    values_1 <- pair_data$value[pair_data$comparison_group == group_1]
    values_2 <- pair_data$value[pair_data$comparison_group == group_2]
    wt <- stats::wilcox.test(values_1, values_2, exact = FALSE)
    delta_result <- effectsize::cliffs_delta(
      value ~ comparison_group,
      data = pair_data,
      ci = 0.95
    )

    tibble::tibble(
      group1 = group_1,
      group2 = group_2,
      n1 = length(values_1),
      n2 = length(values_2),
      median1 = stats::median(values_1),
      median2 = stats::median(values_2),
      wilcoxon_W = unname(wt$statistic),
      p_value = wt$p.value,
      cliff_delta = delta_result$r_rank_biserial,
      delta_ci_low = delta_result$CI_low,
      delta_ci_high = delta_result$CI_high
    )
  }) %>% bind_rows()

  out %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      delta_magnitude = delta_magnitude(cliff_delta),
      meaningful_difference =
        p_adj_BH < alpha_threshold & abs(cliff_delta) >= delta_threshold
    )
}

make_effect_letters <- function(pair_table, level_order) {
  letter_values <- ifelse(pair_table$meaningful_difference, 0, 1)
  names(letter_values) <- paste(pair_table$group1, pair_table$group2, sep = "-")
  raw_letters <- multcompView::multcompLetters(
    letter_values,
    threshold = 0.05
  )$Letters

  letter_df <- tibble::tibble(
    group = names(raw_letters),
    raw_letter = unname(raw_letters),
    display_rank = match(names(raw_letters), level_order)
  )

  tokens <- unique(unlist(strsplit(letter_df$raw_letter, "")))
  token_key <- bind_rows(lapply(tokens, function(token) {
    ranks <- letter_df$display_rank[
      grepl(token, letter_df$raw_letter, fixed = TRUE)
    ]
    tibble::tibble(
      token = token,
      first_rank = min(ranks),
      mean_rank = mean(ranks)
    )
  })) %>%
    arrange(first_rank, mean_rank) %>%
    mutate(display_token = letters[row_number()])

  relabel <- function(raw_label) {
    mapped <- token_key$display_token[
      match(strsplit(raw_label, "")[[1]], token_key$token)
    ]
    paste(sort(mapped), collapse = "")
  }

  letter_df %>%
    transmute(
      group = factor(group, levels = level_order),
      letter = vapply(raw_letter, relabel, character(1))
    ) %>%
    arrange(group)
}

# -----------------------------------------------------------------------------
# Panel a: temporal changes in GH repertoire, transporter repertoire and ratio
# -----------------------------------------------------------------------------
panel_a <- readr::read_csv(panel_a_input, show_col_types = FALSE) %>%
  transmute(
    time = factor(as.character(time), levels = time_levels),
    depolymerisation = suppressWarnings(as.numeric(GH_total_abundance)),
    uptake = suppressWarnings(as.numeric(ABC_total_abundance)),
    ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  mutate(log10_ratio = if_else(is.finite(ratio) & ratio > 0, log10(ratio), NA_real_)) %>%
  filter(!is.na(time))

metric_spec <- tibble::tribble(
  ~metric, ~column, ~axis_label,
  "Depolymerisation", "depolymerisation", "Complex carbohydrate\ndepolymerisation enzymes",
  "Uptake", "uptake", "Carbohydrate\ntransporters",
  "Ratio", "log10_ratio", "log10 uptake-to-\ndepolymerisation ratio"
)

panel_a_pairs <- lapply(seq_len(nrow(metric_spec)), function(i) {
  spec <- metric_spec[i, ]
  pairwise_wilcox_delta(panel_a, spec$column, "time", time_levels) %>%
    mutate(metric = spec$metric, .before = 1)
}) %>% bind_rows()

panel_a_letters <- lapply(seq_len(nrow(metric_spec)), function(i) {
  spec <- metric_spec[i, ]
  metric_pairs <- panel_a_pairs %>% filter(metric == spec$metric)
  make_effect_letters(metric_pairs, time_levels) %>%
    mutate(metric = spec$metric, .before = 1)
}) %>% bind_rows()

panel_a_summary <- lapply(seq_len(nrow(metric_spec)), function(i) {
  spec <- metric_spec[i, ]
  panel_a %>%
    transmute(time, value = .data[[spec$column]]) %>%
    filter(is.finite(value)) %>%
    group_by(time) %>%
    summarise(
      n = n(),
      median = median(value),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      .groups = "drop"
    ) %>%
    mutate(metric = spec$metric, .before = 1) %>%
    left_join(
      panel_a_letters %>% filter(metric == spec$metric) %>% select(time = group, letter),
      by = "time"
    )
}) %>% bind_rows()

readr::write_csv(
  panel_a_pairs,
  file.path(tab_dir, "Figure_litter_panel_a_pairwise_cliffs_delta.csv")
)
readr::write_csv(
  panel_a_summary,
  file.path(tab_dir, "Figure_litter_panel_a_summary_effect_letters.csv")
)

panel_a_plots <- lapply(seq_len(nrow(metric_spec)), function(i) {
  spec <- metric_spec[i, ]
  plot_data <- panel_a %>%
    transmute(time, value = .data[[spec$column]]) %>%
    filter(is.finite(value))
  annotation_data <- panel_a_summary %>%
    filter(metric == spec$metric) %>%
    mutate(
      y = max(plot_data$value) + 0.045 * diff(range(plot_data$value))
    )

  ggplot(plot_data, aes(time, value)) +
    geom_violin(
      fill = "grey88", color = NA, width = 0.82,
      scale = "width", trim = FALSE
    ) +
    geom_jitter(
      aes(color = time), width = 0.16, height = 0,
      size = 2.0, alpha = 0.80
    ) +
    stat_summary(
      fun = median, geom = "crossbar", width = 0.50,
      color = "black", linewidth = 0.55
    ) +
    geom_text(
      data = annotation_data,
      aes(time, y, label = letter),
      inherit.aes = FALSE,
      fontface = "bold", size = 4.8, vjust = 0
    ) +
    scale_color_manual(values = time_cols, guide = "none") +
    scale_x_discrete(labels = time_labels) +
    labs(x = NULL, y = spec$axis_label[[1]]) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      plot.margin = margin(12, 8, 4, 6)
    )
})

panel_a_figure <- wrap_plots(panel_a_plots, ncol = 1)

ggsave(
  file.path(fig_dir, "Figure_litter_panel_a_cliffs_delta.pdf"),
  panel_a_figure, width = 5.3, height = 9.2
)
ggsave(
  file.path(fig_dir, "Figure_litter_panel_a_cliffs_delta.png"),
  panel_a_figure, width = 5.3, height = 9.2, dpi = 400
)

# -----------------------------------------------------------------------------
# Panel c, upper: GH gene number in bacterial versus fungal MAGs
# -----------------------------------------------------------------------------
panel_c <- readr::read_csv(panel_c_input, show_col_types = FALSE) %>%
  transmute(
    tax = recode(as.character(tax), "Bacteria" = "Bacteria", "euk" = "Fungi"),
    GH = suppressWarnings(as.numeric(GH))
  ) %>%
  filter(tax %in% c("Bacteria", "Fungi"), is.finite(GH)) %>%
  mutate(tax = factor(tax, levels = c("Bacteria", "Fungi")))

panel_c_result <- pairwise_wilcox_delta(
  panel_c,
  response = "GH",
  group = "tax",
  level_order = c("Bacteria", "Fungi")
) %>%
  mutate(
    letter_bacteria = "a",
    letter_fungi = if_else(meaningful_difference, "b", "a")
  )

readr::write_csv(
  panel_c_result,
  file.path(tab_dir, "Figure_litter_panel_c_top_cliffs_delta.csv")
)

panel_c_letters <- tibble::tibble(
  tax = factor(c("Bacteria", "Fungi"), levels = c("Bacteria", "Fungi")),
  letter = c(
    panel_c_result$letter_bacteria[[1]],
    panel_c_result$letter_fungi[[1]]
  )
) %>%
  left_join(
    panel_c %>% group_by(tax) %>% summarise(y_max = max(GH), .groups = "drop"),
    by = "tax"
  ) %>%
  mutate(y = y_max + 0.055 * diff(range(panel_c$GH)))

panel_c_figure <- ggplot(panel_c, aes(tax, GH)) +
  geom_violin(
    aes(fill = tax), color = NA, alpha = 0.30,
    width = 0.90, trim = FALSE, scale = "width"
  ) +
  geom_jitter(
    aes(color = tax), width = 0.16,
    size = 2.0, alpha = 0.82
  ) +
  stat_summary(
    fun = median, geom = "crossbar", width = 0.52,
    color = "black", linewidth = 0.60
  ) +
  geom_text(
    data = panel_c_letters,
    aes(tax, y, label = letter),
    inherit.aes = FALSE,
    fontface = "bold", size = 5.0, vjust = 0
  ) +
  scale_fill_manual(values = tax_cols, guide = "none") +
  scale_color_manual(values = tax_cols, guide = "none") +
  labs(x = NULL, y = "Glycoside hydrolase\ngene number") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.75),
    axis.ticks = element_line(linewidth = 0.75),
    plot.margin = margin(14, 8, 6, 8)
  )

ggsave(
  file.path(fig_dir, "Figure_litter_panel_c_top_cliffs_delta.pdf"),
  panel_c_figure, width = 4.5, height = 4.0
)
ggsave(
  file.path(fig_dir, "Figure_litter_panel_c_top_cliffs_delta.png"),
  panel_c_figure, width = 4.5, height = 4.0, dpi = 400
)

combined_preview <- panel_a_figure | panel_c_figure
ggsave(
  file.path(fig_dir, "Figure_litter_panels_a_c_top_cliffs_delta_preview.pdf"),
  combined_preview, width = 10.2, height = 9.2
)
ggsave(
  file.path(fig_dir, "Figure_litter_panels_a_c_top_cliffs_delta_preview.png"),
  combined_preview, width = 10.2, height = 9.2, dpi = 350
)

message("Panel a letters:")
print(panel_a_summary %>% select(metric, time, n, median, letter))
message("Panel c result:")
print(panel_c_result)

# -----------------------------------------------------------------------------
# Panel b, upper: fungal-to-bacterial abundance ratio across sampling times
# Retain the existing ANOVA/Tukey framework and add matched parametric effect
# sizes: omega-squared overall and Hedges' g for pairwise comparisons.
# -----------------------------------------------------------------------------
fb_data <- readr::read_delim(
  file.path(base_dir, "data/raw/LOMA_FBratio.txt"),
  delim = "\t",
  show_col_types = FALSE,
  name_repair = "unique"
) %>%
  transmute(
    time = factor(as.character(samp), levels = time_levels),
    bacteria = suppressWarnings(as.numeric(bac)),
    fungi = suppressWarnings(as.numeric(fun))
  ) %>%
  mutate(fungal_bacterial_ratio = fungi / bacteria) %>%
  filter(!is.na(time), is.finite(fungal_bacterial_ratio))

fb_model <- stats::aov(fungal_bacterial_ratio ~ time, data = fb_data)
fb_anova <- summary(fb_model)[[1]]
fb_omega <- effectsize::omega_squared(fb_model, ci = 0.95)

fb_overall <- tibble::tibble(
  n = nrow(fb_data),
  df_between = fb_anova$Df[[1]],
  df_within = fb_anova$Df[[2]],
  F_value = fb_anova$`F value`[[1]],
  p_value = fb_anova$`Pr(>F)`[[1]],
  omega_squared = fb_omega$Omega2[[1]],
  omega_ci_low = fb_omega$CI_low[[1]],
  omega_ci_high = fb_omega$CI_high[[1]]
)

fb_tukey <- as.data.frame(stats::TukeyHSD(fb_model)$time) %>%
  tibble::rownames_to_column("contrast")

fb_pairwise <- lapply(combn(time_levels, 2, simplify = FALSE), function(pair) {
  pair_data <- fb_data %>%
    filter(time %in% pair) %>%
    mutate(time = droplevels(time))
  hedges <- effectsize::hedges_g(
    fungal_bacterial_ratio ~ time,
    data = pair_data,
    ci = 0.95
  )
  tukey_name <- paste(pair[[2]], pair[[1]], sep = "-")
  tukey_row <- fb_tukey %>% filter(contrast == tukey_name)

  tibble::tibble(
    group1 = pair[[1]],
    group2 = pair[[2]],
    n1 = sum(pair_data$time == pair[[1]]),
    n2 = sum(pair_data$time == pair[[2]]),
    median1 = median(pair_data$fungal_bacterial_ratio[pair_data$time == pair[[1]]]),
    median2 = median(pair_data$fungal_bacterial_ratio[pair_data$time == pair[[2]]]),
    mean_difference_tukey = tukey_row$diff[[1]],
    tukey_ci_low = tukey_row$lwr[[1]],
    tukey_ci_high = tukey_row$upr[[1]],
    p_adj_tukey = tukey_row$`p adj`[[1]],
    hedges_g = hedges$Hedges_g[[1]],
    hedges_g_ci_low = hedges$CI_low[[1]],
    hedges_g_ci_high = hedges$CI_high[[1]]
  )
}) %>% bind_rows()

readr::write_csv(
  fb_overall,
  file.path(tab_dir, "Figure_litter_panel_b_fb_ratio_anova_effect_size.csv")
)
readr::write_csv(
  fb_pairwise,
  file.path(tab_dir, "Figure_litter_panel_b_fb_ratio_tukey_hedges_g.csv")
)

message("Panel b fungal:bacterial ratio overall result:")
print(fb_overall)
message("Saved updated figures and effect-size tables.")
