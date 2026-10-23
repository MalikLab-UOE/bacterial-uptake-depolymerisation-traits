#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(emmeans)
  library(multcompView)
  library(sandwich)
  library(patchwork)
})

set.seed(20260828)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
fig_dir <- file.path(base_dir, "results/figures")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

glo_rds <- file.path(base_dir, "data/processed/glo_input_binary.rds")
all_db_csv <- file.path(base_dir, "data/processed/all_db_merged.csv")
high_env_csv <- file.path(base_dir, "data/raw/high_level_env__value.csv")

out_png <- file.path(fig_dir, "Fig_1_corrected_log10_weighted_habitats.png")
out_pdf <- file.path(fig_dir, "Fig_1_corrected_log10_weighted_habitats.pdf")
out_b_png <- file.path(fig_dir, "Fig_1b_corrected_log10_weighted_habitats.png")
out_b_pdf <- file.path(fig_dir, "Fig_1b_corrected_log10_weighted_habitats.pdf")
out_summary <- file.path(tab_dir, "Fig_1b_corrected_habitat_summary.csv")
out_pairs <- file.path(tab_dir, "Fig_1b_corrected_Tukey_pairwise_comparisons.csv")
out_model <- file.path(tab_dir, "Fig_1b_corrected_model_summary.csv")
group_levels <- c(
  "animal_associated", "soil", "built", "sediment", "aquatic"
)

group_cols <- c(
  "animal_associated" = "#78ADD0",
  "soil" = "#BDBF63",
  "built" = "#A88BC7",
  "sediment" = "#E79862",
  "aquatic" = "#D87378"
)

group_fill_cols <- c(
  "animal_associated" = "#A9CFE8",
  "soil" = "#E5E6A3",
  "built" = "#D9C7EA",
  "sediment" = "#F7C49E",
  "aquatic" = "#F1A9A8"
)

group_labels <- c(
  "animal_associated" = "Animal-\nassociated",
  "soil" = "Soil",
  "built" = "Built",
  "sediment" = "Sediment",
  "aquatic" = "Aquatic"
)

delta_threshold <- 0.147

weighted_quantile <- function(x, weight, probs = c(0.25, 0.5, 0.75)) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  x <- x[keep]
  weight <- weight[keep]
  order_index <- order(x)
  x <- x[order_index]
  weight <- weight[order_index]
  cumulative_weight <- cumsum(weight) / sum(weight)
  vapply(
    probs,
    function(probability) x[which(cumulative_weight >= probability)[1]],
    numeric(1)
  )
}

make_box_summary <- function(x, weight) {
  quartiles <- weighted_quantile(x, weight)
  iqr <- quartiles[3] - quartiles[1]
  lower_fence <- quartiles[1] - 1.5 * iqr
  upper_fence <- quartiles[3] + 1.5 * iqr
  tibble(
    lower = quartiles[1],
    middle = quartiles[2],
    upper = quartiles[3],
    ymin = min(x[x >= lower_fence], na.rm = TRUE),
    ymax = max(x[x <= upper_fence], na.rm = TRUE),
    mean = weighted.mean(x, weight)
  )
}

delta_magnitude <- function(delta) {
  abs_delta <- abs(delta)
  case_when(
    abs_delta < delta_threshold ~ "negligible",
    abs_delta < 0.330 ~ "small",
    abs_delta < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

# Weighted Cliff's delta is the weighted probability that a value from group 1
# exceeds a value from group 2 minus the reverse probability. With unit weights
# this reduces to the standard Cliff's delta.
weighted_cliffs_delta <- function(x, weight_x, y, weight_y) {
  x_df <- tibble(value = x, weight = weight_x) %>%
    group_by(value) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    arrange(value)
  y_df <- tibble(value = y, weight = weight_y) %>%
    group_by(value) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    arrange(value)

  y_cumulative <- cumsum(y_df$weight)
  y_total <- sum(y_df$weight)
  index_le <- findInterval(x_df$value, y_df$value)
  equal_index <- match(x_df$value, y_df$value, nomatch = 0L)
  equal_weight <- ifelse(equal_index > 0, y_df$weight[pmax(equal_index, 1L)], 0)
  less_weight <- ifelse(index_le > 0, y_cumulative[pmax(index_le, 1L)], 0) -
    equal_weight
  greater_weight <- y_total - less_weight - equal_weight

  sum(x_df$weight * (less_weight - greater_weight)) /
    (sum(x_df$weight) * y_total)
}

calculate_pairwise_delta <- function(data, pair_table, weight_column) {
  bind_rows(lapply(seq_len(nrow(pair_table)), function(i) {
    group_1 <- pair_table$group1[i]
    group_2 <- pair_table$group2[i]
    data_1 <- data %>% filter(as.character(habitat) == group_1)
    data_2 <- data %>% filter(as.character(habitat) == group_2)
    delta <- weighted_cliffs_delta(
      data_1$log10_ratio,
      data_1[[weight_column]],
      data_2$log10_ratio,
      data_2[[weight_column]]
    )
    tibble(
      group1 = group_1,
      group2 = group_2,
      cliffs_delta = delta,
      delta_magnitude = delta_magnitude(delta)
    )
  }))
}

make_practical_letters <- function(pair_table, level_order) {
  letter_values <- ifelse(pair_table$practical_difference, 0, 1)
  names(letter_values) <- paste(pair_table$group1, pair_table$group2, sep = "-")
  raw_letters <- multcompView::multcompLetters(
    letter_values,
    threshold = 0.05
  )$Letters

  letter_df <- tibble(
    habitat = names(raw_letters),
    raw_letter = unname(raw_letters),
    display_rank = match(names(raw_letters), level_order)
  )
  tokens <- unique(unlist(strsplit(letter_df$raw_letter, "")))
  token_key <- bind_rows(lapply(tokens, function(token) {
    ranks <- letter_df$display_rank[
      grepl(token, letter_df$raw_letter, fixed = TRUE)
    ]
    tibble(token = token, first_rank = min(ranks), mean_rank = mean(ranks))
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
      habitat = factor(habitat, levels = level_order),
      letter = vapply(raw_letter, relabel, character(1))
    )
}

# Panel a: the same continuous log10 trait used in panel b.
glo <- readRDS(glo_rds) %>%
  transmute(
    GH = suppressWarnings(as.numeric(GH_richness)),
    uptake = suppressWarnings(as.numeric(ABC_richness)),
    ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  filter(
    is.finite(GH),
    is.finite(uptake),
    is.finite(ratio),
    ratio > 0
  ) %>%
  mutate(log10_ratio = log10(ratio))

p_a <- ggplot(glo, aes(GH, uptake, color = log10_ratio)) +
  geom_point(size = 1.55, alpha = 0.52, stroke = 0) +
  scale_color_gradient2(
    low = "#276C58",
    mid = "#F7F7F4",
    high = "#C92C28",
    midpoint = 0,
    limits = c(-2, 2),
    oob = squish,
    breaks = -2:2,
    name = expression("Uptake-to-depolymerisation ratio ("*log[10]*")")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barheight = unit(2.2, "cm"),
      barwidth = unit(0.25, "cm"),
      ticks.colour = "white",
      frame.colour = NA
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    x = "Complex carbohydrate depolymerisation enzymes",
    y = "Carbohydrate transporters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.57, 0.90),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.title = element_text(size = 8.5, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.title = element_text(face = "bold", size = 11, color = "black"),
    axis.text = element_text(size = 9.5, color = "black"),
    plot.margin = margin(8, 14, 8, 8)
  )

# Panel b data: assign each genome to its most frequently detected habitat.
all_db <- read_csv(all_db_csv, show_col_types = FALSE) %>%
  mutate(
    genome_id = as.character(genome_id),
    database = tolower(as.character(database)),
    ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  filter(database == "glo", is.finite(ratio), ratio > 0) %>%
  distinct(genome_id, .keep_all = TRUE) %>%
  dplyr::select(genome_id, ratio)

env_long <- read_csv(high_env_csv, show_col_types = FALSE) %>%
  rename(genome_id = ID) %>%
  mutate(genome_id = as.character(genome_id)) %>%
  pivot_longer(-genome_id, names_to = "habitat", values_to = "detect_count") %>%
  mutate(
    habitat = sub("^L1_", "", habitat),
    detect_count = suppressWarnings(as.numeric(detect_count))
  ) %>%
  filter(
    habitat %in% group_levels,
    is.finite(detect_count),
    detect_count > 0
  )

plot_df <- env_long %>%
  group_by(genome_id) %>%
  slice_max(detect_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(all_db, by = "genome_id") %>%
  mutate(
    habitat = factor(habitat, levels = group_levels),
    log10_ratio = log10(ratio)
  ) %>%
  filter(is.finite(log10_ratio), is.finite(detect_count), detect_count > 0)

# Detection count is used as the analytical weight. Tukey letters therefore
# describe comparisons among weighted estimated means on the log10 scale.
fit <- lm(log10_ratio ~ habitat, data = plot_df, weights = detect_count)
robust_vcov <- sandwich::vcovHC(fit, type = "HC3")
emm <- emmeans(fit, ~ habitat, vcov. = robust_vcov)

pair_df <- summary(pairs(emm, adjust = "tukey"), infer = TRUE) %>%
  as.data.frame() %>%
  mutate(fold_difference = 10^estimate) %>%
  separate_wider_delim(
    contrast,
    delim = " - ",
    names = c("group1", "group2"),
    cols_remove = FALSE
  )

weighted_effects <- calculate_pairwise_delta(
  plot_df,
  pair_df,
  weight_column = "detect_count"
)
pair_df <- pair_df %>%
  left_join(weighted_effects, by = c("group1", "group2")) %>%
  mutate(
    delta_threshold = delta_threshold,
    practical_difference = p.value < 0.05 &
      abs(cliffs_delta) >= delta_threshold
  )
weighted_letters <- make_practical_letters(pair_df, group_levels)

letters_df <- summary(emm, infer = c(TRUE, FALSE), adjust = "none") %>%
  as.data.frame() %>%
  transmute(
    habitat = factor(habitat, levels = group_levels),
    weighted_mean_log10 = emmean,
    SE = SE,
    lower_CL = lower.CL,
    upper_CL = upper.CL
  ) %>%
  left_join(weighted_letters, by = "habitat")

summary_df <- plot_df %>%
  group_by(habitat) %>%
  summarise(
    unique_genomes = n(),
    weighted_n = sum(detect_count),
    median_log10 = median(log10_ratio),
    unweighted_mean_log10 = mean(log10_ratio),
    SD_log10 = sd(log10_ratio),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "habitat") %>%
  mutate(
    geometric_mean_ratio = 10^weighted_mean_log10,
    habitat_label = unname(group_labels[as.character(habitat)])
  )

write_csv(summary_df, out_summary)
write_csv(pair_df, out_pairs)
write_csv(
  tibble(
    unique_genomes = nrow(plot_df),
    model_R2 = summary(fit)$r.squared,
    adjusted_R2 = summary(fit)$adj.r.squared,
    model_df = df.residual(fit),
    covariance_estimator = "HC3 robust"
  ),
  out_model
)

count_df <- summary_df %>%
  mutate(
    count_label = case_when(
      weighted_n >= 1e6 ~ paste0("Weighted n = ", round(weighted_n / 1e6, 1), "M"),
      weighted_n >= 1e3 ~ paste0("Weighted n = ", round(weighted_n / 1e3, 1), "K"),
      TRUE ~ paste0("Weighted n = ", round(weighted_n))
    )
  )

weighted_box_df <- plot_df %>%
  group_by(habitat) %>%
  group_modify(~ make_box_summary(.x$log10_ratio, .x$detect_count)) %>%
  ungroup()

p_b <- ggplot(
  plot_df,
  aes(x = habitat, y = log10_ratio, weight = detect_count)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey60"
  ) +
  geom_violin(
    aes(fill = habitat),
    color = NA,
    alpha = 0.58,
    linewidth = 0,
    width = 0.82,
    trim = TRUE,
    scale = "width",
    quantiles = NULL
  ) +
  geom_jitter(
    aes(color = habitat),
    width = 0.14,
    height = 0,
    size = 0.38,
    alpha = 0.065,
    stroke = 0
  ) +
  geom_boxplot(
    data = weighted_box_df,
    aes(
      x = habitat,
      ymin = ymin,
      lower = lower,
      middle = middle,
      upper = upper,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    stat = "identity",
    width = 0.18,
    linewidth = 0.62,
    color = "grey15",
    fill = "grey30",
    alpha = 0.70
  ) +
  geom_point(
    data = weighted_box_df,
    aes(x = habitat, y = mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.4,
    stroke = 0.55,
    fill = "white",
    color = "grey15"
  ) +
  geom_text(
    data = letters_df,
    aes(x = habitat, y = 1.68, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    data = count_df,
    aes(x = habitat, y = -1.69, label = count_label, color = habitat),
    inherit.aes = FALSE,
    size = 2.45,
    fontface = "bold",
    lineheight = 0.9
  ) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_fill_cols) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(
    breaks = -2:2,
    limits = c(-2.02, 2.02),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = expression("Uptake-to-depolymerisation ratio ("*log[10]*")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(
      size = 9.2,
      face = "bold",
      color = "black",
      lineheight = 0.9
    ),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 9)
  )

combined <- (p_a + p_b + plot_layout(widths = c(1.08, 1))) +
  plot_annotation(tag_levels = "a")

ggsave(out_png, combined, width = 12.1, height = 5.0, dpi = 400, bg = "white")
ggsave(out_pdf, combined, width = 12.1, height = 5.0, bg = "white")
ggsave(out_b_png, p_b, width = 6.1, height = 5.0, dpi = 400, bg = "white")
ggsave(out_b_pdf, p_b, width = 6.1, height = 5.0, bg = "white")
message("Saved: ", out_png)
message("Saved: ", out_pdf)
message(
  "Weighted P + Cliff's delta letters: ",
  paste(summary_df$letter, collapse = ", ")
)
message("Model R2: ", signif(summary(fit)$r.squared, 4))
