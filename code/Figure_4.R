#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(mgcv)
  library(effectsize)
})

set.seed(20260924)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
in_rds <- file.path(base_dir, "data/processed/glo_input_binary.rds")
fig_dir <- file.path(base_dir, "results/figures")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(fig_dir, "Genome_size_GAM_only.png")
out_pdf <- file.path(fig_dir, "Genome_size_GAM_only.pdf")
out_models <- file.path(tab_dir, "Genome_size_GLM_vs_GAM_comparison.csv")
out_extremes <- file.path(tab_dir, "Genome_size_top_bottom10_summary.csv")

dat <- readRDS(in_rds) %>%
  transmute(
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
  mutate(log10_ratio = log10(ratio))

if (nrow(dat) == 0) stop("No valid GLO MAGs available for genome-size analysis.")

# Display every MAG. Low opacity preserves density information in crowded areas.
point_display_n <- nrow(dat)
point_dat <- dat

glm_fit <- glm(
  genome_size_Mb ~ log10_ratio,
  data = dat,
  family = gaussian(link = "identity")
)
gam_fit <- gam(
  genome_size_Mb ~ s(log10_ratio, k = 10),
  data = dat,
  method = "REML"
)

glm_summary <- summary(glm_fit)
gam_summary <- summary(gam_fit)
model_test <- anova(glm_fit, gam_fit, test = "F")

glm_r2 <- 1 - glm_fit$deviance / glm_fit$null.deviance

model_stats <- tibble(
  model = c("Gaussian GLM", "Generalised additive model"),
  n = nrow(dat),
  AIC = c(AIC(glm_fit), AIC(gam_fit)),
  adjusted_R2 = c(NA_real_, gam_summary$r.sq),
  deviance_explained = c(glm_r2, gam_summary$dev.expl),
  effective_degrees_of_freedom = c(1, gam_summary$s.table[1, "edf"]),
  term_p_value = c(
    glm_summary$coefficients["log10_ratio", "Pr(>|t|)"],
    gam_summary$s.table[1, "p-value"]
  ),
  comparison_p_value = c(NA_real_, model_test$`Pr(>F)`[2])
)
write_csv(model_stats, out_models)

prediction_grid <- tibble(
  log10_ratio = seq(
    min(dat$log10_ratio),
    max(dat$log10_ratio),
    length.out = 500
  )
)
gam_pred <- predict(gam_fit, newdata = prediction_grid, se.fit = TRUE)
prediction_grid <- prediction_grid %>%
  mutate(
    fit = as.numeric(gam_pred$fit),
    lower = fit - 1.96 * as.numeric(gam_pred$se.fit),
    upper = fit + 1.96 * as.numeric(gam_pred$se.fit)
  )

gam_label <- paste0(
  "n = ", format(nrow(dat), big.mark = ","), "\n",
  "GAM edf = ", sprintf("%.2f", gam_summary$s.table[1, "edf"]),
  "\nDeviance explained = ", sprintf("%.2f%%", 100 * gam_summary$dev.expl),
  "\n", ifelse(
    gam_summary$s.table[1, "p-value"] < 0.001,
    "P < 0.001",
    paste0("P = ", signif(gam_summary$s.table[1, "p-value"], 3))
  )
)

p_a <- ggplot(point_dat, aes(log10_ratio, genome_size_Mb)) +
  geom_point(
    aes(color = log10_ratio),
    alpha = 0.11,
    size = 0.95,
    stroke = 0
  ) +
  geom_ribbon(
    data = prediction_grid,
    aes(x = log10_ratio, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "grey68",
    alpha = 0.30
  ) +
  geom_line(
    data = prediction_grid,
    aes(x = log10_ratio, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.78
  ) +
  annotate(
    "text",
    x = -1.33,
    y = 9.75,
    label = gam_label,
    hjust = 0,
    vjust = 1,
    size = 3.25,
    fontface = "bold",
    lineheight = 1.02
  ) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5), limits = c(-1.45, 1.58)) +
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(0, 10.2)) +
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
    x = expression("Uptake-to-depolymerisation ratio ("*log[10]*")"),
    y = "Genome size (Mb)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 9.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.45),
    legend.position = "right",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(6, 8, 6, 6)
  )

q10 <- quantile(dat$ratio, 0.10, na.rm = TRUE)
q90 <- quantile(dat$ratio, 0.90, na.rm = TRUE)

extreme_df <- dat %>%
  mutate(
    strategy = case_when(
      ratio <= q10 ~ "Bottom 10%",
      ratio >= q90 ~ "Top 10%",
      TRUE ~ NA_character_
    ),
    strategy = factor(strategy, levels = c("Bottom 10%", "Top 10%"))
  ) %>%
  filter(!is.na(strategy))

extreme_point_df <- extreme_df

wilcox_result <- wilcox.test(genome_size_Mb ~ strategy, data = extreme_df, exact = FALSE)
delta_result <- effectsize::cliffs_delta(
  genome_size_Mb ~ strategy,
  data = extreme_df
)
delta_value <- delta_result$r_rank_biserial
delta_threshold <- 0.147
delta_magnitude <- case_when(
  abs(delta_value) < delta_threshold ~ "negligible",
  abs(delta_value) < 0.330 ~ "small",
  abs(delta_value) < 0.474 ~ "medium",
  TRUE ~ "large"
)
practical_difference <-
  wilcox_result$p.value < 0.05 && abs(delta_value) >= delta_threshold
strategy_letters <- if (practical_difference) c("a", "b") else c("a", "a")

extreme_summary <- extreme_df %>%
  group_by(strategy) %>%
  summarise(
    n = n(),
    median_genome_size_Mb = median(genome_size_Mb),
    q25_genome_size_Mb = quantile(genome_size_Mb, 0.25),
    q75_genome_size_Mb = quantile(genome_size_Mb, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    q10_ratio = q10,
    q90_ratio = q90,
    wilcoxon_p_value = wilcox_result$p.value,
    cliff_delta = delta_value,
    delta_ci_low = delta_result$CI_low,
    delta_ci_high = delta_result$CI_high,
    delta_magnitude = delta_magnitude,
    delta_threshold = delta_threshold,
    practical_difference = practical_difference,
    group_letter = strategy_letters
  )
write_csv(extreme_summary, out_extremes)

strategy_cols <- c("Bottom 10%" = "#8EA9A4", "Top 10%" = "#E48670")
count_df <- extreme_summary %>%
  mutate(label = paste0("n = ", format(n, big.mark = ",")), y = -0.47)
letter_df <- extreme_summary %>%
  transmute(strategy, label = group_letter, y = 8.85)

p_b <- ggplot(extreme_df, aes(strategy, genome_size_Mb, color = strategy)) +
  geom_violin(
    aes(fill = strategy),
    width = 0.82,
    trim = TRUE,
    alpha = 0.48,
    linewidth = 0.65,
    color = "grey25"
  ) +
  geom_jitter(
    data = extreme_point_df,
    width = 0.14,
    size = 0.62,
    alpha = 0.060,
    stroke = 0
  ) +
  geom_boxplot(
    aes(fill = strategy),
    width = 0.15,
    outlier.shape = NA,
    color = "grey20",
    linewidth = 0.58,
    alpha = 0.82
  ) +
  geom_text(
    data = count_df,
    aes(strategy, y, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    color = "grey35"
  ) +
  geom_text(
    data = letter_df,
    aes(strategy, y, label = label),
    inherit.aes = FALSE,
    size = 4.1,
    color = "black",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 9.65,
    label = paste0(
      ifelse(
        wilcox_result$p.value < 0.001,
        "Wilcoxon P < 0.001",
        paste0("Wilcoxon P = ", signif(wilcox_result$p.value, 3))
      ),
      "\nCliff's delta = ", sprintf("%.2f", delta_value),
      " (", delta_magnitude, ")"
    ),
    size = 3.15,
    fontface = "bold",
    lineheight = 1.05
  ) +
  scale_color_manual(values = strategy_cols) +
  scale_fill_manual(values = strategy_cols) +
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(-0.65, 10.2)) +
  labs(x = NULL, y = "Genome size (Mb)") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 9.5),
    axis.text.x = element_text(face = "plain"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.45),
    plot.margin = margin(6, 6, 6, 8)
  )

ggsave(out_png, p_a, width = 6.3, height = 5.1, dpi = 500, bg = "white")
ggsave(out_pdf, p_a, width = 6.3, height = 5.1, bg = "white")

message("Saved: ", out_png)
message("Saved: ", out_pdf)
message("Saved: ", out_models)
message("Saved: ", out_extremes)
