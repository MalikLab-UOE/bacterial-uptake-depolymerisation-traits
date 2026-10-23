#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(multcompView)
  library(plotbiomes)
  library(patchwork)
  library(scales)
  library(readxl)
  library(stringr)
  library(effectsize)
})

set.seed(20260828)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
in_rds <- file.path(base_dir, "data/processed/img_input.rds")
fig_dir <- file.path(base_dir, "results/figures")
tab_dir <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(fig_dir, "Figure_2_ratio_coloured_Whittaker.png")
out_pdf <- file.path(fig_dir, "Figure_2_ratio_coloured_Whittaker.pdf")
out_site_png <- file.path(fig_dir, "Figure_S_bulk_rhizosphere_comparison.png")
out_site_pdf <- file.path(fig_dir, "Figure_S_bulk_rhizosphere_comparison.pdf")
out_summary <- file.path(tab_dir, "Figure_2_ratio_coloured_biome_summary.csv")
out_biome_pairs <- file.path(tab_dir, "Figure_2_biome_pairwise_cliffs_delta.csv")
out_site_summary <- file.path(tab_dir, "Figure_2_bulk_rhizosphere_site_summary.csv")
out_site_tests <- file.path(tab_dir, "Figure_2_bulk_rhizosphere_site_tests.csv")
out_site_exclusions <- file.path(
  tab_dir,
  "Figure_2_bulk_rhizosphere_conflicting_site_metadata_excluded.csv"
)
metadata_xlsx <- file.path(base_dir, "data/raw/Original_eco.xlsx")

min_n <- 100
palette_low <- "#276C58"
palette_mid <- "#F7F7F4"
palette_high <- "#C92C28"

label_wrap <- c(
  "Woodland/shrubland" = "Woodland/shrubland",
  "Temperate grassland/desert" = "Temperate grassland/desert",
  "Tropical seasonal forest/savanna" = "Tropical seasonal forest/savanna",
  "Tundra" = "Tundra",
  "Subtropical desert" = "Subtropical desert",
  "Temperate seasonal forest" = "Temperate seasonal forest",
  "Boreal forest" = "Boreal forest"
)

merged_table <- readRDS(in_rds)

plot_table <- merged_table %>%
  rename(Biome_label_raw = Biome_new) %>%
  filter(
    !is.na(Biome_label_raw),
    is.finite(ABC_GH_ratio),
    ABC_GH_ratio > 0
  ) %>%
  mutate(log10_ratio = log10(ABC_GH_ratio)) %>%
  add_count(Biome_label_raw, name = "n") %>%
  filter(n >= min_n)

summary_df <- plot_table %>%
  group_by(Biome_label_raw) %>%
  summarise(
    n = n(),
    mean_log10 = mean(log10_ratio),
    median_log10 = median(log10_ratio),
    q25_log10 = quantile(log10_ratio, 0.25),
    q75_log10 = quantile(log10_ratio, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(median_log10), desc(mean_log10))

# Order biomes by decreasing median ratio. Mean ratio breaks exact median ties
# so the visual order remains deterministic and interpretable.
ordered_biomes <- summary_df$Biome_label_raw

plot_table <- plot_table %>%
  left_join(
    summary_df %>% select(Biome_label_raw, median_log10),
    by = "Biome_label_raw"
  ) %>%
  mutate(
    Biome_label_raw = factor(Biome_label_raw, levels = rev(ordered_biomes)),
    Biome_label = paste0(
      label_wrap[as.character(Biome_label_raw)],
      " (n=", n, ")"
    ),
    Biome_label = factor(
      Biome_label,
      levels = paste0(
        label_wrap[rev(ordered_biomes)],
        " (n=",
        summary_df$n[match(rev(ordered_biomes), summary_df$Biome_label_raw)],
        ")"
      )
    )
  )

# Require both statistical significance and a non-negligible Cliff's delta
# for the compact-letter display. This prevents the very large MAG sample
# sizes from turning trivial effects into separate visual groups.
kw_res <- plot_table %>% kruskal_test(log10_ratio ~ Biome_label)
dunn_res <- plot_table %>%
  dunn_test(log10_ratio ~ Biome_label, p.adjust.method = "BH")

delta_magnitude <- function(delta) {
  abs_delta <- abs(delta)
  case_when(
    abs_delta < 0.147 ~ "negligible",
    abs_delta < 0.330 ~ "small",
    abs_delta < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

pairwise_effects <- lapply(seq_len(nrow(dunn_res)), function(i) {
  group_1 <- dunn_res$group1[i]
  group_2 <- dunn_res$group2[i]
  pair_data <- plot_table %>%
    filter(as.character(Biome_label) %in% c(group_1, group_2)) %>%
    transmute(
      log10_ratio,
      group = factor(as.character(Biome_label), levels = c(group_1, group_2))
    )
  delta_result <- effectsize::cliffs_delta(log10_ratio ~ group, data = pair_data)

  tibble::tibble(
    group1 = group_1,
    group2 = group_2,
    cliff_delta = delta_result$r_rank_biserial,
    delta_ci_low = delta_result$CI_low,
    delta_ci_high = delta_result$CI_high
  )
}) %>%
  bind_rows()

dunn_res <- dunn_res %>%
  left_join(pairwise_effects, by = c("group1", "group2")) %>%
  mutate(
    delta_magnitude = delta_magnitude(cliff_delta),
    meaningful_difference = p.adj < 0.05 & abs(cliff_delta) >= 0.147
  )

letter_values <- ifelse(dunn_res$meaningful_difference, 0, 1)
names(letter_values) <- paste(dunn_res$group1, dunn_res$group2, sep = "-")
raw_letters <- multcompLetters(letter_values, threshold = 0.05)$Letters

biome_label_key <- plot_table %>%
  distinct(Biome_label, Biome_label_raw) %>%
  mutate(
    Biome_label = as.character(Biome_label),
    Biome_label_raw = as.character(Biome_label_raw)
  )

pairwise_out <- dunn_res %>%
  left_join(
    biome_label_key %>% rename(group1 = Biome_label, biome1 = Biome_label_raw),
    by = "group1"
  ) %>%
  left_join(
    biome_label_key %>% rename(group2 = Biome_label, biome2 = Biome_label_raw),
    by = "group2"
  ) %>%
  select(
    biome1, biome2, n1, n2, statistic, p, p.adj, p.adj.signif,
    cliff_delta, delta_ci_low, delta_ci_high, delta_magnitude,
    meaningful_difference
  )

write.csv(pairwise_out, out_biome_pairs, row.names = FALSE)

cld_df <- data.frame(
  Biome_label = names(raw_letters),
  cld_raw = unname(raw_letters),
  row.names = NULL
) %>%
  left_join(
    plot_table %>%
      distinct(Biome_label, Biome_label_raw) %>%
      mutate(Biome_label_raw = as.character(Biome_label_raw)),
    by = "Biome_label"
  ) %>%
  left_join(summary_df, by = "Biome_label_raw")

cld_ranked <- cld_df %>%
  mutate(display_rank = match(Biome_label_raw, ordered_biomes))

letter_tokens <- unique(unlist(strsplit(cld_ranked$cld_raw, "")))
letter_key <- lapply(letter_tokens, function(token) {
  token_ranks <- cld_ranked$display_rank[
    grepl(token, cld_ranked$cld_raw, fixed = TRUE)
  ]
  tibble::tibble(
    token = token,
    first_rank = min(token_ranks),
    mean_rank = mean(token_ranks)
  )
}) %>%
  bind_rows() %>%
  arrange(first_rank, mean_rank) %>%
  mutate(display_token = letters[row_number()])

relabel_cld <- function(raw_label) {
  mapped <- letter_key$display_token[
    match(strsplit(raw_label, "")[[1]], letter_key$token)
  ]
  paste(sort(mapped), collapse = "")
}

cld_df <- cld_df %>%
  mutate(cld = vapply(cld_raw, relabel_cld, character(1))) %>%
  mutate(
    Biome_label = factor(Biome_label, levels = levels(plot_table$Biome_label)),
    x_pos = 1.52
  )

summary_out <- summary_df %>%
  left_join(
    cld_df %>% select(Biome_label_raw, cld, cld_raw),
    by = "Biome_label_raw"
  ) %>%
  mutate(kruskal_p = kw_res$p)

write.csv(summary_out, out_summary, row.names = FALSE)

# A symmetric range keeps zero white while using the available colour contrast.
colour_limit <- max(abs(summary_df$median_log10), na.rm = TRUE)
colour_limit <- ceiling(colour_limit * 100) / 100

p_box <- ggplot(plot_table, aes(x = log10_ratio, y = Biome_label)) +
  geom_point(
    aes(color = median_log10),
    position = position_jitter(width = 0, height = 0.12),
    alpha = 0.16,
    size = 0.64,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(fill = median_log10),
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.82,
    color = "black",
    linewidth = 0.58,
    show.legend = FALSE
  ) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  geom_text(
    data = cld_df,
    aes(x = x_pos, y = Biome_label, label = cld),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5.2,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = palette_low,
    mid = palette_mid,
    high = palette_high,
    midpoint = 0,
    limits = c(-colour_limit, colour_limit),
    oob = squish
  ) +
  scale_color_gradient2(
    low = palette_low,
    mid = "#BFC2BE",
    high = palette_high,
    midpoint = 0,
    limits = c(-colour_limit, colour_limit),
    oob = squish
  ) +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1, 1.5),
    limits = c(-1.45, 1.68),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression("Uptake-to-depolymerisation ratio ("*log[10]*")"),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10.5, color = "black"),
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(8, 30, 8, 8)
  )

data("Whittaker_biomes", package = "plotbiomes")

whittaker_df <- Whittaker_biomes %>%
  left_join(
    summary_df %>%
      transmute(biome = Biome_label_raw, median_log10),
    by = "biome"
  )

# Manual positions keep labels inside the narrow biome polygons.
label_positions <- tibble::tribble(
  ~biome, ~display_name, ~x, ~y, ~angle,
  "Tundra", "Tundra", -9.8, 17, 24,
  "Boreal forest", "Boreal forest", -2.3, 84, 0,
  "Temperate grassland/desert", "Temperate grassland/\ndesert", 7.3, 13, 0,
  "Woodland/shrubland", "Woodland/\nshrubland", 12.0, 56, 20,
  "Temperate seasonal forest", "Temperate seasonal\nforest", 10.4, 154, 0,
  "Temperate rain forest", "Temperate rain forest", 11.8, 272, 0,
  "Tropical rain forest", "Tropical rain forest", 24.5, 380, 0,
  "Tropical seasonal forest/savanna", "Tropical seasonal\nforest/savanna", 25.0, 170, 0,
  "Subtropical desert", "Subtropical\ndesert", 26.5, 48, 0
) %>%
  left_join(
    summary_df %>% transmute(biome = Biome_label_raw, median_log10),
    by = "biome"
  ) %>%
  filter(!is.na(median_log10)) %>%
  mutate(
    label = paste0(display_name, "\n(", sprintf("%.2f", median_log10), ")")
  )

p_whittaker <- ggplot(
  whittaker_df,
  aes(x = temp_c, y = precp_cm, group = biome, fill = median_log10)
) +
  geom_polygon(color = "white", linewidth = 1.05) +
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label, angle = angle),
    inherit.aes = FALSE,
    size = 2.75,
    lineheight = 0.86,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = palette_low,
    mid = palette_mid,
    high = palette_high,
    midpoint = 0,
    limits = c(-colour_limit, colour_limit),
    oob = squish,
    na.value = "#E3E5E4",
    breaks = c(-colour_limit, 0, colour_limit),
    labels = label_number(accuracy = 0.01),
    name = expression("Median ratio ("*log[10]*")")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.2, "cm"),
      barheight = unit(0.30, "cm"),
      direction = "horizontal"
    )
  ) +
  scale_x_continuous(
    limits = c(-16, 31),
    breaks = c(-10, 0, 10, 20, 30),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(-20, 455),
    breaks = c(0, 100, 200, 300, 400),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression("Temperature ("*degree*C*")"),
    y = "Precipitation (cm)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.25, 0.93),
    legend.justification = c(0.5, 1),
    legend.background = element_blank(),
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 8.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    plot.margin = margin(8, 10, 8, 10)
  )

# Panel c: within-site comparisons for locations containing both explicitly
# labelled bulk-soil and rhizosphere MAG datasets. These are independent MAG
# comparisons within each location, not one-to-one paired observations.
metadata <- read_excel(
  metadata_xlsx,
  sheet = "IMG_bindata_withmeta_norestrict"
) %>%
  transmute(
    Bin.ID = as.character(Bin.ID),
    Project = as.character(Project),
    Geographic.Location = as.character(Geographic.Location),
    project_lower = str_to_lower(coalesce(Project, "")),
    geographic_lower = str_to_lower(coalesce(Geographic.Location, "")),
    habitat_lower = str_to_lower(coalesce(as.character(Habitat), "")),
    ecosystem_lower = str_to_lower(coalesce(as.character(Specific.Ecosystem), "")),
    has_rhizosphere = str_detect(
      paste(habitat_lower, ecosystem_lower, project_lower),
      "rhizosphere"
    ),
    has_bulk_soil = str_detect(
      paste(habitat_lower, ecosystem_lower, project_lower),
      "bulk soil"
    ),
    soil_zone = case_when(
      has_rhizosphere & !has_bulk_soil ~ "Rhizosphere",
      has_bulk_soil & !has_rhizosphere ~ "Bulk soil",
      TRUE ~ NA_character_
    ),
    site_metadata_conflict =
      (str_detect(geographic_lower, "clatskanie") &
         str_detect(project_lower, "corvallis")) |
      (str_detect(geographic_lower, "corvallis") &
         str_detect(project_lower, "clatskanie")),
    site = case_when(
      str_detect(geographic_lower, "clatskanie") ~
        "Poplar - Clatskanie",
      str_detect(geographic_lower, "corvallis") ~
        "Poplar - Corvallis",
      str_detect(geographic_lower, "tennessee") ~
        "Poplar - Tennessee",
      str_detect(project_lower, "poplar common garden site in clatskanie") ~
        "Poplar - Clatskanie",
      str_detect(project_lower, "poplar common garden site in corvallis") ~
        "Poplar - Corvallis",
      str_detect(project_lower, "populus") & str_detect(project_lower, "tennessee") ~
        "Poplar - Tennessee",
      TRUE ~ NA_character_
    )
  )

site_levels <- c(
  "Poplar - Clatskanie",
  "Poplar - Corvallis",
  "Poplar - Tennessee"
)
zone_levels <- c("Bulk soil", "Rhizosphere")
zone_cols <- c("Bulk soil" = "#3B6FB6", "Rhizosphere" = "#D9772B")

site_joined <- merged_table %>%
  transmute(
    Bin.ID = as.character(Bin.ID),
    ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  inner_join(metadata, by = "Bin.ID") %>%
  filter(
    !is.na(site),
    !is.na(soil_zone),
    is.finite(ratio),
    ratio > 0
  )

excluded_conflicting_site <- site_joined %>%
  filter(site_metadata_conflict) %>%
  select(Bin.ID, Project, Geographic.Location, soil_zone, site, ratio)

write.csv(excluded_conflicting_site, out_site_exclusions, row.names = FALSE)

if (nrow(excluded_conflicting_site) != 29) {
  warning(
    "Expected 29 MAGs with conflicting Project and Geographic.Location labels; found ",
    nrow(excluded_conflicting_site),
    "."
  )
}

site_df <- site_joined %>%
  filter(!site_metadata_conflict) %>%
  mutate(
    log10_ratio = log10(ratio),
    site = factor(site, levels = site_levels),
    soil_zone = factor(soil_zone, levels = zone_levels)
  )

run_site_test <- function(df) {
  x <- df$log10_ratio[df$soil_zone == "Rhizosphere"]
  y <- df$log10_ratio[df$soil_zone == "Bulk soil"]
  test <- wilcox.test(x, y, exact = FALSE)
  combined_ranks <- rank(c(x, y), ties.method = "average")
  u_rhizosphere <- sum(combined_ranks[seq_along(x)]) -
    length(x) * (length(x) + 1) / 2

  tibble::tibble(
    n_bulk = length(y),
    n_rhizosphere = length(x),
    median_difference = median(x) - median(y),
    rank_biserial = 2 * u_rhizosphere / (length(x) * length(y)) - 1,
    p_value = test$p.value,
    effect_magnitude = delta_magnitude(
      2 * u_rhizosphere / (length(x) * length(y)) - 1
    )
  )
}

site_tests <- site_df %>%
  group_split(site) %>%
  lapply(function(df) {
    run_site_test(df) %>% mutate(site = as.character(unique(df$site)), .before = 1)
  }) %>%
  bind_rows()

site_summary <- site_df %>%
  group_by(site, soil_zone) %>%
  summarise(
    n = n(),
    mean_log10 = mean(log10_ratio),
    median_log10 = median(log10_ratio),
    q25_log10 = quantile(log10_ratio, 0.25),
    q75_log10 = quantile(log10_ratio, 0.75),
    .groups = "drop"
  )

write.csv(site_summary, out_site_summary, row.names = FALSE)
write.csv(site_tests, out_site_tests, row.names = FALSE)

format_site_p <- function(p) {
  ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, digits = 3, format = "f")))
}

site_count_labels <- site_summary %>%
  mutate(label = paste0("n = ", n), y = -1.27)

site_p_labels <- site_tests %>%
  mutate(
    site = factor(site, levels = site_levels),
    label = paste0(
      format_site_p(p_value),
      "\nCliff's delta = ", sprintf("%.2f", rank_biserial),
      " (", effect_magnitude, ")"
    ),
    x = 1.5,
    y = 1.44
  )

p_sites <- ggplot(site_df, aes(soil_zone, log10_ratio, color = soil_zone)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.42,
    color = "grey65"
  ) +
  geom_violin(
    aes(fill = soil_zone),
    width = 0.82,
    trim = TRUE,
    alpha = 0.13,
    linewidth = 0.75
  ) +
  geom_jitter(width = 0.13, size = 1.02, alpha = 0.48, stroke = 0) +
  geom_boxplot(
    aes(fill = soil_zone),
    width = 0.14,
    outlier.shape = NA,
    color = "grey20",
    linewidth = 0.58,
    alpha = 0.82
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.1,
    stroke = 0.8,
    fill = "white",
    colour = "grey15"
  ) +
  geom_text(
    data = site_count_labels,
    aes(x = soil_zone, y = y, label = label, color = soil_zone),
    inherit.aes = FALSE,
    size = 3.0,
    fontface = "bold"
  ) +
  geom_text(
    data = site_p_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.9,
    lineheight = 1.0
  ) +
  facet_wrap(~site, nrow = 1) +
  scale_color_manual(values = zone_cols) +
  scale_fill_manual(values = zone_cols) +
  scale_y_continuous(
    limits = c(-1.40, 1.55),
    breaks = c(-1, 0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = NULL,
    y = expression("Uptake-to-depolymerisation ratio (" * log[10] * ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.3),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 9.5),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(10, 8, 8, 8)
  )

combined <- (p_box + p_whittaker) +
  plot_layout(widths = c(1.08, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(out_png, combined, width = 11.4, height = 5.35, dpi = 400, bg = "white")
ggsave(out_pdf, combined, width = 11.4, height = 5.35, bg = "white")
ggsave(out_site_png, p_sites, width = 13.4, height = 4.3, dpi = 400, bg = "white")
ggsave(out_site_pdf, p_sites, width = 13.4, height = 4.3, bg = "white")

message("Saved: ", out_png)
message("Saved: ", out_pdf)
message("Saved: ", out_site_png)
message("Saved: ", out_site_pdf)
message("Saved: ", out_summary)
message("Saved: ", out_biome_pairs)
message("Saved: ", out_site_summary)
message("Saved: ", out_site_tests)
message("Ratio colour scale: ", -colour_limit, " to ", colour_limit)
