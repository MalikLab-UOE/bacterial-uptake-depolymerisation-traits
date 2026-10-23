#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

in_rds <- file.path(base_dir, "data/processed/glo_input_binary.rds")
out_pdf <- file.path(base_dir, "results/figures/Supplementary_Figure_S2.pdf")
out_png <- file.path(base_dir, "results/figures/Supplementary_Figure_S2.png")
out_summary <- file.path(base_dir, "results/tables/GLO_phylum_binary_traits_summary.csv")
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

palette_low <- "#276C58"
palette_mid <- "#F7F7F4"
palette_high <- "#C92C28"

phylum_order <- c(
  "Actinomycetota", "Cyanobacteriota", "Spirochaetota", "Bacillota",
  "Omnitrophota", "Pseudomonadota", "Chloroflexota", "Desulfobacterota",
  "Campylobacterota", "Planctomycetota", "Bacteroidota_A", "Nitrospirota",
  "Myxococcota", "Gemmatimonadota", "Verrucomicrobiota", "Acidobacteriota",
  "Bdellovibrionota", "Patescibacteriota", "Bacteroidota"
)

dat <- readRDS(in_rds) %>%
  transmute(
    phylum = sub("^p__", "", as.character(Phylum)),
    depolymerisation = as.numeric(GH_richness),
    uptake = as.numeric(ABC_richness),
    ratio = as.numeric(ABC_GH_ratio),
    genome_size_Mb = as.numeric(genome_size) / 1e6
  ) %>%
  filter(
    phylum %in% phylum_order,
    is.finite(depolymerisation),
    is.finite(uptake),
    is.finite(ratio), ratio > 0,
    is.finite(genome_size_Mb), genome_size_Mb > 0
  ) %>%
  mutate(log10_ratio = log10(ratio))

counts <- dat %>% count(phylum, name = "n")
dat <- dat %>%
  left_join(counts, by = "phylum") %>%
  mutate(
    phylum_label = paste0(phylum, " (n=", n, ")"),
    phylum_label = factor(
      phylum_label,
      levels = paste0(phylum_order, " (n=", counts$n[match(phylum_order, counts$phylum)], ")")
    )
  )

long <- dat %>%
  transmute(
    phylum_label,
    colour_value = log10_ratio,
    `Complex carbohydrate depolymerisation enzymes` = depolymerisation,
    `Carbohydrate transporters` = uptake,
    `Uptake-to-depolymerisation ratio (log10)` = log10_ratio,
    `Genome size (Mb)` = genome_size_Mb
  ) %>%
  pivot_longer(
    cols = -c(phylum_label, colour_value),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c(
    "Complex carbohydrate depolymerisation enzymes",
    "Carbohydrate transporters",
    "Uptake-to-depolymerisation ratio (log10)",
    "Genome size (Mb)"
  )))

# Match the 0-10 Mb display range used in Figure 4 only for the genome-size
# panel. The three functional-trait panels retain all otherwise valid MAGs.
long <- long %>%
  filter(metric != "Genome size (Mb)" | value <= 10)

medians <- long %>%
  group_by(phylum_label, metric) %>%
  summarise(median = median(value), .groups = "drop")

summary_table <- dat %>%
  group_by(phylum) %>%
  summarise(
    n = n(),
    n_genome_size_0_10_Mb = sum(genome_size_Mb <= 10),
    median_ABC_GH_ratio = median(ratio),
    median_ABC_richness = median(uptake),
    median_GH_richness = median(depolymerisation),
    median_genome_size_Mb = median(genome_size_Mb),
    .groups = "drop"
  )
write_csv(summary_table, out_summary)

set.seed(20261005)
p <- ggplot(long, aes(value, 0, colour = colour_value)) +
  geom_jitter(height = 0.23, width = 0, size = 0.42, alpha = 0.34) +
  geom_vline(
    data = medians,
    aes(xintercept = median),
    inherit.aes = FALSE,
    linewidth = 0.7,
    linetype = "longdash"
  ) +
  facet_grid(
    phylum_label ~ metric,
    scales = "free_x",
    switch = "y",
    labeller = labeller(metric = as_labeller(c(
      "Complex carbohydrate depolymerisation enzymes" =
        "bold('Complex carbohydrate depolymerisation enzymes')",
      "Carbohydrate transporters" = "bold('Carbohydrate transporters')",
      "Uptake-to-depolymerisation ratio (log10)" =
        "bold('Uptake-to-depolymerisation ratio (' * log[10] * ')')",
      "Genome size (Mb)" = "bold('Genome size (Mb)')"
    ), default = label_parsed))
  ) +
  scale_colour_gradient2(
    low = palette_low, mid = palette_mid, high = palette_high,
    midpoint = 0, limits = c(-2, 2), oob = scales::squish,
    name = expression(atop(
      "Uptake-to-depolymerisation",
      "ratio (" * log[10] * ")"
    ))
  ) +
  guides(colour = guide_colourbar(
    direction = "horizontal",
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(2.45, "cm"),
    barheight = grid::unit(0.30, "cm")
  )) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y.left = element_text(
      face = "bold",
      angle = 0,
      size = 11,
      hjust = 1,
      margin = margin(r = 6)
    ),
    panel.spacing.y = unit(0.02, "lines"),
    panel.spacing.x = unit(0.65, "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.49, 0.14),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.90),
      colour = "grey55",
      linewidth = 0.25
    ),
    legend.margin = margin(2, 4, 2, 4),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

ggsave(out_pdf, p, width = 18, height = 7.5, bg = "white")
ggsave(out_png, p, width = 18, height = 7.5, dpi = 400, device = ragg::agg_png, bg = "white")
