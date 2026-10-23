#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[1])
base_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

in_csv <- file.path(base_dir, "data/processed/loma_input.csv")
out_pdf <- file.path(base_dir, "results/figures/Figure_5_family_panel.pdf")
out_png <- file.path(base_dir, "results/figures/Figure_5_family_panel.png")
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

time_levels <- c("T1", "T2", "T3", "T4")
time_labels <- c(T1 = "2 months", T2 = "7 months", T3 = "14 months", T4 = "17 months")
time_cols <- c(T1 = "#D52B34", T2 = "#F39A72", T3 = "#8CC5DE", T4 = "#176AA5")

# Text-safe, darker versions of the three hues assigned in Figure 3. The hue
# identity is preserved while improving contrast on a white background.
phylum_cols <- c(
  Pseudomonadota = "#C84F49",
  Bacteroidota = "#9A6F00",
  Actinomycetota = "#707600"
)

dat <- read_csv(in_csv, show_col_types = FALSE) %>%
  transmute(
    phylum = sub("^p__", "", as.character(Phylum)),
    family = sub("^f__", "", as.character(Family)),
    time = factor(as.character(time), levels = time_levels),
    ratio = as.numeric(ABC_GH_ratio)
  ) %>%
  filter(is.finite(ratio), ratio > 0, !is.na(time)) %>%
  mutate(log10_ratio = log10(ratio))

make_panel <- function(
    data,
    taxon,
    order,
    show_legend = FALSE,
    label_cols = NULL,
    show_phylum_key = FALSE
) {
  plot_data <- data %>%
    filter(.data[[taxon]] %in% order) %>%
    mutate(taxon_label = factor(.data[[taxon]], levels = rev(order)))

  p <- ggplot(plot_data, aes(log10_ratio, taxon_label)) +
    geom_violin(fill = "grey91", colour = NA, scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = time), height = 0.10, width = 0, size = 1.05, alpha = 0.78) +
    stat_summary(fun = median, geom = "point", shape = 124, size = 7.8, colour = "grey10") +
    scale_colour_manual(
      values = time_cols,
      breaks = time_levels,
      labels = time_labels,
      name = "Time since start of litter\ndecomposition experiment"
    ) +
    labs(x = expression("Uptake-to-depolymerisation ratio ("*log[10]*")"), y = NULL) +
    theme_classic(base_size = 10.5) +
    theme(
      axis.text.y = if (is.null(label_cols)) {
        element_text(face = "bold", size = 9.5, colour = "grey15")
      } else {
        ggtext::element_markdown(
          face = "bold", size = 9.2, colour = "grey15",
          margin = margin(r = 4)
        )
      },
      axis.text.x = element_text(size = 9, colour = "grey20"),
      axis.title.x = element_text(size = 10.5, margin = margin(t = 5)),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.28),
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.985, 0.245),
      legend.justification = c(1, 0.5),
      legend.title = element_text(size = 7.8, face = "bold", lineheight = 0.9),
      legend.text = element_text(size = 7.6, colour = "grey20"),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.spacing.y = grid::unit(0.02, "cm"),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.90),
        colour = "grey82",
        linewidth = 0.25
      ),
      legend.margin = margin(3, 5, 3, 5),
      plot.margin = margin(7, 9, 7, 7)
    )

  if (!is.null(label_cols)) {
    p <- p + scale_y_discrete(labels = function(x) {
      cols <- unname(label_cols[as.character(x)])
      cols[is.na(cols)] <- "black"
      sprintf("<span style='color:%s'>%s</span>", cols, x)
    })
  }

  if (show_phylum_key) {
    phylum_key <- paste0(
      "<span style='color:#3F3F3F'><b>Phylum</b></span><br>",
      "<span style='color:", phylum_cols[["Pseudomonadota"]],
      "'><b>Pseudomonadota</b></span><br>",
      "<span style='color:", phylum_cols[["Actinomycetota"]],
      "'><b>Actinomycetota</b></span><br>",
      "<span style='color:", phylum_cols[["Bacteroidota"]],
      "'><b>Bacteroidota</b></span>"
    )

    p <- p + annotate(
      "richtext",
      x = -1.05,
      y = Inf,
      label = phylum_key,
      hjust = 0,
      vjust = 1.05,
      size = 2.75,
      lineheight = 1.08,
      fill = scales::alpha("white", 0.94),
      label.colour = "grey82",
      label.size = 0.25,
      label.r = grid::unit(0.08, "lines"),
      label.padding = grid::unit(c(3, 4, 3, 4), "pt")
    )
  }

  p
}

# Orders are specified from top to bottom; make_panel reverses them for the
# underlying factor levels used by ggplot2's y axis.
phylum_order <- c("Actinomycetota", "Pseudomonadota", "Bacteroidota")
family_order <- c(
  "Enterobacteriaceae", "Pseudomonadaceae", "Rhizobiaceae",
  "Rhodobacteraceae", "Geodermatophilaceae", "Micrococcaceae",
  "Microbacteriaceae", "Kineococcaceae", "Beijerinckiaceae",
  "Micromonosporaceae", "Burkholderiaceae", "Cellulomonadaceae",
  "Acetobacteraceae", "Xanthomonadaceae", "Sphingomonadaceae",
  "Weeksellaceae", "Spirosomataceae", "Hymenobacteraceae",
  "Sphingobacteriaceae", "Chitinophagaceae"
)

family_phylum <- dat %>%
  filter(family %in% family_order, phylum %in% names(phylum_cols)) %>%
  count(family, phylum, sort = TRUE) %>%
  group_by(family) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

family_label_cols <- setNames(
  unname(phylum_cols[family_phylum$phylum]),
  family_phylum$family
)

p_a <- make_panel(dat, "phylum", phylum_order) + ggtitle("a") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = -0.08))
p_b <- make_panel(
  dat,
  "family",
  family_order,
  show_legend = TRUE,
  label_cols = family_label_cols,
  show_phylum_key = TRUE
) + ggtitle("b") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = -0.08))

combined <- p_a / p_b + plot_layout(heights = c(1, 2.55))
ggsave(out_pdf, combined, width = 5.2, height = 10.4, bg = "white")
ggsave(out_png, combined, width = 5.2, height = 10.4, dpi = 400, bg = "white")
