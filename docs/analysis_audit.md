# Figure and analysis consistency audit

Audit updated: 21 October 2026

This audit covers only the analyses represented in main Figures 1–5, Supplementary Figures S1–S3 and Data S1.

## Figure 1

- The uptake-to-depolymerisation ratio is calculated from binary presence–absence repertoire sizes.
- Each MAG is assigned to the habitat with its highest detection count.
- Figure 1b uses the selected detection count as the analytical weight in the plotted distribution, weighted summaries, weighted least-squares model and weighted Cliff's delta.
- The archived statistical outputs are `Fig_1b_corrected_habitat_summary.csv`, `Fig_1b_corrected_Tukey_pairwise_comparisons.csv` and `Fig_1b_corrected_model_summary.csv`.

## Figure 2

- The analysis uses the IMG binary uptake-to-depolymerisation repertoire ratio.
- Biome comparisons use a Kruskal–Wallis test followed by BH-adjusted Dunn tests.
- Cliff's delta is used with the statistical comparisons for the displayed group letters.
- Biome colours and parenthetical values represent median log10-transformed ratios.

## Figure 3

- The phylogenetic tree contains 249,310 GlobDB MAGs matched to the binary-repertoire input.
- The ratio ring is based on `glo_input_binary.rds`.
- Pagel's lambda is calculated for the continuous log10-transformed uptake-to-depolymerisation ratio.

## Figure 4

- The genome-size analysis uses 248,230 MAGs after applying the scripted validity filters.
- Genome size is expressed in Mb.
- The archived comparison table reports the linear-model and GAM fits.

## Figure 5

- Panel a compares temporal changes in depolymerisation repertoire, uptake repertoire and their log10 ratio using pairwise Wilcoxon tests, BH correction and Cliff's delta.
- Panel b shows family-level uptake-to-depolymerisation distributions; family-label colours identify phyla and black bars show medians.
- Panel c shows the fungal-to-bacterial abundance analysis using one-way ANOVA followed by Tukey HSD.
- Panel d compares bacterial and fungal glycoside-hydrolase gene numbers using a Wilcoxon rank-sum test and Cliff's delta, and transcript abundance using Welch's t-test.

## Figure S1

- Bulk-soil and rhizosphere MAGs are compared separately within Clatskanie, Corvallis and Tennessee.
- The displayed statistics are MAG-level Wilcoxon rank-sum tests and Cliff's delta.
- Twenty-nine MAGs with conflicting project and geographic-location site assignments are excluded and listed in `Figure_2_bulk_rhizosphere_conflicting_site_metadata_excluded.csv`.

## Figure S2

- The four phylum-level panels are depolymerisation repertoire, uptake repertoire, log10 uptake-to-depolymerisation ratio and genome size.
- Genome size is expressed in Mb on a common 0–10 Mb scale.
- The colour key is inset in the transporter panel.

## Figure S3

- The figure contains the complete GlobDB dataset and the seven most abundant phyla.
- Each dataset has a genome-size GAM and a bottom-versus-top ratio-decile comparison.
- Extreme-group comparisons use two-sided Wilcoxon rank-sum tests and Cliff's delta.

## Data S1

- The workbook contains 41 unique CAZy glycoside-hydrolase families used for the depolymerisation repertoire.
- It contains 211 unique carbohydrate-transporter KEGG Orthology entries used for the uptake repertoire.
- Two KO entries with unavailable source descriptions are explicitly flagged.

## Archive scope

The archive contains only the inputs, code, component outputs and statistical tables needed to reproduce Figures 1–5, Figures S1–S3 and Data S1. Analyses not represented in these numbered outputs are excluded.
