# Bacterial uptake-to-depolymerisation trait analyses

This repository contains the reproducible analyses, analysis-ready data, statistical outputs and final manuscript figures for the bacterial uptake-to-depolymerisation trait study. It covers main Figures 1–5, Supplementary Figures S1–S3 and Data S1.

## Quick start

This repository uses Git LFS for large PDFs and binary data. Clone it with Git LFS enabled:

```bash
git lfs install
git clone <repository-url>
cd <repository-name>
```

Install the R and Python dependencies:

```bash
make setup
```

Validate the archived inputs and deliverables:

```bash
make validate
```

Regenerate the main or supplementary analyses with `make main` or `make supplementary`. Figure 3 is skipped by the default main-figure target because rendering 249,310 tree tips is memory intensive; use `make main-with-tree` on a high-memory system.

## Binary repertoire definition

The uptake and depolymerisation repertoires are presence-absence repertoire sizes, not gene-copy abundances. For each MAG, every matched transporter or depolymerisation feature was converted to 0/1 and then summed within its functional category. The uptake-to-depolymerisation ratio is:

`uptake repertoire / depolymerisation repertoire`

MAGs with zero uptake or zero depolymerisation repertoire were excluded before calculating the ratio. Ratios were log10-transformed for analyses and plotting where stated. `data/processed/glo_input_binary.rds`, `data/processed/img_input.rds` and `data/processed/loma_input.csv` are the analysis-ready binary-repertoire inputs.

## Data S1

`data/supplementary/Data_S1_functional_gene_sets.xlsx` documents the feature sets used to construct the two repertoires. It contains 41 selected CAZy glycoside-hydrolase families for complex-carbohydrate depolymerisation and 211 selected carbohydrate-transporter KEGG Orthology entries for uptake. Two distinct KO entries share the gene name `npr`; their source descriptions were unavailable and are explicitly flagged rather than inferred.

## Directory structure

- `figures/`: final manuscript Figures 1-5 and Supplementary Figures S1-S3.
- `data/supplementary/`: Data S1 functional-gene mapping table.
- `code/`: R scripts for the corresponding analyses and plots.
- `data/processed/`: analysis-ready data read directly by the scripts.
- `data/raw/`: source metadata, detection tables and the GlobDB phylogeny needed by the analyses.
- `data/figure_inputs/`: compact human-readable plotting tables.
- `results/`: statistical outputs and regenerated component plots.
- `docs/`: analysis audit, data dictionary, validation record and software-session information.
- `.github/workflows/`: lightweight syntax checks run by GitHub Actions.
- `requirements.txt` and `code/setup/install_r_dependencies.R`: Python and R dependency specifications.
- `code/validate_binary_inputs.R`: reproducible checks for the binary feature values, repertoire totals and ratios present in the archive.
- `code/input_builders/`: scripts documenting how the GlobDB, IMG and Loma binary repertoire inputs were constructed. The original large annotation matrices are not redistributed in this archive.

## Figure-to-file map

### Figure 1

- Final: `figures/Figure_1.pdf`
- Code: `code/Figure_1.R`
- Inputs: `data/processed/glo_input_binary.rds`, `data/processed/all_db_merged.csv`, `data/raw/high_level_env__value.csv`
- Statistics: `results/tables/Fig_1b_corrected_habitat_summary.csv`, `Fig_1b_corrected_Tukey_pairwise_comparisons.csv`, and `Fig_1b_corrected_model_summary.csv`

Each MAG was assigned to the habitat with its highest detection count. The selected detection count was used as the analytical weight in the violin density, boxplot summaries, weighted least-squares model and weighted Cliff's delta. A weighted `n` is therefore the sum of genome detections, not the number of unique MAGs. In the weighted Cliff's delta, each cross-habitat MAG pair is weighted by the product of the two MAG detection counts.

### Figure 2

- Final: `figures/Figure_2.pdf`
- Code: `code/Figure_2.R`
- Inputs: `data/processed/img_input.rds`, `data/raw/Original_eco.xlsx`
- Statistics: `results/tables/Figure_2_ratio_coloured_biome_summary.csv`, `Figure_2_biome_pairwise_cliffs_delta.csv`

In panel b, biome colours and the values in parentheses both show the median log10-transformed uptake-to-depolymerisation ratio. Biomes without qualifying data are omitted.

### Figure 3

- Final: `figures/Figure_3.pdf`
- Code: `code/Figure_3_tree_visualisation.R`, `code/Figure_3_Pagel_lambda.R`
- Inputs: `data/raw/bacteria_globdb_r226.tree`, `data/processed/glo_input_binary.rds`
- Statistics: `results/tables/continuous_ratio_full_tree_pagel_lambda.csv`

The final PDF has manually positioned phylum labels. The plotting script reproduces the binary-repertoire ratio ring, tree and phylum colours; label placement can vary slightly among graphics devices.

### Figure 4

- Final: `figures/Figure_4.pdf`
- Code: `code/Figure_4.R`
- Input: `data/processed/glo_input_binary.rds`
- Statistics: `results/tables/Genome_size_GLM_vs_GAM_comparison.csv`

### Figure 5

- Final: `figures/Figure_5.pdf`
- Code: `code/Figure_5_panels_a_c.R`, `code/Figure_5_family_panel.R`, `code/Figure_5_panel_b.R`, `code/Figure_5_panel_c_GH_genomes.R`, `code/Figure_5_panel_c_GH_transcripts.R`, `code/assemble_figure5_with_family.py`
- Inputs: `data/processed/loma_input.csv`, `data/raw/LOMA_merged_mags_50_10.csv`, `data/raw/LOMA_FBratio.txt`, `data/raw/GH_tax_contribuction.csv`

The final layout was assembled from the supplied component plots. The family-level resource-acquisition plot is panel b of Figure 5 and is not a separate supplementary figure. Panel a uses pairwise Wilcoxon tests with BH correction and Cliff's delta. The fungal-to-bacterial comparison uses one-way ANOVA and Tukey HSD. The genome and transcript comparisons use a Wilcoxon test and Cliff's delta for GH gene numbers and Welch's t-test for transcript abundances.

## Supplementary figures

### Figure S1: rhizosphere versus bulk soil

- Final: `figures/supplementary/Figure_S1.pdf`
- Code: `code/Figure_2.R`
- Inputs: `data/processed/img_input.rds`, `data/raw/Original_eco.xlsx`
- Statistics: `results/tables/Figure_2_bulk_rhizosphere_site_summary.csv`, `Figure_2_bulk_rhizosphere_site_tests.csv`

The plotted Wilcoxon comparisons treat MAGs as independent observations, as stated in the current Methods. Twenty-nine MAGs with conflicting project and geographic-location site labels were excluded before plotting and are listed in `results/tables/Figure_2_bulk_rhizosphere_conflicting_site_metadata_excluded.csv`.

### Figure S2: phylum-level genomic traits

- Final: `figures/supplementary/Figure_S2.pdf`
- Code: `code/Supplementary_Figure_S2.R`
- Input: `data/processed/glo_input_binary.rds`
- Summary: `results/tables/GLO_phylum_binary_traits_summary.csv`

Genome size is displayed in megabases from 0 to 10 Mb. The uptake-to-depolymerisation ratio is the log10-transformed binary repertoire ratio.
The four trait panels are arranged in one row, with the colour key inset in the transporter panel to reduce unused whitespace.

### Figure S3: phylum-level genome-size GAMs

- Final: `figures/supplementary/Figure_S3.pdf`
- Code: `code/Supplementary_Figure_S3.R`
- Input: `data/processed/glo_input_binary.rds`
- Statistics: `results/tables/Genome_size_GAM_by_phylum_stats.csv`, `results/tables/Genome_size_top_bottom10_by_phylum.csv`

For the complete dataset and each of the seven most abundant phyla, the figure pairs a genome-size GAM with a comparison of genome size between MAGs in the dataset-specific bottom and top 10% of the uptake-to-depolymerisation ratio. All MAGs tied at a decile cutoff are retained. Violin and boxplot colours represent the median log10 ratio of each extreme group on the same absolute colour scale used for the GAM points. Extreme-group comparisons use a two-sided Wilcoxon rank-sum test and Cliff's delta; different letters require both `P < 0.05` and `|delta| >= 0.147`.

## Running the analyses

Run an individual script from any working directory, for example:

```bash
Rscript code/Figure_1.R
```

Run the main-figure analyses with `Rscript run_all.R` or `make main`. Run Supplementary Figures S1–S3 with `Rscript run_supplementary.R` or `make supplementary`. After regenerating the Figure 5 components, assemble the final layout with `make figure5`; this step uses `pypdf` and `reportlab`.

Run all local validation checks with `make validate`. This parses every R and Python script, verifies the binary repertoire inputs, reopens every PDF and checks the Data S1 entry counts.

Figure 3 is skipped by default because rendering 249,310 tips is memory intensive. Include it on a high-memory system with `Rscript run_all.R --include-tree`.

## Large files

The repository is approximately 359 MB and contains several binary files larger than 10 MB. `.gitattributes` assigns PDFs, PNGs, RDS files, Excel files, compressed CSV files, tree files and the largest uncompressed CSV to Git LFS. Install Git LFS before the first `git add` or use GitHub Desktop, which includes Git LFS support.

## Documentation

- Analysis scope and figure mapping: `docs/analysis_audit.md`
- Data dictionary: `docs/data_dictionary.md`
- Validation record: `docs/validation.md`
- Software session: `docs/session_info.txt`

## Citation and license

Add the final manuscript citation and DOI when available. No software or data license has been assigned yet; the authors should select an appropriate license before making the repository public.

## Notes

- `n` denotes unique MAGs unless explicitly labelled as weighted `n` in Figure 1b.
- Compact-letter displays combine adjusted statistical significance with the prespecified absolute Cliff's delta threshold of 0.147 only where the corresponding methods and script explicitly state this rule.
- Source-database redistribution remains subject to the original providers' terms.
