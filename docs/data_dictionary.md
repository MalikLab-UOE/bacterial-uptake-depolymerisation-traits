# Data dictionary

## Binary repertoire variables

- `GH_richness` / `GH_total_abundance`: count of depolymerisation features present after converting each matched feature to 0/1.
- `ABC_richness` / `ABC_total_abundance`: count of uptake features present after converting each matched feature to 0/1.
- `ABC_GH_ratio`: uptake repertoire divided by depolymerisation repertoire.

Some inherited column names contain the word `abundance`; in the analysis-ready files these columns are binary repertoire totals, not gene-copy abundances.

## `data/supplementary/Data_S1_functional_gene_sets.xlsx`

The workbook is the feature-level definition of the two resource-acquisition traits.

- `Depolymerisation`: 41 selected CAZy GH families, their substrate classes and source descriptions.
- `Uptake`: 211 selected carbohydrate-transporter KO entries, gene identifiers and source descriptions.
- `Included_in_analysis`: `Yes` for every feature used by the current binary-input builders.
- `Description_status`: records whether an uptake-feature description was matched to the source annotation table.

The workbook defines trait membership; it does not contain genome-level measurements.

## `GlobDB_MAG_traits.csv.gz`

- `genome_id`: GlobDB genome identifier.
- `phylum`: bacterial phylum assignment.
- `depolymerisation_repertoire`: binary depolymerisation repertoire size.
- `uptake_repertoire`: binary uptake repertoire size.
- `uptake_to_depolymerisation_ratio`: uptake repertoire divided by depolymerisation repertoire.
- `genome_size_bp`: genome size in base pairs.

## `IMG_soil_MAG_biome_traits.csv.gz`

- `Bin_ID`: IMG soil MAG identifier.
- `biome`: Whittaker biome assignment.
- `mean_annual_temperature_C`: mean annual temperature at the sampling location.
- `annual_precipitation_mm`: annual precipitation at the sampling location.
- `uptake_to_depolymerisation_ratio`: binary uptake-to-depolymerisation repertoire ratio.

## `GlobDB_dominant_habitat_assignment.csv`

- `genome_id`: GlobDB genome identifier.
- `biome`: high-level habitat with the greatest detection count for that MAG.
- `detect_count`: number of detections in the selected habitat; this is the Figure 1b analytical weight and its sum is the displayed weighted `n`.

## `Loma_MAG_temporal_traits.csv`

- `id`: MAG identifier.
- `time`: litter-decomposition sampling time (`T1`-`T4`).
- `GH_total_abundance`: binary depolymerisation repertoire size.
- `ABC_total_abundance`: binary uptake repertoire size.
- `ABC_GH_ratio`: uptake-to-depolymerisation repertoire ratio.

## `Loma_MAG_taxonomic_traits.csv`

- `Bin_ID`: Loma Ridge MAG identifier.
- `time`: litter-decomposition sampling time (`T1`-`T4`).
- `phylum`: GTDB phylum assignment.
- `family`: GTDB family assignment.
- `uptake_to_depolymerisation_ratio`: binary repertoire ratio.
