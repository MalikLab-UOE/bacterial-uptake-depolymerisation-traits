# Validation record

Validation was updated on 20 October 2026.

- Binary repertoire status was verified for the GlobDB, IMG and Loma analysis inputs. Ratios reproduce exactly from the corresponding uptake and depolymerisation repertoire totals.
- GlobDB binary input contains 249,311 MAGs; the Figure 3 tree match retains 249,310 tips.
- Figure 1 uses the binary ratio. Figure 1b is detection-count weighted.
- Figure 2 and the bulk-soil/rhizosphere analysis use the IMG binary repertoire ratio.
- The bulk-soil/rhizosphere plotting dataset excludes 29 MAGs with conflicting project and geographic-location site labels.
- Figure 3's tree-visualisation script now reads `glo_input_binary.rds`, consistent with the final figure and Pagel's lambda analysis.
- Figure 4 uses 248,230 MAGs after filtering.
- Figure 5 uses binary repertoire totals for the MAG-level panels.
- The numbered supplementary set contains S1 rhizosphere, S2 phylum traits and S3 phylum GAMs. The family-level Loma plot is panel b of Figure 5.
- Figure S2 displays genome size in Mb from 0 to 10 Mb and uses the compact single-row layout with an inset colour key.
- Figure S3 contains 16 panels arranged in four rows. It includes the complete dataset and the seven most abundant phyla; each dataset has one GAM panel and one dataset-specific bottom-versus-top-decile panel.
- Data S1 contains 41 unique GH-family entries and 211 unique transporter KO entries. The two missing source descriptions are explicitly flagged.
- All R scripts passed a parser syntax check after renaming and path correction.
- `Figure_2.R`, `Supplementary_Figure_S2.R` and `Figure_5_family_panel.R` were rerun from the archive on 20 October 2026. Their regenerated outputs match the final compact layouts and reported statistics.
- `Figure_4.R` was rerun to replace a truncated PDF copy; the regenerated one-page PDF passed strict parsing and visual inspection without changing the analysis.
- The final set of eight numbered PDFs (Figures 1–5 and S1–S3) was reopened successfully after the updates; the revised Figures 2, 4, 5, S1 and S2 were also rendered and showed no clipping, overlap or missing labels.

`docs/analysis_audit.md` records only the analyses represented in Figures 1–5, Figures S1–S3 and Data S1.
