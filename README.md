# Mucin-type O-glycosylation pan-cancer analysis: code and source tables

Materials for the manuscript **“Integrative pan-cancer and PAAD-focused multi-omics analysis of mucin-type O-glycosylation enzyme transcriptional landscapes.”** The analysis covers 22 mucin-type O-glycosylation enzymes and a frozen two-gene PAAD index: GIScore = -0.344 × z(GALNT16) + 0.202 × z(GALNT14).

## Contents

- `supporting_materials/Additional_file_1/`: analysis scripts, dated specifications, model-lock records, data-source inventory, numerical tables, PAAD KEGG heatmap, and audit/provenance files.
- `supporting_materials/Additional_file_3/`: single-cell source tables, annotations, inferCNV values, and the patient-proportion source table.
- `supporting_materials/Additional_file_4/`: Visium spot-level source records, Figure 6 sensitivity reconstruction, and audit records.
- `supporting_materials/Additional_file_5/`: computational pathology predictions, AUC/permutation analyses, site-confounding analyses, and source tables.
- `archive/Additional_file_2_original.xlsx`: the original 17-worksheet GIScore model-results workbook, retained unchanged.
- `Supplementary_Tables.xlsx`: 21-worksheet journal file with the original 17 model-results worksheets and four curated, figure-linked source-summary worksheets (55 selected tables; 5,931 rows).
- `Supplementary_Figures_S1-S4.pdf`: the four cited supplementary figures in order (five pages total).
- `selected_table_manifest.csv`: the source files selected for the journal workbook.
- `SHA256SUMS.txt`: release-wide checksums (excludes itself and `repository_inventory.csv` to avoid a recursive manifest). `repository_inventory.csv` records byte size and SHA-256 for every other release file, including `SHA256SUMS.txt`. Each source package also has its own checksum manifest.

The two files intended for journal upload are `Supplementary_Figures_S1-S4.pdf` and `Supplementary_Tables.xlsx`. The workbook contains selected summary tables for readability; all 123 individual CSV/CSV.GZ source tables (620,339 rows), including the detailed cell- and spot-level records, remain in this repository.

## Public datasets

The study uses public resources identified in the manuscript, including TCGA, CPTAC, GEO accessions GSE62452, GSE57495, GSE85916, GSE202051, GSE205013, GSE254829, and CRA001160. Elhossiny et al.'s processed pancreatic single-cell and spatial objects are available from [Zenodo record 21867195, version 3](https://doi.org/10.5281/zenodo.21867195). UNI2-h embeddings are not redistributed; users must obtain access from the provider under its stated terms.

## Reproducibility notes

The archived analysis outputs and their specifications are provided as recorded; no whole-study environment lock was available. Archived scripts may require upstream datasets and local path configuration. See `PUBLIC_RELEASE_NOTES.md` and the package README/data dictionaries before reuse. The locked GIScore coefficients are provided for transparency; sensitivity refits do not replace the frozen model.

## Checksums and reuse

Use the root and package-level `SHA256SUMS.txt` files to verify released files. Public accessibility does not override upstream dataset/provider terms. No additional reuse license was supplied with this release.
