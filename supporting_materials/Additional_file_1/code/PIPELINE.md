# Primary pipeline (ignore backups)

Run order is numbered. Publication scripts are frozen copies; live copies live in `03_analysis/` and should match `*_FROZEN_source.*` in the figure folders.

## Figure 1 — pan-cancer survey

| Step | Script | Output |
|---|---|---|
| Data | `03_analysis/16_figure1_galnt_data.R` | `04_results/2026-10-12_Figure1_GALNT_run01/` |
| Plot | `03_analysis/17_figure1_galnt_publication.R` | `05_figures/2026-10-12_Figure1_GALNT_run01/` |

## Figure 2 — Hallmark pathways

| Step | Script | Output |
|---|---|---|
| Scores | `03_analysis/18_figure2_galnt_scores.R` (reuses `02_intermediate/figure2`) | `04_results/2026-10-12_Figure2_GALNT_run01/` |
| Plot | `03_analysis/19_figure2_galnt_publication.R` | `05_figures/2026-10-12_Figure2_GALNT_run01/` |

## Figure 3 — consensus subtypes

| Step | Script | Output |
|---|---|---|
| Six-system clustering | `03_analysis/13_figure3_multisystem.R` | `04_results/2026-10-12_Figure3_multisystem_run01/` |
| Confounders | `03_analysis/14_galnt_confounders.R` | same `/GALNT/` |
| Plot | `03_analysis/15_figure3_galnt_publication.R` | `05_figures/2026-10-12_Figure3_GALNT_run01/` |

## Figure 4 — GIScore (locked)

| Step | Script | Output |
|---|---|---|
| Model library | `03_analysis/20a_models_lib.R` | — |
| Train + hash | `03_analysis/20_figure4_train.R` | `FROZEN_model.rds`, `FROZEN_model_sha256.txt` |
| Validate | `03_analysis/21_figure4_validate.R` | `validation_*.csv` |
| Downstream | `03_analysis/22_figure4_downstream.R` | clinical / Reactome / immune proxies |
| Plot | `03_analysis/23_figure4_publication.R` | Figure 4 + S4 |
| Supp tables | `03_analysis/24_figure4_supp.R` | — |
| Sensitivity (post-freeze) | `03_analysis/65_figure4_uni_threshold_sensitivity.R` | `04_results/2026-10-13_Figure4_uni_threshold_sensitivity/` |

Lock check: `21` must verify `FROZEN_model_sha256.txt` before reading GEO.

GIScore formula in the frozen model:

```
−0.344 · z(GALNT16) + 0.202 · z(GALNT14)
```

GALNT14 entered via a top-2 fallback of the univariate prefilter (univariate P = 0.0506). Do not treat this as a p < 0.05 gene unless the manuscript states the fallback.

## Figure 5 — single-cell

| Step | Script | Output |
|---|---|---|
| Extract | `03_analysis/25_figure5_extract.py` | `02_intermediate/fig5_sc/` |
| Analysis | `03_analysis/26_figure5_analysis.py` | `04_results/2026-10-12_Figure5_GALNT_run01/` |
| TCGA composition | `03_analysis/27_figure5_tcga_composition.R` | `H_tcga_composition_*.csv` |
| Template data | `03_analysis/29_figure5_template_data.py` | — |
| Fig1B dots | `03_analysis/63_figure5_fig1B_galnt_dots.py` | `B_fig1B_GALNT14_16_dots.csv` |
| Malignant states | `03_analysis/64_figure5_malignant_states.py` | C–E tables |
| Plot | `03_analysis/30_figure5_publication.R` | Figure 5 + S5 |

**Provenance (must stay in Methods somewhere):**

- A, B, S5D: RefA TME (`ann_fig1B`; 200,731 cells; 100,000 shown).
- C–E: `adata_malignant_r5_clean.h5ad` (41,073 cells, 6 states). CTGF, CYR61, GALNTL5 absent.
- F, S5A: extraction from 25/26.
- G–I and S5B GSE205013 GIScore: **earlier PanIN-object extraction, not rebased** onto RefA / clean-epithelium.

## Figure 6 — GSE254829 visium + Elhossiny GIScore

| Step | Script | Output |
|---|---|---|
| Prep | `03_analysis/31_figure6_visium_prep.py` | `02_intermediate/fig6_visium/` |
| Analysis | `03_analysis/32_figure6_analysis.R` | `04_results/2026-10-12_Figure6_GALNT_run01/` |
| Plot | `03_analysis/33_figure6_publication.R` | Figure 6 + S6 |

## Figure 7 — visium domains / co-localisation

| Step | Script | Output |
|---|---|---|
| Signature | `03_analysis/34_figure7_signature.py` | — |
| Spatial | `03_analysis/35_figure7_spatial_analysis.R` | — |
| Co-localisation | `03_analysis/36_figure7_colocalization.R` | — |
| Epithelial niche | `03_analysis/38_figure7_epithelial_niche.R` | — |
| Detectability | `03_analysis/39_figure7_visium_detectable.R` | `H7f_detectability.csv` |
| Export | `03_analysis/40_figure7_export.R` | — |
| Plot | `03_analysis/41_figure7_publication.R` | Figure 7 |

Skip `37_figure7_publication_v1_GIScore_superseded.R`.

The epithelial-detectable gene list (GALNT1/2/3/4/7/10/11/12) and the PanIN-associated set (GALNT3/4/5/6/10, C1GALT1) were defined after looking at detection/slopes. Specs call this post hoc / exploratory.

## Figure 8 — UNI2-h H&E

| Step | Script | Output |
|---|---|---|
| Cluster | `03_analysis/42_figure8_cluster.py` | — |
| Stats | `03_analysis/43_figure8_stats.py` | `H8a_cluster_vs_GIScore.csv` |
| MIL | `03_analysis/44_figure8_mil.py` | `H8b_*.json` |
| Slide images | `03_analysis/45_figure8_slide_images.py` | — |
| Assign tiles | `03_analysis/46_figure8_assign_all_tiles.py` | — |
| CPTAC transfer | `03_analysis/48_figure8_cptac_external.py` | `H8e_CPTAC_transfer.json` |
| Attention | `03_analysis/49_figure8_attention.py` | — |
| Zoom crops | `03_analysis/50_figure8_zoom_crops.py` | — |
| Plot | `03_analysis/51_figure8_publication_v2.py` | Figure 8 |
| Regression / ridge / ink QC | `52`–`62` as needed | `04_results/2026-10-12_Figure8_*` |
| S8 | same publication script / `05_figures/2026-10-12_FigureS8_pathology_run01/` | Figure S8 |

Skip `47` (older publication), `54_*_REJECTED.py`.

## Seeds

Default seed 666 unless a panel specifies otherwise (Fig 8 tiles: 20261012).
