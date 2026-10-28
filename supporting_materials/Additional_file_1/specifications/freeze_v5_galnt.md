# Freeze v5 - GALNT / mucin-type O-glycosylation becomes the primary system (2026-10-12)

## Decision and disclosure
- Selected by the user after the six-system comparison (freeze v4). It satisfied the pre-written priority rule (1): global BH<0.05 in >=1 digestive cancer (PAAD: p=8.3e-4, BH_global=0.014, HR 2.25, n=177). No other system satisfied it.
- This is still a post hoc system choice (made after seeing results of 6 systems). All six systems are reported (ALL_systems_cluster_survival.csv). Neddylation and SUMO are parked, not deleted (NEDDYLATION_FROZEN.md).
- Focus digestive cancer: PAAD (v4 rule). The earlier n>=200 eligibility line was a v2 (Neddylation) criterion; PAAD has 177 tumors (<200). Disclosed, not hidden.

## Gene set (unchanged from freeze v4): 22 genes
- Initiation (GalNAc-T family, 20): GALNT1-GALNT18, GALNTL5, GALNTL6
- Extension (2): C1GALT1, C1GALT1C1 (Cosmc)
- No eraser exists for secretory-pathway mucin-type O-GalNAc (see chat 2026-10-12); no writer/eraser grouping. Functional grouping = initiation vs extension. Capping enzymes (ST6GALNAC family) not added (would be a post hoc expansion).
- Per-cancer expression filter (median log2(norm+1) < 1 dropped; >=6 genes) as in v4. Aliases: GALNT15=GALNTL2, GALNT16=GALNTL1, GALNT17=WBSCR17, GALNT18=GALNTL4.

## Figure 3 (GALNT) display rule
- k-means consensus clustering results from 04_results/2026-10-12_Figure3_multisystem_run01/GALNT.
- Panels A-C: cancers with BH_within<0.05 for GALNT (7: LGG, KIRP, KIRC, PAAD, BLCA, UCEC, HNSC) + the two largest digestive cancers not otherwise included (COAD n=432, STAD n=389), flagged if not significant. LIHC/CHOL/ESCA/READ in tables.
- Panel D Venn: PAAD + 3 other cancers with smallest p.
- Confounder checks (supplement, same cancers): Cox with cluster + age + sex + stage; PAAD additionally acinar-contamination score (PRSS1, CPA1, CELA3A, CTRB1, PNLIP, CEL) and stromal score (COL1A1, COL1A2, COL3A1, DCN, LUM, FAP); LGG additionally IDH1/2 mutation (MC3) and grade. Expression scores are proxies for purity (no ABSOLUTE/ESTIMATE purity downloaded).
- Results of confounder checks must be reported even if they remove the PAAD/LGG signal.

## Figure 1 and 2 (GALNT) - scripts and outputs (2026-10-12)
- Fig 1: 03_analysis/16_figure1_galnt_data.R (tables), 17_figure1_galnt_publication.R; 04_results/2026-10-12_Figure1_GALNT_run01, 05_figures/2026-10-12_Figure1_GALNT_run01
  - Raw added: 01_raw/uniprot/hpa_subcellular_galnt22.tsv (HPA; GALNT1/13/L5/L6 have no annotation), 01_raw/toil/toil_galnt22_log2norm.csv (GTEx/TOIL)
  - Same thresholds as figure1_replication_spec_v1.md; extra rule: gene x cancer with median log2(norm+1) < 1 is not tested (v4). Aliases: GALNT15-18 = GALNTL2/L1/WBSCR17/L4.
- Fig 2: 18_figure2_galnt_scores.R (reuses ssGSEA scores in 02_intermediate/figure2), 19_figure2_galnt_publication.R; 04_results/2026-10-12_Figure2_GALNT_run01, 05_figures/2026-10-12_Figure2_GALNT_run01. Focus cancer = PAAD.
- Caveats: PAAD tumor-vs-normal uses GTEx pancreas (TCGA PAAD has 4 normals); long genes (GALNT17/13) have passenger-like mutation counts; CPTAC mRNA-protein pairs exist for only 12/22 genes; purity not controlled in Fig 2.

## Figure 4 status (2026-10-12)
GIScore (Uni+CoxBoost; GALNT16 and GALNT14) built in TCGA-PAAD ductal n = 146 and validated in 3 GEO cohorts; the pooled criterion was met. Details and caveats in figure4_spec_v1.md (addendum).

## Figure 1 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 1 (panels A-F, 7.2 x 9.5 in; 05_figures/2026-10-12_Figure1_GALNT_run01/Figure1.{png,svg,tiff}) is frozen after the layout fix (title/tag theme, height 9.5 in). Source: 03_analysis/17_figure1_galnt_publication.R (read-only; copy Figure1_FROZEN_source.R in the figure folder); data script 16_figure1_galnt_data.R; SHA-256 of figure files, scripts and 04_results/2026-10-12_Figure1_GALNT_run01 in FROZEN_SHA256.txt. Any further change needs a dated amendment.

### Figure 1 amendment 1 (2026-10-12, user request; cosmetic)
Panel C tag overlapped the heatmap corner; added plot.margin = margin(13, 2, 2, 7) to pC in 17_figure1_galnt_publication.R. No data, statistics or other panels changed (result CSVs verified identical to the freeze checksums; seed 666 unchanged). Figure re-frozen; new SHA-256 in FROZEN_SHA256.txt. Previous version: 05_figures/2026-10-12_Figure1_GALNT_run01/archive_before_C_tag_fix/.

## Figure 2 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 2 (panels A-F, 7.2 x 9.4 in; 05_figures/2026-10-12_Figure2_GALNT_run01/Figure2.{png,svg,tiff}) is frozen after the layout fix and the redesign of panel A (four-card schematic). Source: 03_analysis/19_figure2_galnt_publication.R (read-only; copy Figure2_FROZEN_source.R); score script 18_figure2_galnt_scores.R; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure2_GALNT_run01 in FROZEN_SHA256.txt. Legend must state: Hallmark only is shown (KEGG in supplement), PAAD chosen after a digestive-system restriction (post hoc, disclosed; was KIRC under the pan-cancer ranking), panel A icons are schematic. Fig 1 was re-frozen after amendment 1 (panel C tag margin). Any further change needs a dated amendment.

## Figure 3 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 3 (7.2 x 9.5 in; 05_figures/2026-10-12_Figure3_GALNT_run01/Figure3.{png,svg,tiff}) is frozen after the layout fix (height 9.5 in; panel C title "KEGG enrichment (C2 vs C1 DEGs)"). Source: 03_analysis/15_figure3_galnt_publication.R (read-only; copy Figure3_FROZEN_source.R); upstream 13_figure3_multisystem.R and 14_galnt_confounders.R; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure3_multisystem_run01/GALNT in FROZEN_SHA256.txt. Any further change needs a dated amendment.

## Figure 4 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 4 (7.2 x 9.7 in; 05_figures/2026-10-12_Figure4_GALNT_run01/Figure4.{png,svg,tiff}) and the accompanying FigureS4_validation_multivariable are frozen after the layout fix (K labels above thin bars, M annotation, row heights). Source: 03_analysis/23_figure4_publication.R (read-only; copy Figure4_FROZEN_source.R); upstream 20_figure4_train.R, 21_figure4_validate.R, 22_figure4_downstream.R; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure4_GALNT_run01 in FROZEN_SHA256.txt. Any further change needs a dated amendment.

## Figure 5 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 5 (panels A-K, 7.2 x 10.0 in; 05_figures/2026-10-12_Figure5_GALNT_run01/Figure5.{png,svg,tiff}) and Supplementary Fig S5 (A three-cohort GALNT dot plots, B treatment/site, C inferCNV GSE202051, D per-patient proportions; 7.2 x 7.1 in; FigureS5_GALNT_cohorts_treatment.{png,svg}) are frozen. Source: 03_analysis/30_figure5_publication.R (read-only; copy Figure5_FROZEN_source.R); upstream 25_figure5_extract.py, 26_figure5_analysis.py, 27_figure5_tcga_composition.R, 29_figure5_template_data.py, 63_figure5_fig1B_galnt_dots.py, 64_figure5_malignant_states.py; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure5_GALNT_run01 in FROZEN_SHA256.txt. Superseded layouts: 03_analysis/30_figure5_publication_before_*.R and 28_figure5_publication_v1_superseded.R (not part of the freeze).
Data provenance (record in Methods):
- A, B: RefA object runs/2026-03-19_refA_tme (ann_fig1B, 18 groups plotted; 100,000 of 200,731 cells shown, seed 666).
- C-E: clean malignant epithelium runs/2026-05-10_Epithelial_scvi/artifacts/adata_malignant_r5_clean.h5ad (41,073 cells, six states); CTGF, CYR61 and GALNTL5 absent from the object.
- F, S5A: GALNT dot plots from 25/26 extraction; G (paired GIScore), H (subtype), I (composition), S5B: GSE202051/CRA001160 and the earlier PanIN-object GSE205013 extraction (standardized/GSE205013.h5ad). The GSE205013 GIScore values in G-I and S5A/B were NOT rebased onto the RefA / clean-epithelium objects; state this in the Methods and re-check if the GSE205013 object is ever replaced.
- S5D: per-patient proportions from the RefA TME_umap_table (ann_fig1B, PDAC patients with >=200 cells).
- J: TCGA-PAAD compartment correlations; K: adjusted Cox (GIScore per SD, OS).
Unfreeze only for user-requested amendments; verify result CSVs against FROZEN_SHA256.txt, then re-freeze.

## Figure 6 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 6 (panels A-H, 7.2 x 9.9 in; 05_figures/2026-10-12_Figure6_GALNT_run01/Figure6.{png,svg,tiff}) and Supplementary Fig S6 (spot-level GALNT16 / GALNT14 / GIScore by manual grade; FigureS6_spot_level_by_grade.{png,svg}) are frozen. Source: 03_analysis/33_figure6_publication.R (read-only; copy Figure6_FROZEN_source.R); upstream 31_figure6_visium_prep.py, 32_figure6_analysis.R; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure6_GALNT_run01 in FROZEN_SHA256.txt.
Content: A-C GSE254829 Visium slides HG2 and HG4 (chosen by a rule fixed before analysis); D all-gene grade slopes (exploratory); E six enzymes (post hoc list, exploratory); F Elhossiny epithelial pseudobulk trend (GIScore p = 0.90); G GL vs PD (p = 0.35); H depth caveat.
Interpretation to state in Results/Discussion: the frozen GIScore does not change along the lesion axis (null in both datasets; GALNT14 slope is opposite to the hypothesis, p = 8.8e-4, S6B). F and G are low-information because GALNT16 counts track pseudobulk library size (Spearman rho 0.95, n = 20; library sizes differ ~200-fold; this rho is computed from H2_epithelial_pseudobulk_scores.csv and is not printed on the figure). The user did not request adding it to panel H; if desired, put it in the legend.

## Figure 7 (GALNT) FREEZE (2026-10-12, user-approved)
Figure 7 (panels A-F, 7.2 x 8.4 in; 05_figures/2026-10-12_Figure7_GALNT_run01/Figure7.{png,svg,tiff}) is frozen after two cosmetic fixes: rows wrapped in wrap_elements (flush to page width) and a single tag size (10 pt bold) for all panels (panel B, theme_void, had fallen back to the default tag size). Source: 03_analysis/41_figure7_publication.R (read-only; copy Figure7_FROZEN_source.R); upstream 34-40_figure7_*; SHA-256 of figures, scripts and 04_results/2026-10-12_Figure7_GALNT_run01 in FROZEN_SHA256.txt. The superseded GIScore-spatial version (Figure7_v1_GIScore_spatial_superseded.*) is archived in the same folder and is a candidate for Supplementary Fig S7 (its hash is recorded).
Items that must be in the legend (not drawn on the figure): the dashed line in A separates the post hoc epithelial-detectable set (GALNT1, 2, 3, 4, 7, 10, 11, 12) from the rest; dagger in D-E = same sign in all 7 PDAC sections; C p-values are section-level exact tests (GL vs normal Mann-Whitney, GL vs PD paired signed-rank; smallest possible paired p with n = 7 is 0.016), uncorrected for 8 genes x 2 contrasts; the gene list was defined post hoc (A5b) after the frozen criterion kept only GALNT1/2/11; GALNT14/16 detection is 5-13%, so GIScore cannot be evaluated in space; "PanIN/GL tumor" is a mixed domain label. Panel F title is posed as a question ("same direction?"); GALNT12 is ~0 in GSE254829.

## Figure 4 sensitivity note (2026-10-13, user question: "can the univariate prefilter be set to p < 0.2 to keep more genes?")
The frozen GIScore (prefilter p < 0.05, 2 genes) is NOT changed. Because the threshold was fixed before validation (figure4_spec_v1.md), picking a new threshold after seeing external results would be post hoc; instead a grid declared in advance (0.05, 0.10, 0.20, 0.30) was run as a sensitivity analysis: 03_analysis/65_figure4_uni_threshold_sensitivity.R -> 04_results/2026-10-13_Figure4_uni_threshold_sensitivity/ (new folder; frozen Fig 4 files untouched). Training-only 10x5 CV C-index (prefilter inside each fold): 0.602 / 0.596 / 0.594 / 0.588. Genes kept: 0.05 -> GALNT16, GALNT14; 0.10 -> + C1GALT1C1; 0.20 and 0.30 -> GALNT10, GALNT14, GALNT16, C1GALT1C1 (GALNT17 passes the prefilter but is dropped by CoxBoost). External pooled HR (3 GEO cohorts, random effects): 1.65 / 1.56 / 1.75 / 1.69, all three cohorts same direction at every threshold; GSE85916 stays non-significant (C-index 0.52-0.53). Conclusion: a looser prefilter adds genes but does not improve CV or external performance; the result is not sensitive to the threshold. Suggested use: Supplementary table (not a replacement of GIScore). Small differences from the frozen validation numbers at 0.05 (e.g. GSE57495 HR 2.06 vs 1.88) come from the random folds inside cv.CoxBoost when refitting.
