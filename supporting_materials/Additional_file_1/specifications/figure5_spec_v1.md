# Figure 5 spec v1 (frozen before any GALNT expression is read from the single-cell cohorts)

Date: 2026-10-12. Purpose: place the GALNT system and the frozen GIScore (GALNT16 -0.344, GALNT14 +0.202, z-scored) at single-cell resolution in PDAC. This is an exploratory/mechanistic-context figure, not a validation of Fig 4.

## Data (read-only, ${PANIN_ROOT}/02_intermediate/2026-09-16_six_cohort_input_bundle/standardized)
| Cohort | Tech | Cells | Used for |
|---|---|---|---|
| GSE202051 (Hwang 2022) | snRNA | 224,988 | lineages, malignant vs non-malignant epithelium, basal/classical, chemo vs naive |
| GSE205013 (Chen 2022) | scRNA | 251,080 | lineages, tumor states (Basal_Squamous vs Ductal etc.), primary vs liver metastasis, immune/CAF composition |
| CRA001160 (Peng 2019) | scRNA | 57,530 | lineages, normal ductal vs tumor ductal |
Raw counts (`layers["counts"]`) only. Per cell: log1p(counts / library size x 1e4), library size = sum over ALL genes.
Not used here: GSE242230, GSE291124, GSE253429 (lineage not annotated), spatial/Xenium/Visium HD (reserved for a later PanIN spatial figure).

## Unit of analysis
Patient (not cell). Patient x lineage pseudobulk = mean of per-cell log-normalized values; require >= 20 cells in the unit. Cells are only used for dot plots. All p-values use patient-level units.

## Gene sets
- GALNT system (22): GALNT1-18, GALNTL5, GALNTL6, C1GALT1, C1GALT1C1 (aliases GALNTL2/L1/WBSCR17/L4 mapped to GALNT15-18).
- GIScore (frozen, Fig 4): -0.344*z(GALNT16) + 0.202*z(GALNT14); z computed across patient pseudobulks of the SAME lineage within the cohort (the frozen TCGA center/scale cannot be transferred across technologies).
- Moffitt basal-like / classical gene lists (originally written from memory; CORRECTED 2026-10-12 against WikiPathways WP5390 = Moffitt 2015 Nat Genet, 25+25 genes; the lists are in `fig5_genes.json`; GSE202051 also carries the authors' precomputed Moffitt_basal/Moffitt_classical, used as the primary there). Subtype score = mean z(basal genes) - mean z(classical genes), z across malignant patient pseudobulks.
- Broad signatures for composition: T/NK, myeloid, CAF, endothelial fractions among non-epithelial cells (annotation-derived, not expression-derived).

## Pre-specified analyses (direction fixed from TCGA, Figs 3-4)
A. Lineage specificity: dot plot (mean expression, fraction expressing) of the 22 genes x lineage in all 3 cohorts. Descriptive.
B. Malignant vs non-malignant epithelium (GSE202051 malignant vs non-malignant epithelial; CRA001160 tumor ductal vs normal ductal; GSE205013 has no normal reference, skipped). Patient-paired Wilcoxon signed-rank when both compartments exist per patient, otherwise Mann-Whitney. Genes: GALNT16, GALNT14, GIScore. No direction pre-specified (TCGA/GTEx contrast had a batch caveat).
C. GIScore vs subtype (malignant cells, patient level): Spearman between GIScore and (basal - classical) score. Hypothesis: positive (TCGA high-score tumors are enriched in keratinization). One-sided in direction, report two-sided p. Replicate in each cohort that has >= 8 eligible patients.
D. Tumor states (GSE205013, ann_L3 within epithelium: Ductal, Secretory_Ductal, Cycling_Tumor, Basal_Squamous, Boundary): GIScore by state, patient x state pseudobulk, paired by patient, Friedman test then Basal_Squamous vs Ductal paired Wilcoxon.
E. Microenvironment composition: Spearman between malignant GIScore (patient) and the T/NK fraction and the myeloid fraction among non-epithelial cells (primary pancreas only in GSE205013). Hypothesis (from TCGA immune-cold finding): negative for T/NK. BH across the 2 compositions x cohorts.
F. Treatment: GIScore in malignant cells, naive vs chemo (Mann-Whitney, unpaired patients), GSE202051 and GSE205013 (primary pancreas). Exploratory.
G. Primary vs liver metastasis (GSE205013): GIScore and GALNT16/GALNT14 in tumor epithelium, Mann-Whitney. Exploratory.

## Multiplicity and reporting
BH within each panel family. Everything reported, including null results. No cell-level p-values.

## Known limitations to disclose
Cohorts differ by technology (snRNA vs scRNA); patients per cohort are modest; malignant/non-malignant calls are the authors'; GSE205013 contains chemo-treated and metastatic samples; GIScore z-scoring is within-cohort; dropout makes low-expressed GALNT16 sparse, so patient pseudobulks (not cells) are used.

---
## Addendum (2026-10-12, after execution)
Scripts: 25 (extract), 26 (analysis), 27 (post hoc TCGA composition), 28 (figure). Results: 04_results/2026-10-12_Figure5_GALNT_run01. Figure: 05_figures/2026-10-12_Figure5_GALNT_run01/Figure5.{svg,png,tiff}.
Outcomes against the pre-specified hypotheses:
- A: GALNT16 is detected in <1% of malignant cells in all 3 cohorts (0.1-0.7%); highest in CAF/fibroblast (3-6%) and endocrine (1-11%). GALNT14 6-11% in malignant cells. GIScore at single-cell level is therefore driven mostly by GALNT14 dropout-limited signal.
- B: GALNT16 is LOWER in malignant than non-malignant epithelium (paired Wilcoxon p = 0.0019 GSE202051, n=23; p = 0.0024 CRA001160, n=14); GALNT14 unchanged; GIScore higher in malignant (p = 0.025, 8.5e-4), driven by GALNT16 loss.
- C: GIScore vs basal-classical: rho 0.38 (p 0.020, n=38) in GSE202051 only; 0.16 (p 0.44) GSE205013; 0.05 (p 0.81) CRA001160. NOT replicated (hypothesis supported in 1/3, nominal only; BH across 9 tests 0.18).
- D: In GSE205013 GIScore differs across states (Friedman p 4.4e-4, n=14) but Basal_Squamous is LOWER than Ductal (p 0.0012), i.e. opposite to the TCGA keratinization hypothesis.
- E: T/NK fraction vs GIScore: rho -0.02 (GSE202051, lymphoid), +0.29, +0.21: the TCGA immune-cold finding is NOT reproduced in single-cell composition (all BH > 0.1).
- F/G: treatment and primary-vs-liver comparisons all ns.
- Post hoc (27): in TCGA-PAAD GIScore correlates with bulk composition (immune rho -0.61, neural -0.56, endocrine -0.40, acinar -0.33, tumor-epithelial +0.64; GALNT16 itself +0.65 with immune, +0.67 with neural). The GIScore HR survives adjustment (1.50 -> 1.64-1.66; non-tumor PC1 alone p = 0.11), so GIScore keeps prognostic information beyond these proxies, but GALNT16 behaves as a non-tumor-compartment marker in bulk. The tumor-intrinsic enzyme interpretation is NOT supported by the single-cell data and the paper should be framed accordingly.
- Gene lists for Moffitt subtype (GSE205013, CRA001160) were originally written from memory and were partly wrong (8 of 25 classical genes incorrect, e.g. S100P/GATA6/CLDN18 were not in the list; 2 basal symbols outdated). VERIFIED/CORRECTED against WikiPathways WP5390 (Moffitt 2015): basal CTSL2->CTSV, ANXA8L2->ANXA8L1; classical now incl. PRR15L, CEACAM6, CYP3A7, MYO1A, KRT20, CDH17, SPINK4, LOC400573. Present: GSE205013 25/24, CRA001160 24/24 (LOC400573 absent from all matrices). Re-run: GIScore vs subtype rho 0.16/0.07 (was 0.155/0.051), still non-significant; conclusions unchanged. Only the supplier (WikiPathways) was checked, not the paper supplement table.

## Layout revision (2026-10-12)
Fig 5 was re-laid out to follow the template Fig 5 (Fig5_template.png, downloaded from the paper): A UMAP of main cell types, B marker dot plot (+GALNT14/GALNT16), C per-patient proportions, D inferCNV score (authors' score, GSE202051), E malignant-state UMAP, F state x GALNT dot plot, then our results G-K (paired malignant vs non-malignant, subtype, TME composition, TCGA composition, adjusted HR). Not reproduced from the template: CellChat panels, myeloid sub-clustering and ICI-response panels (not run / no data). UMAPs show 80,000 / 40,000 randomly sampled cells (seed 666). The 3-cohort GALNT dot plots and treatment/site contrasts moved to FigureS5. Scripts: 29 (export), 30 (figure); 28 superseded.

### Layout amendment (2026-10-12, user request)
Former Fig 5C (per-patient cell proportions, GSE205013) moved to Supplementary Fig S5C (horizontal stacked bars, 27 patients). Main Fig 5 re-lettered: A UMAP, B markers, C inferCNV, D malignant states, E GALNT by state, F paired malignant vs non-malignant, G subtype, H TME composition, I TCGA composition, J adjusted HR (7.2 x 8.55 in). Other changes: G/H annotation moved to subtitle (no overlap with points), E legend on the right, UMAP1 axis labels, larger C-axis fonts removed with panel, stars in H enlarged, S5 dot-plots with collected legend. No data or statistics changed. Backup: 03_analysis/30_figure5_publication_before_v2_review.R.
Layout amendment 2 (2026-10-12, user: distribution odd): UMAP axis titles replaced by corner arrows (no floating "UMAP1" labels); B legends moved to the right so A and B share top/bottom edges; panel titles shortened (C "inferCNV (GSE202051)", H "Composition vs GIScore", I "TCGA-PAAD bulk") to stop titles running into neighbouring tags; F strips shortened (units in y label); row widths re-balanced. Backup: 03_analysis/30_figure5_publication_before_v3_layout.R. No data changed.

### Layout amendment 3 (2026-10-12, user-specified source): Fig 5A, B, S5C redrawn on the GSE205013 Figure 1B annotation
Source (read directly, nothing copied): ${PDAC_CAF_CHEMO_ROOT}/runs/2026-03-19_refA_tme (object artifacts/RefA_TME_RefClean_Fig1B.h5ad; tables/TME_umap_table.tsv.gz with `ann_fig1B`; figures/ad4_dotplot_ann_fig1B.tsv). 19 groups (Background excluded; 200,731 cells), marker_panels_fig1B (3 markers per group, 57 genes), expression layer log1p_counts, dot colour = min-max of group means per gene, dot size = fraction expressing (as in the source dotplot). Group colours were sampled from the source Fig1A UMAP. Fig 5A: 100,000 randomly sampled cells (seed 666). Fig 5B additionally shows GALNT14 and GALNT16 (03_analysis/63_figure5_fig1B_galnt_dots.py -> 04_results/.../B_fig1B_GALNT14_16_dots.csv; method verified by reproducing EPCAM from the source table, max abs diff 6e-8). S5C patient proportions now use ann_fig1B (PDAC patients with >= 200 cells).
Earlier version of A/B/S5C used the PanIN standardized GSE205013.h5ad (ann_L2, 8 groups, no mast markers); it is superseded. NOT yet rebased on the RefA object: Fig 5D/E (malignant states, ann_L3 of the PanIN h5ad), 5F-5H GSE205013 pseudobulk and composition numbers; consistency check pending user decision.

### Layout amendment 4 (2026-10-12, user-specified source): Fig 5C-E rebased on the GSE205013 clean malignant epithelium
Source: ${PDAC_CAF_CHEMO_ROOT}/runs/2026-05-10_Epithelial_scvi/artifacts/adata_malignant_r5_clean.h5ad (41,073 cells; 6 states = ann_L2/ann_celltype: Gastric_classical, Classical_mucinous, Hypoxic_inflammatory_intermediate, Inflammatory_basal, Wound_squamous, Basal_plastic; X_umap, layer log1p_counts). Extraction: 03_analysis/64_figure5_malignant_states.py -> 04_results/2026-10-12_Figure5_GALNT_run01/{C_malignant_clean_umap.csv.gz (UMAP coordinates exported from the h5ad), D_malignant_clean_marker_dots.csv, E_malignant_clean_GALNT_dots.csv, CDE_malignant_clean_meta.json}.
C: UMAP by state. D: marker dotplot with the panels of notebooks/Epi/Tumor_tme_scvi.ipynb (cell 60; scanpy dotplot, layer log1p_counts, standard_scale='var', dot_max 0.6; CTGF and CYR61 are absent from the object and were absent in the source figure too); state order from sc.tl.dendrogram(groupby='ann_celltype') on X_pca. E: 21 of 22 GALNT-system genes (GALNTL5 absent from the object; GALNT15/16/17/18 via aliases GALNTL2/GALNTL1/WBSCR17/GALNTL4) per state, dot size = fraction expressing (cap 60%), colour = z of group mean (-1..2).
Moved: inferCNV score (GSE202051) -> Supplementary Fig S5C; per-patient proportions -> S5D. Fig 5 re-lettered: A global UMAP, B global marker dotplot, C malignant states, D state markers, E GALNT by state, F paired malignant vs non-malignant, G subtype, H composition, I TCGA composition, J adjusted HR (7.2 x 9.75 in).
Still on the PanIN GSE205013 extraction (not rebased): the GSE205013 GIScore pseudobulk in 5G/5H (malignant GIScore per patient), S5A dot plot of GSE205013, S5B treatment/site. Backup script: 03_analysis/30_figure5_publication_before_CDE_rebase.R.

### Layout amendment 5 (2026-10-12): row gutters
Cause of the blank right margin (up to 1.3 in on rows 3-4): nested patchwork rows shared gutters, so the widest legend / text margin of any row (legends in rows 1-2, 58-pt label margin in row 5) shrank every other row. Fix: each row is wrapped in wrap_elements(full = row) so rows are flush to the page width (content now spans 52-2800 px of 2880). Titles of A/C indented 13 pt to clear the tag letters; F y-axis title set on two lines. Backup: 03_analysis/30_figure5_publication_before_rowwrap.R.

### Layout amendment 6 (2026-10-12, user request): panels re-organised
Fig 5 (7.2 x 10.0 in): A/C UMAPs drawn with coord_fixed (natural aspect, no stretching); F (new) = GSE202051 GALNT-system dot plot by lineage (same data/code as S5A first panel), next to E; G = paired malignant vs non-malignant pseudobulk, GIScore only (single-gene GALNT16/GALNT14 panels removed; their statistics remain in 04_results/.../B_malignant_vs_nonmalignant.csv); H subtype; I composition; J TCGA; K adjusted HR.
Supplementary Fig S5 (7.2 x 9.2 in): A dot plots of three cohorts; B/C earlier PanIN-object UMAPs (main cell types; malignant states Ductal/Secretory_Ductal/Boundary/Cycling_Tumor/Basal_Squamous), titled as such; D treatment/site; E inferCNV (GSE202051); F per-patient proportions. Backup: 03_analysis/30_figure5_publication_before_v6_reorg.R.
Caveat: the S5B/C annotation scheme differs from main-figure C-E (clean-epithelium six states).

### Layout amendment 7 (2026-10-12, user request)
S5 earlier-PanIN UMAPs (former S5B, C) removed; S5 = A dot plots (3 cohorts), B treatment/site, C inferCNV, D per-patient proportions (7.2 x 7.1 in). Main A/C UMAP points enlarged (0.08 -> 0.3 pt, C cells shuffled) and labels given a white halo, because the colours were too pale. Backup: 03_analysis/30_figure5_publication_before_v6_reorg.R (pre-amendment 6).
