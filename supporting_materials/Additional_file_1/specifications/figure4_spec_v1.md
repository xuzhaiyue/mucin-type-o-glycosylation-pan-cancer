# Figure 4 spec v1 (GALNT / PAAD prognostic score) - frozen BEFORE external expression data were opened (2026-10-12)

Template: Zhu et al. 2026 Fig. 4 (AIScore: univariate forest, 101 ML combinations by C-index, LASSO/RSF diagnostics, KM in training + 2 validation cohorts, clinical association donuts, Reactome enrichment, TIDE, MDSC).
Our analogue: GALNT imbalance score (GIScore) in PAAD.

## Cohorts
- Training: TCGA-PAAD, histological type "Pancreas-Adenocarcinoma Ductal Type" only (n=154; 153 with OS; 92 events). Rationale: external cohorts are PDAC; non-ductal PAAD (pNET etc.) differ biologically. (Fig. 3 used all 177 PAAD; ductal-only result reported as sensitivity.)
- Validation (opened only after the model is frozen and its hash saved): GSE62452 (tumors with OS; GPL6244), GSE57495 (GPL15048; n=63), GSE85916 (GPL13667; n=80). GSE28735 excluded (same patients as GSE62452). GSE71729 and GSE21501 not used (reference-based/clinical annotation incomplete in series matrix).
- Probe to gene: max-mean probe per gene; gene z-scored within each cohort; aliases GALNT15-18 = GALNTL2/L1/WBSCR17/L4.

## Candidate genes
22 GALNT genes (freeze v5) that (i) have median log2(norm+1) >= 1 in training tumors and (ii) are present on all three validation platforms (GALNT4 is absent on GPL15048 -> excluded). Decided without looking at validation expression.

## Model building (training data only)
- 16 pre-specified models, all producing a risk score: Lasso, Ridge, ElasticNet(0.5), CoxBoost, RSF, GBM, StepCox[both]; each also with a univariate prefilter (p<0.05 within the training fold); plus Lasso+RSF and Lasso+GBM.
- Model selection: highest mean C-index in 10 x repeated 5-fold CV inside TCGA training cohort. Ties within 0.005 -> simpler model. Validation cohorts are NOT used for selection (deviation from the template, which ranked by training + validation average C-index and thereby leaked validation information).
- The selected model is refit on all training samples; coefficients/model saved with sha256 before validation data are read.
- Score applied in each cohort to z-scored genes; high/low split at the cohort median (training: training median).

## Primary validation criteria (fixed now)
- Success = pooled (inverse-variance random-effects) validation HR(high vs low) > 1 with p < 0.05 AND the same direction in >= 2/3 cohorts. Anything else is reported as "not validated".
- All C-index values (also for all 16 models in validation cohorts, shown only for transparency) are reported.

## Panels
A univariate Cox forest (continuous z) in training; B C-index (training CV + 3 validation cohorts) for all 16 models; C model diagnostics (LASSO lambda CV, RSF error, variable importance, as applicable to the selected model); D-G KM in training and 3 validation cohorts; H-J clinical association donuts + chi-square (TCGA, GSE62452/85916/57495 as available); K Reactome enrichment of genes differentially expressed between score-high and low in TCGA-PAAD (signature has too few genes for ORA; deviation); L immune-evasion surrogate (TIDE web tool cannot be run offline: we export a TIDE-ready matrix and use signature proxies); M correlation of score with MDSC signature (ssGSEA proxy).
- Multivariable Cox (age, sex, stage, grade where available; + immune/stromal signature proxies) in every cohort.

---
## Addendum (2026-10-12, after execution)
- Actual training set: TCGA-PAAD ductal, n = 146 (85 deaths) after RNA matching (the n = 154/153 above is superseded). The sensitivity set, all TCGA-PAAD (n = 177), was run afterwards.
- Selected model (training 10x5 CV only): Uni+CoxBoost, CV C-index 0.602. The model has 2 genes: GALNT16 (beta -0.344) and GALNT14 (beta +0.202). GALNT14 entered through the top-2 fallback of the univariate prefilter (only GALNT16 had p < 0.05). The model was frozen (FROZEN_model.rds, sha256 file) BEFORE any external expression was read; the script verifies the hash before validation.
- Validation (scripts 21): GSE62452 n = 65, HR 2.01, log-rank p = 0.021, C = 0.606; GSE57495 n = 63, HR 1.88, p = 0.039, C = 0.617; GSE85916 n = 79, HR 1.19, p = 0.52, C = 0.527. Pooled random-effects HR 1.60 (1.14-2.24), p = 0.0065, I2 = 3%. The pre-specified success criterion was MET (pooled HR > 1 and p < 0.05, same direction in 3/3 cohorts; nominally significant in 2/3).
- Caveats to disclose: (1) GALNT16 alone drives the signal; (2) modest discrimination (C about 0.60); (3) GSE85916 not validated; (4) GSE62452/GSE57495 stage and grade are coarse and few; the score HR in GSE62452 after stage + grade adjustment is 1.28 (0.95-1.73, p = 0.099), underpowered; (5) the other 15 models were refit on the full training set and evaluated on validation cohorts for transparency only; they did not influence selection; (6) the immune signature panel uses expression proxies, TIDE itself was not run (TIDE-ready file exported); (7) low-score tumors are enriched for acinar/neuronal/immune genes, which may partly reflect tumor purity, not controlled in the DEG/Reactome panel.
- Scripts: 20, 20a (model library), 21, 22, 23. Results: 04_results/2026-10-12_Figure4_GALNT_run01. Figure: 05_figures/2026-10-12_Figure4_GALNT_run01/Figure4.{svg,png,tiff}.
