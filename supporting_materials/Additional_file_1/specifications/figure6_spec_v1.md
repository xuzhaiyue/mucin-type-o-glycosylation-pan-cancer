# Figure 6 spec v1 (FROZEN before any GALNT value was read from the spatial datasets)

Date: 2026-10-12. User chose: not to copy the template Fig 6 (Visium tumor/stroma/interface + SPP1 macrophage co-localisation) but to ask
"how does the GALNT system / GIScore change along the lesion axis normal -> PanIN -> PDAC, in space?".
(Correction: the template paper has 7 main figures + a workflow Fig 8; Fig 6 = spatial, Fig 7 = digital pathology. Fig6/7/8 templates saved in 01_raw/template_figs.)

## Data (what was inspected so far: metadata and counts of annotation classes only, no GALNT expression)
1. GSE254829 (Visium, human PanIN): 7 slides (PanIN-LG1-3, PanIN-HG1-4), 16,162 spots in codatable. Epithelial spots carry a manual histology
   label `paningrade` in {normal, low_grade, high_grade} (354 / 110 / 159 spots); all other spots = "not_manually_annotated_as_epithelial".
   Per-spot manual % of islets / normal epithelium / smooth muscle / fat / acini / collagen is also given. Spot = mixture of cells (55 um).
2. Elhossiny et al. 2026 `pseudobulk_epi.rds` (Visium PDAC, BayesPrism-deconvolved EPITHELIAL expression, 20 region pseudobulks):
   Donor_PanIN n=4 (healthy donors), AdjNormal_PanIN n=2, GLTumor (gland-forming) n=7, PDTumor (poorly differentiated) n=7 (GL/PD paired within patient;
   16 donors/patients). Deconvolution gives epithelial-specific expression, which removes the stromal/acinar mixing problem of raw spots.

## Frozen quantities
- Genes: GALNT16, GALNT14, GIScore = -0.344*z(GALNT16) + 0.202*z(GALNT14) (frozen Fig 4 model; z-score within the analysed unit set of each dataset).
- Exploratory: all 22 GALNT-system genes (BH within the family).
- Expression scale: log1p(counts / library x 1e4) for spots; for EpithelialExp use the stored values as is (log scale checked, not rescaled).

## Hypotheses (directions from TCGA: high GIScore = worse OS = low GALNT16 + high GALNT14)
- H1 (Visium PanIN, within-slide): GIScore increases normal -> low grade -> high grade epithelium. Test: slide-level mean per grade; slides with >=2 grades;
  one-sided-direction reported but two-sided p; Spearman trend over grade (0,1,2) at spot level with slide as a random intercept (lme4), and slide-paired
  Wilcoxon of HG vs normal if >=5 slides have both. Sensitivity: add acinar % and collagen % as covariates (spots are mixtures).
- H2 (Elhossiny epithelial pseudobulk): GIScore increases Donor < AdjNormal < GLTumor < PDTumor. Test: Jonckheere-Terpstra / Spearman trend over ordered classes;
  GL vs PD paired signed-rank (n = number of patients with both).
- H3 (null-honest check): the effect of GALNT16 is not explained by tumor-epithelial fraction (BayesPrism "Tumor-Epithelial" fraction) -> report correlation.
- Each gene-level test reported regardless of direction. No re-definition of grades, no spot or sample exclusion after looking at results.
  Only pre-listed exclusions: spots with `nontissue`==100 and library size < 500 UMI.

## Planned panels (template-inspired grid, 7.2 in wide)
A H&E (hires image) of 2 representative slides with annotated normal/LG/HG epithelial spots; B same slides coloured by GIScore (spot level);
C boxplots of GALNT16 / GALNT14 / GIScore by grade (slide-paired lines); D forest/heat of all 22 GALNT genes (grade trend rho);
E Elhossiny epithelial pseudobulk by region (donor, adj normal, GL, PD) with paired GL/PD lines; F tumor-fraction control; G summary of the
lesion-axis (normal -> LG -> HG -> GL -> PD) slope for GIScore across both datasets.

## Caveats written in advance
- GSE254829: n = 7 slides, few annotated spots per grade, spots are cell mixtures; grades are cross-sectional, not tracked lineages.
- Elhossiny: deconvolved expression, n is small (4/2/7/7), and "AdjNormal" n=2.
- The two datasets use different platforms; no pooling of values, only agreement of directions.
- GALNT16 is detected in <1% of malignant cells in scRNA-seq (Fig 5); any spot/pseudobulk GALNT16 signal may reflect non-malignant cells. GIScore remains
  a composition-sensitive marker.

## Amendment A1 (2026-10-12, made after inspecting the Elhossiny object structure but BEFORE computing any GALNT statistic)
- EpithelialExp `data` layer has a floor artefact (zero counts map to 2.721) and non-standard scale -> use `counts` and compute log1p(CPM) (counts / colSums x 1e6) instead.
- Assay `BayesPrismFraction` turned out to be a 0/1 indicator (not a continuous fraction) -> H3 as written (correlation with tumor-epithelial fraction) is NOT testable.
  H3 is replaced by: report raw GALNT16/GALNT14 counts per pseudobulk (detection) and flag that GALNT16 has 0-13 counts per sample (very low information).
- H1: slides with both normal and LG: 4 (HG2, LG1, LG2, LG3); normal and HG: 4 (HG1-4); LG and HG: 1 (HG2) -> the slide-paired Wilcoxon condition (>=5 slides) is not met;
  the mixed model (slide random intercept) is the primary analysis. Spot GIScore z-scored over all annotated epithelial spots (normal+LG+HG) used in the analysis.

## Results addendum (2026-10-12; analysis script 32_figure6_analysis.R; outputs 04_results/2026-10-12_Figure6_GALNT_run01/)
- H1 (Visium PanIN, 7 slides, 623 annotated epithelial spots: 354 normal / 110 LG / 159 HG): NOT SUPPORTED for GIScore. Mixed-model slope per grade: GIScore -0.011 (p=0.57; with acini/collagen covariates -0.030, p=0.16);
  GALNT16 -0.014 (p=0.26); GALNT14 -0.037 (p=0.0009, direction OPPOSITE to the hypothesis, persists with covariates p=0.0015). LG and HG vs normal: GIScore p=0.67 / 0.46.
- Exploratory (all 22 genes, BH): GALNT3, GALNT4, GALNT5, GALNT6, GALNT10, C1GALT1 increase with grade (BH<0.001, largest GALNT6 beta 0.18/grade, GALNT3 0.20); not part of the frozen score; spots are cell mixtures; some random-intercept fits were singular.
- H2 (Elhossiny epithelial pseudobulk, n=20): NOT SUPPORTED. Trend over Donor<AdjNormal<GL<PD: GIScore rho -0.05 (JT p=0.90); GALNT16 rho 0.16 (p=0.52); GALNT14 rho 0.02 (p=1.0). GL vs PD paired (n=7): GIScore p=0.35. Normal(n=6) vs tumor(n=14): p=0.77.
- Elhossiny caveat: library sizes differ by up to ~200x (donors 0.15-0.9M counts vs tumors 1-38M); GALNT16 counts per sample 0-71 (mostly <20): very low information. The exploratory all-gene trends in this dataset (e.g. GALNT6 rho -0.77, opposite to Visium) are most likely depth-driven and are NOT interpretable.
- H3 untestable (see A1).
- Conclusion: along the lesion axis the frozen GIScore does not change; the Visium signal is in other family members (GALNT3/4/5/6/10, C1GALT1), not replicated in the second dataset.

## Figure built (33_figure6_publication.R; 05_figures/2026-10-12_Figure6_GALNT_run01/Figure6.{svg,png,tiff})
A-C: slides HG2 (all three grades) and HG4 (most HG spots) chosen by rule, not by outcome; D-F spot points + slide means; G exploratory 22-gene trend (Visium); H-J Elhossiny pseudobulk (trend, GL/PD pairs, depth caveat).
Template Fig 6 (tumor/stroma/interface, SPP1+ macrophage co-localisation, Reactome, spatial signalling) NOT reproduced; Fig 7 (digital pathology) not attempted.
Caption must state: lesion labels are manual annotations on 7 slides, spots are mixtures, GALNT16 is detected in ~9% of annotated spots, Elhossiny counts differ ~200x in depth.

## Layout amendment 1 (2026-10-12, user review; not frozen yet)
D-F single-gene/GIScore spot panels replaced by one panel D (GIScore only, violins + slide means; GALNT16/GALNT14 spot statistics remain in H1_trend_mixed_model.csv and are quoted in the legend); slide block compressed (A-C titles shortened: H&E / Manual grade / GIScore (spots)); rows wrapped in wrap_elements (no shared gutters); re-lettered: E all-gene trend, F six enzymes, G Elhossiny classes, H GL vs PD, I depth caveat. H title kept as "PanIN-associated enzymes (one line per slide)" at the user's choice. Backup: 03_analysis/33_figure6_publication_before_v2_review.R. No data or statistics changed.

## Layout amendment 2 (2026-10-12, user request)
Spot-level grade panel moved to Supplementary Fig S6 (FigureS6_spot_level_by_grade.{png,svg}; A GALNT16, B GALNT14, C GIScore; violins + slide means; slopes from H1_trend_mixed_model.csv). Main Fig 6 (7.2 x 9.9 in): the whole top row = slide panels A-C (full width, shared legend row so the two slides stay column-aligned); D all-gene trend; E six enzymes; F Elhossiny classes; G GL vs PD; H depth caveat. No data or statistics changed.
