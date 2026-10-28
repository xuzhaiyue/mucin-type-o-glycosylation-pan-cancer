# Figure 7 spec v1 (FROZEN before any GALNT/SPP1 spot-level statistic was computed)

Date: 2026-10-12. Narrative requested by user: PanIN stage = GALNT3/4/5/6/10 + C1GALT1 (Fig 6H, exploratory, Visium only, robust to acini/collagen/depth adjustment, direction consistent in 3-4/4 slides);
tumor stage = GIScore (Fig 4, TCGA prognostic). Fig 7 asks the template's Fig 6 question for our score: in tumor tissue, does GIScore co-localise with SPP1+ macrophages and which
pathways / signalling directions mark co-high regions?
NOTE (stated before running): single-cell (Fig 5) and the region analyses in Fig 6 did NOT show a tumor-stage rise of GIScore; Fig 7 tests a spatial association, it does not prove a tumor-stage role.

## Data (inspected: structure only)
`Elhossiny_CancerDiscovery_2026/spatial.rds` (Seurat, Visium v2): 164,198 spots, 14 sections (5 healthy donors, 2 adjacent-normal, 7 PDAC patients / 8 sections incl. Patient_22 two sections), 18,085 genes,
manual `SpatialDomain` (17 domains incl. PanIN/GLTumor, PDTumor, Macrophages), `RCTD_sub_wts` (20 cell types, incl. Tumor-Epithelial, Macrophages, Fibroblasts). SPP1, CD44, MDK, NCL, GALNT3/6/14/16 present.

## Frozen definitions
- Tumor spots: PDAC sections, SpatialDomain in {PanIN/GLTumor, PDTumor}. (Sections with <100 such spots are excluded and listed.)
- Expression: log1p(counts/nCount x 1e4) per spot.
- GIScore(spot) = -0.344*z(GALNT16) + 0.202*z(GALNT14), z over tumor spots of all sections (frozen Fig 4 coefficients).
- SPP1+ macrophage signature (Macro_SPP1): derived in scRNA GSE205013 (myeloid cells; macrophages with SPP1 in top quartile vs the rest of myeloid; top 20 genes by Wilcoxon log2FC with
  pct_in>=0.3, SPP1 included). Spot score = mean z of signature genes. Derivation is run once and saved; no tuning on spatial data.
- Macrophage abundance covariate = RCTD "Macrophages" weight; tumor-epithelial covariate = RCTD "Tumor-Epithelial" weight; depth covariate = log10 nCount.
- "Co-high" spot = GIScore and Macro_SPP1 both in the top quartile within the section; "GIScore-high only", "Macro-high only" analogous; others = rest.

## Hypotheses and tests
- H7a: within-section Spearman(GIScore, Macro_SPP1) > 0 in tumor spots; combine sections by random-effects meta of Fisher-z (DerSimonian-Laird); report count of sections with rho>0.
  Adjusted version: partial correlation controlling Macrophage weight, Tumor-Epithelial weight, log10 nCount (primary for interpretation; unadjusted is the template analogue).
- H7b: Reactome ORA on genes up in co-high vs rest spots (within-section Wilcoxon, combined by sum of ranks / pooled pseudobulk edgeR-free approach: mean log fold change per section, genes with
  consistent sign in >= 6/8 sections and combined FDR < 0.05, excluding the signature genes and GALNT16/GALNT14 themselves). Background = expressed genes.
- H7c: ligand-receptor spatial coupling for SPP1-CD44 and MDK-NCL: for each spot, neighbour-mean ligand (6 nearest spots) x receptor in the spot; compare mean coupling in co-high vs rest spots against a
  within-section permutation null (1000 label shuffles). Directionality arrows (template J/K) are replaced by this score because no signalling-direction tool is installed (state this in caption).
- Null results are reported. Sections/spots are not re-selected after seeing results.

## Panels planned
A-D four PDAC sections (H&E, tumor/stroma domains, co-high map), E-H spot scatter GIScore vs Macro_SPP1 per section, I pooled forest (unadjusted vs adjusted rho), J Reactome bars,
K-L SPP1-CD44 and MDK-NCL coupling vs permutation null.

## Caveats written in advance
GALNT16 is sparse in Visium (few counts per spot); GIScore spot values are noisy; macrophage-rich regions inflate any macrophage-signature association (hence partial correlation);
n = 7 patients; correlations are spatial association only.

## Results addendum (2026-10-12; scripts 34 (signature), 35 (analysis, F7MODE=raw|smooth); outputs 04_results/2026-10-12_Figure7_GALNT_run01/)
Signature derived once in GSE205013 macrophages (n=36,970; SPP1-high = top quartile, n see json): SPP1, CTSD, APOE, GPNMB, CSTB, CTSB, APOC1, CTSL, RNASE1, IFI6, CD63, FABP5, MIF, NUPR1, LGALS1, ALOX5AP, FTL, FBP1, CD68, PLIN2.
Sections: 7 PDAC sections had >=100 tumor spots (Patient_22_SP_2 has 291, too few for DE; the spec expected 8, there are 7 PDAC sections incl. two for Patient_22).
- H7a (frozen spec, raw spot GIScore): NOT SUPPORTED. Pooled within-section Spearman rho = 0.008 (95% CI -0.012-0.028, p=0.42); adjusted for RCTD macrophage/tumor-epithelial weight and depth rho = 0.032 (p=0.066); 4/7 (unadjusted) and 5/7 (adjusted) sections positive.
- SPEC DEFECT FOUND: 89% of tumor spots have GALNT16 = 0 and 83% share one single GIScore value (both genes zero), so "top quartile" captured the tie (74% of spots "GIScore-high"). Consequently the frozen co-high group ~ Macro_SPP1-high group, and the frozen-spec H7b (Reactome) and H7c (coupling) outputs are NOT interpretable as GIScore effects (files without suffix).
- Amendment A2 (post hoc, made after seeing the above; reported as exploratory): GIScore, GALNT16/14 and Macro_SPP1 averaged over the spot and its 6 nearest tumor-region neighbours (same section); same tests (files with suffix _smooth). Still 36% of spots share the modal value.
  H7a smoothed: pooled rho 0.035 (p=0.32; adjusted 0.051, p=0.11), I2 = 96%.
- Circularity: SPP1 is part of the Macro_SPP1 signature, so the SPP1-CD44 coupling (z = 42.9, 6/6 sections) in "Both" spots is guaranteed by construction and must NOT be presented as evidence. MDK-NCL coupling: Z = -2.1 (2/6 sections positive), no support.
- Reactome ORA of "Both vs rest" (1294-1561 consistently up genes): neutrophil degranulation, ECM organization, collagen formation/degradation, integrin interactions, interleukin signalling (p.adjust down to 1e-32) = stromal / macrophage programs of macrophage-rich regions; not attributable to GIScore (control "Both vs Macro_SPP1-high only" not run).
- Conclusion: no spatial evidence that GIScore co-localises with SPP1+ macrophages; the template Fig 6 storyline cannot be reproduced for GIScore.

## Amendment A3 (2026-10-12, user question: "does GIScore co-localise with myCAF / iCAF or some other cell type?") - FROZEN before running
Systematic co-localisation scan in the PDAC Visium sections (7 sections incl. Patient_22 x2), instead of the single SPP1-macrophage question.
- Targets (22): RCTD weights of all 20 cell types + myCAF score + iCAF score. CAF gene lists (text-reported markers of Elyada 2019 Cancer Discov, checked online, not the full supplementary table):
  myCAF = TAGLN, MYL9, TPM1, TPM2, MMP11, POSTN, HOPX, ACTA2; iCAF = IL6, CXCL8 (IL8), CXCL1, CXCL2, CCL2, CXCL12, CFD, LMNA, DPT, HAS1, HAS2, AGTR1. Score = mean z of log1p-CP10k over genes present.
- Metrics: GIScore (frozen coefficients), GALNT16, GALNT14, each as raw spot value (primary). Sensitivity: the same metrics averaged over the spot + 6 nearest neighbours (same section, all spots).
- Scopes: (A) all PDAC spots; (B) tumor-domain spots only (SpatialDomain in PanIN/GLTumor, PDTumor).
- Statistic: per-section Spearman rho adjusted for log10 nCount (rank-residual partial correlation); sections combined by random-effects DL on Fisher z; BH across the 22 targets within each metric x scope. Null results reported.
- Predictions written BEFORE running (from Fig 5 scRNA): P1 GALNT16 positively associated with Fibroblast/Pericyte weights and CAF scores, hence GIScore negatively (GIScore = -0.344 z(GALNT16) + 0.202 z(GALNT14));
  P2 GALNT14 positively associated with Tumor-Epithelial/Ducts weight; P3 GIScore not associated with Macrophages (consistent with H7a).
- Interpretation rule: association of GIScore with a stromal weight indicates a composition/neighbourhood marker, not a tumor-intrinsic program.

## Results addendum A3 (36_figure7_colocalization.R -> H7d_colocalization_scan.csv; 7 PDAC sections; pooled per-section Spearman adjusted for depth)
All effects are SMALL (|rho| <= 0.15) although many are significant because of ~60k spots.
- P1 (GALNT16 with fibroblast/pericyte neighbourhood, GIScore opposite): CONFIRMED in direction. Smoothed, all PDAC spots: GIScore vs Fibroblasts rho -0.148 (BH 7e-7), iCAF -0.100 (6e-4), Endothelial -0.079, Pericytes -0.058, T-cells -0.077, Macrophages -0.081; vs Tumor-Epithelial +0.161.
  GALNT16 vs Fibroblasts +0.133, iCAF +0.107, Tumor-Epithelial -0.149. Raw (unsmoothed) values have the same direction with roughly half the size. Tumor-domain spots only: GIScore vs Fibroblasts -0.112 (smoothed; raw -0.061), Tumor-Epithelial +0.119.
- myCAF vs iCAF: GIScore vs myCAF score -0.048 (smoothed, all spots; not BH<0.05 in the raw scope) vs iCAF -0.100; GALNT16 vs myCAF +0.029 vs iCAF +0.107 -> GALNT16 sits in the iCAF-like (inflammatory) CAF neighbourhood, not the myCAF one.
- P2 (GALNT14 with tumor-epithelial / ducts): NOT SUPPORTED (GALNT14 vs Tumor-Epithelial +0.017/+0.033; vs Ducts +0.04/+0.08; vs myCAF -0.03/-0.04).
- P3 (GIScore vs macrophages null): only partly: rho -0.037 raw / -0.081 smoothed (BH<0.05), consistent with a stromal-neighbourhood gradient, not a macrophage-specific link.
- Interpretation: GIScore-high spots are tumor-epithelium-rich and fibroblast/iCAF/immune-poor; GIScore behaves as an inverse stromal-neighbourhood marker (via GALNT16), in line with the Fig 5 composition result. Not evidence of a tumor-intrinsic spatial program.

## Figure built (37_figure7_publication.R; 05_figures/2026-10-12_Figure7_GALNT_run01/Figure7.{svg,png,tiff})
A-C representative section (Patient_24_SP: most tumor spots, rule fixed in advance): domains, iCAF score, neighbour-smoothed GIScore; D-E heatmaps of pooled rho (GALNT16, GALNT14, GIScore vs 22 targets; smoothed; all PDAC spots / tumor-domain spots; raw version in H7d CSV);
F GIScore vs CAF subtypes, raw vs smoothed; G GIScore vs SPP1+ macrophage signature per section (frozen raw = primary; smoothed = post hoc).
NOT in the figure (circular or uninformative, kept in CSV): H7b Reactome of "co-high" spots, H7c SPP1-CD44 coupling (SPP1 is inside the signature).
Template Fig 6 panels (H&E per section, SPP1+ macro maps, Reactome, signalling arrows) not reproduced. Fig 8 (digital pathology) not started.

## Amendment A4 (2026-10-12, user question: "what environment does the GIScore epithelium sit in?") - FROZEN before running
Niche analysis restricted to epithelial-dominant spots so that the question is about the epithelium and not about composition.
- Epithelial spots: PDAC spots with RCTD Tumor-Epithelial weight >= 0.5.
- Niche features (computed from the 6 nearest spots in the same section, EXCLUDING the spot itself, so GALNT values of the spot and its niche are not mixed): mean RCTD weight of each of the 20 cell types, mean myCAF and iCAF scores, tumor-domain fraction of neighbours (purity),
  distance (in median nearest-neighbour spacings) to the nearest Fibro1-4 domain spot.
- Metrics (raw spot values only, no smoothing): GIScore (continuous), GALNT16 detected vs not, GALNT14 detected vs not (the latter two as per-section Spearman with the 0/1 indicator).
- Statistic: per-section Spearman partial for log10 nCount and the spot's own Tumor-Epithelial weight; random-effects DL meta across the sections; BH across features within each metric.
- Predictions: GIScore-high epithelium has purer (more epithelial) neighbourhoods, fewer fibroblast/iCAF/immune neighbours and lies farther from fibroblast domains; GALNT16+ epithelial spots lie closer to fibroblasts/iCAF.
- Caveat: Visium spots (55 um) mix cells and transcripts diffuse between neighbours, so a GALNT16 signal in an "epithelial" spot may be stromal bleed.

## Results addendum A4 (38_figure7_epithelial_niche.R -> H7e_epithelial_niche.csv; 13,574 epithelial-dominant spots, 7 sections; niche = 6 NN excluding self)
- Prediction (purer, fibroblast/iCAF/immune-poor, farther from fibroblasts) NOT SUPPORTED in size: GIScore vs niche |rho| <= 0.035 for all features. Fibroblast neighbours rho -0.026 (BH 0.027, 0/7 sections positive, I2=5%), Endocrine -0.035 (BH 0.045, 1/7 positive); purity +0.033 (BH 0.52); iCAF -0.008 (BH 0.77); myCAF +0.017; distance to fibroblast domain +0.019 (BH 0.87).
- GALNT16-detected epithelial spots (10.9% of spots): neighbours richer in Endocrine (rho +0.054, BH 4e-4, 7/7 sections) and poorer in myCAF (-0.046, BH 0.004, 0/7), lower purity (-0.043, BH 0.036); no iCAF link. GALNT14-detected (9.1%): more Ducts neighbours (+0.031, BH 0.002, 6/7); fibroblast/myCAF neighbours slightly fewer (n.s.).
- Interpretation: the spot-level GIScore-stroma association of A3 is carried by the composition of the spot itself; within epithelial-dominant spots the niche shows only tiny effects. GALNT16 signal in epithelial spots is more compatible with islet/endocrine proximity (or bleed) than with a fibroblast niche.

## Amendment A5 (2026-10-12, user: "analyse with the Visium data" after the GALNT14/16 detection problem) - FROZEN before running
Use only well-detected GALNT-system genes, in the two Visium datasets (Elhossiny PDAC/normal; GSE254829 PanIN).
- D0 detectability: fraction of spots with >=1 UMI per gene in (a) all Elhossiny spots (tissue-wise) and (b) all GSE254829 spots (lib>=500). DETECTABLE = detected in >=30% of spots in BOTH datasets (threshold fixed now; list produced by the script before any other statistic).
- D1 lesion axis replication in Elhossiny spots (replicates Fig 6H, which failed there only because of pseudobulk depth): groups (ordinal g): 0 = SpatialDomain "Ducts" in Healthy or AdjNormal sections; 1 = "PanIN/GLTumor" in PDAC; 2 = "PDTumor" in PDAC.
  Model per detectable gene: lmer(logCP10k ~ g + log10 nCount + (1|section)); also categorical contrasts g1 vs g0, g2 vs g0; BH across detectable genes. Prediction: the six PanIN-associated enzymes of Fig 6H (GALNT3/4/5/6/10, C1GALT1) that are detectable increase with g.
  Caveat: "Ducts" is a spatial domain label of mixed spots; healthy duct spots are few (1,044 + 664).
- D2 co-localisation (PDAC sections, all spots): per-section partial Spearman (log nCount) of each detectable gene vs 20 RCTD weights + myCAF + iCAF scores; DL meta over 7 sections; BH within gene. No prediction other than: GALNT3/6 (known epithelial-lineage enzymes) positive with Tumor-Epithelial/Ducts.
- D3 niche of epithelial-dominant spots (as A4): detectable-gene expression vs neighbour composition (6 NN excluding self), partial for depth and own epithelial weight; meta; BH within gene.
- Null results reported; the genes list is not edited after seeing results.

### A5b (written after the detectability table was printed: the frozen all-spot criterion returned only GALNT1, GALNT2, GALNT11; the six Fig 6H enzymes fail it because their detection in GSE254829 over ALL spots is 3-26%)
Post hoc second list (explicitly exploratory): detected in >=30% of EPITHELIAL spots (Elhossiny: PDAC spots with Tumor-Epithelial weight >=0.5; GSE254829: manually annotated normal/LG/HG epithelial spots) in BOTH datasets. D1-D3 run on both lists; the frozen list is the primary.

## Results addendum A5 (39_figure7_visium_detectable.R -> H7f_*.csv; Elhossiny Visium, 7 PDAC sections + 7 normal/adjacent sections)
IMPORTANT: spot-level mixed models (23,895 spots, section random intercept only) ignore spatial autocorrelation, so p/BH values (down to 1e-282) are NOT meaningful; interpret effect sizes and the number of sections (k=7) with the same sign.
- D0 detectability: frozen all-spot criterion -> only GALNT1, GALNT2, GALNT11 (>=30% in both datasets). Post hoc epithelial-spot criterion (A5b) -> GALNT1, 2, 3, 4, 7, 10, 11, 12. GALNT14 (5-9%) and GALNT16 (9-13%) fail in both; GALNT5/6/C1GALT1 (Fig 6H set) fail both criteria (GALNT6 31%/22% epithelial; C1GALT1 52%/23%).
- D1 lesion axis (normal duct n=1,708 -> PanIN/GL tumor n=8,991 -> PD tumor n=13,196; log1p CP10k, depth-adjusted):
  frozen list: GALNT2 0.43 -> 0.61 -> 0.65 (slope +0.11/step), GALNT11 0.31 -> 0.47 -> 0.48, GALNT1 0.35 -> 0.50 -> 0.46 (rise at PanIN/GL, plateau in PD).
  post hoc list: GALNT3 0.47 -> 0.93 -> 0.44, GALNT4 0.40 -> 0.70 -> 0.42, GALNT7 0.25 -> 0.61 -> 0.34, GALNT12 0.28 -> 0.62 -> 0.38 = PEAK at PanIN/GL tumor, back to ~normal in poorly differentiated tumor (inverted U; the linear slope is not the right summary); GALNT10 0.39 -> 0.73 -> 0.75 = rise and stay.
- D2 co-localisation (7 PDAC sections, depth-adjusted): GALNT3/4/7/12 with Tumor-Epithelial weight rho +0.30/+0.23/+0.19/+0.27 (7/7 sections) and with Fibroblasts -0.19/-0.16/-0.09/-0.15 (0-1/7 positive); GALNT10 with myCAF +0.18 (7/7); GALNT1/2/11 weakly with myCAF/fibroblasts (+0.06 to +0.12).
- D3 niche of epithelial-dominant spots: GALNT3/4/7 epithelial spots have fewer myCAF neighbours (rho -0.14/-0.13/-0.09; 0-1/7 sections positive), GALNT10 epithelial spots sit in fibroblast/myCAF/endothelial-rich neighbourhoods (+0.09/+0.09/+0.07; 7/7); GALNT2 +0.08 myCAF.
- Selection caveat: the post hoc list was defined after seeing that the frozen criterion kept only GALNT1/2/11; it is exploratory. "PanIN/GLTumor" is a mixed domain label (PanIN + gland-forming tumor); "Ducts" healthy/adjacent domain n small.
- Fig 6H (GSE254829) rises with PanIN grade for GALNT3/4/5/6/10/C1GALT1; Elhossiny supports GALNT3/4/10 (PanIN/GL domain up vs normal duct) but shows decline for GALNT3/4 in PD tumor; GALNT5/6/C1GALT1 not testable (low detection).

## Figure 7 REPLACED (2026-10-12): new Fig 7 = PanIN/gland-forming-stage enzymes in Visium (41_figure7_publication.R after 40_figure7_export.R; 05_figures/2026-10-12_Figure7_GALNT_run01/Figure7.{svg,png,tiff}).
The GIScore-spatial version is kept as Figure7_v1_GIScore_spatial_superseded.* (candidate for Supplement S7; script 37_figure7_publication_v1_GIScore_superseded.R).
Panels: A detection (% spots) of the 22 genes in four Visium settings; B one PDAC section (Patient_24_SP: most tumor spots, rule fixed in advance): domains, tumor-epithelial weight, GALNT3, GALNT10;
C lesion axis at SECTION level (7 normal-duct sections, 7 PanIN/GL and 7 PD sections; GL vs normal Mann-Whitney exact, GL vs PD paired signed-rank exact, uncorrected): GALNT3 p=0.001/0.016, GALNT4 0.053/0.016, GALNT7 0.002/0.016, GALNT12 0.011/0.031, GALNT10 0.002/0.219,
GALNT1 0.128/0.688, GALNT2 0.165/0.016 (PD higher), GALNT11 0.128/0.047. (With 7 pairs the smallest possible exact two-sided paired p is 0.016.) These section-level p-values replace the unreliable spot-level p-values; no multiple-testing correction across 8 genes x 2 contrasts.
D-E co-localisation / epithelial niche (post hoc list; dagger = same sign in all 7 sections); F concordance with Fig 6H (GSE254829 slope per grade vs Elhossiny PanIN/GL minus normal duct): GALNT3, 4, 10 positive in both, GALNT7 small, GALNT12 ~0 in GSE254829.
Status: exploratory; gene list defined post hoc (A5b); PanIN/GLTumor is a mixed domain; GIScore (GALNT14/16) cannot be assessed in space (detection 5-13%).
Fig 8 (digital pathology) still not started.
