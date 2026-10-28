# Freeze v4 - six-system comparison for Figure 3 (2026-10-12), written BEFORE running OST/GALNT/GPI/UFM1 and BEFORE the k-means re-run

## Why
Neddylation (v2) and SUMO (v3) were chosen/switched after seeing Fig. 3 results (post hoc). To avoid further outcome-driven switching, all candidate systems are now run with ONE pre-specified pipeline and ALL are reported. No system is dropped on the basis of its survival result.

## Gene sets (outcome-blind; MSigDB 2026.1 via msigdbr)
- NEDD8 (16): freeze_v2_enzyme16.md
- SUMO (22): freeze_v3_sumo22.md
- OST / N-glycosylation (12): GOCC_OLIGOSACCHARYLTRANSFERASE_COMPLEX (13 genes) minus TREX1 (DNase, not an OST subunit) = DAD1, DDOST, KRTCAP2, MAGT1, OST4, OSTC, RPN1, RPN2, STT3A, STT3B, TMEM258, TUSC3
- GALNT / mucin-type O-glycosylation (22): GOMF_POLYPEPTIDE_N_ACETYLGALACTOSAMINYLTRANSFERASE_ACTIVITY (GALNT1-18, GALNTL6) + from REACTOME_O_LINKED_GLYCOSYLATION_OF_MUCINS: GALNTL5, C1GALT1, C1GALT1C1
- GPI anchor biosynthesis (24): KEGG_GLYCOSYLPHOSPHATIDYLINOSITOL_GPI_ANCHOR_BIOSYNTHESIS (25) minus GPLD1 (GPI-degrading phospholipase, not synthesis)
- UFM1 (8): GOBP_PROTEIN_UFMYLATION (UFM1, UBA5, UFC1, UFL1, DDRGK1, CDK5RAP3, UFSP1) + UFSP2 (main de-UFMylase; in the GO ubiquitin-like peptidase set)

## Pipeline (same for all systems; supersedes the hierarchical-clustering version used for the first two runs, which is kept as a sensitivity result)
- TCGA EB++ log2(norm+1), primary tumors, 33 cancers (>=30 tumors), OS from TCGA-CDR truncated at 120 months.
- Per cancer, drop genes with median expression < 1 (log2 scale; essentially not expressed); require >=6 genes remaining, else the cancer is NA for that system. Then gene-wise z-score.
- ConsensusClusterPlus, k = 2 fixed, clusterAlg = "km" (k-means, Euclidean), reps = 1000, pItem = 0.8, pFeature = 1, seed 666. (Reason for km: hierarchical Ward clustering often isolated tiny clusters, e.g. COAD 16 vs 416.)
- C1 = cluster with better OS, C2 = worse. Log-rank p; flag "unbalanced" if minority cluster < 10% of patients (reported, excluded from the headline count).
- Multiplicity: BH across 33 cancers within each system AND BH across all systems x cancers (global). Headline = global BH < 0.05; nominal p<0.05 reported as descriptive.
- Digestive cancers (LIHC, STAD, COAD, READ, ESCA, PAAD, CHOL) always listed.
- The system carried forward to Figures 1, 2, 4+ will be chosen by the STATED priority: (1) global BH<0.05 in >=1 digestive cancer, else (2) biological rationale (digestive-system relevance, gene count) decided with the user and disclosed as such. A system chosen after seeing these results must be labelled post hoc.
