# Figure 3 spec v1 (frozen before running)

Template: Zhu et al. 2026 Fig. 3 (substrate-expression subtypes with different OS). Enzyme-only analogue.

## Analysis rules
- Data: TCGA PanCan Atlas EB++ log2(norm+1), primary tumors, one sample/patient, 33 cancers with >=30 tumors (same as Fig. 2). OS from TCGA-CDR, truncated at 120 months (as in template axis).
- Features: 16 Neddylation enzymes, per-cancer gene-wise z-score.
- Clustering: ConsensusClusterPlus, k = 2 fixed (template shows C1/C2), hierarchical (ward.D2, Euclidean on z-scores), reps = 1000, pItem = 0.8, pFeature = 1, seed 666. Cluster with better OS (lower HR) named C1, other C2 (poor prognosis).
- Prognosis test: log-rank p (unadjusted, as template) + BH across 33 cancers reported. Cancers with log-rank p < 0.05 are "prognostically distinct". Panels A-D show those (max 9 by smallest p; LIHC always shown and flagged if not significant).
- Panel B: limma (C2 vs C1) on the 16 enzymes; "enrichment score" = logFC; direction annotation = sign of median logFC.
- Panel C: genome-wide limma C2 vs C1, DEG = BH<0.05 & |logFC|>=0.5; over-representation (clusterProfiler::enricher, universe = tested genes) against MSigDB KEGG_LEGACY; shown = top 20 pathways by number of cancers with FDR<0.05, then by mean -log10 FDR. Pathway category assigned manually from KEGG BRITE.
  Deviation: template enriched pathways of substrates; with 16 enzymes ORA is not meaningful, so DEGs between clusters are used.
- Panel D: Venn of significant KEGG pathways (FDR<0.05) in LIHC + 3 other cancers with smallest log-rank p.
- Caveats to disclose: purity/proliferation not adjusted; k=2 fixed, not chosen by data; unadjusted p used for selection (BH reported); clustering on 16 genes largely captures one dominant axis.
