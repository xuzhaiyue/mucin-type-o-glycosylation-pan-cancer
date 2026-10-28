# Figure 2 spec v1 (2026-10-12) - enzyme-only design

Template Fig. 2 (npj Precis Oncol 2026, PMID 42665670) = enzyme-substrate SVM (ASI) models. Under the user-chosen "enzyme_only" design
(UniProt has only 8 curated NEDD8 substrates) that figure cannot be replicated. Our Figure 2 is the enzyme-only analogue of the template's
pathway analysis (template Fig. S1/S2: enzyme-pathway correlations across cancers), promoted to a main figure.

Gene set: 16 enzymes (freeze_v2_enzyme16.md).
Data: TCGA PanCan Atlas EB++ log2(norm+1), primary tumors (LAML: 03), one sample/patient, cancers with >= 30 tumors.
Pathways: MSigDB Hallmark (50) + KEGG legacy (186) via msigdbr 26.1.0 (human), sizes 10-500 after intersect with expressed genes.
Pathway activity: ssGSEA (GSVA 2.4.8) per cancer (Gaussian-kernel GSVA skipped for speed; disclosed).
Association: Spearman rho enzyme expression vs pathway score per cancer; BH within cancer across enzyme x pathway pairs.
Pairs where the enzyme is a member of the pathway are set to NA (circularity).
Proliferation control (mandatory): rank-based partial correlation adjusting the 8-gene proliferation score (MKI67,TOP2A,PCNA,MCM2,MCM6,CCNB1,CDK1,AURKA).
Sig: BH < 0.05 (also |rho| >= 0.3 for "strong").
Planned panels: A workflow; B pan-cancer Hallmark landscape (median rho, #cancers sig); C top pathways by #cancers for conjugation vs
deneddylation groups (raw vs proliferation-adjusted); D proliferation dependence per enzyme; E KIRC enzyme-pathway heatmap.
No gene/pathway chosen using survival or response.

### Layout amendment (2026-10-12, user request): Fig 2A redesigned
Panel A redrawn as a four-card schematic (numbered header, icon, one-line method; initiation/extension colour key at top right). Content unchanged (ssGSEA, Spearman, rank partial correlation controlling proliferation, strong = BH<0.05 and |rho|>=0.3, PAAD focus). Row heights 1.15/5.15/2.6 (total height still 9.4 in). Panels B-F unchanged. Backup script: 03_analysis/19_figure2_galnt_publication_before_A_redesign.R.
