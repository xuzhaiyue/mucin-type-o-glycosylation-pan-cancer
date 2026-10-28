# Figure 1 replication spec v1 (2026-10-11)

Template: Zhu et al., npj Precis Oncol 2026 (PMID 42665670), Fig. 1 "Multi-omics investigations on acetylation-related enzymes in 33 cancers".
Gene set: 16 Neddylation enzymes (freeze_v2_enzyme16.md): 12 conjugation + 4 deneddylation (template: 15 acetyltransferases + 18 deacetylases).

| Panel | Template | Replication (this project) | Known deviation |
|---|---|---|---|
| A | Overview of enzymes (writers/erasers) | Tile schematic: 16 enzymes by role/family | Roles from Reactome/UniProt |
| B | Top-10 mutated enzymes | Nonsilent mutation frequency (TCGA MC3 gene-level matrix, UCSC Xena), top 10 of 16 | Template used GDC MAF + maftools; MC3 matrix is the PanCan Atlas call set |
| C | Mutation frequency by cancer (+ mutation load) | 16 genes x 33 cancers heatmap; top bar = mean number of genes with nonsilent mutation per sample (load proxy) | Load proxy = genes mutated per sample, not mutations per sample |
| D | mRNA tumor vs normal + prognosis across 29 cancers; asterisk if up/down ratio > 25% | limma (all genes, unpaired) on Xena EB++ log2(norm+1); cancers with >= 5 TCGA normals; OS log-rank on median split (all cancers, >= 10 events), nominal p < 0.05 | Number of cancers differs from template (29); data = PanCan Atlas EB++, not TCGAbiolinks |
| E | BRCA CPTAC mRNA-protein Spearman vs random pairs | CPTAC-BRCA 2018 (LinkedOmics) mRNA and protein, 16 matched pairs vs 5000 random gene pairs | Source LinkedOmics |
| F | CPTAC protein tumor vs normal; asterisk where mRNA aligned | 7 CPTAC cohorts with normals (CCRCC, COAD, HNSCC, LSCC, LUAD, PDAC, UCEC), limma; asterisk = protein sig and TCGA mRNA sig in same direction | PDAC values are MD-abundance log2, others ratios |

Thresholds (fixed before running): expression sig = BH < 0.05 and |log2FC| >= 0.5; protein sig = BH < 0.05 (no FC cut); survival sig = log-rank p < 0.05 (nominal, as in template);
primary tumors (01; LAML 03), one sample per patient; normals code 11.
No gene was added or removed after seeing results.
