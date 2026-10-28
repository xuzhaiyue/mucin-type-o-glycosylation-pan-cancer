"""Fig 5C-E data: GSE205013 clean malignant epithelium (adata_malignant_r5_clean.h5ad; 6 states, 41,073 cells).
C: UMAP (obsm X_umap) coloured by state (ann_L2 == ann_celltype, 1:1).  D: marker dotplot with the panels of the source notebook
(notebooks/Epi/Tumor_tme_scvi.ipynb cell 60: marker_panels_raw; layer log1p_counts, groupby ann_celltype, standard_scale='var' = min-max of group means per gene,
dot size = fraction expressing, state order = scanpy dendrogram on X_pca as in the notebook).  E: 22 GALNT-system genes per state (mean log1p, fraction expressing)."""
from __future__ import annotations
import json
from pathlib import Path
import h5py, numpy as np, pandas as pd

SRC = Path("${PDAC_CAF_CHEMO_ROOT}/runs/2026-05-10_Epithelial_scvi/artifacts/adata_malignant_r5_clean.h5ad")
OUT = Path("${PROJECT_ROOT}/04_results/2026-10-12_Figure5_GALNT_run01"); OUT.mkdir(parents=True, exist_ok=True)
PANELS = {
    "Epithelial_Core": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "Classical_PDAC": ["GATA6", "TFF1", "TFF2", "MUC1", "AGR2", "CEACAM5"],
    "Basal_Squamoid": ["KRT5", "KRT6A", "KRT17", "S100A2", "TP63", "CAV1"],
    "YAP_TAZ_Adhesion": ["CTGF", "CYR61", "ANKRD1", "AXL", "ITGB5", "ITGAV"],
    "EMT_Mesenchymal": ["VIM", "ZEB1", "SNAI1", "SNAI2", "CDH2"],
    "Hypoxia_Acidosis": ["HIF1A", "CA9", "LDHA", "SLC16A3", "ENO1"],
    "Oxidative_Stress": ["NQO1", "HMOX1", "SQSTM1", "FTL", "FTH1"],
    "Cycling_Malignant": ["MKI67", "TOP2A", "BIRC5", "CENPF", "UBE2C"],
}
GALNT = [f"GALNT{i}" for i in range(1, 19)] + ["GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"]
ALIAS = {"GALNT15": ["GALNTL2"], "GALNT16": ["GALNTL1"], "GALNT17": ["WBSCR17"], "GALNT18": ["GALNTL4"]}


def cat(f: h5py.File, key: str) -> np.ndarray:
    c = [x.decode() if isinstance(x, bytes) else x for x in f[f"obs/{key}/categories"][:]]; return np.array(c, dtype=object)[f[f"obs/{key}/codes"][:]]


def main() -> None:
    f = h5py.File(SRC, "r")
    names = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/_index"][:]]); pos = {g: i for i, g in enumerate(names)}
    state = cat(f, "ann_L2"); state2 = cat(f, "ann_celltype"); patient = cat(f, "Patient_ID"); um = f["obsm/X_umap"][:]
    assert len(set(zip(state, state2))) == len(set(state)), "ann_L2 and ann_celltype are not 1:1"
    pd.DataFrame({"UMAP1": um[:, 0], "UMAP2": um[:, 1], "state": state, "patient": patient}).to_csv(OUT / "C_malignant_clean_umap.csv.gz", index=False)
    # dendrogram order exactly as scanpy in the notebook (groupby ann_celltype on X_pca)
    import anndata as ad, scanpy as sc
    a = ad.AnnData(X=np.zeros((len(state2), 1), dtype=np.float32)); a.obs["ann_celltype"] = pd.Categorical(state2); a.obsm["X_pca"] = f["obsm/X_pca"][:]
    sc.tl.dendrogram(a, groupby="ann_celltype", use_rep="X_pca"); order_ct = list(a.uns["dendrogram_ann_celltype"]["categories_ordered"])
    ct2state = dict(zip(state2, state)); order = [ct2state[c] for c in order_ct]; print("dendrogram order:", order)
    # collect genes
    want = {}
    for mod, gl in PANELS.items():
        for g in gl:
            if g in pos: want.setdefault(g, mod)
            else: print("missing marker", mod, g)
    gal = {}
    for g in GALNT:
        for cand in [g] + ALIAS.get(g, []):
            if cand in pos: gal[g] = cand; break
        else: print("missing GALNT gene", g)
    genes = sorted(set(want) | set(gal.values()), key=lambda x: pos[x]); cols = np.array([pos[g] for g in genes]); cidx = {g: i for i, g in enumerate(genes)}
    L = f["layers/log1p_counts"]; n, _ = [int(x) for x in L.attrs["shape"]]; ip = L["indptr"][:]; M = np.zeros((n, len(genes)), np.float32); colmap = {c: j for j, c in enumerate(cols)}
    for s in range(0, n, 10000):
        e = min(s + 10000, n); a_, b_ = ip[s], ip[e]; ind = L["indices"][a_:b_]; dat = L["data"][a_:b_]; rows = np.repeat(np.arange(e - s), np.diff(ip[s:e + 1]))
        for c, j in colmap.items():
            m = ind == c; M[s + rows[m], j] = dat[m]
    df = pd.DataFrame(M, columns=genes); df["state"] = state
    mean = df.groupby("state").mean(); frac = df.groupby("state").agg(lambda x: float((x > 0).mean())); nn = df.groupby("state").size()
    rows = []
    for g, mod in want.items():
        m = mean[g].loc[order]; disp = (m - m.min()) / (m.max() - m.min()) if m.max() > m.min() else m * 0
        for st in order: rows.append(dict(module=mod, gene=g, state=st, mean_expression=m[st], display_expression=disp[st], fraction_expressing=frac.loc[st, g], n_cells=int(nn[st])))
    pd.DataFrame(rows).to_csv(OUT / "D_malignant_clean_marker_dots.csv", index=False)
    rows = []
    for g, src in gal.items():
        m = mean[src].loc[order]; z = (m - m.mean()) / m.std(ddof=0) if m.std(ddof=0) > 0 else m * 0
        for st in order: rows.append(dict(gene=g, source_gene=src, state=st, mean_expression=m[st], scaled_mean=z[st], fraction_expressing=frac.loc[st, src], n_cells=int(nn[st])))
    pd.DataFrame(rows).to_csv(OUT / "E_malignant_clean_GALNT_dots.csv", index=False)
    json.dump({"state_order": order, "n_cells": int(len(state)), "n_cells_per_state": {k: int(v) for k, v in nn.items()}, "module_of_gene": want}, open(OUT / "CDE_malignant_clean_meta.json", "w"), indent=1)
    print(nn.loc[order]); print(frac.loc[order, [gal["GALNT14"], gal["GALNT16"]]].round(3))


if __name__ == "__main__":
    main()
