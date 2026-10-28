"""Figure 5 (template-style layout) data export for GSE205013: UMAP, marker/GALNT dot data per cell type and per malignant state, patient composition; GSE202051 inferCNV score by lineage."""
import json
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

ROOT = Path("${PROJECT_ROOT}"); SC = ROOT / "02_intermediate/fig5_sc"; RES = ROOT / "04_results/2026-10-12_Figure5_GALNT_run01"
BASE = Path("${PANIN_ROOT}/02_intermediate/2026-09-16_six_cohort_input_bundle/standardized")
ALIAS = {"GALNTL2": "GALNT15", "GALNTL1": "GALNT16", "WBSCR17": "GALNT17", "GALNTL4": "GALNT18"}
G = json.load(open(SC / "fig5_genes.json")); GAL = G["GALNT"]

c = "GSE205013"; M = sp.load_npz(SC / f"{c}_counts.npz").tocsr(); genes = [ALIAS.get(g, g) for g in (SC / f"{c}_genes.txt").read_text().split("\n")]
ob = pd.read_csv(SC / f"{c}_obs.csv.gz", low_memory=False); lib = ob["lib"].values.astype(float); L = sp.diags(1e4 / lib) @ M; L.data = np.log1p(L.data); X = pd.DataFrame(L.toarray(), columns=genes)
X = X.loc[:, ~X.columns.duplicated()]
main = {"Tumor_Epithelial": "Epithelial (malignant)", "Tumor_Basal": "Epithelial (malignant)", "T_NK": "T/NK", "Myeloid": "Myeloid", "CAF": "Fibroblast/CAF", "B_Plasma": "B/Plasma", "Endothelial": "Endothelial", "Pericyte_SMC": "Pericyte/SMC", "Acinar": "Acinar", "Mast": "Mast"}
ob["main"] = ob["ann_L2"].map(main); keep = ob["keep_ref_clean"].astype(str).str.lower().isin(["true", "1"]); ob.loc[~keep, "main"] = np.nan
with h5py.File(BASE / f"{c}.h5ad", "r") as h: um = h["obsm/X_umap"][:]
ob["UMAP1"], ob["UMAP2"] = um[:, 0], um[:, 1]
rng = np.random.default_rng(666)
a = ob[ob.main.notna()]; a.iloc[rng.choice(len(a), min(80000, len(a)), replace=False)][["UMAP1", "UMAP2", "main"]].to_csv(RES / "A_umap_main_GSE205013.csv", index=False)
mal = ob[(ob.main == "Epithelial (malignant)") & ob.ann_L3.isin(["Ductal", "Secretory_Ductal", "Cycling_Tumor", "Basal_Squamous", "Boundary"])]
mal.iloc[rng.choice(len(mal), min(40000, len(mal)), replace=False)][["UMAP1", "UMAP2", "ann_L3"]].to_csv(RES / "E_umap_malignant_states_GSE205013.csv", index=False)
mk = ["EPCAM", "KRT19", "KRT8", "KRT18", "CD3E", "NKG7", "MS4A1", "JCHAIN", "CD68", "LYZ", "CD14", "COL1A1", "DCN", "LUM", "PECAM1", "VWF", "RGS5", "ACTA2", "PRSS1", "CPA1"]
want = [g for g in mk + GAL if g in X.columns]
def dots(key: pd.Series, name: str) -> pd.DataFrame:
    out = []
    for l, idx in ob.groupby(key).groups.items():
        sub = X.loc[idx, want]; out.append(pd.DataFrame({name: l, "gene": want, "mean_log": sub.mean().values, "pct_expr": (sub > 0).mean().values * 100, "n_cells": len(idx)}))
    return pd.concat(out)
dots(ob["main"], "celltype").to_csv(RES / "B_dot_markers_GSE205013.csv", index=False)
st = ob["ann_L3"].where(ob.main.eq("Epithelial (malignant)") & ob.ann_L3.isin(["Ductal", "Secretory_Ductal", "Cycling_Tumor", "Basal_Squamous", "Boundary"])); dots(st, "state").to_csv(RES / "F_dot_states_GSE205013.csv", index=False)
comp = ob.dropna(subset=["main"]).groupby(["patient_id", "main"]).size().unstack(fill_value=0); comp = comp[comp.sum(1) >= 200]; comp.div(comp.sum(1), axis=0).reset_index().to_csv(RES / "C_patient_proportions_GSE205013.csv", index=False); print("patients", len(comp))
o2 = pd.read_csv(SC / "GSE202051_obs.csv.gz", usecols=["patient_id", "Level 1 Annotation", "cnv_score"]); o2.to_csv(RES / "D_cnv_score_GSE202051.csv.gz", index=False)
print(ob.main.value_counts().to_dict(), "cnv non-null", o2.cnv_score.notna().sum())
