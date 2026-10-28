"""Figure 5 step 2: single-cell analyses of the GALNT system / frozen GIScore (spec: 00_meta/figure5_spec_v1.md).
Patient = unit of analysis; pseudobulk = mean of per-cell log1p(CP10k); >=20 cells per unit."""
import json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path("${PROJECT_ROOT}"); SC = ROOT / "02_intermediate/fig5_sc"
RES = ROOT / "04_results/2026-10-12_Figure5_GALNT_run01"; RES.mkdir(parents=True, exist_ok=True)
G = json.load(open(SC / "fig5_genes.json")); GAL = G["GALNT"]; ALIAS = {"GALNTL2": "GALNT15", "GALNTL1": "GALNT16", "WBSCR17": "GALNT17", "GALNTL4": "GALNT18"}
BETA = {"GALNT16": -0.3436912, "GALNT14": 0.2020544}  # frozen in FROZEN_model.rds (Fig 4)
MINC = 20

def load(c: str):
    M = sp.load_npz(SC / f"{c}_counts.npz").tocsr(); genes = (SC / f"{c}_genes.txt").read_text().split("\n"); ob = pd.read_csv(SC / f"{c}_obs.csv.gz", low_memory=False)
    genes = [ALIAS.get(g, g) for g in genes]; lib = ob["lib"].values.astype(float); lib[lib == 0] = np.nan
    L = sp.diags(1e4 / lib) @ M; L.data = np.log1p(L.data)
    X = pd.DataFrame(L.toarray() if L.shape[1] < 200 else L.toarray(), columns=genes, index=ob.index)
    X = X.loc[:, ~X.columns.duplicated()]; return ob, X

def harmonize(c: str, ob: pd.DataFrame) -> pd.Series:
    if c == "CRA001160":
        m = {"Malignant_epithelial": "Malignant epithelial", "Nonmalignant_ductal": "Non-malignant epithelial", "Acinar": "Acinar", "Endocrine": "Endocrine", "Fibroblast": "CAF/fibroblast", "Stellate": "CAF/fibroblast", "Endothelial": "Endothelial", "Myeloid": "Myeloid", "T_cell": "T/NK", "B_cell": "B/Plasma"}
        return ob["cell_type"].map(m)
    if c == "GSE205013":
        m = {"Tumor_Epithelial": "Malignant epithelial", "Tumor_Basal": "Malignant epithelial", "T_NK": "T/NK", "Myeloid": "Myeloid", "CAF": "CAF/fibroblast", "B_Plasma": "B/Plasma", "Endothelial": "Endothelial", "Pericyte_SMC": "Pericyte/SMC", "Acinar": "Acinar", "Mast": "Myeloid"}
        h = ob["ann_L2"].map(m); h[~ob["keep_ref_clean"].astype(str).str.lower().isin(["true", "1"])] = np.nan; return h
    m = {"Epithelial (malignant)": "Malignant epithelial", "Epithelial (non-malignant)": "Non-malignant epithelial", "Cancer-associated fibroblast": "CAF/fibroblast", "Endothelial": "Endothelial", "Myeloid": "Myeloid", "Pericyte": "Pericyte/SMC", "Vascular smooth muscle": "Pericyte/SMC", "Endocrine": "Endocrine", "Lymphoid": "T/NK/B (lymphoid)"}
    return ob["Level 1 Annotation"].map(m)

def pseudobulk(X, ob, key, cols):
    d = X[cols].copy(); d["_p"] = ob["patient_id"].astype(str).values; d["_l"] = key.values; d = d.dropna(subset=["_l"])
    n = d.groupby(["_p", "_l"]).size().rename("n_cells"); pb = d.groupby(["_p", "_l"])[cols].mean(); pb = pb.join(n); return pb[pb.n_cells >= MINC].reset_index().rename(columns={"_p": "patient", "_l": "lineage"})

def zs(v): v = np.asarray(v, float); return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)
def giscore(df): return BETA["GALNT16"] * zs(df["GALNT16"]) + BETA["GALNT14"] * zs(df["GALNT14"])

dot, pbs, outB, outC, outD, outE, outF, outG, cells_summary = [], [], [], [], [], [], [], [], []
for c in ["GSE202051", "GSE205013", "CRA001160"]:
    ob, X = load(c); lin = harmonize(c, ob); ob["lineage_h"] = lin.values; gal = [g for g in GAL if g in X.columns]
    print(c, "cells", len(ob), "genes GALNT", len(gal), "lineages", ob.lineage_h.value_counts().to_dict())
    cells_summary.append(pd.DataFrame({"cohort": c, "lineage": ob.lineage_h.value_counts().index, "n_cells": ob.lineage_h.value_counts().values}))
    # A: dot plot data (cell level, descriptive)
    for l, idx in ob.groupby("lineage_h").groups.items():
        sub = X.loc[idx, gal]; dot.append(pd.DataFrame({"cohort": c, "lineage": l, "gene": gal, "mean_log": sub.mean().values, "pct_expr": (sub > 0).mean().values * 100, "n_cells": len(idx)}))
    need = list(dict.fromkeys(gal + [g for g in G["moffitt_basal"] + G["moffitt_classical"] if g in X.columns]))
    pb = pseudobulk(X, ob, ob["lineage_h"], need); pb["cohort"] = c
    # patient-level metadata
    meta = ob.groupby(ob["patient_id"].astype(str)).agg(**{k: (k, lambda s: s.dropna().astype(str).mode().iloc[0] if s.notna().any() else np.nan) for k in ["Tissue", "Treatment", "Disease_State"] if k in ob.columns} , **({"condition": ("condition", lambda s: s.dropna().astype(str).mode().iloc[0])} if "condition" in ob.columns else {}))
    # malignant GIScore (z across eligible malignant patients in the cohort)
    mal = pb[pb.lineage == "Malignant epithelial"].copy(); mal["GIScore"] = giscore(mal)
    mal = mal.set_index("patient")
    # B: malignant vs non-malignant epithelium
    nm = pb[pb.lineage == "Non-malignant epithelial"].set_index("patient")
    if len(nm) >= 5:
        both = pb[pb.lineage.isin(["Malignant epithelial", "Non-malignant epithelial"])].copy(); both["GIScore"] = giscore(both)
        w = both.pivot(index="patient", columns="lineage", values=["GALNT16", "GALNT14", "GIScore"]); paired = w.dropna()
        for g in ["GALNT16", "GALNT14", "GIScore"]:
            x = w[g]["Malignant epithelial"].dropna(); y = w[g]["Non-malignant epithelial"].dropna(); pr = paired[g]
            if len(pr) >= 5: st = stats.wilcoxon(pr["Malignant epithelial"], pr["Non-malignant epithelial"]); test = "paired Wilcoxon"; n = len(pr); dmed = float(np.median(pr["Malignant epithelial"] - pr["Non-malignant epithelial"]))
            else: st = stats.mannwhitneyu(x, y); test = "Mann-Whitney"; n = f"{len(x)}/{len(y)}"; dmed = float(x.median() - y.median())
            outB.append({"cohort": c, "gene": g, "test": test, "n_patients": n, "median_diff_mal_minus_nonmal": dmed, "p": st.pvalue, "n_paired": len(pr)})
        both.assign(cohort=c).to_csv(RES / f"B_pseudobulk_{c}.csv", index=False)
    # C: subtype
    if c == "GSE202051":
        d = ob[ob.lineage_h == "Malignant epithelial"].copy(); d["patient"] = d["patient_id"].astype(str)
        s = d.groupby("patient").agg(basal=("Moffitt_basal", "mean"), classical=("Moffitt_classical", "mean"), n=("cell", "size")); s = s[s.n >= MINC]; mal["subtype_score"] = (zs(s.basal.reindex(mal.index)) - zs(s.classical.reindex(mal.index))); mal["subtype_source"] = "authors' Moffitt scores"
    else:
        bg = [g for g in G["moffitt_basal"] if g in mal.columns]; cg = [g for g in G["moffitt_classical"] if g in mal.columns]
        mal["subtype_score"] = np.mean([zs(mal[g]) for g in bg], axis=0) - np.mean([zs(mal[g]) for g in cg], axis=0); mal["subtype_source"] = f"Moffitt lists ({len(bg)} basal / {len(cg)} classical genes present)"
    for g in ["GIScore", "GALNT16", "GALNT14"]:
        ok = mal[[g, "subtype_score"]].dropna()
        if len(ok) >= 8: r, p = stats.spearmanr(ok[g], ok["subtype_score"]); outC.append({"cohort": c, "x": g, "n_patients": len(ok), "spearman_rho": r, "p": p, "subtype_source": mal.subtype_source.iloc[0]})
    mal.reset_index().assign(cohort=c).to_csv(RES / f"C_malignant_pseudobulk_{c}.csv", index=False)
    # E: composition among non-epithelial, non-endocrine, non-acinar cells
    d = ob[["patient_id", "lineage_h"]].dropna().copy(); d["patient"] = d["patient_id"].astype(str)
    d = d[~d.lineage_h.isin(["Malignant epithelial", "Non-malignant epithelial", "Endocrine", "Acinar"])]
    comp = d.groupby(["patient", "lineage_h"]).size().unstack(fill_value=0); comp = comp[comp.sum(1) >= 200]; frac = comp.div(comp.sum(1), axis=0)
    if c == "GSE205013":
        keep = meta[(meta.Tissue == "Pancreas")].index  # primary pancreas only
        frac = frac.loc[frac.index.intersection(keep)]
    cols_c = {"GSE202051": {"Lymphoid (T/NK/B)": "T/NK/B (lymphoid)", "Myeloid": "Myeloid", "CAF/fibroblast": "CAF/fibroblast"}, "GSE205013": {"T/NK": "T/NK", "Myeloid": "Myeloid", "CAF/fibroblast": "CAF/fibroblast"}, "CRA001160": {"T/NK": "T/NK", "Myeloid": "Myeloid", "CAF/fibroblast": "CAF/fibroblast"}}[c]
    for lab, col in cols_c.items():
        if col not in frac.columns: continue
        j = pd.concat([mal["GIScore"], frac[col].rename("frac")], axis=1, join="inner").dropna()
        if len(j) >= 8: r, p = stats.spearmanr(j.GIScore, j.frac); outE.append({"cohort": c, "compartment": lab, "n_patients": len(j), "spearman_rho": r, "p": p})
    # F: treatment; G: tissue
    if c in ("GSE202051", "GSE205013"):
        tcol = "Treatment" if c == "GSE205013" else "condition"; m2 = mal.join(meta[[tcol] + (["Tissue"] if c == "GSE205013" else [])], how="left")
        if c == "GSE205013": f_ = m2[m2.Tissue == "Pancreas"]
        else: f_ = m2
        for g in ["GIScore", "GALNT16", "GALNT14"]:
            a = f_[f_[tcol].isin(["Naive"])][g].dropna(); b = f_[f_[tcol].isin(["Chemo"])][g].dropna()
            if len(a) >= 4 and len(b) >= 4: outF.append({"cohort": c, "gene": g, "n_naive": len(a), "n_chemo": len(b), "median_naive": a.median(), "median_chemo": b.median(), "p": stats.mannwhitneyu(a, b).pvalue})
    if c == "GSE205013":
        m2 = mal.join(meta[["Tissue"]], how="left")
        for g in ["GIScore", "GALNT16", "GALNT14"]:
            a = m2[m2.Tissue == "Pancreas"][g].dropna(); b = m2[m2.Tissue == "Liver"][g].dropna()
            if len(a) >= 4 and len(b) >= 4: outG.append({"cohort": c, "gene": g, "n_pancreas": len(a), "n_liver": len(b), "median_pancreas": a.median(), "median_liver": b.median(), "p": stats.mannwhitneyu(a, b).pvalue})
        # D: tumor states (ann_L3 within epithelium), patient x state pseudobulk
        st = ob["ann_L3"].where(ob["ann_L3"].isin(["Ductal", "Secretory_Ductal", "Cycling_Tumor", "Basal_Squamous", "Boundary"]) & ob.lineage_h.eq("Malignant epithelial"))
        pbs_ = pseudobulk(X, ob, st, ["GALNT16", "GALNT14"]); pbs_["GIScore"] = giscore(pbs_)  # z across all patient x state units
        pbs_.to_csv(RES / "D_state_pseudobulk_GSE205013.csv", index=False)
        pv = pbs_.pivot(index="patient", columns="lineage", values="GIScore")
        for a_, b_ in [("Basal_Squamous", "Ductal"), ("Basal_Squamous", "Secretory_Ductal"), ("Cycling_Tumor", "Ductal")]:
            pr = pv[[a_, b_]].dropna()
            if len(pr) >= 5: outD.append({"contrast": f"{a_} vs {b_}", "n_patients": len(pr), "median_diff": float((pr[a_] - pr[b_]).median()), "p": stats.wilcoxon(pr[a_], pr[b_]).pvalue})
        full = pv.dropna()
        if len(full) >= 5: outD.append({"contrast": "Friedman, 5 states", "n_patients": len(full), "median_diff": np.nan, "p": stats.friedmanchisquare(*[full[k] for k in full.columns]).pvalue})
    pbs.append(pb)

def tab(rows, name, pcol="p", by=None):
    d = pd.DataFrame(rows)
    if len(d): d["BH"] = multipletests(d[pcol], method="fdr_bh")[1]
    d.to_csv(RES / name, index=False); print("\n==", name); print(d.to_string(index=False, float_format=lambda x: f"{x:.3g}")); return d
pd.concat(dot).to_csv(RES / "A_dotplot_data.csv", index=False); pd.concat(cells_summary).to_csv(RES / "A_cells_per_lineage.csv", index=False)
tab(outB, "B_malignant_vs_nonmalignant.csv"); tab(outC, "C_giscore_vs_subtype.csv"); tab(outD, "D_tumor_states_GSE205013.csv"); tab(outE, "E_composition_vs_giscore.csv"); tab(outF, "F_treatment.csv"); tab(outG, "G_primary_vs_liver.csv")
