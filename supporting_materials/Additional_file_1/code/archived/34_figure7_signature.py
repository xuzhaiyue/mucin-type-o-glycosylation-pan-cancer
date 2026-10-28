"""Fig 7 step 1: derive the SPP1+ macrophage signature in scRNA GSE205013 (frozen definition in 00_meta/figure7_spec_v1.md).
Macrophages = ann_L3 in {TREM2_C1QC_TAM, Inflammatory_Myeloid}; SPP1-high = top quartile of SPP1 (log1p CP10k) among macrophages (cells with SPP1>0 if the quartile is 0, noted);
signature = top 20 genes by log2FC (SPP1-high vs other macrophages) among genes with pct_in>=0.3 and Wilcoxon BH<0.05, SPP1 included."""
import json
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path("${PROJECT_ROOT}"); OUT = ROOT / "04_results/2026-10-12_Figure7_GALNT_run01"; OUT.mkdir(parents=True, exist_ok=True)
H5 = Path("${PANIN_ROOT}/02_intermediate/2026-09-16_six_cohort_input_bundle/standardized/GSE205013.h5ad")
o = pd.read_csv(ROOT / "02_intermediate/fig5_sc/GSE205013_obs.csv.gz"); mac = np.where(o.ann_L3.isin(["TREM2_C1QC_TAM", "Inflammatory_Myeloid"]).values)[0]; lib = o["lib"].values
f = h5py.File(H5, "r"); var = f["var"]; idx = var.attrs["_index"] if "_index" in var.attrs else "_index"; names = np.array([x.decode() if isinstance(x, bytes) else x for x in var[idx][:]])
g = f["layers/counts"]; n = int(g.attrs["shape"][0]); assert g.attrs["encoding-type"] == "csc_matrix"; indptr = g["indptr"][:]
inmac = np.zeros(n, bool); inmac[mac] = True; pos_in_mac = -np.ones(n, int); pos_in_mac[mac] = np.arange(len(mac))
def col(j):  # normalised log1p CP10k values of macrophage cells for gene column j
    a, b = indptr[j], indptr[j + 1]; ii = g["indices"][a:b]; dd = g["data"][a:b].astype(float); k = inmac[ii]; v = np.zeros(len(mac)); v[pos_in_mac[ii[k]]] = np.log1p(dd[k] / lib[ii[k]] * 1e4); return v
spp1 = col(int(np.where(names == "SPP1")[0][0])); q = np.quantile(spp1, 0.75); thr_note = "top quartile (q75=%.3f)" % q
hi = spp1 >= q if q > 0 else spp1 > 0
print("macrophages", len(mac), "SPP1-high", hi.sum(), thr_note)
rows = []
for j in range(len(names)):
    v = col(j); pin = (v[hi] > 0).mean(); mi, mo = v[hi].mean(), v[~hi].mean()
    rows.append((names[j], mi, mo, pin, (v[~hi] > 0).mean(), np.log2((np.expm1(mi) + 1) / (np.expm1(mo) + 1))))
    if j % 3000 == 0: print(j, flush=True)
D = pd.DataFrame(rows, columns=["gene", "mean_hi", "mean_lo", "pct_hi", "pct_lo", "log2fc"])
cand = D[(D.pct_hi >= 0.3) & (D.log2fc > 0)].sort_values("log2fc", ascending=False).head(300).copy()
pv = []
for gname in cand.gene:
    v = col(int(np.where(names == gname)[0][0])); pv.append(mannwhitneyu(v[hi], v[~hi], alternative="greater").pvalue)
cand["p"] = pv; cand["BH"] = multipletests(cand.p, method="fdr_bh")[1]; sig = cand[cand.BH < 0.05].sort_values("log2fc", ascending=False).head(20)
genes = list(dict.fromkeys(["SPP1"] + sig.gene.tolist()))[:20]
cand.to_csv(OUT / "S7_macro_spp1_candidates.csv", index=False); D.to_csv(OUT / "S7_macro_spp1_all_genes.csv.gz", index=False)
json.dump({"definition": "SPP1-high macrophages (GSE205013) vs other macrophages", "threshold": thr_note, "n_macrophages": int(len(mac)), "n_hi": int(hi.sum()), "genes": genes}, open(OUT / "macro_spp1_signature.json", "w"), indent=1)
print(genes)
