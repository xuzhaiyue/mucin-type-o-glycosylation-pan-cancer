"""Fig 5B support: mean log1p_counts and fraction expressing of GALNT14 / GALNT16 per ann_fig1B group, computed on the
RefA_TME_RefClean_Fig1B.h5ad object exactly as the existing marker dotplot (ad4_dotplot_ann_fig1B.tsv: Background excluded, layer log1p_counts,
display = min-max of group means per gene). The UMAP table and the 57-gene dotplot table are read directly from the source run in the R script (not copied)."""
from __future__ import annotations
from pathlib import Path
import h5py, numpy as np, pandas as pd

SRC = Path("${PDAC_CAF_CHEMO_ROOT}/runs/2026-03-19_refA_tme")
OUT = Path("${PROJECT_ROOT}/04_results/2026-10-12_Figure5_GALNT_run01"); OUT.mkdir(parents=True, exist_ok=True)
GENES = ["GALNT14", "GALNT16"]


def main() -> None:
    f = h5py.File(SRC / "artifacts/RefA_TME_RefClean_Fig1B.h5ad", "r")
    var = f["var"]; names = np.array([x.decode() if isinstance(x, bytes) else x for x in var["_index"][:]]); pos = {g: i for i, g in enumerate(names)}
    cat = f["obs/ann_fig1B/categories"][:]; cat = [x.decode() if isinstance(x, bytes) else x for x in cat]; codes = f["obs/ann_fig1B/codes"][:]
    lab = np.array(cat, dtype=object)[codes]
    g = f["layers/log1p_counts"]; n, ncol = [int(x) for x in g.attrs["shape"]]; indptr = g["indptr"][:]
    cols = {gene: pos[gene] for gene in GENES}; vals = {gene: np.zeros(n, np.float32) for gene in GENES}
    chunk = 20000
    for s in range(0, n, chunk):
        e = min(s + chunk, n); a, b = indptr[s], indptr[e]; ind = g["indices"][a:b]; dat = g["data"][a:b]; ip = indptr[s:e + 1] - a
        rows = np.repeat(np.arange(e - s), np.diff(ip))
        for gene, c in cols.items():
            m = ind == c; vals[gene][s + rows[m]] = dat[m]
        print(e, "/", n, flush=True)
    df = pd.DataFrame({"group": lab, **vals}); df = df[df.group != "Background"]
    rows = []
    for gene in GENES:
        s = df.groupby("group")[gene].agg(mean_expression="mean", fraction_expressing=lambda x: float((x > 0).mean()), n_cells_group="size").reset_index(); s["gene"] = gene
        z = (s.mean_expression - s.mean_expression.mean()) / s.mean_expression.std(ddof=0); s["scaled_expression"] = z
        s["display_expression"] = (s.mean_expression - s.mean_expression.min()) / (s.mean_expression.max() - s.mean_expression.min()); rows.append(s)
    out = pd.concat(rows); out.to_csv(OUT / "B_fig1B_GALNT14_16_dots.csv", index=False)
    print(out.pivot(index="group", columns="gene", values="fraction_expressing").round(3))
    print("cells (Background excluded):", len(df), "| source table totals check:", pd.read_csv(SRC / "figures/ad4_dotplot_ann_fig1B.tsv", sep="\t").query("gene=='EPCAM'").n_cells_group.sum())


if __name__ == "__main__":
    main()
