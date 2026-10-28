"""Fig 8 A10-Q4 step 3: before/after ink filtering. Frozen clusters (no re-clustering); same sampled tiles minus ink tiles; H8a model as in 43_figure8_stats.py.
POST HOC QC sensitivity analysis; does not replace the primary analysis."""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

R = "${PROJECT_ROOT}"
RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; Q = f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01"
CS = [f"C{i}" for i in range(1, 7)]


def fractions(tc: pd.DataFrame, n: int) -> pd.DataFrame:
    ct = pd.crosstab(tc.slide, tc.cluster).reindex(index=range(n), columns=range(1, 7), fill_value=0); f = ct.div(ct.sum(1).replace(0, np.nan), axis=0); f.columns = CS; f["n"] = ct.sum(1); return f


def h8a(d: pd.DataFrame, fr: pd.DataFrame) -> pd.DataFrame:
    d = d.copy(); d[CS] = fr[CS].values; d = d.dropna(subset=CS).copy(); d["gi"] = (d.risk - d.risk.mean()) / d.risk.std(); d["grade"] = d.gradeg.fillna("G1-2"); rows = []
    clr = np.log(d[CS].clip(lower=1e-3)); clr = clr.sub(clr.mean(1), axis=0)
    for c in CS:
        d["y"] = clr[c]; m1 = smf.ols("y ~ gi + C(grade) + C(tss_grp)", d).fit(cov_type="HC3"); rho, prho = spearmanr(d.risk, d[c])
        rows.append(dict(cluster=c, mean_frac=d[c].mean(), rho=rho, beta_adj=m1.params["gi"], lo=m1.conf_int().loc["gi", 0], hi=m1.conf_int().loc["gi", 1], p_adj=m1.pvalues["gi"], n=len(d)))
    t = pd.DataFrame(rows); t["bh_adj"] = multipletests(t.p_adj, method="fdr_bh")[1]; return t


def main() -> None:
    d = pd.read_csv(f"{RES}/H8_slide_table.csv"); n = len(d); tc = pd.read_csv(f"{R}/02_intermediate/fig8/tile_clusters.csv.gz")
    have = [p for p in d.patient if os.path.exists(f"{Q}/ink_tiles/{p}.csv.gz")]; print("slides with ink QC:", len(have), "of", n)
    parts = []
    for i, p in enumerate(d.patient):
        if p not in have: continue
        t = pd.read_csv(f"{Q}/ink_tiles/{p}.csv.gz", usecols=["tile_idx", "core", "edge"]); s = tc[tc.slide == i].merge(t, on="tile_idx", how="left"); parts.append(s)
    tcq = pd.concat(parts); tcq[["core", "edge"]] = tcq[["core", "edge"]].fillna(False).astype(bool); tcq["ink"] = tcq.core | tcq.edge
    print("sampled tiles:", len(tcq), "core %.4f  core+edge %.4f" % (tcq.core.mean(), tcq.ink.mean()))
    cl = tcq.groupby("cluster").agg(n=("ink", "size"), core=("core", "sum"), core_or_edge=("ink", "sum")); cl["frac_of_cluster_core"] = cl.core / cl.n; cl["frac_of_cluster_core_or_edge"] = cl.core_or_edge / cl.n
    cl["share_of_all_ink"] = cl.core_or_edge / cl.core_or_edge.sum(); print(cl.round(4)); cl.to_csv(f"{Q}/ink_by_cluster.csv")
    idx = d.index[d.patient.isin(have)]
    out = {}
    for name, flt in [("original", tcq), ("minus_core", tcq[~tcq.core]), ("minus_core_edge", tcq[~tcq.ink])]:
        fr = fractions(flt, n); t = h8a(d.loc[idx], fr.loc[idx]); out[name] = (fr, t); t.to_csv(f"{Q}/H8a_{name}.csv", index=False)
        print("\n==", name, "(n slides %d)" % len(idx)); print(t.round(4).to_string())
    f0, f1 = out["original"][0], out["minus_core_edge"][0]
    pat = pd.DataFrame({"patient": d.patient, "tss": d.tss, "n_sampled": f0.n, "n_after": f1.n, "ink_core_edge_frac": 1 - f1.n / f0.n, **{f"{c}_before": f0[c] for c in CS}, **{f"{c}_after": f1[c] for c in CS}}).loc[idx]
    pat["dC1"] = pat.C1_after - pat.C1_before; pat["dC5"] = pat.C5_after - pat.C5_before; pat.to_csv(f"{Q}/patient_level_before_after.csv", index=False)
    print("\nslides by ink fraction:"); print(pat.sort_values("ink_core_edge_frac", ascending=False).head(12)[["patient", "tss", "ink_core_edge_frac", "C5_before", "C5_after", "C1_before", "C1_after"]].round(3).to_string())
    print("slides with ink_core_edge_frac > 1%:", int((pat.ink_core_edge_frac > 0.01).sum()), " > 5%:", int((pat.ink_core_edge_frac > 0.05).sum()))
    json.dump(dict(n_slides=len(idx), sampled=len(tcq), core=float(tcq.core.mean()), core_or_edge=float(tcq.ink.mean())), open(f"{Q}/ink_overall.json", "w"))


if __name__ == "__main__":
    main()
