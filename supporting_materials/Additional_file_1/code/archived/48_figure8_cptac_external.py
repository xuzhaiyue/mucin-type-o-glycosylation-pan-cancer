"""Fig 8 external check in CPTAC-PDA (spec A6)."""
from __future__ import annotations
import glob, json, os, sys
import h5py, numpy as np, pandas as pd, statsmodels.formula.api as smf, torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(__file__))
import importlib; mil = importlib.import_module("44_figure8_mil")

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; SEED = 20261012


def main() -> None:
    rna = pd.read_csv(f"{R}/01_raw/CPTAC_PDAC_LinkedOmics/mRNA_RSEM_UQ_log2_Tumor.cct", sep="\t", index_col=0)
    g = rna.loc[["GALNT14", "GALNT16"]].T
    z = (g - g.mean()) / g.std()
    gi = (-0.344 * z.GALNT16 + 0.202 * z.GALNT14).rename("giscore")
    print("CPTAC GIScore n", len(gi), "corr(G14,G16)", g.corr().iloc[0, 1])
    W = np.load(f"{OUT}/map_W.npy"); C = np.load(f"{OUT}/centroids.npy")
    rng = np.random.default_rng(SEED)
    # TCGA standardisation stats (same as 44) on cached sample
    X = np.load(f"{OUT}/feats_f16.npy").astype(np.float32); mu = X.mean(0); sd = X.std(0) + 1e-6; del X
    bag, frac, ntile, rows = {}, {}, {}, []
    for pid in gi.index:
        fs = sorted(glob.glob(f"{R}/01_raw/uni2h_features/cptac_pda/{pid}-*.h5")); F = []
        for f in fs:
            with h5py.File(f, "r") as h:
                n = h["features"].shape[1]; idx = np.sort(rng.choice(n, min(2000, n), replace=False)); F.append(h["features"][0][idx].astype(np.float32)); ntile.setdefault(pid, 0); ntile[pid] += n
        F = np.vstack(F); Zh = np.hstack([F, np.ones((len(F), 1), np.float32)]) @ W
        cl = ((Zh[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
        frac[pid] = np.bincount(cl, minlength=6) / len(cl)
        sel = rng.choice(len(F), min(2000, len(F)), replace=False); bag[pid] = torch.from_numpy((F[sel] - mu) / sd)
        rows.append(dict(patient=pid, n_slides=len(fs), n_tiles=ntile[pid]))
    d = pd.DataFrame(rows).set_index("patient").join(gi); cs = [f"C{i}" for i in range(1, 7)]
    d[cs] = pd.DataFrame(frac).T.values; d["log_tiles"] = np.log10(d.n_tiles)
    d["gi"] = (d.giscore - d.giscore.mean()) / d.giscore.std()
    out = []
    for c in cs:
        d["y"] = np.log(d[c].clip(lower=1e-3)) - np.log(d[cs].clip(lower=1e-3)).mean(1)
        m = smf.ols("y ~ gi + log_tiles", d).fit(cov_type="HC3"); rho, p = spearmanr(d.giscore, d[c])
        out.append(dict(cluster=c, mean_frac=d[c].mean(), rho=rho, p_rho=p, beta=m.params["gi"], lo=m.conf_int().loc["gi", 0], hi=m.conf_int().loc["gi", 1], p_two=m.pvalues["gi"],
                        p_onesided_neg=(m.pvalues["gi"] / 2 if m.params["gi"] < 0 else 1 - m.pvalues["gi"] / 2)))
    t = pd.DataFrame(out); t["bh_two"] = multipletests(t.p_two, method="fdr_bh")[1]; t.to_csv(f"{RES}/H8d_CPTAC_cluster_vs_GIScore.csv", index=False)
    pd.set_option("display.width", 200); print(t.round(4).to_string())
    d.to_csv(f"{RES}/H8d_CPTAC_patient_table.csv")
    # H8e: train on all TCGA, apply to CPTAC
    tc, bags = mil.load_feats(); y = (tc.group == "High").astype(int).values; allidx = np.arange(len(tc))
    pids = list(d.index); yc = (d.giscore > d.giscore.median()).astype(int).values
    P = np.zeros((5, len(pids)))
    for s in range(5):
        m = mil.ABMIL().to(mil.dev); torch.manual_seed(SEED + s); r = np.random.default_rng(SEED + s)
        opt = torch.optim.Adam(m.parameters(), lr=2e-4, weight_decay=1e-4); lf = torch.nn.BCEWithLogitsLoss()
        for _ in range(30):
            m.train()
            for i in r.permutation(allidx):
                o = m(bags[i].to(mil.dev)); loss = lf(o, torch.tensor([float(y[i])], device=mil.dev)); opt.zero_grad(); loss.backward(); opt.step()
        m.eval()
        with torch.no_grad(): P[s] = [torch.sigmoid(m(bag[p].to(mil.dev))).item() for p in pids]
        print("seed", s, roc_auc_score(yc, P[s]), flush=True)
    pm = P.mean(0); auc = roc_auc_score(yc, pm)
    M = np.vstack([b.mean(0).numpy() for b in bags]); Mc = np.vstack([bag[p].mean(0).numpy() for p in pids])
    sc = StandardScaler().fit(M); lr = LogisticRegression(C=0.01, max_iter=2000).fit(sc.transform(M), y); pb = lr.predict_proba(sc.transform(Mc))[:, 1]; aucb = roc_auc_score(yc, pb)
    pr = np.random.default_rng(SEED); null = np.array([roc_auc_score(pr.permutation(yc), pm) for _ in range(1000)]); nullb = np.array([roc_auc_score(pr.permutation(yc), pb) for _ in range(1000)])
    res = dict(n=len(pids), mil_auc=float(auc), mil_perm_p=float((np.sum(null >= auc) + 1) / 1001), base_auc=float(aucb), base_perm_p=float((np.sum(nullb >= aucb) + 1) / 1001),
               spearman_mil_vs_giscore=float(spearmanr(pm, d.giscore)[0]))
    json.dump(res, open(f"{RES}/H8e_CPTAC_transfer.json", "w"), indent=1); print(res)
    pd.DataFrame({"patient": pids, "giscore": d.giscore.values, "y": yc, "mil": pm, "base": pb}).to_csv(f"{RES}/H8e_CPTAC_predictions.csv", index=False)


if __name__ == "__main__":
    main()
