"""Fig 8 A8: Ridge baseline re-run with the previously validated nested-alpha recipe."""
from __future__ import annotations
import glob, importlib, json, os, sys
import h5py, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(__file__)); mil = importlib.import_module("44_figure8_mil")
R = "${PROJECT_ROOT}"; OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; SEED = 20261012
GRID = {"model__alpha": np.logspace(-2, 4, 13)}


def pipe() -> GridSearchCV:
    return GridSearchCV(Pipeline([("scale", StandardScaler()), ("model", Ridge())]), GRID, scoring="neg_mean_squared_error", cv=KFold(4, shuffle=True, random_state=SEED), refit=True)


def main() -> None:
    d, bags = mil.load_feats(); n = len(d); t = ((d.risk - d.risk.mean()) / d.risk.std()).values; y = (d.group == "High").astype(int).values
    M = np.vstack([b.mean(0).numpy() for b in bags]); oof = np.zeros((5, n)); alphas = []
    for r in range(5):
        for f, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(M, y)):
            g = pipe(); g.cv = KFold(4, shuffle=True, random_state=SEED + r * 100 + f); g.fit(M[tr], t[tr]); oof[r, te] = g.predict(M[te]); alphas.append(g.best_params_["model__alpha"])
    res = dict(alpha_median=float(np.median(alphas)), alpha_values=[float(a) for a in sorted(set(alphas))],
               pearson_by_rep=[float(pearsonr(o, t)[0]) for o in oof], spearman_by_rep=[float(spearmanr(o, t)[0]) for o in oof], auc_by_rep=[float(roc_auc_score(y, o)) for o in oof])
    for k in ("pearson", "spearman", "auc"): res[k + "_mean"] = float(np.mean(res[k + "_by_rep"]))
    # paired bootstrap vs binary ABMIL (rep-mean OOF)
    ob = pd.read_csv(f"{RES}/H8b_oof_predictions.csv"); rm = oof.mean(0); rng = np.random.default_rng(SEED); diff = []
    for _ in range(2000):
        i = rng.integers(0, n, n)
        if y[i].min() != y[i].max(): diff.append(roc_auc_score(y[i], rm[i]) - roc_auc_score(y[i], ob.mil_oof_mean.values[i]))
    res["auc_repmean_ridge"] = float(roc_auc_score(y, rm)); res["auc_repmean_ABMIL_binary"] = float(roc_auc_score(y, ob.mil_oof_mean)); res["paired_diff_ridge_minus_ABMIL_CI95"] = np.percentile(diff, [2.5, 97.5]).tolist()
    # CPTAC transfer
    rna = pd.read_csv(f"{R}/01_raw/CPTAC_PDAC_LinkedOmics/mRNA_RSEM_UQ_log2_Tumor.cct", sep="\t", index_col=0); g2 = rna.loc[["GALNT14", "GALNT16"]].T; z = (g2 - g2.mean()) / g2.std()
    gi = -0.344 * z.GALNT16 + 0.202 * z.GALNT14; tc = ((gi - gi.mean()) / gi.std()).values; pids = list(gi.index)
    X = np.load(f"{OUT}/feats_f16.npy").astype(np.float32); mu = X.mean(0); sd = X.std(0) + 1e-6; del X
    r_ = np.random.default_rng(SEED); Mc = []
    for pid in pids:
        F = []
        for f in sorted(glob.glob(f"{R}/01_raw/uni2h_features/cptac_pda/{pid}-*.h5")):
            with h5py.File(f, "r") as h:
                k = h["features"].shape[1]; idx = np.sort(r_.choice(k, min(2000, k), replace=False)); F.append(h["features"][0][idx].astype(np.float32))
        F = np.vstack(F); sel = r_.choice(len(F), min(2000, len(F)), replace=False); Mc.append(((F[sel] - mu) / sd).mean(0))
    Mc = np.vstack(Mc); g = pipe().fit(M, t); pc = g.predict(Mc); yc = (tc > np.median(tc)).astype(int)
    pr = np.random.default_rng(SEED); null = np.array([spearmanr(pr.permutation(tc), pc)[0] for _ in range(2000)])
    res["CPTAC"] = dict(alpha=float(g.best_params_["model__alpha"]), pearson=float(pearsonr(pc, tc)[0]), spearman=float(spearmanr(pc, tc)[0]), auc=float(roc_auc_score(yc, pc)),
                        perm_p=float((np.sum(null >= spearmanr(pc, tc)[0]) + 1) / 2001))
    pd.DataFrame({"patient": d.patient, "y": y, "giscore": d.risk, "ridge_oof_mean": rm}).to_csv(f"{RES}/H8g_ridge_oof_TCGA.csv", index=False)
    pd.DataFrame({"patient": pids, "giscore": gi.values, "y": yc, "ridge": pc}).to_csv(f"{RES}/H8g_ridge_pred_CPTAC.csv", index=False)
    # exploratory diagnostic (no refit): tuned ridge applied to ONE slide per CPTAC patient (lowest slide index) vs all slides
    r2 = np.random.default_rng(SEED); M1 = []
    for pid in pids:
        f = sorted(glob.glob(f"{R}/01_raw/uni2h_features/cptac_pda/{pid}-*.h5"))[0]
        with h5py.File(f, "r") as h:
            k = h["features"].shape[1]; idx = np.sort(r2.choice(k, min(2000, k), replace=False)); F = h["features"][0][idx].astype(np.float32)
        M1.append(((F - mu) / sd).mean(0))
    p1 = g.predict(np.vstack(M1)); res["CPTAC_one_slide"] = dict(pearson=float(pearsonr(p1, tc)[0]), spearman=float(spearmanr(p1, tc)[0]), auc=float(roc_auc_score(yc, p1)))
    json.dump(res, open(f"{RES}/H8g_ridge_audit.json", "w"), indent=1); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
