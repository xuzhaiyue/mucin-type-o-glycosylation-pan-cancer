"""Fig 8 A7: continuous-GIScore regression (ABMIL-MSE and ridge) in TCGA CV and CPTAC transfer."""
from __future__ import annotations
import glob, importlib, json, os, sys
import h5py, numpy as np, pandas as pd, torch, torch.nn as nn
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(__file__)); mil = importlib.import_module("44_figure8_mil")
R = "${PROJECT_ROOT}"; OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; SEED = 20261012
dev = mil.dev


def train_reg(bags, t, tr, seed, epochs=30):
    torch.manual_seed(seed); rng = np.random.default_rng(seed); m = mil.ABMIL().to(dev)
    opt = torch.optim.Adam(m.parameters(), lr=2e-4, weight_decay=1e-4); lf = nn.MSELoss()
    for _ in range(epochs):
        m.train()
        for i in rng.permutation(tr):
            loss = lf(m(bags[i].to(dev)), torch.tensor([float(t[i])], device=dev)); opt.zero_grad(); loss.backward(); opt.step()
    m.eval(); return m


def pred(m, b):
    with torch.no_grad(): return m(b.to(dev)).item()


def main() -> None:
    d, bags = mil.load_feats(); n = len(d); t = ((d.risk - d.risk.mean()) / d.risk.std()).values; y = (d.group == "High").astype(int).values
    M = np.vstack([b.mean(0).numpy() for b in bags]); reps = 5
    oofA = np.zeros((reps, n)); oofR = np.zeros((reps, n))
    for r in range(reps):
        for f, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(M, y)):
            m = train_reg(bags, t, tr, SEED + r * 10 + f); oofA[r, te] = [pred(m, bags[i]) for i in te]
            sc = StandardScaler().fit(M[tr]); oofR[r, te] = Ridge(alpha=1e5).fit(sc.transform(M[tr]), t[tr]).predict(sc.transform(M[te]))
        print("rep", r, flush=True)
    res = {}
    for nm, o in (("ABMIL_reg", oofA), ("ridge_meanpool", oofR)):
        res[nm] = dict(rho_by_rep=[float(spearmanr(x, t)[0]) for x in o], rho_mean=float(np.mean([spearmanr(x, t)[0] for x in o])),
                       auc_by_rep=[float(roc_auc_score(y, x)) for x in o], auc_mean=float(np.mean([roc_auc_score(y, x) for x in o])),
                       rho_of_repmean=float(spearmanr(o.mean(0), t)[0]), auc_of_repmean=float(roc_auc_score(y, o.mean(0))))
    print(json.dumps(res, indent=1), flush=True)
    # CPTAC transfer
    rna = pd.read_csv(f"{R}/01_raw/CPTAC_PDAC_LinkedOmics/mRNA_RSEM_UQ_log2_Tumor.cct", sep="\t", index_col=0); g = rna.loc[["GALNT14", "GALNT16"]].T; z = (g - g.mean()) / g.std()
    gi = (-0.344 * z.GALNT16 + 0.202 * z.GALNT14); tc = ((gi - gi.mean()) / gi.std()).values; pids = list(gi.index)
    X = np.load(f"{OUT}/feats_f16.npy").astype(np.float32); mu = X.mean(0); sd = X.std(0) + 1e-6; del X
    rng = np.random.default_rng(SEED); cb = []
    for pid in pids:
        F = []
        for f in sorted(glob.glob(f"{R}/01_raw/uni2h_features/cptac_pda/{pid}-*.h5")):
            with h5py.File(f, "r") as h:
                nn_ = h["features"].shape[1]; idx = np.sort(rng.choice(nn_, min(2000, nn_), replace=False)); F.append(h["features"][0][idx].astype(np.float32))
        F = np.vstack(F); sel = rng.choice(len(F), min(2000, len(F)), replace=False); cb.append(torch.from_numpy((F[sel] - mu) / sd))
    P = np.zeros((5, len(pids))); allidx = np.arange(n)
    for s in range(5):
        m = train_reg(bags, t, allidx, SEED + s); P[s] = [pred(m, b) for b in cb]
    pm = P.mean(0); yc = (tc > np.median(tc)).astype(int)
    sc = StandardScaler().fit(M); rr = Ridge(alpha=1e5).fit(sc.transform(M), t).predict(sc.transform(np.vstack([b.mean(0).numpy() for b in cb])))
    pr = np.random.default_rng(SEED); null = np.array([spearmanr(pr.permutation(tc), pm)[0] for _ in range(2000)]); nullr = np.array([spearmanr(pr.permutation(tc), rr)[0] for _ in range(2000)])
    ext = dict(ABMIL_reg_rho=float(spearmanr(pm, tc)[0]), ABMIL_reg_auc=float(roc_auc_score(yc, pm)), ABMIL_reg_perm_p=float((np.sum(null >= spearmanr(pm, tc)[0]) + 1) / 2001),
               ridge_rho=float(spearmanr(rr, tc)[0]), ridge_auc=float(roc_auc_score(yc, rr)), ridge_perm_p=float((np.sum(nullr >= spearmanr(rr, tc)[0]) + 1) / 2001))
    res["CPTAC"] = ext; print(ext); json.dump(res, open(f"{RES}/H8f_regression.json", "w"), indent=1)


if __name__ == "__main__":
    main()
