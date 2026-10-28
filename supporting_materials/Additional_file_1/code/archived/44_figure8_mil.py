"""Fig 8 step 3: H8b/H8c ABMIL on frozen UNI2-h features (spec figure8_spec_v1.md).
Amendment A1 (before any result): the MIL uses the fixed 2,000-tile set per slide from step 1 (not re-sampled per epoch)."""
from __future__ import annotations
import os, sys, time, json
import h5py, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"
SEED = 20261012
dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
torch.set_num_threads(4)


def load_feats() -> tuple[pd.DataFrame, list[torch.Tensor]]:
    d = pd.read_csv(f"{RES}/H8_slide_table.csv")
    cache = f"{OUT}/feats_f16.npy"
    sid = np.load(f"{OUT}/slide_id.npy"); tix = np.load(f"{OUT}/tile_idx.npy")
    if os.path.exists(cache):
        X = np.load(cache)
    else:
        X = np.zeros((len(sid), 1536), np.float16)
        for i, r in d.iterrows():
            m = np.where(sid == i)[0]
            with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{r.file}", "r") as h:
                X[m] = h["features"][0][tix[m]].astype(np.float16)
        np.save(cache, X)
    mu = X.astype(np.float32).mean(0); sd = X.astype(np.float32).std(0) + 1e-6
    bags = [torch.from_numpy(((X[sid == i].astype(np.float32) - mu) / sd)) for i in range(len(d))]
    return d, bags


class ABMIL(nn.Module):
    def __init__(self, k: int = 1, h: int = 256, p: float = 0.25) -> None:
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(1536, h), nn.ReLU(), nn.Dropout(p))
        self.a = nn.Sequential(nn.Linear(h, 128), nn.Tanh(), nn.Linear(128, 1))
        self.out = nn.Linear(h, k)

    def forward(self, x: torch.Tensor, ret_att: bool = False):
        z = self.fc(x); w = torch.softmax(self.a(z).squeeze(-1), 0)
        o = self.out((w[:, None] * z).sum(0))
        return (o, w) if ret_att else o


def fit_predict(bags, y, tr, te, k=1, epochs=30, seed=0, lr=2e-4):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    m = ABMIL(k).to(dev); opt = torch.optim.Adam(m.parameters(), lr=lr, weight_decay=1e-4)
    lossf = nn.BCEWithLogitsLoss() if k == 1 else nn.CrossEntropyLoss()
    for _ in range(epochs):
        m.train()
        for i in rng.permutation(tr):
            x = bags[i].to(dev); o = m(x)
            t = torch.tensor([float(y[i])], device=dev) if k == 1 else torch.tensor(int(y[i]), device=dev)
            loss = lossf(o, t) if k == 1 else lossf(o[None], t[None])
            opt.zero_grad(); loss.backward(); opt.step()
    m.eval(); pr = []
    with torch.no_grad():
        for i in te:
            o = m(bags[i].to(dev)); pr.append(torch.sigmoid(o).item() if k == 1 else torch.softmax(o, 0).cpu().numpy())
    return np.array(pr), m


def cv(bags, y, strat, reps=5, k=1, epochs=30, tag=""):
    n = len(y); oof = np.zeros((reps, n) if k == 1 else (reps, n, k)); fa = []
    for r in range(reps):
        for f, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(np.zeros(n), strat)):
            pr, _ = fit_predict(bags, y, tr, te, k, epochs, seed=SEED + r * 10 + f); oof[r, te] = pr
            fa.append(roc_auc_score(y[te], pr) if k == 1 else roc_auc_score(np.eye(k)[y[te]][:, np.unique(y[te])], pr[:, np.unique(y[te])], average="macro"))
        print(tag, "rep", r, time.strftime("%H:%M:%S"), flush=True)
    return oof, np.array(fa)


def meanpool_cv(M, y, strat, reps=5, seed_off=0):
    n = len(y); oof = np.zeros((reps, n)); fa = []
    for r in range(reps):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r + seed_off).split(M, strat):
            sc = StandardScaler().fit(M[tr]); lr = LogisticRegression(C=0.01, max_iter=2000).fit(sc.transform(M[tr]), y[tr])
            p = lr.predict_proba(sc.transform(M[te]))[:, 1]; oof[r, te] = p; fa.append(roc_auc_score(y[te], p))
    return oof, np.array(fa)


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    d, bags = load_feats(); n = len(d)
    y = (d.group == "High").astype(int).values
    if mode == "main":
        M = np.vstack([b.mean(0).numpy() for b in bags])
        oofm, fa2 = cv(bags, y, y, reps=5, epochs=30, tag="MIL")
        pooled = [roc_auc_score(y, o) for o in oofm]
        bo, bf = meanpool_cv(M, y, y)
        res = dict(MIL_fold_auc_mean=float(fa2.mean()), MIL_fold_auc_sd=float(fa2.std()), MIL_pooled_auc_by_rep=pooled, MIL_pooled_auc_mean=float(np.mean(pooled)),
                   base_fold_auc_mean=float(bf.mean()), base_pooled_auc_by_rep=[roc_auc_score(y, o) for o in bo])
        json.dump(res, open(f"{RES}/H8b_cv_auc.json", "w"), indent=1); print(res)
        pd.DataFrame({"patient": d.patient, "y": y, "mil_oof_mean": oofm.mean(0), "base_oof_mean": bo.mean(0), "risk": d.risk}).to_csv(f"{RES}/H8b_oof_predictions.csv", index=False)
        # H8c: TSS / grade
        top = d.tss_grp.map({"IB": 0, "HZ": 1, "FB": 2}).fillna(3).astype(int).values
        o3, f3 = cv(bags, top, top, reps=1, k=4, tag="TSS")
        g = (d.gradeg == "G3-4").astype(int).values
        og, fg = cv(bags, g, g, reps=1, tag="grade")
        # TSS-stratified folds for GIScore
        os_, fs_ = cv(bags, y, d.tss_grp.values, reps=1, tag="GI_TSSfolds")
        # leave-one-site-out (3 biggest sites)
        lso = {}
        for s in ("IB", "HZ", "FB"):
            te = np.where(d.tss_grp == s)[0]; tr = np.where(d.tss_grp != s)[0]
            pr, _ = fit_predict(bags, y, tr, te, seed=SEED); lso[s] = dict(n=len(te), auc=float(roc_auc_score(y[te], pr)))
        res2 = dict(TSS_macroAUC_fold_mean=float(f3.mean()), grade_fold_auc_mean=float(fg.mean()), GI_TSSstratified_fold_auc_mean=float(fs_.mean()),
                    GI_TSSstratified_pooled_auc=float(roc_auc_score(y, os_[0])), leave_site_out=lso)
        json.dump(res2, open(f"{RES}/H8c_confound_auc.json", "w"), indent=1); print(res2)
    elif mode == "perm":
        nperm = int(sys.argv[2]); out = []
        for p in range(nperm):
            yp = np.random.default_rng(SEED + 1000 + p).permutation(y)
            o, _ = cv(bags, yp, yp, reps=1, epochs=30, tag=f"perm{p}")
            out.append(roc_auc_score(yp, o[0])); pd.Series(out).to_csv(f"{RES}/H8b_perm_auc.csv", index=False, header=["pooled_auc"])
    elif mode == "baseline_perm":
        M = np.vstack([b.mean(0).numpy() for b in bags]); out = []
        for p in range(200):
            yp = np.random.default_rng(SEED + 1000 + p).permutation(y)
            o, _ = meanpool_cv(M, yp, yp, reps=1); out.append(roc_auc_score(yp, o[0]))
        pd.Series(out).to_csv(f"{RES}/H8b_baseline_perm_auc.csv", index=False, header=["pooled_auc"])


if __name__ == "__main__":
    main()
