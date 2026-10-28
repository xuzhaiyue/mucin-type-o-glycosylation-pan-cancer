"""Fig 8 A10-Q4 step 4: re-run the frozen ABMIL/mean-pool analyses on bags WITHOUT ink (core+edge) tiles. Same seeds, same standardisation constants, same folds.
POST HOC QC sensitivity analysis; primary results unchanged."""
from __future__ import annotations
import os, sys, json, importlib
import numpy as np, pandas as pd, torch
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(__file__)); mil = importlib.import_module("44_figure8_mil")
R = mil.R; OUT = mil.OUT; RES = mil.RES; Q = f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01"


def main() -> None:
    d = pd.read_csv(f"{RES}/H8_slide_table.csv"); n = len(d); y = (d.group == "High").astype(int).values
    X = np.load(f"{OUT}/feats_f16.npy"); sid = np.load(f"{OUT}/slide_id.npy"); tix = np.load(f"{OUT}/tile_idx.npy")
    mu = X.astype(np.float32).mean(0); sd = X.astype(np.float32).std(0) + 1e-6
    keep = np.ones(len(sid), bool)
    for i, p in enumerate(d.patient):
        t = pd.read_csv(f"{Q}/ink_tiles/{p}.csv.gz", usecols=["tile_idx", "core", "edge"]).set_index("tile_idx"); m = np.where(sid == i)[0]
        ink = (t.core | t.edge).reindex(tix[m]).fillna(False).values.astype(bool); keep[m[ink]] = False
    print("tiles kept", keep.sum(), "of", len(keep), "removed %.4f" % (1 - keep.mean()))
    bags = [torch.from_numpy((X[(sid == i) & keep].astype(np.float32) - mu) / sd) for i in range(n)]
    res = dict(tiles_removed_frac=float(1 - keep.mean()))
    top = d.tss_grp.map({"IB": 0, "HZ": 1, "FB": 2}).fillna(3).astype(int).values
    o3, f3 = mil.cv(bags, top, top, reps=1, k=4, tag="TSS"); res["TSS_macroAUC_fold_mean"] = float(f3.mean())
    M = np.vstack([b.mean(0).numpy() for b in bags]); bo, bf = mil.meanpool_cv(M, y, y); res["base_fold_auc_mean"] = float(bf.mean()); res["base_pooled_auc_by_rep"] = [roc_auc_score(y, o) for o in bo]
    oofm, fa = mil.cv(bags, y, y, reps=5, epochs=30, tag="MIL"); pooled = [roc_auc_score(y, o) for o in oofm]
    res.update(MIL_fold_auc_mean=float(fa.mean()), MIL_pooled_auc_by_rep=pooled, MIL_pooled_auc_mean=float(np.mean(pooled)))
    g = (d.gradeg == "G3-4").astype(int).values; og, fg = mil.cv(bags, g, g, reps=1, tag="grade"); res["grade_fold_auc_mean"] = float(fg.mean())
    os_, fs_ = mil.cv(bags, y, d.tss_grp.values, reps=1, tag="GI_TSSfolds"); res["GI_TSSstratified_fold_auc_mean"] = float(fs_.mean())
    json.dump(res, open(f"{Q}/H8_ink_filtered_confound_auc.json", "w"), indent=1); print(res)
    print("ORIGINAL:", json.load(open(f"{RES}/H8c_confound_auc.json")), {k: v for k, v in json.load(open(f"{RES}/H8b_cv_auc.json")).items() if "mean" in k})


if __name__ == "__main__":
    main()
