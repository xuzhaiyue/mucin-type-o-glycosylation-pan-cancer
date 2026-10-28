"""Fig 8 step 1 (spec figure8_spec_v1.md): sample UNI2-h tiles, PCA(32), k-means(6), t-SNE, per-slide fractions."""
from __future__ import annotations
import glob, os, re, json
import h5py, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, adjusted_rand_score

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; os.makedirs(OUT, exist_ok=True)
RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; os.makedirs(RES, exist_ok=True)
SEED, NT = 20261012, 2000


def pick_slides() -> pd.DataFrame:
    sc = pd.read_csv(f"{R}/04_results/2026-10-12_Figure4_GALNT_run01/training_scores.csv")
    sc["patient"] = sc["sample"].str[:12]
    cl = pd.read_csv(f"{R}/04_results/2026-10-12_Figure4_GALNT_run01/H_clinical_table_TCGA_samples.csv")[["patient", "gradeg", "stageg"]]
    fs = pd.DataFrame({"file": sorted(os.path.basename(f) for f in glob.glob(f"{R}/01_raw/uni2h_features/tcga_paad/*.h5"))})
    fs["patient"] = fs.file.str[:12]
    fs["dx"] = fs.file.str.extract(r"-DX(\d)")[0].astype(int)
    fs = fs.sort_values(["patient", "dx"]).drop_duplicates("patient")
    d = sc.merge(fs, on="patient", how="left").merge(cl, on="patient", how="left")
    d["tss"] = d.patient.str[5:7]
    miss = d[d.file.isna()].patient.tolist()
    json.dump(miss, open(f"{RES}/excluded_no_slide.json", "w"))
    d = d.dropna(subset=["file"]).reset_index(drop=True)
    vc = d.tss.value_counts()
    d["tss_grp"] = np.where(d.tss.map(vc) >= 5, d.tss, "other")
    return d


def main() -> None:
    d = pick_slides()
    print("slides:", len(d), "excluded:", open(f"{RES}/excluded_no_slide.json").read())
    rng = np.random.default_rng(SEED)
    X, sid, tix, ntile = [], [], [], []
    for i, r in d.iterrows():
        with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{r.file}", "r") as h:
            f = h["features"]; n = f.shape[1]
            idx = np.sort(rng.choice(n, min(NT, n), replace=False))
            X.append(f[0][idx].astype(np.float32)); sid += [i] * len(idx); tix += idx.tolist(); ntile.append(n)
    X = np.vstack(X); sid = np.array(sid); tix = np.array(tix)
    d["n_tiles_total"] = ntile; d["n_tiles_used"] = np.bincount(sid)
    print("tiles", X.shape)
    X = (X - X.mean(0)) / (X.std(0) + 1e-6)
    pca = PCA(64, random_state=SEED).fit(X)
    Z64 = pca.transform(X).astype(np.float32); Z = Z64[:, :32]
    print("PCA var 32/64:", pca.explained_variance_ratio_[:32].sum(), pca.explained_variance_ratio_.sum())
    km = KMeans(6, n_init=10, random_state=SEED).fit(Z); lab = km.labels_
    # sensitivity
    sens = {}
    km64 = KMeans(6, n_init=10, random_state=SEED).fit(Z64); sens["ARI_k6_pca32_vs_pca64"] = adjusted_rand_score(lab, km64.labels_)
    for k in (5, 7):
        sens[f"ARI_k6_vs_k{k}"] = adjusted_rand_score(lab, KMeans(k, n_init=10, random_state=SEED).fit_predict(Z))
    sub = rng.choice(len(Z), 50000, replace=False)
    sil = {}
    for k in range(3, 11):
        sil[k] = float(silhouette_score(Z[sub], KMeans(k, n_init=5, random_state=SEED).fit_predict(Z[sub])))
    pd.Series(sil, name="silhouette").to_csv(f"{RES}/S8_silhouette_k3_10.csv")
    json.dump(sens, open(f"{RES}/S8_sensitivity.json", "w"), indent=1); print(sens, sil)
    # fractions
    t = pd.DataFrame({"slide": sid, "cl": lab})
    fr = pd.crosstab(t.slide, t.cl).reindex(range(len(d)), fill_value=0)
    frac = fr.div(fr.sum(1), axis=0); frac.columns = [f"C{c+1}" for c in frac.columns]
    pd.concat([d, frac], axis=1).to_csv(f"{RES}/H8_slide_table.csv", index=False)
    # tile table + t-SNE (stratified: equal number per cluster, up to 5000)
    tl = pd.DataFrame({"slide": sid, "tile_idx": tix, "cluster": lab + 1})
    tl.to_csv(f"{OUT}/tile_clusters.csv.gz", index=False)
    ts = np.concatenate([rng.choice(np.where(lab == c)[0], min(5000, (lab == c).sum()), replace=False) for c in range(6)])
    E = TSNE(2, perplexity=30, random_state=SEED, init="pca", n_jobs=4).fit_transform(Z[ts])
    pd.DataFrame({"tsne1": E[:, 0], "tsne2": E[:, 1], "cluster": lab[ts] + 1, "slide": sid[ts], "tile_idx": tix[ts]}).to_csv(f"{RES}/H8_tsne.csv", index=False)
    # centroid-nearest tiles from distinct slides
    rep = []
    for c in range(6):
        ii = np.where(lab == c)[0]
        dd = ((Z[ii] - km.cluster_centers_[c]) ** 2).sum(1); o = ii[np.argsort(dd)]
        seen, k = set(), 0
        for j in o:
            if sid[j] in seen: continue
            seen.add(sid[j]); rep.append(dict(cluster=c + 1, slide=int(sid[j]), file=d.file[sid[j]], tile_idx=int(tix[j]), patient=d.patient[sid[j]]))
            k += 1
            if k == 12: break
    pd.DataFrame(rep).to_csv(f"{RES}/H8_representative_tiles.csv", index=False)
    np.save(f"{OUT}/Z32.npy", Z); np.save(f"{OUT}/slide_id.npy", sid); np.save(f"{OUT}/tile_idx.npy", tix)
    print("cluster sizes", np.bincount(lab), "done")


if __name__ == "__main__":
    main()
