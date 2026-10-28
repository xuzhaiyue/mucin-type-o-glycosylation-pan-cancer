"""Fig 8 step 5: assign ALL tiles of the 4 display slides (+ attention) to the frozen k-means clusters.
PCA is linear: recover the standardise+PCA map by least squares from the cached sampled tiles (tiles x 1536 -> Z32), check agreement, then assign nearest centroid."""
from __future__ import annotations
import h5py, numpy as np, pandas as pd

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"


def main() -> None:
    X = np.load(f"{OUT}/feats_f16.npy").astype(np.float32); Z = np.load(f"{OUT}/Z32.npy")
    lab = pd.read_csv(f"{OUT}/tile_clusters.csv.gz").cluster.values - 1
    A = np.hstack([X, np.ones((len(X), 1), np.float32)])
    W, *_ = np.linalg.lstsq(A[::3], Z[::3], rcond=None)
    Zh = A @ W; print("map R2", 1 - ((Zh - Z) ** 2).sum() / ((Z - Z.mean(0)) ** 2).sum())
    C = np.vstack([Z[lab == c].mean(0) for c in range(6)])
    pred = ((Zh[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
    print("agreement with frozen labels", (pred == lab).mean())
    np.save(f"{OUT}/map_W.npy", W); np.save(f"{OUT}/centroids.npy", C)
    d = pd.read_csv(f"{RES}/H8_slide_table.csv"); disp = pd.read_csv(f"{RES}/H8_display_slides.csv")
    for k in disp[disp.role == "display"].key:
        f = k + ".h5"
        with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{f}", "r") as h:
            F = h["features"][0][:].astype(np.float32); co = h["coords"][0][:]
        # use the same global mean/sd as cached sample
        z = np.hstack([F, np.ones((len(F), 1), np.float32)]) @ W
        c = ((z[:, None, :] - C[None]) ** 2).sum(-1).argmin(1) + 1
        pd.DataFrame({"x": co[:, 0], "y": co[:, 1], "cluster": c}).to_csv(f"{OUT}/alltiles_{k[:12]}.csv.gz", index=False)
        print(k[:12], np.bincount(c, minlength=7)[1:] / len(c))


if __name__ == "__main__":
    main()
