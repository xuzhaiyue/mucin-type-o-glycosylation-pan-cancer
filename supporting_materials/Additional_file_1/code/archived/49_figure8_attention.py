"""Fig 8 attention maps (replaces Grad-CAM; spec A5 revised): leave-one-slide-out ABMIL, attention for all tiles of the display slides.
Zoom window rule (fixed before viewing): 6x6-tile window (3072 px level-0) centred on the tile with the highest 5x5-neighbourhood-mean attention."""
from __future__ import annotations
import importlib, os, sys
import h5py, numpy as np, pandas as pd, torch, tiffslide
from scipy.ndimage import uniform_filter
from PIL import Image, PngImagePlugin
PngImagePlugin.MAX_TEXT_CHUNK = 200 * 1024 * 1024
sys.path.insert(0, os.path.dirname(__file__)); mil = importlib.import_module("44_figure8_mil")

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"; SEED = 20261012


def main() -> None:
    d, bags = mil.load_feats(); y = (d.group == "High").astype(int).values
    d["key"] = d.file.str.replace(".h5", "", regex=False)
    X = np.load(f"{OUT}/feats_f16.npy").astype(np.float32); mu = X.mean(0); sd = X.std(0) + 1e-6; del X
    disp = pd.read_csv(f"{RES}/H8_display_slides.csv"); dk = disp[disp.role == "display"].key.tolist()
    gdc = pd.read_csv(f"{OUT}/gdc_paad_dx_slides.csv"); gdc["key"] = gdc.file_name.str.replace(".svs", "", regex=False)
    rows = []
    for k in dk:
        i = int(d.index[d.key == k][0]); tr = np.array([j for j in range(len(d)) if j != i])
        pr, m = mil.fit_predict(bags, y, tr, np.array([i]), seed=SEED + i)
        with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{k}.h5", "r") as h:
            F = h["features"][0][:].astype(np.float32); co = h["coords"][0][:]
        with torch.no_grad():
            _, w = m(torch.from_numpy((F - mu) / sd).to(mil.dev), ret_att=True)
        w = w.cpu().numpy(); rk = pd.Series(w).rank(pct=True).values
        pd.DataFrame({"x": co[:, 0], "y": co[:, 1], "attn": w, "rank": rk}).to_csv(f"{OUT}/attn_{k[:12]}.csv.gz", index=False)
        rows.append(dict(key=k, patient=d.patient[i], risk=d.risk[i], group=d.group[i], loo_prob_high=float(pr[0]), correct=int((pr[0] > 0.5) == (y[i] == 1))))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{RES}/H8_attention_loo_predictions.csv", index=False)


if __name__ == "__main__":
    main()
