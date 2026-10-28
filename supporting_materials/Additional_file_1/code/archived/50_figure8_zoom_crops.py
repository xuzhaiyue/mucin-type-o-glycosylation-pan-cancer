"""Fig 8: zoom window per display slide (rule in 49_figure8_attention.py) read remotely at level 1; saves crop + window geometry."""
from __future__ import annotations
import importlib, json, os, sys
import numpy as np, pandas as pd, tiffslide
from scipy.ndimage import uniform_filter
sys.path.insert(0, os.path.dirname(__file__)); s45 = importlib.import_module("45_figure8_slide_images")

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"
T, NW = 512, 10


def main() -> None:
    gdc = pd.read_csv(f"{OUT}/gdc_paad_dx_slides.csv"); gdc["key"] = gdc.file_name.str.replace(".svs", "", regex=False)
    pred = pd.read_csv(f"{RES}/H8_attention_loo_predictions.csv")
    # one High + one Low patient: most extreme GIScore among LOO-correctly-classified display slides
    sel = [pred[(pred.group == g) & (pred.correct == 1)].sort_values("risk", ascending=(g == "Low")).iloc[0] for g in ("High", "Low")]
    info = {}
    for r in sel:
        k = r.key; a = pd.read_csv(f"{OUT}/attn_{k[:12]}.csv.gz")
        gx = (a.x // T).astype(int); gy = (a.y // T).astype(int)
        G = np.full((gy.max() + 1, gx.max() + 1), np.nan); G[gy, gx] = a["rank"].values
        cov = (~np.isnan(G)).astype(float); Gf = np.nan_to_num(G, nan=0.0)
        covw = uniform_filter(cov, NW, mode="constant"); meanw = uniform_filter(Gf, NW, mode="constant") / np.maximum(covw, 1e-6)
        meanw[covw < 0.80] = -1   # rule A9: window must be >= 80% covered by released tiles; then highest mean attention
        cy, cx = np.unravel_index(np.argmax(meanw), meanw.shape)
        x0 = int(max(0, cx - NW // 2)) * T; y0 = int(max(0, cy - NW // 2)) * T
        print("coverage of chosen window", round(float(covw[cy, cx]), 2), "mean rank", round(float(meanw[cy, cx]), 2), flush=True)
        s = s45.retry(lambda: tiffslide.TiffSlide(s45.src(k, gdc)))
        lvl = 1 if len(s.level_dimensions) > 1 else 0; ds = float(s.level_downsamples[lvl]); sz = int(NW * T / ds)
        im = s45.retry(lambda: s.read_region((x0, y0), lvl, (sz, sz))).convert("RGB"); im.save(f"{OUT}/zoom_{k[:12]}.png")
        W = np.full((NW, NW), np.nan)
        for j in range(NW):
            for i in range(NW):
                m = (gx == x0 // T + i) & (gy == y0 // T + j)
                if m.any(): W[j, i] = a["rank"].values[m][0]
        np.save(f"{OUT}/zoomattn_{k[:12]}.npy", W)
        info[k[:12]] = dict(x0=x0, y0=y0, size=NW * T, group=r.group, risk=float(r.risk), prob=float(r.loo_prob_high), level=lvl, ds=ds)
        print(k[:12], info[k[:12]], flush=True)
    json.dump(info, open(f"{OUT}/zoom_info.json", "w"), indent=1)
    print("ZOOM_DONE")


if __name__ == "__main__":
    main()
