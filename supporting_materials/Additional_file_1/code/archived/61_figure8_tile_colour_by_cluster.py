"""Fig 8 A10-Q4 item 5 helper: thumbnail-derived mean colour of every sampled tile (image + cluster labels only; no outcome data)."""
from __future__ import annotations
import os, sys, json, importlib
import numpy as np, pandas as pd, h5py
from PIL import Image, PngImagePlugin
PngImagePlugin.MAX_TEXT_CHUNK = 200 * 1024 * 1024; Image.MAX_IMAGE_PIXELS = None
sys.path.insert(0, os.path.dirname(__file__)); rule = importlib.import_module("57_figure8_ink_rule")
R = "${PROJECT_ROOT}"; Q = f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01"; TH = f"{R}/02_intermediate/fig8_ink"
d = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_pathology_run01/H8_slide_table.csv", usecols=["patient", "file", "tss"]); tc = pd.read_csv(f"{R}/02_intermediate/fig8/tile_clusters.csv.gz"); rows = []
for i, r in enumerate(d.itertuples()):
    meta = json.load(open(f"{TH}/thumb_{r.patient}.json")); ds = meta["ds"]; img = np.asarray(Image.open(f"{TH}/thumb_{r.patient}.png").convert("RGB")); lab = rule.rgb2lab(img)
    with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{r.file}", "r") as h: co = h["coords"][0][:]
    s = tc[tc.slide == i]; p = max(1, int(round(512 / ds))); H, W = lab.shape[:2]
    for _, t in s.iterrows():
        x, y = co[int(t.tile_idx)]; x0, y0 = int(x / ds), int(y / ds); patch = lab[y0:min(H, y0 + p), x0:min(W, x0 + p)].reshape(-1, 3)
        if len(patch) == 0: continue
        tis = patch[patch[:, 0] < 85]; m = tis.mean(0) if len(tis) else patch.mean(0)
        rows.append((i, int(t.tile_idx), int(t.cluster), m[0], m[1], m[2], len(tis) / len(patch)))
    print(r.patient, flush=True)
pd.DataFrame(rows, columns=["slide", "tile_idx", "cluster", "L", "a", "b", "tissue_frac"]).to_csv(f"{Q}/tile_thumb_colour.csv.gz", index=False)
