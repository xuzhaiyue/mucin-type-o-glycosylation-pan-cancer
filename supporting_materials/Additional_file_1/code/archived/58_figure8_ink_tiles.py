"""Fig 8 A10-Q4 step 2: flag ink tiles for every slide that has a thumbnail (image + coordinates only; NO outcome data are loaded).
Outputs: 04_results/2026-10-12_Figure8_ink_QC_run01/ink_tiles/<patient>.csv.gz, qc_overlay/<patient>.png, ink_slide_summary.csv."""
from __future__ import annotations
import os, sys, json, importlib
import numpy as np, pandas as pd, h5py
from scipy import ndimage as ndi
from PIL import Image, PngImagePlugin
PngImagePlugin.MAX_TEXT_CHUNK = 200 * 1024 * 1024; Image.MAX_IMAGE_PIXELS = None
sys.path.insert(0, os.path.dirname(__file__)); rule = importlib.import_module("57_figure8_ink_rule")

R = "${PROJECT_ROOT}"
TH = f"{R}/02_intermediate/fig8_ink"; OUT = f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01"
for s in ("ink_tiles", "qc_overlay"): os.makedirs(f"{OUT}/{s}", exist_ok=True)
CORE = 0.10


def main() -> None:
    d = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_pathology_run01/H8_slide_table.csv", usecols=["patient", "file", "tss"])   # no outcome columns
    tc = pd.read_csv(f"{R}/02_intermediate/fig8/tile_clusters.csv.gz"); rows = []
    for i, r in enumerate(d.itertuples()):
        tf = f"{TH}/thumb_{r.patient}.png"
        if not os.path.exists(tf) or not os.path.exists(f"{TH}/thumb_{r.patient}.json"): continue
        meta = json.load(open(f"{TH}/thumb_{r.patient}.json")); ds = meta["ds"]; img = np.asarray(Image.open(tf).convert("RGB")); mk = rule.ink_mask(img)
        with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{r.file}", "r") as h: co = h["coords"][0][:]
        p = max(1, int(round(512 / ds))); ii = np.pad(mk.astype(np.int64).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
        x0 = np.clip((co[:, 0] / ds).astype(int), 0, mk.shape[1]); y0 = np.clip((co[:, 1] / ds).astype(int), 0, mk.shape[0])
        x1 = np.clip(x0 + p, 0, mk.shape[1]); y1 = np.clip(y0 + p, 0, mk.shape[0]); area = np.maximum(1, (x1 - x0) * (y1 - y0))
        frac = (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]) / area; core = frac >= CORE
        g = {(a // 512, b // 512) for a, b in co[core]}; edge = np.array([(not c) and any(((a // 512 + dx, b // 512 + dy) in g) for dx in (-1, 0, 1) for dy in (-1, 0, 1)) for (a, b), c in zip(co, core)])
        t = pd.DataFrame(dict(tile_idx=np.arange(len(co)), x=co[:, 0], y=co[:, 1], ink_frac=frac, core=core, edge=edge)); t.to_csv(f"{OUT}/ink_tiles/{r.patient}.csv.gz", index=False)
        s = tc[tc.slide == i]; sm = s.merge(t[["tile_idx", "core", "edge"]], on="tile_idx", how="left")
        ov = img.copy(); ov[ndi.binary_dilation(mk, iterations=2)] = [255, 255, 0]; O = Image.fromarray(ov); sc = 1000 / max(O.size); O.resize((int(O.size[0] * sc), int(O.size[1] * sc))).save(f"{OUT}/qc_overlay/{r.patient}.png")
        rows.append(dict(slide=i, patient=r.patient, tss=r.tss, ds=ds, n_tiles=len(co), core_frac_all=float(core.mean()), edge_frac_all=float(edge.mean()), n_sampled=len(sm), core_frac_sampled=float(sm.core.mean()),
                         edge_frac_sampled=float(sm.edge.mean()), ink_pixel_frac=float(mk.mean()), **{f"core_in_C{c}": int(((sm.cluster == c) & sm.core).sum()) for c in range(1, 7)}, **{f"edge_in_C{c}": int(((sm.cluster == c) & sm.edge).sum()) for c in range(1, 7)}))
        print(r.patient, "core %.4f edge %.4f" % (core.mean(), edge.mean()), flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/ink_slide_summary.csv", index=False)


if __name__ == "__main__":
    main()
