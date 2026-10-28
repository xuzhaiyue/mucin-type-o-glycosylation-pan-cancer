"""A10-Q4 item 5: look at non-ink C5 tiles of TCGA-F2-A44G (and 2 other high-C5 slides): morphology vs stain density vs ink. Image + cluster only."""
from __future__ import annotations
import os, h5py, numpy as np, pandas as pd, tiffslide
from concurrent.futures import ThreadPoolExecutor
from PIL import Image
R = "${PROJECT_ROOT}"; Q = f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01"
d = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_pathology_run01/H8_slide_table.csv", usecols=["patient", "file"]); tc = pd.read_csv(f"{R}/02_intermediate/fig8/tile_clusters.csv.gz")
g = pd.read_csv(f"{R}/02_intermediate/fig8/gdc_paad_dx_slides.csv"); g["key"] = g.file_name.str.replace(".svs", "", regex=False)
jobs = []
for p in ["TCGA-F2-A44G", "TCGA-2J-AABU", "TCGA-US-A77G"]:
    i = int(d.index[d.patient == p][0]); ink = pd.read_csv(f"{Q}/ink_tiles/{p}.csv.gz").set_index("tile_idx"); s = tc[(tc.slide == i) & (tc.cluster == 5)]
    s = s[~(ink.loc[s.tile_idx, "core"] | ink.loc[s.tile_idx, "edge"]).values].sample(4, random_state=20261012)
    for t in s.tile_idx: jobs.append((p, i, int(t)))
def one(j):
    p, i, t = j; r = d.iloc[i]; key = r.file.replace(".h5", "")
    with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{r.file}", "r") as h: x, y = h["coords"][0][t]
    s = tiffslide.TiffSlide(f"https://api.gdc.cancer.gov/data/{g.loc[g.key == key, 'file_id'].iloc[0]}"); im = s.read_region((int(x), int(y)), 0, (512, 512)).convert("RGB").resize((256, 256), Image.LANCZOS); print(p, t, flush=True); return im
with ThreadPoolExecutor(12) as ex: ims = list(ex.map(one, jobs))
M = Image.new("RGB", (4 * 258, 3 * 258), "white")
for k, im in enumerate(ims): M.paste(im, ((k % 4) * 258, (k // 4) * 258))
M.save(f"{Q}/C5_nonink_tiles_F2A44G_2JAABU_USA77G.png"); print("done")
