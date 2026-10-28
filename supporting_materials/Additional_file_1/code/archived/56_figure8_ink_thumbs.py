"""Fig 8 A10-Q4 step 1: low-resolution thumbnails (max side ~2400 px) of all 145 TCGA-PAAD diagnostic slides for ink QC (remote range reads).
No outcome information is used or loaded."""
from __future__ import annotations
import os, json, time
import numpy as np, pandas as pd, tiffslide
from concurrent.futures import ThreadPoolExecutor

R = "${PROJECT_ROOT}"
F8 = f"{R}/02_intermediate/fig8"; OUT = f"{R}/02_intermediate/fig8_ink"; os.makedirs(OUT, exist_ok=True)
LOCAL = os.environ.get("TCGA_PAAD_SVS_ROOT")


def retry(fn, n=8):
    for i in range(n):
        try: return fn()
        except Exception as e:
            print("retry", i, type(e).__name__, flush=True); time.sleep(5 * (i + 1))
    raise RuntimeError("failed")


def one(args):
    key, patient, gdc = args
    fn = f"{OUT}/thumb_{patient}.png"
    if os.path.exists(fn): return patient
    loc = f"{LOCAL}/{key}.svs" if LOCAL else ""; u = loc if loc and os.path.exists(loc) else f"https://api.gdc.cancer.gov/data/{gdc.loc[gdc.key == key, 'file_id'].iloc[0]}"
    s = retry(lambda: tiffslide.TiffSlide(u)); lvl = min(range(len(s.level_dimensions)), key=lambda i: abs(max(s.level_dimensions[i]) - 2400))
    im = retry(lambda: s.read_region((0, 0), lvl, s.level_dimensions[lvl])).convert("RGB")
    im.save(fn); json.dump(dict(level=lvl, ds=float(s.level_downsamples[lvl]), w=im.size[0], h=im.size[1], dims=list(s.dimensions)), open(f"{OUT}/thumb_{patient}.json", "w"))
    print("thumb", patient, lvl, im.size, flush=True); return patient


if __name__ == "__main__":
    gdc = pd.read_csv(f"{F8}/gdc_paad_dx_slides.csv"); gdc["key"] = gdc.file_name.str.replace(".svs", "", regex=False)
    d = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_pathology_run01/H8_slide_table.csv"); d["key"] = d.file.str.replace(".h5", "", regex=False)
    # order: pilot / flagged slides first so QC can start early
    first = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_stain_pilot_run01/pilot_slides.csv").patient.tolist()
    d["o"] = d.patient.map(lambda p: 0 if p in first else 1); d = d.sort_values("o")
    if os.environ.get("REV"): d = d.iloc[::-1]
    with ThreadPoolExecutor(8) as ex: list(ex.map(one, [(r.key, r.patient, gdc) for r in d.itertuples()]))
    print("THUMBS_DONE")
