"""Fig 8 step 4: representative tiles + whole-slide cluster maps (remote range-reads from GDC via tiffslide; no full SVS download)."""
from __future__ import annotations
import os, h5py, numpy as np, pandas as pd, tiffslide
from concurrent.futures import ThreadPoolExecutor
from PIL import Image

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"
TD = f"{OUT}/tiles"; os.makedirs(TD, exist_ok=True)
LOCAL = os.environ.get("TCGA_PAAD_SVS_ROOT")


def src(key: str, gdc: pd.DataFrame) -> str:
    p = f"{LOCAL}/{key}.svs" if LOCAL else ""
    return p if p and os.path.exists(p) else f"https://api.gdc.cancer.gov/data/{gdc.loc[gdc.key == key, 'file_id'].iloc[0]}"


def coords(file: str) -> np.ndarray:
    with h5py.File(f"{R}/01_raw/uni2h_features/tcga_paad/{file}", "r") as h:
        return h["coords"][0][:]


def retry(fn, n=6):
    import time
    for i in range(n):
        try: return fn()
        except Exception as e:
            print("retry", i, type(e).__name__, flush=True); time.sleep(5 * (i + 1))
    raise RuntimeError("failed after retries")


def tiles_for_slide(args):
    key, rows, gdc = args
    s = retry(lambda: tiffslide.TiffSlide(src(key, gdc))); co = coords(key + ".h5")
    for _, r in rows.iterrows():
        fn = f"{TD}/C{r.cluster}_{r.patient}_{r.tile_idx}.png"
        if os.path.exists(fn): continue
        x, y = co[r.tile_idx]
        im = retry(lambda: s.read_region((int(x), int(y)), 0, (512, 512)))
        im.convert("RGB").resize((256, 256), Image.LANCZOS).save(fn)
    return key


def thumb(args):
    key, gdc = args
    fn = f"{OUT}/thumb_{key[:12]}.png"
    if os.path.exists(fn): return key
    s = retry(lambda: tiffslide.TiffSlide(src(key, gdc)))
    lvl = min(range(len(s.level_dimensions)), key=lambda i: abs(max(s.level_dimensions[i]) - 2500))
    im = retry(lambda: s.read_region((0, 0), lvl, s.level_dimensions[lvl])).convert("RGB"); im.save(fn)
    pd.Series(dict(level=lvl, ds=float(s.level_downsamples[lvl]), w=im.size[0], h=im.size[1])).to_json(f"{OUT}/thumb_{key[:12]}.json")
    return key


def main() -> None:
    gdc = pd.read_csv(f"{OUT}/gdc_paad_dx_slides.csv"); gdc["key"] = gdc.file_name.str.replace(".svs", "", regex=False)
    rep = pd.read_csv(f"{RES}/H8_representative_tiles.csv").groupby("cluster").head(5)
    rep["key"] = rep.file.str.replace(".h5", "", regex=False)
    disp = pd.read_csv(f"{RES}/H8_display_slides.csv"); dk = disp[disp.role == "display"].key.tolist()
    jobs = [(k, g, gdc) for k, g in rep.groupby("key")]
    with ThreadPoolExecutor(4) as ex:
        for k in ex.map(tiles_for_slide, jobs): print("tiles", k[:12], flush=True)
        for k in ex.map(thumb, [(k, gdc) for k in dk]): print("thumb", k[:12], flush=True)
    print("IMAGES_DONE")


if __name__ == "__main__":
    main()
