"""Fig 6 data prep: GSE254829 Visium PanIN slides -> per-spot GALNT expression + manual annotation. Spec: 00_meta/figure6_spec_v1.md"""
import gzip, io, json, tarfile
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp

ROOT = Path("${PROJECT_ROOT}")
RAW = Path("${PANIN_ROOT}/01_raw/GSE254829"); OUT = ROOT / "02_intermediate/fig6_visium"; OUT.mkdir(parents=True, exist_ok=True)
GENES = [f"GALNT{i}" for i in range(1, 19)] + ["GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"]
ALIAS = {"GALNT15": "GALNTL2", "GALNT16": "GALNTL1", "GALNT17": "WBSCR17", "GALNT18": "GALNTL4"}
codes = pd.read_csv(RAW / "GSE254829_codatable_may202024.csv.gz"); codes["bc"] = codes.iloc[:, 0].str.extract(r"^([ACGT]+-\d)")[0]  # suffixes after '-1' differ per slide (_1_1, _2_1, none...)
codes = codes.rename(columns={codes.columns[0]: "spot_id"})
rows = []
outer = tarfile.open(RAW / "GSE254829_RAW.tar")
for m in outer.getmembers():
    name = m.name.split("_", 1)[1].replace(".tar.gz", "")
    if name not in set(codes["sample"]): continue
    inner = tarfile.open(fileobj=io.BytesIO(outer.extractfile(m).read()), mode="r:gz")
    names = {x.name: x for x in inner.getmembers()}
    get = lambda suffix: inner.extractfile(next(v for k, v in names.items() if k.endswith(suffix))).read()
    feats = pd.read_csv(io.BytesIO(get("features.tsv.gz")), sep="\t", header=None, compression="gzip"); bcs = pd.read_csv(io.BytesIO(get("barcodes.tsv.gz")), header=None, compression="gzip")[0].values
    M = sio.mmread(io.BytesIO(gzip.decompress(get("matrix.mtx.gz")))).tocsr()
    sym = feats[1].values; lib = np.asarray(M.sum(0)).ravel(); idx = {g: i for i, g in enumerate(sym)}
    cols = {}
    for g in GENES:
        gg = g if g in idx else ALIAS.get(g)
        if gg in idx: cols[g] = np.log1p(np.asarray(M[idx[gg]].todense()).ravel() / np.maximum(lib, 1) * 1e4)
    df = pd.DataFrame(cols); df["bc"] = bcs; df["lib"] = lib; df["sample"] = name
    pos = pd.read_csv(io.BytesIO(get("tissue_positions_list.csv")), header=None, names=["bc", "in_tissue", "arow", "acol", "px_row", "px_col"])
    df = df.merge(pos[["bc", "px_row", "px_col"]], on="bc", how="left")
    rows.append(df)
    for fn in ["tissue_hires_image.png", "scalefactors_json.json"]:
        d = OUT / name; d.mkdir(exist_ok=True); (d / fn).write_bytes(get(fn))
    print(name, M.shape, "genes found:", len(cols), "missing:", [g for g in GENES if g not in cols])
A = pd.concat(rows, ignore_index=True)
A = A.merge(codes.drop(columns=["pixel_row", "pixel_col"]), on=["sample", "bc"], how="inner")
A.to_csv(OUT / "spots_galnt_annot.csv.gz", index=False)
print(A.shape, A["paningrade"].value_counts().to_dict()); print(A.groupby(["sample", "paningrade"]).size().unstack(fill_value=0))
