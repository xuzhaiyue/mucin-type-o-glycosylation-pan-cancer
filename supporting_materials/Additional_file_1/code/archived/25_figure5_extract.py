"""Figure 5 step 1: extract GALNT-system, subtype and marker genes + metadata from the three annotated PanIN-project cohorts.
Spec: 00_meta/figure5_spec_v1.md. Output: 02_intermediate/fig5_sc/<cohort>_{obs.csv.gz, counts.npz, genes.txt}.
Counts are raw (layers['counts']); library size = sum over all genes."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

ROOT = Path("${PROJECT_ROOT}")
OUT = ROOT / "02_intermediate/fig5_sc"; OUT.mkdir(parents=True, exist_ok=True)
BASE = Path("${PANIN_ROOT}/02_intermediate/2026-09-16_six_cohort_input_bundle/standardized")
GALNT = [f"GALNT{i}" for i in range(1, 19)] + ["GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"]
ALIAS = {"GALNT15": ["GALNTL2"], "GALNT16": ["GALNTL1"], "GALNT17": ["WBSCR17"], "GALNT18": ["GALNTL4"]}
MOF_BASAL = "VGLL1 UCA1 S100A2 LY6D SPRR3 SPRR1B LEMD1 KRT15 CTSV DHRS9 AREG CST6 SERPINB3 KRT6C KRT6A SERPINB4 FAM83A SCEL FGFBP1 KRT7 KRT17 GPR87 TNS4 SLC2A1 ANXA8L1".split()
MOF_CLASS = "BTNL8 FAM3D PRR15L AGR3 CTSE LOC400573 LYZ TFF2 TFF1 ANXA10 LGALS4 PLA2G10 CEACAM6 VSIG2 TSPAN8 ST6GALNAC1 AGR2 TFF3 CYP3A7 MYO1A CLRN3 KRT20 CDH17 SPINK4 REG4".split()  # verified vs WikiPathways WP5390 (Moffitt 2015), 25 each
MARKERS = "EPCAM KRT19 KRT8 KRT18 PTPRC CD3E CD8A CD4 NKG7 CD68 LYZ CD14 CD19 MS4A1 JCHAIN COL1A1 COL1A2 DCN LUM PECAM1 VWF RGS5 ACTA2 PRSS1 CPA1 INS GCG CD274 FOXP3 GZMB MKI67".split()
json.dump({"GALNT": GALNT, "alias": ALIAS, "moffitt_basal": MOF_BASAL, "moffitt_classical": MOF_CLASS, "markers": MARKERS}, open(OUT / "fig5_genes.json", "w"), indent=1)
WANT = list(dict.fromkeys(['CTSL2', 'ANXA8L2', 'ANXA8'] + GALNT + sum(ALIAS.values(), []) + MOF_BASAL + MOF_CLASS + MARKERS))
KEEP = {
    "CRA001160": ["sample_id", "condition", "patient_id", "author_celltype", "cell_type", "lineage"],
    "GSE205013": ["sample_id", "patient_id", "Tissue", "Disease_State", "Treatment", "Group", "lineage", "ann_L1", "ann_L2", "ann_L3", "keep_ref_clean", "doublet_final"],
    "GSE202051": ["sample_id", "patient_id", "condition", "condition_detail", "response", "lineage", "Level 1 Annotation", "Moffitt_basal", "Moffitt_classical", "Score_Basal", "Score_Classical", "cnv_score", "new_celltypes", "Level 2 Annotation", "Level 3 Annotation"],
}

def run(cohort: str, chunk: int = 20000) -> None:
    f = h5py.File(BASE / f"{cohort}.h5ad", "r")
    var = f["var"]; idx = var.attrs["_index"] if "_index" in var.attrs else "_index"
    names = np.array([x.decode() if isinstance(x, bytes) else x for x in var[idx][:]])
    pos = {g: i for i, g in enumerate(names)}
    genes = [g for g in WANT if g in pos]; cols = np.array([pos[g] for g in genes])
    g = f["layers/counts"]; n, ncol = [int(x) for x in g.attrs["shape"]]; enc = g.attrs["encoding-type"]
    if enc == "csc_matrix":
        indptr = g["indptr"][:]; lib = np.zeros(n); nnz = int(indptr[-1]); step = 50_000_000
        for s in range(0, nnz, step):
            e = min(s + step, nnz); lib += np.bincount(g["indices"][s:e], weights=g["data"][s:e].astype(float), minlength=n); print(cohort, "lib", e, "/", nnz, flush=True)
        cs = []
        for c in cols:
            a, b = indptr[c], indptr[c + 1]; cs.append(sp.csc_matrix((g["data"][a:b], g["indices"][a:b], np.array([0, b - a])), shape=(n, 1)))
        M = sp.hstack(cs).tocsr()
    else:
        indptr = g["indptr"][:]; blocks = []; lib = np.zeros(n)
        for s in range(0, n, chunk):
            e = min(s + chunk, n); a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((g["data"][a:b], g["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, ncol))
            lib[s:e] = np.asarray(m.sum(1)).ravel(); blocks.append(m[:, cols].tocsr()); print(cohort, e, "/", n, flush=True)
        M = sp.vstack(blocks).tocsr()
    sp.save_npz(OUT / f"{cohort}_counts.npz", M); (OUT / f"{cohort}_genes.txt").write_text("\n".join(genes))
    import anndata as ad
    ob = ad.read_h5ad(BASE / f"{cohort}.h5ad", backed="r").obs
    ob = ob[[c for c in KEEP[cohort] if c in ob.columns]].copy(); ob.insert(0, "cell", ob.index.astype(str)); ob["lib"] = lib
    ob.to_csv(OUT / f"{cohort}_obs.csv.gz", index=False); print(cohort, "done; genes kept", len(genes), "missing", [w for w in WANT if w not in pos])

if __name__ == "__main__":
    for c in (sys.argv[1:] or ["CRA001160", "GSE202051", "GSE205013"]): run(c)
