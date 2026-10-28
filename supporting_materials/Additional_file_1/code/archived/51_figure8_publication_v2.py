"""Fig 8 publication figure v2, laid out after the template digital-pathology figure (workflow + t-SNE, cluster tile panels,
whole-slide cluster maps for a GIScore-High and a GIScore-Low patient with H&E + attention zoom), plus statistics rows. 7.2 in wide."""
from __future__ import annotations
import glob, json, os
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Patch
from matplotlib.colors import ListedColormap
from scipy.ndimage import zoom, gaussian_filter
from sklearn.metrics import roc_curve, roc_auc_score
from PIL import Image, PngImagePlugin
PngImagePlugin.MAX_TEXT_CHUNK = 200 * 1024 * 1024; Image.MAX_IMAGE_PIXELS = None

R = "${PROJECT_ROOT}"
OUT = f"{R}/02_intermediate/fig8"; RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"
FD = f"{R}/05_figures/2026-10-12_Figure8_pathology_run01"; os.makedirs(FD, exist_ok=True)
mpl.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6.5, "axes.linewidth": 0.5, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2, "ytick.major.size": 2})
CC = ["#D62728", "#1F77B4", "#8FA31E", "#9467BD", "#FF9F1C", "#7F7F7F"]
NAMES = {1: "tumour cells in inflamed stroma", 2: "loose stroma, fat, oedema", 3: "fibrous stroma", 4: "pale sparse stroma", 5: "dense, deeply stained epithelium", 6: "pink fibro-fatty tissue"}
W, H = 7.2, 9.95
fig = plt.figure(figsize=(W, H))


def axi(x: float, y: float, w: float, h: float, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def tag(x: float, y: float, t: str) -> None:
    fig.text(x / W, 1 - y / H, t, fontsize=10, fontweight="bold", va="top", ha="left")


def tiles_of(c: int) -> list[str]:
    return sorted(f for f in glob.glob(f"{OUT}/tiles/C{c}_*.png") if not os.path.basename(f).startswith("._"))


# A10-Q4: C5 examples = nearest-to-centroid tiles of distinct slides with ink-flagged tiles excluded (rank from H8_representative_tiles.csv)
def c5_examples() -> list[str]:
    rep = pd.read_csv(f"{RES}/H8_representative_tiles.csv"); rep = rep[rep.cluster == 5]; out = []
    for x in rep.itertuples():
        t = pd.read_csv(f"{R}/04_results/2026-10-12_Figure8_ink_QC_run01/ink_tiles/{x.patient}.csv.gz").set_index("tile_idx").loc[x.tile_idx]
        f = f"{OUT}/tiles/C5_{x.patient}_{int(x.tile_idx)}.png"
        if not (t.core or t.edge) and os.path.exists(f): out.append(f)
    return out


def rimg(f: str) -> Image.Image:
    return Image.open(f).convert("RGB")


def cluster_canvas(k: str, th: Image.Image, ds: float) -> np.ndarray:
    at = pd.read_csv(f"{OUT}/alltiles_{k[:12]}.csv.gz"); g = np.asarray(th.convert("L"), float)[..., None] / 255.0; cv = np.repeat(0.80 + 0.20 * g, 3, axis=2)
    p = 512 / ds
    for c in range(1, 7):
        col = np.array(mpl.colors.to_rgb(CC[c - 1]))
        for x, y in zip(at.x[at.cluster == c], at.y[at.cluster == c]):
            cv[int(y / ds):int(np.ceil(y / ds + p)), int(x / ds):int(np.ceil(x / ds + p))] = col
    im = Image.fromarray((cv * 255).astype(np.uint8)); sc = 900 / max(im.size)
    return np.asarray(im.resize((int(im.size[0] * sc), int(im.size[1] * sc)), Image.BOX))


# ---------------- data
d = pd.read_csv(f"{RES}/H8_slide_table.csv"); a = pd.read_csv(f"{RES}/H8a_cluster_vs_GIScore.csv"); ac = pd.read_csv(f"{RES}/H8d_CPTAC_cluster_vs_GIScore.csv")
ts = pd.read_csv(f"{RES}/H8_tsne.csv"); auc = json.load(open(f"{RES}/H8b_cv_auc.json")); conf = json.load(open(f"{RES}/H8c_confound_auc.json"))
oof = pd.read_csv(f"{RES}/H8b_oof_predictions.csv"); ext = json.load(open(f"{RES}/H8e_CPTAC_transfer.json")); extp = pd.read_csv(f"{RES}/H8e_CPTAC_predictions.csv")
pm = pd.read_csv(f"{RES}/H8b_perm_auc.csv").pooled_auc; pb = pd.read_csv(f"{RES}/H8b_baseline_perm_auc.csv").pooled_auc
zinfo = json.load(open(f"{OUT}/zoom_info.json"))
mil_obs = float(np.mean(auc["MIL_pooled_auc_by_rep"])); base_obs = float(np.mean(auc["base_pooled_auc_by_rep"]))
p_mil = (np.sum(pm >= mil_obs) + 1) / (len(pm) + 1); p_base = (np.sum(pb >= base_obs) + 1) / (len(pb) + 1)

# ---------------- A workflow (template style: tiles -> dashed boxes -> output)
tag(0.10, 0.08, "A")
ax = axi(0.30, 0.22, 4.25, 2.35); ax.axis("off"); ax.set_xlim(0, 4.25); ax.set_ylim(0, 2.35)
ax.text(2.12, 2.30, "In silico pathological region annotation", ha="center", va="top", fontsize=7.2, style="italic")
mos = [tiles_of(1)[0], tiles_of(2)[0], tiles_of(3)[0], tiles_of(4)[0], tiles_of(5)[2], tiles_of(6)[0], tiles_of(1)[1], tiles_of(3)[1], tiles_of(5)[4]]
for i, f in enumerate(mos):
    ax.imshow(rimg(f), extent=[0.02 + (i % 3) * 0.43, 0.02 + (i % 3) * 0.43 + 0.41, 1.78 - (i // 3) * 0.43 - 0.41 + 0.0, 1.78 - (i // 3) * 0.43], zorder=2)
ax.set_xlim(0, 4.25); ax.set_ylim(0, 2.35)
ax.annotate("", xy=(1.52, 1.15), xytext=(1.36, 1.15), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#3A6EA5"))
ax.add_patch(FancyBboxPatch((1.55, 0.35), 1.45, 1.75, boxstyle="round,pad=0.02,rounding_size=0.06", fc="none", ec="#3A6EA5", lw=0.8, ls=(0, (3, 2))))
ax.text(2.27, 2.0, "Unsupervised learning", ha="center", va="center", fontsize=6.4, fontweight="bold")
for j, s in enumerate(["UNI2-h tile features\n(ViT-H/14, 1,536-d)", "PCA (32 components)", "k-means clustering\n(k = 6)"]):
    y = 1.50 - j * 0.46
    ax.add_patch(FancyBboxPatch((1.63, y - 0.15), 1.29, 0.32, boxstyle="round,pad=0.01,rounding_size=0.05", fc="white", ec="#3A6EA5", lw=0.8)); ax.text(2.27, y + 0.01, s, ha="center", va="center", fontsize=5.8)
    if j < 2: ax.annotate("", xy=(2.27, y - 0.30), xytext=(2.27, y - 0.17), arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#3A6EA5"))
ax.annotate("", xy=(3.2, 1.15), xytext=(3.03, 1.15), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#3A6EA5"))
ax.add_patch(FancyBboxPatch((3.22, 0.35), 0.98, 1.75, boxstyle="round,pad=0.02,rounding_size=0.06", fc="none", ec="#3A6EA5", lw=0.8, ls=(0, (3, 2))))
ax.text(3.71, 2.0, "Region labelling", ha="center", va="center", fontsize=6.4, fontweight="bold")
for j, lab_ in enumerate(["tumour, inflamed", "loose, fat, oedema", "fibrous stroma", "pale sparse stroma", "dense epithelium", "pink fibro-fatty"]):
    yy = 1.72 - j * 0.15
    ax.add_patch(Rectangle((3.29, yy - 0.04), 0.08, 0.08, fc=CC[j], ec="none")); ax.text(3.41, yy, f"C{j+1} {lab_}", fontsize=4.7, va="center")
ax.plot([3.29, 4.15], [0.80, 0.80], color="#3A6EA5", lw=0.4)
ax.text(3.71, 0.62, "per-slide fractions\nvs frozen GIScore\n(High 73 / Low 72)", ha="center", va="center", fontsize=5.0)
ax.text(3.71, 0.42, "145 slides, 269,259 tiles", ha="center", va="center", fontsize=4.4, color="#555")
ax.text(0.02, 0.28, "tiles: 256 px at 20x, TCGA-PAAD DX1 slides", fontsize=5.0, color="#555", va="center")

# ---------------- B t-SNE
tag(4.70, 0.08, "B")
ax = axi(4.95, 0.33, 2.15, 2.0); t2 = ts.sample(frac=1, random_state=1)
ax.scatter(t2.tsne1, t2.tsne2, c=[CC[c - 1] for c in t2.cluster], s=0.45, lw=0, rasterized=True)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1", labelpad=1); ax.set_ylabel("t-SNE 2", labelpad=1); ax.set_title("t-SNE of UNI2-h tile features", fontsize=6.8, pad=3)
for s in ("left", "bottom"): ax.spines[s].set_linewidth(0.5)
ax.legend(handles=[Patch(color=CC[i], label=f"C{i+1}") for i in range(6)], ncol=6, frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.045), handlelength=0.9, columnspacing=0.8, handletextpad=0.3)

# ---------------- C cluster tiles
tag(0.10, 2.78, "C")
gw = 2.3; tw = 0.70
for c in range(1, 7):
    gx = 0.30 + ((c - 1) % 3) * (gw + 0.06); gy = 2.80 + ((c - 1) // 3) * 0.93
    fig.text((gx + gw / 2) / W, 1 - (gy + 0.01) / H, f"Cluster {c}: {NAMES[c]}", fontsize=6.0, ha="center", va="top", color=CC[c - 1], fontweight="bold")
    for j, f in enumerate((c5_examples() if c == 5 else tiles_of(c))[:3]):
        ax = axi(gx + j * (tw + 0.03), gy + 0.17, tw, tw); ax.imshow(rimg(f)); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_color(CC[c - 1]); s.set_linewidth(1.0)

fig.text(0.30 / W, 1 - 4.63 / H, "C5 examples: nearest-to-centroid tiles of distinct slides; ink-marked tiles excluded.", fontsize=4.6, color="#444", va="top")
# ---------------- D/E slides
disp = {k: v for k, v in zinfo.items()}
for row, (grp, lab) in enumerate([("High", "D"), ("Low", "E")]):
    k12 = [k for k, v in disp.items() if v["group"] == grp][0]; v = disp[k12]
    y0 = 4.84 + row * 1.5
    tag(0.10, y0 + 0.05, lab)
    fig.text(0.34 / W, 1 - (y0 + 0.82) / H, f"GIScore-{grp} patient\n{k12} ({v['risk']:+.2f})", fontsize=6.0, rotation=90, va="center", ha="center", fontweight="bold")
    thf = f"{OUT}/thumb_{k12}.png"; th = rimg(thf); meta = json.load(open(f"{OUT}/thumb_{k12}.json")); ds = meta["ds"]
    key = [kk for kk in glob.glob(f"{OUT}/alltiles_{k12}.csv.gz")][0]
    cv = cluster_canvas(k12 + "_", th, ds)
    ax = axi(0.55, y0 + 0.05, 1.35, 1.35); ax.imshow(cv); ax.axis("off"); ax.set_anchor("C")
    at = pd.read_csv(key); fr = at.cluster.value_counts(normalize=True)
    ax.set_title("cluster map", fontsize=6.2, pad=1.5)
    ax = axi(2.0, y0 + 0.05, 1.35, 1.35); ax.imshow(th); ax.axis("off"); ax.set_anchor("C"); ax.set_title("H&E (diagnostic slide)", fontsize=6.2, pad=1.5)
    x0, y0w, sz = v["x0"], v["y0"], v["size"]
    ax.add_patch(Rectangle((x0 / ds, y0w / ds), sz / ds, sz / ds, fill=False, ec="k", lw=0.8))
    zi = rimg(f"{OUT}/zoom_{k12}.png"); A0 = np.load(f"{OUT}/zoomattn_{k12}.npy"); msk = zoom((~np.isnan(A0)).astype(float), zi.size[0] / A0.shape[0], order=0) > 0.5
    A = np.where(np.isnan(A0), np.nanmean(A0), A0)
    Az = np.ma.masked_where(~msk, gaussian_filter(zoom(A, zi.size[0] / A0.shape[0], order=1), zi.size[0] / A0.shape[0] * 0.5))
    for j, (ttl, kind) in enumerate([("original", 0), ("attention", 1), ("overlay", 2)]):
        ax = axi(3.62 + j * 1.07, y0 + 0.225, 1.0, 1.0); ax.axis("off"); ax.set_title(ttl, fontsize=6.2, pad=1.5)
        if kind == 0: ax.imshow(zi)
        elif kind == 1: im = ax.imshow(Az, cmap="jet", vmin=0, vmax=1)
        else: ax.imshow(zi); ax.imshow(Az, cmap="jet", vmin=0, vmax=1, alpha=0.42)
    cax = axi(6.88, y0 + 0.225, 0.07, 1.0); cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5); cb.outline.set_linewidth(0.4)
    fig.text(3.62 / W, 1 - (y0 + 1.27) / H, f"attention: leave-one-slide-out ABMIL, P(High) = {v['prob']:.2f}; window centred on peak attention", fontsize=5.0, color="#444", va="top")
    fig.text(0.55 / W, 1 - (y0 + 1.42) / H, "C1 " + f"{fr.get(1, 0)*100:.0f}%" + "   C2 " + f"{fr.get(2, 0)*100:.0f}%" + "   C3 " + f"{fr.get(3, 0)*100:.0f}%", fontsize=5.4, va="top", color="#444")
# (cluster colour key: see panels A and B)

fig.text(0.55 / W, 1 - 7.88 / H, "Released UNI2-h embeddings cover bands of ~46-53% of the tissue (faded H&E = not covered); zoom windows show covered tiles only.", fontsize=5.0, color="#444", va="top")
# ---------------- statistics
cs = [f"C{i}" for i in range(1, 7)]
tag(0.10, 8.12, "F")
ax = axi(0.55, 8.36, 3.2, 1.2)
for i, c in enumerate(cs):
    for k, (g, off, col) in enumerate([("Low", -0.19, "#4C78A8"), ("High", 0.19, "#E4572E")]):
        vv = d.loc[d.group == g, c].values
        ax.boxplot(vv, positions=[i + off], widths=0.3, showfliers=False, patch_artist=True, medianprops=dict(color="k", lw=0.6),
                   boxprops=dict(lw=0.4, facecolor=col, alpha=0.8), whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
        ax.scatter(np.random.default_rng(i * 2 + k).normal(i + off, 0.045, len(vv)), vv, s=0.9, c="k", alpha=0.25, lw=0)
ax.set_xticks(range(6)); ax.set_xticklabels(cs); ax.set_ylabel("fraction of tiles per slide"); ax.set_ylim(0, 1)
r1 = a.iloc[0]; ax.text(0, 0.80, f"β = {r1.beta_adj:.2f}\nBH = {r1.bh_adj:.3f}", ha="center", fontsize=5.6, color="#C0392B")
ax.legend(handles=[Patch(color="#4C78A8", label="GIScore Low (n = 72)"), Patch(color="#E4572E", label="GIScore High (n = 73)")], frameon=False, fontsize=5.5, loc="upper right")
ax.set_title("TCGA-PAAD cluster composition", fontsize=6.6, pad=3)

tag(3.95, 8.12, "G")
ax = axi(4.45, 8.36, 2.65, 1.2); y = np.arange(6)[::-1]; ax.axvline(0, color="#999", lw=0.4)
ax.errorbar(a.beta_adj, y + 0.14, xerr=[a.beta_adj - a.lo, a.hi - a.beta_adj], fmt="o", ms=2.6, lw=0.8, color="k", capsize=0, label="TCGA-PAAD (n = 145)")
ax.errorbar(ac.beta, y - 0.14, xerr=[ac.beta - ac.lo, ac.hi - ac.beta], fmt="s", ms=2.4, lw=0.8, color="#2E86AB", capsize=0, label="CPTAC-PDA (n = 140)")
ax.set_yticks(y); ax.set_yticklabels(cs); ax.set_xlabel("β per 1 SD GIScore (CLR fraction)", labelpad=1); ax.set_xlim(-0.8, 0.8)
ax.legend(frameon=False, fontsize=5.4, loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=1, borderaxespad=0.1)
ax.text(0.99, 0.99, f"C1: TCGA BH = {a.bh_adj[0]:.3f}\nCPTAC one-sided p = {ac.p_onesided_neg[0]:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=5.2, color="#444")

fig.savefig(f"{FD}/Figure8.png", dpi=400); fig.savefig(f"{FD}/Figure8.svg"); fig.savefig(f"{FD}/Figure8.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

# ---------------- Supplementary Fig S8 (former Fig 8H, I): prediction and site confounding
FD = f"{R}/05_figures/2026-10-12_FigureS8_pathology_run01"; os.makedirs(FD, exist_ok=True)
W, H = 7.2, 2.3; fig = plt.figure(figsize=(W, H))
tag(0.10, 0.08, "A")
ax = axi(0.55, 0.33, 2.55, 1.3)
rg_t = pd.read_csv(f"{RES}/H8g_ridge_oof_TCGA.csv"); rg_c = pd.read_csv(f"{RES}/H8g_ridge_pred_CPTAC.csv"); rr = json.load(open(f"{RES}/H8g_ridge_audit.json"))
assert (rg_c.patient.values == extp.patient.values).all() and (rg_t.patient.values == oof.patient.values).all()
rngb = np.random.default_rng(20261012)


def boot(y_, s1, s2, n=2000):
    y_ = np.asarray(y_); s1 = np.asarray(s1); s2 = np.asarray(s2); dl, a1, a2 = [], [], []
    for _ in range(n):
        i = rngb.integers(0, len(y_), len(y_))
        if y_[i].min() == y_[i].max(): continue
        u, v = roc_auc_score(y_[i], s1[i]), roc_auc_score(y_[i], s2[i]); a1.append(u); a2.append(v); dl.append(u - v)
    q = lambda x: np.percentile(x, [2.5, 97.5])
    return q(a1), q(a2), q(dl)


ci_rt, ci_at, d_t = boot(oof.y, rg_t.ridge_oof_mean, oof.mil_oof_mean)
ci_rc, ci_ac, d_c = boot(extp.y, rg_c.ridge, extp.mil)
au = lambda y_, s_: roc_auc_score(y_, s_)
for (yv, sv, col, ls, lw, lab) in [(oof.y, rg_t.ridge_oof_mean, "#E4572E", "-", 1.1, f"TCGA ridge {au(oof.y, rg_t.ridge_oof_mean):.2f} ({ci_rt[0]:.2f}-{ci_rt[1]:.2f})"),
                                    (oof.y, oof.mil_oof_mean, "#E4572E", "--", 0.8, f"TCGA ABMIL {au(oof.y, oof.mil_oof_mean):.2f} ({ci_at[0]:.2f}-{ci_at[1]:.2f})"),
                                    (extp.y, rg_c.ridge, "#2E86AB", "-", 1.1, f"CPTAC ridge {au(extp.y, rg_c.ridge):.2f} ({ci_rc[0]:.2f}-{ci_rc[1]:.2f})"),
                                    (extp.y, extp.mil, "#2E86AB", "--", 0.8, f"CPTAC ABMIL {au(extp.y, extp.mil):.2f} ({ci_ac[0]:.2f}-{ci_ac[1]:.2f})")]:
    fpr, tpr, _ = roc_curve(yv, sv); ax.plot(fpr, tpr, ls, lw=lw, color=col, label=lab)
ax.plot([0, 1], [0, 1], color="#aaa", lw=0.4); ax.set_xlabel("1 - specificity", labelpad=1); ax.set_ylabel("sensitivity", labelpad=1)
ax.legend(frameon=False, fontsize=5.0, loc="lower right", title="AUC (95% CI)", title_fontsize=5.0, borderaxespad=0.2, handlelength=1.6)
fig.text(0.55 / W, 1 - 1.93 / H, f"Paired AUC difference, ridge - ABMIL: TCGA {au(oof.y, rg_t.ridge_oof_mean) - au(oof.y, oof.mil_oof_mean):+.2f} ({d_t[0]:+.2f} to {d_t[1]:+.2f}); CPTAC {au(extp.y, rg_c.ridge) - au(extp.y, extp.mil):+.2f} ({d_c[0]:+.2f} to {d_c[1]:+.2f})", fontsize=4.6, color="#444", va="top")
fig.text(0.55 / W, 1 - 2.04 / H, "ABMIL pre-specified; ridge post hoc (nested alpha). CPTAC: models trained on all TCGA slides.", fontsize=4.6, color="#444", va="top")
ax.set_title("Predicting GIScore High vs Low from H&E", fontsize=6.4, pad=3)
json.dump(dict(tcga=dict(ridge_ci=ci_rt.tolist(), abmil_ci=ci_at.tolist(), diff_ci=d_t.tolist()), cptac=dict(ridge_ci=ci_rc.tolist(), abmil_ci=ci_ac.tolist(), diff_ci=d_c.tolist())), open(f"{RES}/H8h_paired_auc_bootstrap.json", "w"), indent=1)

tag(3.95, 0.08, "B")
ax = axi(4.45, 0.33, 2.65, 1.3)
items = [("GIScore\n(ABMIL)", auc["MIL_fold_auc_mean"], "#E4572E"), ("GIScore\n(site-strat.)", conf["GI_TSSstratified_fold_auc_mean"], "#E4572E"),
         ("Grade", conf["grade_fold_auc_mean"], "#888888"), ("Tissue\nsource site", conf["TSS_macroAUC_fold_mean"], "#2E86AB")]
for i, (l, v_, c) in enumerate(items):
    ax.bar(i, v_, color=c, width=0.62); ax.text(i, v_ + 0.015, f"{v_:.2f}", ha="center", fontsize=5.6)
for j, (sn, r) in enumerate(conf["leave_site_out"].items()):
    ax.scatter(4.4 + j * 0.45, r["auc"], s=9, c="#E4572E", zorder=3); ax.text(4.4 + j * 0.45, r["auc"] + 0.035, f"{sn}\n{r['n']}", ha="center", fontsize=4.6)
ax.axhline(0.5, color="#aaa", lw=0.4, ls="--"); ax.set_xticks(list(range(4)) + [4.85]); ax.set_xticklabels([i[0] for i in items] + ["Leave-site-\nout"], fontsize=5.0)
ax.set_ylim(0.4, 1.05); ax.set_ylabel("AUC (TCGA CV)", labelpad=1); ax.set_title("Site confounding of H&E features", fontsize=6.4, pad=3)

fig.savefig(f"{FD}/FigureS8.png", dpi=400); fig.savefig(f"{FD}/FigureS8.svg"); fig.savefig(f"{FD}/FigureS8.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print("saved", len(pm), "MIL perms; p =", p_mil)
