# Fig 7 A3: GIScore / GALNT16 / GALNT14 co-localisation scan vs RCTD cell-type weights and myCAF/iCAF scores (spec A3).
suppressPackageStartupMessages({library(Seurat); library(data.table); library(Matrix); library(FNN); library(metafor)})
root <- "${PROJECT_ROOT}"; out <- file.path(root, "04_results/2026-10-12_Figure7_GALNT_run01"); b16 <- -0.344; b14 <- 0.202; zs <- function(x) (x - mean(x)) / sd(x)
x <- readRDS("${PANIN_ROOT}/01_raw/Elhossiny_CancerDiscovery_2026/spatial.rds"); if (length(Layers(x[["Spatial"]], search = "counts")) > 1) x[["Spatial"]] <- JoinLayers(x[["Spatial"]])
cn <- GetAssayData(x, assay = "Spatial", layer = "counts"); md <- as.data.table(x@meta.data, keep.rownames = "spot")
W <- as.data.table(t(as.matrix(GetAssayData(x, assay = "RCTD_sub_wts", layer = "data"))), keep.rownames = "spot"); md <- merge(md, W, by = "spot"); md <- md[match(colnames(cn), spot)]
pd <- md[tissue == "PDAC"]; m <- cn[, pd$spot]; m@x <- m@x / rep(colSums(m), diff(m@p)) * 1e4; m@x <- log1p(m@x); rownames(m) <- rownames(cn)
myc <- c("TAGLN", "MYL9", "TPM1", "TPM2", "MMP11", "POSTN", "HOPX", "ACTA2"); ic <- c("IL6", "CXCL8", "CXCL1", "CXCL2", "CCL2", "CXCL12", "CFD", "LMNA", "DPT", "HAS1", "HAS2", "AGTR1")
cat("myCAF present:", sum(myc %in% rownames(m)), "/", length(myc), "; iCAF present:", sum(ic %in% rownames(m)), "/", length(ic), "\n")
sc <- function(g) { g <- intersect(g, rownames(m)); colMeans(t(scale(t(as.matrix(m[g, ]))))) }
pd[, `:=`(myCAF = sc(myc), iCAF = sc(ic), L16 = as.numeric(m["GALNT16", ]), L14 = as.numeric(m["GALNT14", ]), lc = log10(nCount_Spatial))]
tcols <- c(colnames(W)[-1], "myCAF", "iCAF"); pd[, tumor := SpatialDomain %in% c("PanIN/GLTumor", "PDTumor")]
## sensitivity smoothing over 6 nearest neighbours in the same section (all spots)
pd[, `:=`(S16 = { nn <- get.knn(cbind(pxl_col_in_fullres, pxl_row_in_fullres), k = 6)$nn.index; (L16 + rowSums(matrix(L16[nn], ncol = 6))) / 7 }, S14 = { nn <- get.knn(cbind(pxl_col_in_fullres, pxl_row_in_fullres), k = 6)$nn.index; (L14 + rowSums(matrix(L14[nn], ncol = 6))) / 7 }), by = sample_Id]
run <- function(d, metrics, tag) { rbindlist(lapply(names(metrics), function(mn) rbindlist(lapply(tcols, function(tg) { per <- d[, { a <- resid(lm(rank(get(metrics[[mn]])) ~ lc)); b <- resid(lm(rank(get(tg)) ~ lc)); .(rho = cor(a, b), n = .N) }, by = sample_Id]; per[, z := atanh(rho)]; per[, se := 1 / sqrt(n - 4)]
  r <- rma(yi = per$z, sei = per$se, method = "DL"); data.table(scope = tag, metric = mn, target = tg, k = nrow(per), pooled_rho = tanh(r$b[1]), lo = tanh(r$ci.lb), hi = tanh(r$ci.ub), p = r$pval, I2 = r$I2, sections_pos = sum(per$rho > 0)) })))) }
mk <- function(d, sm) { d[, GIs := b16 * zs(L16) + b14 * zs(L14)]; d[, GIs_sm := b16 * zs(S16) + b14 * zs(S14)]; d }
A <- mk(copy(pd), TRUE); B <- mk(copy(pd[tumor == TRUE]), TRUE)
raw <- c(GIScore = "GIs", GALNT16 = "L16", GALNT14 = "L14"); smo <- c(GIScore = "GIs_sm", GALNT16 = "S16", GALNT14 = "S14")
res <- rbind(run(A, raw, "A_all_PDAC_spots"), run(B, raw, "B_tumor_spots"), run(A, smo, "A_all_PDAC_spots_smoothed"), run(B, smo, "B_tumor_spots_smoothed"))
res[, BH := p.adjust(p, "BH"), by = .(scope, metric)]; fwrite(res, file.path(out, "H7d_colocalization_scan.csv"))
for (sc_ in unique(res$scope)) { cat("\n==", sc_, "\n"); print(dcast(res[scope == sc_], target ~ metric, value.var = "pooled_rho")[order(-GIScore)][, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]) }
cat("\nGIScore BH<0.05:\n"); print(res[metric == "GIScore" & BH < 0.05][order(scope, p)][, .(scope, target, pooled_rho = round(pooled_rho, 3), p = signif(p, 2), BH = signif(BH, 2), sections_pos)])
fwrite(A[, .(spot, sample_Id, SpatialDomain, tumor, x = pxl_col_in_fullres, y = pxl_row_in_fullres, GIScore = GIs, GIScore_sm = GIs_sm, L16, L14, myCAF, iCAF, Fibroblasts, `Tumor-Epithelial`, Macrophages)], file.path(out, "H7d_spots_for_figure.csv.gz"))
