# Fig 7 A4: niche of GIScore / GALNT16 / GALNT14 epithelial spots (epithelial-dominant spots; neighbours exclude self).
suppressPackageStartupMessages({library(Seurat); library(data.table); library(Matrix); library(FNN); library(metafor)})
root <- "${PROJECT_ROOT}"; out <- file.path(root, "04_results/2026-10-12_Figure7_GALNT_run01"); b16 <- -0.344; b14 <- 0.202; zs <- function(x) (x - mean(x)) / sd(x)
x <- readRDS("${PANIN_ROOT}/01_raw/Elhossiny_CancerDiscovery_2026/spatial.rds"); if (length(Layers(x[["Spatial"]], search = "counts")) > 1) x[["Spatial"]] <- JoinLayers(x[["Spatial"]])
cn <- GetAssayData(x, assay = "Spatial", layer = "counts"); md <- as.data.table(x@meta.data, keep.rownames = "spot")
W <- as.data.table(t(as.matrix(GetAssayData(x, assay = "RCTD_sub_wts", layer = "data"))), keep.rownames = "spot"); md <- merge(md, W, by = "spot"); md <- md[match(colnames(cn), spot)]
pd <- md[tissue == "PDAC"]; m <- cn[, pd$spot]; m@x <- m@x / rep(colSums(m), diff(m@p)) * 1e4; m@x <- log1p(m@x); rownames(m) <- rownames(cn)
myc <- c("TAGLN", "MYL9", "TPM1", "TPM2", "MMP11", "POSTN", "HOPX", "ACTA2"); ic <- c("IL6", "CXCL8", "CXCL1", "CXCL2", "CCL2", "CXCL12", "CFD", "LMNA", "DPT", "HAS1", "HAS2", "AGTR1")
sc <- function(g) colMeans(t(scale(t(as.matrix(m[intersect(g, rownames(m)), ]))))); pd[, `:=`(myCAF = sc(myc), iCAF = sc(ic), L16 = as.numeric(m["GALNT16", ]), L14 = as.numeric(m["GALNT14", ]), lc = log10(nCount_Spatial))]
pd[, epi_w := `Tumor-Epithelial`]; pd[, GIs := b16 * zs(L16) + b14 * zs(L14)]; ct <- setdiff(colnames(W), "spot"); pd[, isT := SpatialDomain %in% c("PanIN/GLTumor", "PDTumor")]; pd[, isF := grepl("^Fibro", SpatialDomain)]
feats <- c(ct, "myCAF", "iCAF")
nich <- function(d) { xy <- cbind(d$pxl_col_in_fullres, d$pxl_row_in_fullres); k <- get.knn(xy, k = 7); nn <- k$nn.index[, 1:6]; sp1 <- median(k$nn.dist[, 1])
  F <- sapply(feats, function(f) rowMeans(matrix(d[[f]][nn], ncol = 6))); F <- as.data.table(F); F[, purity := rowMeans(matrix(d$isT[nn] * 1, ncol = 6))]
  fi <- which(d$isF); F[, dist_fibro := if (length(fi)) { dd <- get.knnx(xy[fi, , drop = FALSE], xy, k = 1)$nn.dist[, 1]; dd / sp1 } else NA_real_]; cbind(d[, .(spot, sample_Id, GIs, L16, L14, lc, epi_w)], F) }
N <- rbindlist(lapply(unique(pd$sample_Id), function(s) nich(pd[sample_Id == s]))); N <- N[epi_w >= 0.5]; cat("epithelial-dominant spots per section:\n"); print(N[, .N, by = sample_Id])
N[, `:=`(d16 = as.numeric(L16 > 0), d14 = as.numeric(L14 > 0))]; fe <- c(feats, "purity", "dist_fibro"); mets <- c(GIScore = "GIs", `GALNT16 detected` = "d16", `GALNT14 detected` = "d14")
res <- rbindlist(lapply(names(mets), function(mn) rbindlist(lapply(fe, function(f) { per <- N[, { a <- resid(lm(rank(get(mets[[mn]])) ~ lc + epi_w)); b <- resid(lm(rank(get(f)) ~ lc + epi_w)); .(rho = suppressWarnings(cor(a, b)), n = .N) }, by = sample_Id][is.finite(rho)]
  per[, `:=`(z = atanh(rho), se = 1 / sqrt(n - 5))]; r <- rma(yi = per$z, sei = per$se, method = "DL"); data.table(metric = mn, feature = f, k = nrow(per), pooled_rho = tanh(r$b[1]), lo = tanh(r$ci.lb), hi = tanh(r$ci.ub), p = r$pval, I2 = r$I2, sections_pos = sum(per$rho > 0)) }))))
res[, BH := p.adjust(p, "BH"), by = metric]; fwrite(res, file.path(out, "H7e_epithelial_niche.csv")); fwrite(N, file.path(out, "H7e_epithelial_spots_niche.csv.gz"))
cat("\nGIScore niche (top by |rho|):\n"); print(res[metric == "GIScore"][order(-abs(pooled_rho))][1:14, .(feature, pooled_rho = round(pooled_rho, 3), lo = round(lo, 3), hi = round(hi, 3), BH = signif(BH, 2), sections_pos, I2 = round(I2))])
cat("\nGALNT16 detected niche:\n"); print(res[metric == "GALNT16 detected"][order(-abs(pooled_rho))][1:10, .(feature, pooled_rho = round(pooled_rho, 3), BH = signif(BH, 2), sections_pos)])
cat("\nGALNT14 detected niche:\n"); print(res[metric == "GALNT14 detected"][order(-abs(pooled_rho))][1:8, .(feature, pooled_rho = round(pooled_rho, 3), BH = signif(BH, 2), sections_pos)])
cat("\nfraction of epithelial spots with GALNT16>0:", mean(N$d16), " GALNT14>0:", mean(N$d14), "\n")
