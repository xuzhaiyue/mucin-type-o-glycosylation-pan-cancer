# Figure 7 analysis (spec: 00_meta/figure7_spec_v1.md): GIScore x SPP1+ macrophage signature in PDAC Visium tumor spots (Elhossiny 2026).
suppressPackageStartupMessages({library(Seurat); library(data.table); library(Matrix); library(presto); library(FNN); library(metafor); library(clusterProfiler); library(msigdbr)})
root <- "${PROJECT_ROOT}"; out <- file.path(root, "04_results/2026-10-12_Figure7_GALNT_run01"); set.seed(666)
MODE <- Sys.getenv("F7MODE", "raw"); sfx <- if (MODE == "smooth") "_smooth" else ""  # raw = frozen spec; smooth = post hoc amendment A2 (see spec)
sig <- jsonlite::fromJSON(file.path(out, "macro_spp1_signature.json"))$genes; b16 <- -0.344; b14 <- 0.202; zs <- function(x) (x - mean(x)) / sd(x)
x <- readRDS("${PANIN_ROOT}/01_raw/Elhossiny_CancerDiscovery_2026/spatial.rds")
ly <- Layers(x[["Spatial"]], search = "counts"); cat("count layers:", length(ly), "\n"); if (length(ly) > 1) x[["Spatial"]] <- JoinLayers(x[["Spatial"]])
cn <- GetAssayData(x, assay = "Spatial", layer = "counts"); md <- as.data.table(x@meta.data, keep.rownames = "spot")
w <- t(as.matrix(GetAssayData(x, assay = "RCTD_sub_wts", layer = "data"))); w <- data.table(spot = rownames(w), mac_w = w[, "Macrophages"], epi_w = w[, "Tumor-Epithelial"])
md <- merge(md, w, by = "spot"); md <- md[match(colnames(cn), spot)]
pd <- md[tissue == "PDAC"]; tum <- pd[SpatialDomain %in% c("PanIN/GLTumor", "PDTumor")]
cat("PDAC spots per section; tumor spots per section\n"); print(merge(pd[, .N, by = sample_Id], tum[, .N, by = sample_Id], by = "sample_Id", suffixes = c("_all", "_tumor"), all = TRUE))
secs <- tum[, .N, by = sample_Id][N >= 100]$sample_Id; tum <- tum[sample_Id %in% secs]; cat("sections kept:", secs, "\n")
lognorm <- function(spots) { m <- cn[, spots, drop = FALSE]; m@x <- m@x / rep(colSums(m), diff(m@p)) * 1e4; m@x <- log1p(m@x); m }
Lt <- lognorm(tum$spot); rownames(Lt) <- rownames(cn)
tum[, GIScore := b16 * zs(as.numeric(Lt["GALNT16", ])) + b14 * zs(as.numeric(Lt["GALNT14", ]))]
sg <- intersect(sig, rownames(Lt)); cat("signature genes present:", length(sg), "/", length(sig), "\n"); Z <- t(scale(t(as.matrix(Lt[sg, ])))); tum[, Macro_SPP1 := colMeans(Z)]
tum[, c16 := as.numeric(cn["GALNT16", tum$spot])]; tum[, lc := log10(nCount_Spatial)]

if (MODE == "smooth") {  # A2: average over the spot and its 6 nearest tumor-region neighbours (same section) to counter zero-inflation
  tum[, `:=`(L16 = as.numeric(Lt["GALNT16", ]), L14 = as.numeric(Lt["GALNT14", ]))]; tum[, `:=`(Macro_raw = Macro_SPP1)]
  sm <- function(v, xy) { nn <- get.knn(xy, k = 6)$nn.index; (v + rowSums(matrix(v[nn], ncol = 6))) / 7 }
  tum[, `:=`(s16 = sm(L16, cbind(pxl_col_in_fullres, pxl_row_in_fullres)), s14 = sm(L14, cbind(pxl_col_in_fullres, pxl_row_in_fullres)), Macro_SPP1 = sm(Macro_raw, cbind(pxl_col_in_fullres, pxl_row_in_fullres))), by = sample_Id]
  tum[, GIScore := b16 * zs(s16) + b14 * zs(s14)]; cat("smoothed: distinct GIScore values", uniqueN(round(tum$GIScore, 6)), " share of modal value", max(table(round(tum$GIScore, 6))) / nrow(tum), "\n") }
## H7a
pc <- function(d, adj) { if (!adj) { r <- cor(d$GIScore, d$Macro_SPP1, method = "spearman") } else { a <- resid(lm(rank(GIScore) ~ mac_w + epi_w + lc, d)); b <- resid(lm(rank(Macro_SPP1) ~ mac_w + epi_w + lc, d)); r <- cor(a, b) }; n <- nrow(d); data.table(rho = r, n = n, se = 1 / sqrt(n - 3 - ifelse(adj, 3, 0))) }
H7a <- rbindlist(lapply(c(FALSE, TRUE), function(adj) tum[, pc(.SD, adj), by = sample_Id][, adjusted := adj]))
H7a[, z := atanh(rho)]; fwrite(H7a, file.path(out, paste0("H7a_per_section", sfx, ".csv"))); print(H7a)
meta <- rbindlist(lapply(c(FALSE, TRUE), function(adj) { d <- H7a[adjusted == adj]; m <- rma(yi = d$z, sei = d$se, method = "DL"); data.table(adjusted = adj, k = nrow(d), pooled_rho = tanh(m$b[1]), lo = tanh(m$ci.lb), hi = tanh(m$ci.ub), p = m$pval, I2 = m$I2, sections_rho_gt0 = sum(d$rho > 0)) })); fwrite(meta, file.path(out, paste0("H7a_meta", sfx, ".csv"))); print(meta)
## co-high labels (per section quartiles)
tum[, `:=`(gh = GIScore >= quantile(GIScore, 0.75), mh = Macro_SPP1 >= quantile(Macro_SPP1, 0.75)), by = sample_Id]; tum[, grp := fifelse(gh & mh, "Both", fifelse(gh, "GIScore_high", fifelse(mh, "Macro_SPP1_high", "Other")))]
print(tum[, .N, by = .(sample_Id, grp)][order(sample_Id, grp)][, dcast(.SD, sample_Id ~ grp, value.var = "N")])
## H7b: DE co-high vs rest (presto Wilcoxon), Stouffer across sections, sign consistency, then Reactome ORA
excl <- unique(c(sg, "GALNT16", "GALNT14")); de <- list()
for (s in secs) { d <- tum[sample_Id == s]; nb <- sum(d$grp == "Both"); if (nb < 20) { cat("skip DE", s, "n_both =", nb, "\n"); next }
  m <- Lt[, d$spot]; keep <- Matrix::rowSums(m > 0) >= 0.05 * ncol(m); r <- as.data.table(wilcoxauc(m[keep, ], ifelse(d$grp == "Both", "Both", "Rest")))[group == "Both"]; r[, `:=`(sample_Id = s, z = sign(auc - 0.5) * qnorm(pmin(pval, 1 - 1e-12) / 2, lower.tail = FALSE))]; de[[s]] <- r }
de <- rbindlist(de); k <- uniqueN(de$sample_Id); cat("DE sections:", k, "\n")
agg <- de[!feature %in% excl, .(n_sec = .N, Z = sum(z) / sqrt(.N), n_up = sum(logFC > 0), mean_logFC = mean(logFC)), by = feature]; agg[, p := 2 * pnorm(-abs(Z))]; agg[, BH := p.adjust(p, "BH")]
need <- ceiling(0.75 * k); agg[, consistent_up := n_up >= need & Z > 0 & BH < 0.05]; fwrite(agg[order(p)], file.path(out, paste0("H7b_DE_cohigh_vs_rest", sfx, ".csv"))); cat("up genes:", sum(agg$consistent_up), " (consistency >=", need, "of", k, "sections)\n")
rg <- unique(as.data.frame(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))[, c("gs_name", "gene_symbol")])
up <- agg[consistent_up == TRUE]$feature; o <- as.data.frame(enricher(up, TERM2GENE = rg, universe = agg$feature, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1))
if (nrow(o)) { o$Description <- sub("^REACTOME_", "", o$ID); fwrite(o, file.path(out, paste0("H7b_reactome_ora", sfx, ".csv"))); print(head(o[, c("Description", "Count", "pvalue", "p.adjust")], 15)) } else cat("no ORA terms\n")
## H7c: neighbour ligand x receptor coupling, permutation within section
LR <- list(`SPP1-CD44` = c("SPP1", "CD44"), `MDK-NCL` = c("MDK", "NCL")); res <- list()
for (s in secs) { d_all <- pd[sample_Id == s]; M <- lognorm(d_all$spot); rownames(M) <- rownames(cn); xy <- as.matrix(d_all[, .(pxl_col_in_fullres, pxl_row_in_fullres)]); nn <- get.knn(xy, k = 6)$nn.index
  isT <- d_all$spot %in% tum$spot; lab <- tum[match(d_all$spot, spot), grp]; for (nm in names(LR)) { L <- as.numeric(M[LR[[nm]][1], ]); R <- as.numeric(M[LR[[nm]][2], ]); cp <- rowMeans(matrix(L[nn], ncol = 6)) * R
    ti <- which(isT); both <- which(lab == "Both"); if (length(both) < 20) next; obs <- mean(cp[both]) - mean(cp[setdiff(ti, both)]); nb <- length(both)
    perm <- replicate(1000, { b <- sample(ti, nb); mean(cp[b]) - mean(cp[setdiff(ti, b)]) }); res[[paste(s, nm)]] <- data.table(sample_Id = s, pair = nm, obs = obs, perm_mean = mean(perm), perm_sd = sd(perm), z = (obs - mean(perm)) / sd(perm), p_perm = (1 + sum(perm >= obs)) / 1001) } }
H7c <- rbindlist(res); fwrite(H7c, file.path(out, paste0("H7c_coupling_per_section", sfx, ".csv"))); print(H7c)
H7c_m <- H7c[, .(k = .N, Z = sum(z) / sqrt(.N), p = pnorm(-sum(z) / sqrt(.N)), sections_z_gt0 = sum(z > 0)), by = pair]; fwrite(H7c_m, file.path(out, paste0("H7c_coupling_meta", sfx, ".csv"))); print(H7c_m)
## spot-level table for figure panels (all PDAC spots with domain + scores for tumor)
sp_out <- merge(pd[, .(spot, sample_Id, SpatialDomain, mac_w, epi_w, x = pxl_col_in_fullres, y = pxl_row_in_fullres)], tum[, .(spot, GIScore, Macro_SPP1, grp, c16)], by = "spot", all.x = TRUE); fwrite(sp_out, file.path(out, paste0("H7_spots_for_figure", sfx, ".csv.gz")))
cat("tumor spots by section (for figure selection):\n"); print(tum[, .N, by = sample_Id][order(-N)])
