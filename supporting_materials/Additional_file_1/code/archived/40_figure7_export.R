# Fig 7 (new): data export for panels (section-level group means for the lesion axis; one-section maps). Spec A5.
suppressPackageStartupMessages({library(Seurat); library(data.table); library(Matrix)})
root <- "${PROJECT_ROOT}"; out <- file.path(root, "04_results/2026-10-12_Figure7_GALNT_run01")
G <- c("GALNT1", "GALNT2", "GALNT3", "GALNT4", "GALNT7", "GALNT10", "GALNT11", "GALNT12")
x <- readRDS("${PANIN_ROOT}/01_raw/Elhossiny_CancerDiscovery_2026/spatial.rds"); if (length(Layers(x[["Spatial"]], search = "counts")) > 1) x[["Spatial"]] <- JoinLayers(x[["Spatial"]])
cn <- GetAssayData(x, assay = "Spatial", layer = "counts"); md <- as.data.table(x@meta.data, keep.rownames = "spot"); W <- as.data.table(t(as.matrix(GetAssayData(x, assay = "RCTD_sub_wts", layer = "data"))), keep.rownames = "spot"); md <- merge(md, W, by = "spot"); md <- md[match(colnames(cn), spot)]
m <- cn; m@x <- m@x / rep(colSums(m), diff(m@p)) * 1e4; m@x <- log1p(m@x); E <- as.data.table(t(as.matrix(m[G, ]))); E[, spot := colnames(cn)]; md <- merge(md, E, by = "spot")
d1 <- md[(SpatialDomain == "Ducts" & tissue %in% c("Healthy", "AdjNormal")) | (tissue == "PDAC" & SpatialDomain %in% c("PanIN/GLTumor", "PDTumor"))]; d1[, grp := fifelse(SpatialDomain == "Ducts", "Normal duct", fifelse(SpatialDomain == "PanIN/GLTumor", "PanIN/GL tumor", "PD tumor"))]
sm <- d1[, c(lapply(.SD, mean), .(n = .N)), by = .(sample_Id, tissue, grp), .SDcols = G]; fwrite(sm, file.path(out, "H7f_section_group_means.csv")); print(sm[, .N, by = grp]); print(sm[, .(n_spots_min = min(n), median = median(n)), by = grp])
myc <- c("TAGLN", "MYL9", "TPM1", "TPM2", "MMP11", "POSTN", "HOPX", "ACTA2"); s <- md[sample_Id == "Patient_24_SP"]; mm <- m[, s$spot]; s[, myCAF := colMeans(t(scale(t(as.matrix(mm[myc, ])))))]
fwrite(s[, .(spot, SpatialDomain, x = pxl_col_in_fullres, y = pxl_row_in_fullres, epi_w = `Tumor-Epithelial`, fib_w = Fibroblasts, myCAF, GALNT3, GALNT4, GALNT7, GALNT10, GALNT12)], file.path(out, "H7f_section_map_Patient_24_SP.csv.gz"))
