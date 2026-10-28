# Figure 6 analysis (spec: 00_meta/figure6_spec_v1.md). H1 Visium PanIN grade; H2/H3 Elhossiny epithelial pseudobulk.
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest); library(clinfun); library(Seurat)})
root <- "${PROJECT_ROOT}"; out <- file.path(root, "04_results/2026-10-12_Figure6_GALNT_run01"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"); genes <- genes[genes %in% c(paste0("GALNT", 1:18), "C1GALT1", "C1GALT1C1", "GALNTL5", "GALNTL6")]
b16 <- -0.344; b14 <- 0.202; zs <- function(x) (x - mean(x)) / sd(x)
## ---------- H1: GSE254829 Visium
sp <- fread(file.path(root, "02_intermediate/fig6_visium/spots_galnt_annot.csv.gz"))
sp <- sp[paningrade %in% c("normal", "low_grade", "high_grade") & lib >= 500 & nontissue < 100]
sp[, grade := factor(paningrade, levels = c("normal", "low_grade", "high_grade"))]; sp[, g := as.integer(grade) - 1L]
sp[, GIScore := b16 * zs(GALNT16) + b14 * zs(GALNT14)]
fwrite(sp[, .(sample, bc, px_row, px_col, grade, GALNT16, GALNT14, GIScore, lib, acini, collagen, islets, `normal epithelium`)], file.path(out, "H1_spots_used.csv"))
cat("spots used:\n"); print(sp[, .N, by = .(sample, grade)][order(sample, grade)])
fitone <- function(d, y, cov = character()) {
  f <- as.formula(paste(y, "~ g", if (length(cov)) paste("+", paste0("`", cov, "`", collapse = "+")), "+ (1|sample)")); m <- tryCatch(lmer(f, data = d), error = function(e) NULL)
  if (is.null(m)) return(data.table(beta_per_grade = NA, se = NA, p = NA)); s <- coef(summary(m)); data.table(beta_per_grade = s["g", "Estimate"], se = s["g", "Std. Error"], p = s["g", "Pr(>|t|)"])
}
fitcat <- function(d, y) { m <- lmer(as.formula(paste(y, "~ grade + (1|sample)")), data = d); s <- coef(summary(m)); data.table(contrast = c("low_grade vs normal", "high_grade vs normal"), beta = s[2:3, "Estimate"], se = s[2:3, "Std. Error"], p = s[2:3, "Pr(>|t|)"]) }
res <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(y) cbind(metric = y, fitone(sp, y), model = "primary: ~grade + (1|slide)")))
res2 <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(y) cbind(metric = y, fitone(sp, y, c("acini", "collagen")), model = "sens: + acini% + collagen%")))
resc <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(y) cbind(metric = y, fitcat(sp, y))))
fwrite(rbind(res, res2), file.path(out, "H1_trend_mixed_model.csv")); fwrite(resc, file.path(out, "H1_grade_contrasts.csv")); print(rbind(res, res2)); print(resc)
## slide-level means (for paired lines)
sl <- sp[, .(GALNT16 = mean(GALNT16), GALNT14 = mean(GALNT14), GIScore = mean(GIScore), n = .N), by = .(sample, grade)]; fwrite(sl, file.path(out, "H1_slide_means.csv"))
## all-gene trend (exploratory, BH)
gt <- rbindlist(lapply(genes, function(gn) { if (!gn %in% names(sp)) return(NULL); d <- copy(sp); d[, y := get(gn)]; cbind(gene = gn, fitone(d, "y"), frac_detect = mean(d$y > 0)) })); gt[, BH := p.adjust(p, "BH")]; fwrite(gt, file.path(out, "H1_all_gene_trend.csv")); print(gt)
## spot-level GIScore across all spots for the map panel (z with the same frozen means/sds from annotated epithelial spots)
allsp <- fread(file.path(root, "02_intermediate/fig6_visium/spots_galnt_annot.csv.gz")); mu16 <- mean(sp$GALNT16); sd16 <- sd(sp$GALNT16); mu14 <- mean(sp$GALNT14); sd14 <- sd(sp$GALNT14)
allsp[, GIScore := b16 * (GALNT16 - mu16) / sd16 + b14 * (GALNT14 - mu14) / sd14]; fwrite(allsp[, .(sample, bc, px_row, px_col, paningrade, GALNT16, GALNT14, GIScore, lib)], file.path(out, "H1_all_spots_map.csv"))
## ---------- H2/H3: Elhossiny epithelial pseudobulk
x <- readRDS("${PANIN_ROOT}/01_raw/Elhossiny2026_processed/pseudobulk_epi.rds")
cn <- as.matrix(GetAssayData(x, assay = "EpithelialExp", layer = "counts")); cpm <- log1p(sweep(cn, 2, colSums(cn), "/") * 1e6)
md <- as.data.table(x@meta.data, keep.rownames = "unit"); md[, cls := factor(annotation, levels = c("Donor_PanIN", "AdjNormal_PanIN", "GLTumor", "PDTumor"))]; md[, ord := as.integer(cls)]
gg <- intersect(genes, rownames(cpm)); E <- data.table(unit = colnames(cpm), t(cpm[gg, ])); E <- merge(md, E, by = "unit")
E[, GIScore := b16 * zs(GALNT16) + b14 * zs(GALNT14)]; E[, c16 := cn["GALNT16", unit]]; E[, c14 := cn["GALNT14", unit]]; E[, libsize := colSums(cn)[unit]]
fwrite(E[, .(unit, patient_or_donor_Id, sample_Id, cls, GALNT16, GALNT14, GIScore, c16, c14, libsize)], file.path(out, "H2_epithelial_pseudobulk_scores.csv")); print(E[, .(unit, cls, GALNT16 = round(GALNT16, 2), GALNT14 = round(GALNT14, 2), GIScore = round(GIScore, 2), c16, c14, libsize)])
jt <- function(y) { r <- jonckheere.test(y, E$ord, alternative = "two.sided", nperm = 20000); data.table(JT_stat = unname(r$statistic), p_jt = r$p.value, rho = suppressWarnings(cor(y, E$ord, method = "spearman")), p_rho = suppressWarnings(cor.test(y, E$ord, method = "spearman", exact = FALSE)$p.value)) }
h2 <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(m) cbind(metric = m, jt(E[[m]])))); fwrite(h2, file.path(out, "H2_trend.csv")); print(h2)
## per-class medians and tumor-only (GL vs PD paired within patient)
ep <- E[cls %in% c("GLTumor", "PDTumor")]; ep[, pat := patient_or_donor_Id]
pr <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(m) { w <- dcast(ep[, .(pat, sample_Id, cls, v = get(m))], pat + sample_Id ~ cls, value.var = "v"); w <- w[complete.cases(w)]; t <- wilcox.test(w$PDTumor, w$GLTumor, paired = TRUE, exact = FALSE); data.table(metric = m, n_pairs = nrow(w), median_diff_PD_minus_GL = median(w$PDTumor - w$GLTumor), p_signed_rank = t$p.value) })); fwrite(pr, file.path(out, "H2_GL_vs_PD_paired.csv")); print(pr)
# normal-like (donor + adj normal) vs tumor
E[, grp := ifelse(cls %in% c("Donor_PanIN", "AdjNormal_PanIN"), "Normal", "Tumor")]
nt <- rbindlist(lapply(c("GALNT16", "GALNT14", "GIScore"), function(m) { t <- wilcox.test(E[grp == "Tumor"][[m]], E[grp == "Normal"][[m]], exact = FALSE); data.table(metric = m, n_normal = E[grp == "Normal", .N], n_tumor = E[grp == "Tumor", .N], median_normal = median(E[grp == "Normal"][[m]]), median_tumor = median(E[grp == "Tumor"][[m]]), p_wilcox = t$p.value) })); fwrite(nt, file.path(out, "H2_normal_vs_tumor.csv")); print(nt)
## all-gene trend (exploratory, BH)
ag <- rbindlist(lapply(gg, function(g) { y <- E[[g]]; data.table(gene = g, rho = suppressWarnings(cor(y, E$ord, method = "spearman")), p = suppressWarnings(cor.test(y, E$ord, method = "spearman", exact = FALSE)$p.value)) })); ag[, BH := p.adjust(p, "BH")]; fwrite(ag, file.path(out, "H2_all_gene_trend.csv")); print(ag)
