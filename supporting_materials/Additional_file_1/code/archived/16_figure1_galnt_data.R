# Figure 1 data for GALNT / mucin-type O-glycosylation (freeze v5; 22 genes). Same rules as figure1_replication_spec_v1.md:
# expression sig = BH<0.05 & |log2FC|>=0.5; OS log-rank on median split (nominal p<0.05); normals code 11 (>=5), GTEx (TOIL) for cancers without TCGA normals.
# Extra rule from freeze v4: gene x cancer is skipped for expression/survival if median log2(norm+1) < 1 (not expressed).
suppressPackageStartupMessages({library(data.table); library(limma); library(survival); library(readxl)})
set.seed(666); root <- "${PROJECT_ROOT}"; raw <- file.path(root, "01_raw")
res <- file.path(root, "04_results/2026-10-12_Figure1_GALNT_run01"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1")
al <- c(GALNTL2 = "GALNT15", GALNTL1 = "GALNT16", WBSCR17 = "GALNT17", GALNTL4 = "GALNT18")
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
cdr <- cdr[, c("bcr_patient_barcode", "type", "OS", "OS.time")]; cdr$OS <- suppressWarnings(as.numeric(cdr$OS)); cdr$OS.time <- suppressWarnings(as.numeric(cdr$OS.time)); ptype <- setNames(cdr$type, cdr$bcr_patient_barcode)

## ---- B, C: mutation (MC3 non-silent gene matrix)
mut <- fread(cmd = paste("gzip -dc", file.path(raw, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz")), data.table = FALSE, na.strings = ""); rownames(mut) <- mut[[1]]; mut[[1]] <- NULL; rownames(mut)[rownames(mut) == "WBSCR17"] <- "GALNT17"
cat("genes missing from MC3 matrix:", setdiff(genes, rownames(mut)), "\n"); gm <- intersect(genes, rownames(mut))
load_ps <- colSums(mut > 0, na.rm = TRUE); M <- matrix(FALSE, length(genes), ncol(mut), dimnames = list(genes, colnames(mut))); M[gm, ] <- as.matrix(mut[gm, ]) > 0; rm(mut)
msamp <- data.frame(sample = colnames(M), patient = substr(colnames(M), 1, 12)); msamp$type <- ptype[msamp$patient]; msamp$load <- load_ps[msamp$sample]
keep <- !is.na(msamp$type) & !duplicated(msamp$patient); msamp <- msamp[keep, ]; M <- M[, msamp$sample]
freq_all <- data.frame(gene = genes, n_mut = rowSums(M), n = ncol(M), pct = 100 * rowMeans(M)); fwrite(freq_all[order(-freq_all$pct), ], file.path(res, "B_mutation_frequency_pancancer.csv"))
ct_levels <- sort(unique(msamp$type)); freq_ct <- do.call(rbind, lapply(ct_levels, function(ct) { j <- msamp$type == ct; data.frame(type = ct, gene = genes, n_samples = sum(j), n_mut = rowSums(M[, j, drop = FALSE]), pct = 100 * rowMeans(M[, j, drop = FALSE])) }))
fwrite(freq_ct, file.path(res, "C_mutation_frequency_by_cancer.csv")); fwrite(aggregate(load ~ type, msamp, mean), file.path(res, "C_mutation_load_proxy_by_cancer.csv")); cat("mutation max pct:", max(freq_all$pct), "\n")

## ---- D: TCGA mRNA tumor vs normal (+GTEx) and OS
xena <- file.path(raw, "TCGA_EBpp_RNASeqV2.geneExp.xena.gz"); hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0)); smp <- hdr[-1]
si <- data.frame(sample = smp, patient = substr(smp, 1, 12), code = substr(smp, 14, 15)); si$type <- ptype[si$patient]; si <- si[!is.na(si$type), ]
prim <- si[si$code == "01" | (si$type == "LAML" & si$code == "03"), ]; prim <- prim[!duplicated(prim$patient), ]; norm <- si[si$code == "11", ]; norm <- norm[!duplicated(norm$patient), ]
nn <- table(norm$type); cancers_expr <- names(nn)[nn >= 5]; cat("cancers with >=5 normals:", length(cancers_expr), "\n")
use <- c(prim$sample, norm$sample[norm$type %in% cancers_expr]); E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(use, hdr)), data.table = FALSE)
gE <- E[[1]]; E <- E[!duplicated(gE) & !is.na(gE), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; E <- as.matrix(E); rownames(E)[rownames(E) %in% names(al)] <- al[rownames(E)[rownames(E) %in% names(al)]]
cat("absent from EB++:", setdiff(genes, rownames(E)), "\n")
tn_list <- list()
for (ct in cancers_expr) { t_s <- prim$sample[prim$type == ct]; n_s <- norm$sample[norm$type == ct]; X <- E[, c(n_s, t_s)]; X <- X[rowSums(!is.finite(X)) == 0, ]; X <- X[apply(X, 1, var) > 0, ]
  g <- factor(c(rep("N", length(n_s)), rep("T", length(t_s))), levels = c("N", "T")); fit <- eBayes(lmFit(X, model.matrix(~ g)), trend = TRUE); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none"); tt$BH <- p.adjust(tt$P.Value, "BH")
  i <- match(genes, rownames(tt)); ex <- sapply(genes, function(gn) if (gn %in% rownames(X)) max(median(X[gn, t_s]), median(X[gn, n_s])) else NA)
  tn_list[[ct]] <- data.frame(type = ct, gene = genes, log2FC = tt$logFC[i], p = tt$P.Value[i], BH = tt$BH[i], n_tumor = length(t_s), n_normal = length(n_s), expressed = ex >= 1) }
tn <- rbindlist(tn_list); tn$sig <- ifelse(tn$expressed %in% TRUE & tn$BH < 0.05 & tn$log2FC >= 0.5, "up", ifelse(tn$expressed %in% TRUE & tn$BH < 0.05 & tn$log2FC <= -0.5, "down", "ns")); tn$ref <- "TCGA normal"
fwrite(tn, file.path(res, "D_tumor_vs_normal_limma.csv"))
## GTEx (TOIL) for cancers without >=5 TCGA normals; mapping fixed as in 05_panelD_gtex_toil.R
ph <- fread(file.path(raw, "toil/TcgaTargetGTEX_phenotype.txt.gz"), data.table = FALSE, encoding = "Latin-1"); XT <- as.matrix(read.csv(file.path(raw, "toil/toil_galnt22_log2norm.csv"), row.names = 1, check.names = FALSE)); ph <- ph[match(rownames(XT), ph$sample), ]
tcga_name <- c(ACC = "Adrenocortical Cancer", CESC = "Cervical & Endocervical Cancer", GBM = "Glioblastoma Multiforme", LAML = "Acute Myeloid Leukemia", LGG = "Brain Lower Grade Glioma", OV = "Ovarian Serous Cystadenocarcinoma", PAAD = "Pancreatic Adenocarcinoma", SKCM = "Skin Cutaneous Melanoma", UCS = "Uterine Carcinosarcoma")
gtex_cat <- list(ACC = "Adrenal Gland", CESC = "Cervix - E", GBM = "Brain - ", LAML = "Whole Blood", LGG = "Brain - ", OV = "^Ovary$", PAAD = "^Pancreas$", SKCM = "Skin - ", UCS = "^Uterus$")
runT <- function(ct) { tum <- ph$sample[ph$`_study` == "TCGA" & ph$detailed_category == tcga_name[ct] & ph$`_sample_type` %in% c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood")]; nor <- ph$sample[ph$`_study` == "GTEX" & grepl(gtex_cat[[ct]], ph$detailed_category)]
  g <- factor(c(rep("N", length(nor)), rep("T", length(tum))), levels = c("N", "T")); Y <- t(XT[c(nor, tum), genes]); Y[!is.finite(Y)] <- NA
  fit <- eBayes(lmFit(Y, model.matrix(~ g))); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none"); tt$BH <- p.adjust(tt$P.Value, "BH"); i <- match(genes, rownames(tt))
  ex <- sapply(genes, function(gn) max(median(XT[tum, gn], na.rm = TRUE), median(XT[nor, gn], na.rm = TRUE)))
  data.frame(type = ct, gene = genes, log2FC = tt$logFC[i], p = tt$P.Value[i], BH = tt$BH[i], n_tumor = length(tum), n_normal = length(nor), expressed = ex >= 1) }
new_ct <- c("ACC", "CESC", "GBM", "LAML", "LGG", "OV", "PAAD", "SKCM", "UCS"); gt <- rbindlist(lapply(new_ct, runT))
gt$sig <- ifelse(gt$expressed & gt$BH < 0.05 & gt$log2FC >= 0.5, "up", ifelse(gt$expressed & gt$BH < 0.05 & gt$log2FC <= -0.5, "down", "ns")); gt$ref <- "GTEx normal"
fwrite(gt, file.path(res, "D_tumor_vs_normal_GTEx_TOIL.csv")); print(unique(gt[, .(type, n_tumor, n_normal)]))
tn <- rbind(tn, gt[, names(tn), with = FALSE]); fwrite(tn, file.path(res, "D_tumor_vs_normal_combined.csv"))
ratio <- tn[, .(up_ratio = mean(sig == "up"), down_ratio = mean(sig == "down"), n_cancers = .N), by = gene]; ratio$star <- ifelse(ratio$up_ratio > 0.25, "up", ifelse(ratio$down_ratio > 0.25, "down", "none")); fwrite(ratio, file.path(res, "D_up_down_ratio.csv")); print(ratio[star != "none"])
## survival (median split among expressed gene-cancer pairs)
sv_list <- list(); prim$time <- cdr$OS.time[match(prim$patient, cdr$bcr_patient_barcode)]; prim$event <- cdr$OS[match(prim$patient, cdr$bcr_patient_barcode)]
for (ct in sort(unique(prim$type))) { p <- prim[prim$type == ct & is.finite(prim$time) & prim$time > 0 & prim$event %in% c(0, 1) & prim$sample %in% colnames(E), ]; if (nrow(p) < 30 || sum(p$event) < 10) next
  for (g in genes) { if (!(g %in% rownames(E))) next; x <- E[g, p$sample]; x <- x[is.finite(x)]; p2 <- p[p$sample %in% names(x), ]; if (length(x) < 30 || median(x) < 1) next; p <- p2; hi <- x > median(x); if (!any(hi) || all(hi)) next
    lr <- survdiff(Surv(p$time, p$event) ~ hi); pv <- pchisq(lr$chisq, 1, lower.tail = FALSE); hr <- exp(coef(coxph(Surv(p$time, p$event) ~ hi)))
    sv_list[[paste(ct, g)]] <- data.frame(type = ct, gene = g, n = nrow(p), events = sum(p$event), HR_high_vs_low = hr, p_logrank = pv) } }
sv <- rbindlist(sv_list); sv$call <- ifelse(sv$p_logrank < 0.05 & sv$HR_high_vs_low > 1, "Risk", ifelse(sv$p_logrank < 0.05, "Protective", "n.s.")); fwrite(sv, file.path(res, "D_os_logrank_median_split.csv"))
cat("survival tests:", nrow(sv), " nominal sig:", sum(sv$p_logrank < 0.05), "\n"); print(sv[type == "PAAD" & p_logrank < 0.05])
