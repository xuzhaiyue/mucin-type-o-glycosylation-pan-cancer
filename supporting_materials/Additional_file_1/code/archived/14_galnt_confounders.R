# GALNT cluster confounder checks (freeze v5). Cox: cluster + clinical (age, sex, stage) [M1]; + expression proxies of stroma/immune/acinar and LGG IDH/grade [M2].
suppressPackageStartupMessages({library(data.table); library(readxl); library(survival)})
root <- "${PROJECT_ROOT}"; res <- file.path(root, "04_results/2026-10-12_Figure3_multisystem_run01/GALNT")
cls <- fread(file.path(res, "cluster_assignments.csv.gz")); show <- c("LGG","KIRP","KIRC","PAAD","BLCA","UCEC","HNSC","COAD","STAD","LIHC")
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
cdr$age <- suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis)); cdr$stage <- suppressWarnings(as.numeric(as.roman(sub("^Stage ([IV]+).*", "\\1", cdr$ajcc_pathologic_tumor_stage))))
cdr$stage[!grepl("^Stage [IV]+", cdr$ajcc_pathologic_tumor_stage)] <- NA; cdr$grade <- ifelse(cdr$histological_grade %in% c("G2", "G3"), cdr$histological_grade, NA)
xena <- file.path(root, "01_raw/TCGA_EBpp_RNASeqV2.geneExp.xena.gz"); hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0))
sg <- list(stromal = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP"), immune = c("PTPRC","CD3E","CD8A","CD68","CD4","CD19"), acinar = c("PRSS1","CPA1","CELA3A","CTRB1","PNLIP","CEL"))
need <- unique(unlist(sg)); keep <- cls[type %in% show]; E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(keep$sample, hdr)), data.table = FALSE); E <- E[!duplicated(E[[1]]), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; E <- as.matrix(E[need, , drop = FALSE])
mc3 <- fread(cmd = paste("gzip -dc", file.path(root, "01_raw/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz")), select = 1:9105, data.table = FALSE, na.strings = ""); rownames(mc3) <- mc3[[1]]; mc3[[1]] <- NULL
idh <- (unlist(mc3["IDH1", ]) + unlist(mc3["IDH2", ])) > 0; names(idh) <- substr(colnames(mc3), 1, 12)
out <- list(); assoc <- list()
for (ct in show) {
  d <- as.data.frame(cls[type == ct]); d$C2 <- as.numeric(d$C == "C2"); m <- merge(d, cdr[, c("bcr_patient_barcode", "age", "gender", "stage", "grade")], by.x = "patient", by.y = "bcr_patient_barcode")
  Zs <- t(scale(t(E[, m$sample]))); for (s in names(sg)) m[[s]] <- as.numeric(scale(colMeans(Zs[sg[[s]], , drop = FALSE], na.rm = TRUE)))
  m$idh <- as.numeric(idh[m$patient]); m$male <- as.numeric(m$gender == "MALE")
  fm <- function(extra) { v <- c("C2", extra); v <- v[sapply(v, function(x) x == "C2" || (mean(is.na(m[[x]])) < 0.2 && length(unique(na.omit(m[[x]]))) > 1))]; f <- as.formula(paste("Surv(time, event) ~", paste(v, collapse = " + "))); fit <- coxph(f, m); s <- summary(fit)$coefficients["C2", ]; list(n = fit$n, vars = paste(v[-1], collapse = "+"), HR = unname(s["exp(coef)"]), p = unname(s["Pr(>|z|)"])) }
  b <- fm(character(0)); m1 <- fm(c("age", "male", "stage")); ex <- c("age", "male", "stage", "stromal", "immune"); if (ct == "PAAD") ex <- c(ex, "acinar"); if (ct == "LGG") ex <- c(ex, "idh", "grade"); m2 <- fm(ex)
  out[[ct]] <- data.table(type = ct, n_unadj = b$n, HR_unadj = b$HR, p_unadj = b$p, M1_vars = m1$vars, n_M1 = m1$n, HR_M1 = m1$HR, p_M1 = m1$p, M2_vars = m2$vars, n_M2 = m2$n, HR_M2 = m2$HR, p_M2 = m2$p)
  a <- sapply(c("stromal", "immune", if (ct == "PAAD") "acinar", if (ct == "LGG") "idh"), function(s) { x <- m[[s]]; if (s == "idh") { t <- table(m$C2, x); c(stat = NA, p = fisher.test(t)$p.value, extra = paste0("IDHmut in C1 ", round(mean(x[m$C2 == 0], na.rm = TRUE), 2), " vs C2 ", round(mean(x[m$C2 == 1], na.rm = TRUE), 2))) } else c(stat = mean(x[m$C2 == 1], na.rm = TRUE) - mean(x[m$C2 == 0], na.rm = TRUE), p = wilcox.test(x ~ m$C2)$p.value, extra = "") })
  assoc[[ct]] <- data.table(type = ct, covariate = colnames(a), diff_C2_minus_C1_z = suppressWarnings(as.numeric(a["stat", ])), wilcox_or_fisher_p = as.numeric(a["p", ]), note = a["extra", ])
}
out <- rbindlist(out); assoc <- rbindlist(assoc); fwrite(out, file.path(res, "confounder_cox.csv")); fwrite(assoc, file.path(res, "confounder_cluster_association.csv"))
options(width = 200); print(out[, .(type, n_unadj, HR_unadj = round(HR_unadj, 2), p_unadj = signif(p_unadj, 2), HR_M1 = round(HR_M1, 2), p_M1 = signif(p_M1, 2), HR_M2 = round(HR_M2, 2), p_M2 = signif(p_M2, 2), M2_vars)]); print(assoc[, .(type, covariate, diff = round(diff_C2_minus_C1_z, 2), p = signif(wilcox_or_fisher_p, 2), note)])
