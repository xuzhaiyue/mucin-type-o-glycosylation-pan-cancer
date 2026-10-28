# Figure 4 step 2: apply the FROZEN model to external PAAD cohorts (opened only now). Spec: 00_meta/figure4_spec_v1.md
suppressPackageStartupMessages({library(data.table); library(survival)})
source("${PROJECT_ROOT}/03_analysis/20a_models_lib.R"); set.seed(666)
root <- "${PROJECT_ROOT}"; res <- file.path(root, "04_results/2026-10-12_Figure4_GALNT_run01"); geo <- file.path(root, "01_raw/geo_paad")
fz <- readRDS(file.path(res, "FROZEN_model.rds")); chk <- system(paste("shasum -a 256", shQuote(file.path(res, "FROZEN_model.rds"))), intern = TRUE); stopifnot(substr(chk, 1, 64) == substr(readLines(file.path(res, "FROZEN_model_sha256.txt"))[1], 1, 64)); cat("frozen model hash verified;", fz$selected, "genes:", fz$model$genes, "\n")
genes <- fz$genes; al <- c(GALNTL2 = "GALNT15", GALNTL1 = "GALNT16", WBSCR17 = "GALNT17", GALNTL4 = "GALNT18")
parse_gse <- function(gse, plat) {
  f <- file.path(geo, paste0(gse, "_series_matrix.txt.gz")); lines <- readLines(gzfile(f)); ch <- lines[grepl("^!Sample_characteristics_ch1", lines)]; acc <- gsub("\"", "", strsplit(lines[grepl("^!Sample_geo_accession", lines)], "\t")[[1]][-1])
  cl <- as.data.frame(lapply(ch, function(x) gsub("\"", "", strsplit(x, "\t")[[1]][-1])), stringsAsFactors = FALSE); cl <- as.data.frame(lapply(cl, function(v) v)); names(cl) <- sapply(ch, function(x) sub(":.*", "", gsub("\"", "", strsplit(x, "\t")[[1]][2])))
  for (j in seq_along(cl)) cl[[j]] <- sub("^[^:]*: ?", "", cl[[j]]); cl$gsm <- acc
  b <- which(grepl("^!series_matrix_table_begin", lines)); e <- which(grepl("^!series_matrix_table_end", lines)); M <- read.delim(text = lines[(b + 1):(e - 1)], check.names = FALSE, row.names = 1); colnames(M) <- gsub("\"", "", colnames(M)); rownames(M) <- gsub("\"", "", rownames(M))
  ann <- fread(file.path(geo, paste0(plat, "_probe2symbol.tsv")), colClasses = "character"); ann <- ann[probe %in% rownames(M)]; ann <- ann[, .(symbol = trimws(unlist(strsplit(gsub("///", "//", symbol), "//")))), by = probe]; ann$symbol[ann$symbol %in% names(al)] <- al[ann$symbol[ann$symbol %in% names(al)]]; ann <- ann[symbol %in% genes]
  mm <- rowMeans(M[ann$probe, , drop = FALSE], na.rm = TRUE); ann$mm <- mm[ann$probe]; ann <- ann[order(-mm)][!duplicated(symbol)]; G <- t(M[ann$probe, , drop = FALSE]); colnames(G) <- ann$symbol; list(clin = cl, G = G, range = range(as.matrix(M), na.rm = TRUE)) }
cohorts <- list()
x <- parse_gse("GSE62452", "GPL6244"); cl <- x$clin; cl <- cl[cl$tissue == "Pancreatic tumor", ]; cl$time <- suppressWarnings(as.numeric(cl[["survival months"]])); cl$event <- suppressWarnings(as.numeric(cl[["survival status"]])); cl$stage <- cl$Stage; cl$grade <- cl$grading; cohorts$GSE62452 <- list(cl = cl[is.finite(cl$time) & cl$time > 0 & cl$event %in% c(0, 1), ], G = x$G, range = x$range)
x <- parse_gse("GSE57495", "GPL15048"); cl <- x$clin; cl$time <- suppressWarnings(as.numeric(cl[["overall survival (month)"]])); cl$event <- as.numeric(toupper(cl[["vital.status"]]) == "DEAD"); cl$stage <- cl$Stage; cohorts$GSE57495 <- list(cl = cl[is.finite(cl$time) & cl$time > 0, ], G = x$G, range = x$range); print(names(cl))
x <- parse_gse("GSE85916", "GPL13667"); cl <- x$clin; cl$time <- suppressWarnings(as.numeric(cl[["os.year"]])) * 12; cl$event <- suppressWarnings(as.numeric(cl[["death"]])); cohorts$GSE85916 <- list(cl = cl[is.finite(cl$time) & cl$time > 0 & cl$event %in% c(0, 1), ], G = x$G, range = x$range)
for (n in names(cohorts)) cat(n, "n =", nrow(cohorts[[n]]$cl), "events =", sum(cohorts[[n]]$cl$event), "value range", round(cohorts[[n]]$range, 1), "genes on platform:", ncol(cohorts[[n]]$G), "\n")
## frozen score on each cohort
score_cohort <- function(n) { co <- cohorts[[n]]; G <- co$G[co$cl$gsm, , drop = FALSE]; miss <- setdiff(genes, colnames(G)); for (m in miss) G <- cbind(G, setNames(rep(0, nrow(G)), m)); Zc <- scale(G[, genes, drop = FALSE]); Zc[is.na(Zc)] <- 0; nd <- data.frame(Zc, check.names = FALSE); co$cl$risk <- fz$model$pred(nd); co$cl$cohort <- n; co$cl }
sc <- rbindlist(lapply(names(cohorts), function(n) { d <- score_cohort(n); d$group <- ifelse(d$risk > median(d$risk), "High", "Low"); as.data.table(d)[, .(cohort, gsm, time, event, risk, group, stage = if ("stage" %in% names(d)) stage else NA_character_, grade = if ("grade" %in% names(d)) grade else NA_character_)] }), fill = TRUE)
tr <- as.data.table(fz$train); tr$cohort <- "TCGA-PAAD"; tr$group <- ifelse(tr$risk > median(tr$risk), "High", "Low")
cat("risk summaries:\n"); print(sc[, .(n = .N, nNA = sum(is.na(risk)), sd = sd(risk, na.rm = TRUE), nHigh = sum(group == "High")), by = cohort])
fwrite(sc, file.path(res, "validation_scores.csv")); fwrite(tr[, .(cohort, sample, time, event, risk, group)], file.path(res, "training_scores.csv"))
summ <- rbindlist(lapply(c("TCGA-PAAD", names(cohorts)), function(n) { d <- as.data.frame(if (n == "TCGA-PAAD") tr else sc[cohort == n]); d$hi <- as.numeric(d$group == "High"); lr <- survdiff(Surv(time, event) ~ hi, d); p <- 1 - pchisq(lr$chisq, 1); cx <- summary(coxph(Surv(time, event) ~ hi, d))$coefficients; cz <- summary(coxph(Surv(time, event) ~ I(scale(risk)[, 1]), d))$coefficients
  data.table(cohort = n, n = nrow(d), events = sum(d$event), logrank_p = p, logHR = cx[1], se = cx[3], HR_high_vs_low = exp(cx[1]), HR_per_SD = cz[2], p_per_SD = cz[5], cindex = cidx(d$time, d$event, d$risk)) }))
v <- summ[cohort != "TCGA-PAAD"]; w <- 1 / v$se^2; fe <- sum(w * v$logHR) / sum(w); Q <- sum(w * (v$logHR - fe)^2); tau2 <- max(0, (Q - (nrow(v) - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (v$se^2 + tau2); re <- sum(wr * v$logHR) / sum(wr); sre <- sqrt(1 / sum(wr))
pool <- data.table(pooled_HR = exp(re), lo = exp(re - 1.96 * sre), hi = exp(re + 1.96 * sre), p = 2 * pnorm(-abs(re / sre)), I2 = max(0, (Q - (nrow(v) - 1)) / Q), same_direction = sum(v$logHR > 0), n_cohorts = nrow(v), nominal_sig_same_dir = sum(v$logHR > 0 & v$logrank_p < 0.05))
fwrite(summ, file.path(res, "validation_summary.csv")); fwrite(pool, file.path(res, "validation_pooled.csv")); options(width = 200); print(summ); print(pool)
cat("PRIMARY CRITERION (pooled HR>1, p<0.05, >=2/3 same direction):", pool$pooled_HR > 1 & pool$p < 0.05 & pool$same_direction >= 2, "\n")
## all 16 models: refit on full training (transparency only; selection was done on training CV) and evaluate in the validation cohorts
D <- data.frame(scale(t(readRDS(file.path(res, "tcga_paad_ductal_expr_allgenes.rds"))[genes, fz$train$sample])), check.names = FALSE); D$time <- fz$train$time; D$event <- fz$train$event
allm <- rbindlist(lapply(names(models), function(m) { o <- tryCatch(models[[m]](D, genes), error = function(e) NULL); if (is.null(o)) return(NULL)
  rbindlist(lapply(names(cohorts), function(n) { d <- cohorts[[n]]$cl; G <- cohorts[[n]]$G[d$gsm, , drop = FALSE]; for (mm in setdiff(genes, colnames(G))) G <- cbind(G, setNames(rep(0, nrow(G)), mm)); Zc <- scale(G[, genes, drop = FALSE]); Zc[is.na(Zc)] <- 0; r <- tryCatch(o$pred(data.frame(Zc, check.names = FALSE)), error = function(e) rep(NA, nrow(d))); data.table(model = m, cohort = n, cindex = cidx(d$time, d$event, as.numeric(r))) })) }))
allm <- rbind(allm, fz$cvs[, .(model, cohort = "TCGA-CV", cindex = mean)]); fwrite(allm, file.path(res, "B_all_models_cindex.csv")); print(dcast(allm, model ~ cohort, value.var = "cindex"))
## export clinical/validation covariates for panels H-J
fwrite(sc, file.path(res, "validation_scores.csv"))
