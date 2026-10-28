# Figure 4 step 1: build and select the GALNT score in TCGA-PAAD (ductal) ONLY. Spec: 00_meta/figure4_spec_v1.md. No external cohort is read here.
suppressPackageStartupMessages({library(data.table); library(readxl); library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(CoxBoost); library(parallel)})
set.seed(666); root <- "${PROJECT_ROOT}"; res <- file.path(root, "04_results/2026-10-12_Figure4_GALNT_run01"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
genes0 <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"); al <- c(GALNTL2 = "GALNT15", GALNTL1 = "GALNT16", WBSCR17 = "GALNT17", GALNTL4 = "GALNT18")
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
cl <- cdr[cdr$type == "PAAD" & cdr$histological_type == "Pancreas-Adenocarcinoma Ductal Type", ]; cl$time <- suppressWarnings(as.numeric(cl$OS.time)) / 30.4375; cl$event <- suppressWarnings(as.numeric(cl$OS)); cl <- cl[is.finite(cl$time) & cl$time > 0 & cl$event %in% c(0, 1), ]
xena <- file.path(root, "01_raw/TCGA_EBpp_RNASeqV2.geneExp.xena.gz"); hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0)); smp <- hdr[-1]
si <- data.frame(sample = smp, patient = substr(smp, 1, 12), code = substr(smp, 14, 15)); si <- si[si$code == "01" & si$patient %in% cl$bcr_patient_barcode & !duplicated(si$patient), ]
E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(si$sample, hdr)), data.table = FALSE); g0 <- E[[1]]; E <- E[!duplicated(g0) & !is.na(g0), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; E <- as.matrix(E)
rownames(E)[rownames(E) %in% names(al)] <- al[rownames(E)[rownames(E) %in% names(al)]]
saveRDS(E[rowSums(!is.finite(E)) == 0, ], file.path(res, "tcga_paad_ductal_expr_allgenes.rds"))
dat <- data.frame(patient = si$patient, sample = si$sample); dat <- merge(dat, cl[, c("bcr_patient_barcode", "time", "event")], by.x = "patient", by.y = "bcr_patient_barcode"); cat("training n =", nrow(dat), "events =", sum(dat$event), "\n")
gx <- genes0[genes0 %in% rownames(E)]; gx <- gx[apply(E[gx, dat$sample], 1, median) >= 1]; genes <- setdiff(gx, "GALNT4"); cat("candidate genes (", length(genes), "):", genes, "\n")
Xr <- t(E[genes, dat$sample]); Z <- scale(Xr); zc <- attr(Z, "scaled:center"); zs <- attr(Z, "scaled:scale"); D <- data.frame(Z, check.names = FALSE); D$time <- dat$time; D$event <- dat$event
cidx <- function(time, event, risk) unname(concordance(Surv(time, event) ~ risk, reverse = TRUE)$concordance)
## univariate Cox (panel A)
uni <- rbindlist(lapply(genes, function(g) { f <- summary(coxph(Surv(time, event) ~ D[[g]], D))$coefficients; ci <- exp(confint(coxph(Surv(time, event) ~ D[[g]], D))); data.table(gene = g, HR = f[2], lo = ci[1], hi = ci[2], p = f[5]) })); uni$BH <- p.adjust(uni$p, "BH"); fwrite(uni, file.path(res, "A_univariate_cox_training.csv")); print(uni[order(p)][1:6])
## model library: each returns a function(newdata) -> risk
pre_uni <- function(d, gs) { p <- sapply(gs, function(g) summary(coxph(Surv(time, event) ~ d[[g]], d))$coefficients[5]); k <- names(p)[p < 0.05]; if (length(k) < 2) k <- names(sort(p))[1:2]; k }
f_glmnet <- function(a) function(d, gs) { x <- as.matrix(d[, gs, drop = FALSE]); cv <- cv.glmnet(x, Surv(d$time, d$event), family = "cox", alpha = a, nfolds = 10); m <- glmnet(x, Surv(d$time, d$event), family = "cox", alpha = a, lambda = cv$lambda.min); list(pred = function(nd) as.numeric(predict(m, as.matrix(nd[, gs, drop = FALSE]), type = "link")), genes = gs, cv = cv, fit = m) }
f_rsf <- function(d, gs) { m <- rfsrc(Surv(time, event) ~ ., data = d[, c(gs, "time", "event")], ntree = 500, nodesize = 5, seed = -666, importance = "permute"); list(pred = function(nd) predict(m, nd[, gs, drop = FALSE])$predicted, genes = gs, fit = m) }
f_gbm <- function(d, gs) { m <- gbm(Surv(time, event) ~ ., data = d[, c(gs, "time", "event")], distribution = "coxph", n.trees = 500, interaction.depth = 2, shrinkage = 0.02, cv.folds = 3, n.minobsinnode = 8, n.cores = 1, verbose = FALSE); bi <- max(10, gbm.perf(m, method = "cv", plot.it = FALSE)); list(pred = function(nd) predict(m, nd[, gs, drop = FALSE], n.trees = bi), genes = gs, fit = m) }
f_cb <- function(d, gs) { x <- as.matrix(d[, gs, drop = FALSE]); m <- CoxBoost(d$time, d$event, x, stepno = 100, penalty = 100, criterion = "pscore"); cvm <- cv.CoxBoost(d$time, d$event, x, maxstepno = 100, K = 5, penalty = 100, type = "verweij"); m <- CoxBoost(d$time, d$event, x, stepno = max(1, cvm$optimal.step), penalty = 100); list(pred = function(nd) as.numeric(predict(m, as.matrix(nd[, gs, drop = FALSE]), type = "lp")), genes = gs, fit = m) }
f_step <- function(d, gs) { f0 <- coxph(as.formula(paste("Surv(time, event) ~", paste(sprintf("`%s`", gs), collapse = "+"))), d); m <- suppressWarnings(step(f0, direction = "both", trace = 0)); list(pred = function(nd) as.numeric(predict(m, nd, type = "lp")), genes = names(coef(m)), fit = m) }
lasso_sel <- function(d, gs) { o <- f_glmnet(1)(d, gs); cf <- as.matrix(coef(o$fit)); k <- rownames(cf)[cf[, 1] != 0]; if (length(k) < 2) gs else k }
base <- list(Lasso = f_glmnet(1), Ridge = f_glmnet(0), ElasticNet = f_glmnet(0.5), CoxBoost = f_cb, RSF = f_rsf, GBM = f_gbm, StepCox = f_step)
models <- list(); for (n in names(base)) models[[n]] <- local({ f <- base[[n]]; function(d, gs) f(d, gs) }); for (n in names(base)) models[[paste0("Uni+", n)]] <- local({ f <- base[[n]]; function(d, gs) f(d, pre_uni(d, gs)) })
models[["Lasso+RSF"]] <- function(d, gs) f_rsf(d, lasso_sel(d, gs)); models[["Lasso+GBM"]] <- function(d, gs) f_gbm(d, lasso_sel(d, gs))
simp <- c(Lasso = 1, Ridge = 2, ElasticNet = 3, StepCox = 4, "Uni+Lasso" = 5, "Uni+Ridge" = 6, "Uni+ElasticNet" = 7, "Uni+StepCox" = 8, CoxBoost = 9, "Uni+CoxBoost" = 10, RSF = 11, "Uni+RSF" = 12, GBM = 13, "Uni+GBM" = 14, "Lasso+RSF" = 15, "Lasso+GBM" = 16)
## repeated 5-fold CV inside training only
REP <- 10; K <- 5; job <- expand.grid(model = names(models), rep = 1:REP, stringsAsFactors = FALSE)
run_job <- function(i) { set.seed(1000 + job$rep[i]); fold <- sample(rep(1:K, length.out = nrow(D))); pr <- rep(NA_real_, nrow(D)); cs <- c()
  for (k in 1:K) { tr <- D[fold != k, ]; te <- D[fold == k, ]; o <- tryCatch(models[[job$model[i]]](tr, genes), error = function(e) NULL); if (is.null(o)) next; r <- tryCatch(o$pred(te), error = function(e) NULL); if (!is.null(r)) { pr[fold == k] <- as.numeric(r); cs <- c(cs, cidx(te$time, te$event, as.numeric(r))) } }
  data.frame(model = job$model[i], rep = job$rep[i], cindex = mean(cs, na.rm = TRUE)) }
t0 <- Sys.time(); out <- mclapply(seq_len(nrow(job)), run_job, mc.cores = 4, mc.set.seed = FALSE); cv <- rbindlist(out); cat("CV done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cvs <- cv[, .(mean = mean(cindex, na.rm = TRUE), sd = sd(cindex, na.rm = TRUE)), by = model][order(-mean)]; fwrite(cv, file.path(res, "B_cv_cindex_repeats.csv")); fwrite(cvs, file.path(res, "B_cv_cindex_summary.csv")); print(cvs)
cand <- cvs[mean >= max(mean) - 0.005]$model; sel <- cand[which.min(simp[cand])]; cat("selected model:", sel, " (candidates within 0.005:", paste(cand, collapse = ", "), ")\n")
final <- models[[sel]](D, genes); risk_tr <- final$pred(D); cat("training apparent C-index:", cidx(D$time, D$event, risk_tr), "\n")
diag <- list(lasso = f_glmnet(1)(D, genes), rsf = f_rsf(D, genes))
saveRDS(list(selected = sel, genes = genes, center = zc, scale = zs, model = final, diag = diag, train = data.frame(dat, risk = risk_tr), cvs = cvs, models_all = names(models), uni = uni), file.path(res, "FROZEN_model.rds"))
fwrite(data.frame(gene = final$genes), file.path(res, "FROZEN_model_genes.csv")); if (!is.null(final$fit) && inherits(final$fit, "glmnet")) fwrite(data.table(gene = rownames(coef(final$fit)), coef = as.numeric(coef(final$fit))), file.path(res, "FROZEN_model_coef.csv"))
sh <- system(paste("shasum -a 256", shQuote(file.path(res, "FROZEN_model.rds"))), intern = TRUE); writeLines(c(sh, paste("frozen at", format(Sys.time())), paste("selected:", sel), paste("genes:", paste(final$genes, collapse = ","))), file.path(res, "FROZEN_model_sha256.txt")); print(sh)
