# Fig 4 sensitivity (declared BEFORE running): how does the univariate prefilter threshold change the CoxBoost model?
# The frozen GIScore (threshold p < 0.05; spec figure4_spec_v1.md) is NOT replaced. Grid fixed in advance: 0.05 (frozen), 0.10, 0.20, 0.30.
# Selection of any threshold by external C-index is not allowed; all four are reported. Training-only repeated 5-fold CV (prefilter re-done inside each fold).
suppressPackageStartupMessages({library(data.table); library(survival); library(parallel)})
root <- "${PROJECT_ROOT}"; res <- file.path(root, "04_results/2026-10-12_Figure4_GALNT_run01"); out <- file.path(root, "04_results/2026-10-13_Figure4_uni_threshold_sensitivity"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
src21 <- readLines(file.path(root, "03_analysis/21_figure4_validate.R")); k <- grep("^## frozen score on each cohort", src21)
eval(parse(text = src21[1:(k - 1)]))                      # loads 20a_models_lib.R, verifies the frozen model hash, parses the 3 GEO cohorts (read only)
set.seed(666)
E <- readRDS(file.path(res, "tcga_paad_ductal_expr_allgenes.rds")); tr <- fz$train
Xr <- t(E[genes, tr$sample]); Z <- scale(Xr); D <- data.frame(Z, check.names = FALSE); D$time <- tr$time; D$event <- tr$event
thr <- c(0.05, 0.10, 0.20, 0.30)
pre_t <- function(d, gs, t) { p <- sapply(gs, function(g) summary(coxph(Surv(time, event) ~ d[[g]], d))$coefficients[5]); kk <- names(p)[p < t]; if (length(kk) < 2) kk <- names(sort(p))[1:2]; kk }
mods <- lapply(thr, function(t) local({ tt <- t; function(d, gs) f_cb(d, pre_t(d, gs, tt)) })); names(mods) <- paste0("Uni", thr, "+CoxBoost")
## genes selected on the full training set
selg <- rbindlist(lapply(seq_along(thr), function(i) { gs <- pre_t(D, genes, thr[i]); m <- mods[[i]](D, genes); data.table(threshold = thr[i], n_prefilter = length(gs), prefilter_genes = paste(gs, collapse = ","), n_final = length(m$fit$xnames[m$fit$coefficients[nrow(m$fit$coefficients), ] != 0]), final_genes = paste(m$fit$xnames[m$fit$coefficients[nrow(m$fit$coefficients), ] != 0], collapse = ","), train_apparent_C = cidx(D$time, D$event, m$pred(D))) }))
print(selg); fwrite(selg, file.path(out, "selected_genes_by_threshold.csv"))
## training-only repeated CV
REP <- 10; K <- 5; job <- expand.grid(model = names(mods), rep = 1:REP, stringsAsFactors = FALSE)
run_job <- function(i) { set.seed(1000 + job$rep[i]); fold <- sample(rep(1:K, length.out = nrow(D))); cs <- c()
  for (kf in 1:K) { trn <- D[fold != kf, ]; te <- D[fold == kf, ]; o <- tryCatch(mods[[job$model[i]]](trn, genes), error = function(e) NULL); if (is.null(o)) next; r <- tryCatch(o$pred(te), error = function(e) NULL); if (!is.null(r)) cs <- c(cs, cidx(te$time, te$event, as.numeric(r))) }
  data.frame(model = job$model[i], rep = job$rep[i], cindex = mean(cs, na.rm = TRUE)) }
cv <- rbindlist(mclapply(seq_len(nrow(job)), run_job, mc.cores = 4, mc.set.seed = FALSE)); cvs <- cv[, .(mean = mean(cindex, na.rm = TRUE), sd = sd(cindex, na.rm = TRUE)), by = model]; print(cvs); fwrite(cvs, file.path(out, "cv_cindex_by_threshold.csv"))
## external cohorts (secondary; reported for all thresholds, not used to choose)
ext <- rbindlist(lapply(seq_along(thr), function(i) { m <- mods[[i]](D, genes); gg <- m$fit$xnames
  rbindlist(lapply(names(cohorts), function(n) { co <- cohorts[[n]]; G <- co$G[co$cl$gsm, , drop = FALSE]; for (mm in setdiff(genes, colnames(G))) G <- cbind(G, setNames(rep(0, nrow(G)), mm)); Zc <- scale(G[, genes, drop = FALSE]); Zc[is.na(Zc)] <- 0; nd <- data.frame(Zc, check.names = FALSE); d <- co$cl; d$risk <- m$pred(nd)
    d$hi <- as.numeric(d$risk > median(d$risk)); cx <- summary(coxph(Surv(time, event) ~ hi, d))$coefficients; cz <- summary(coxph(Surv(time, event) ~ I(scale(risk)[, 1]), d))$coefficients
    data.table(threshold = thr[i], cohort = n, n = nrow(d), events = sum(d$event), cindex = cidx(d$time, d$event, d$risk), HR_high_vs_low = exp(cx[1]), logHR = cx[1], se = cx[3], p = cx[5], HR_per_SD = cz[2], p_per_SD = cz[5]) })) }))
pool <- ext[, { w <- 1 / se^2; fe <- sum(w * logHR) / sum(w); Q <- sum(w * (logHR - fe)^2); tau2 <- max(0, (Q - (.N - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + tau2); re <- sum(wr * logHR) / sum(wr); s <- sqrt(1 / sum(wr)); .(pooled_HR = exp(re), lo = exp(re - 1.96 * s), hi = exp(re + 1.96 * s), p = 2 * pnorm(-abs(re / s)), same_dir = sum(logHR > 0)) }, by = threshold]
options(width = 200); print(ext[, .(threshold, cohort, n, cindex = round(cindex, 3), HR = round(HR_high_vs_low, 2), p = round(p, 3))]); print(pool)
fwrite(ext, file.path(out, "external_by_threshold.csv")); fwrite(pool, file.path(out, "external_pooled_by_threshold.csv")); cat("sensitivity done\n")
