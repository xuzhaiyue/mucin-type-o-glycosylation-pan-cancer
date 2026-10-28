# Figure 3 multi-system clustering (freeze v4). Usage: Rscript 13_figure3_multisystem.R NEDD8 SUMO OST GALNT GPI UFM1
suppressPackageStartupMessages({library(data.table); library(readxl); library(ConsensusClusterPlus); library(limma); library(survival); library(clusterProfiler); library(msigdbr)})
set.seed(666); root <- "${PROJECT_ROOT}"; args <- commandArgs(trailingOnly = TRUE)
sets <- list(
  NEDD8 = c("NEDD8","NAE1","UBA3","UBE2M","UBE2F","RBX1","RNF7","DCUN1D1","DCUN1D2","DCUN1D3","DCUN1D4","DCUN1D5","COPS5","COPS6","SENP8","UCHL3"),
  SUMO = c("SUMO1","SUMO2","SUMO3","SAE1","UBA2","UBE2I","CBX4","NSMCE2","PIAS1","PIAS2","PIAS3","PIAS4","RANBP2","ZBED1","SENP1","SENP2","SENP3","SENP5","SENP6","SENP7","DESI1","USPL1"),
  OST = c("DAD1","DDOST","KRTCAP2","MAGT1","OST4","OSTC","RPN1","RPN2","STT3A","STT3B","TMEM258","TUSC3"),
  GALNT = c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"),
  GPI = c("DPM2","GPAA1","PGAP1","PIGA","PIGB","PIGC","PIGF","PIGG","PIGH","PIGK","PIGL","PIGM","PIGN","PIGO","PIGP","PIGQ","PIGS","PIGT","PIGU","PIGV","PIGW","PIGX","PIGY","PIGZ"),
  UFM1 = c("UFM1","UBA5","UFC1","UFL1","DDRGK1","CDK5RAP3","UFSP1","UFSP2"))
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
ptype <- setNames(cdr$type, cdr$bcr_patient_barcode)
os <- data.frame(patient = cdr$bcr_patient_barcode, time = suppressWarnings(as.numeric(cdr$OS.time)) / 30.4375, event = suppressWarnings(as.numeric(cdr$OS)))
os$event[os$time > 120 & !is.na(os$time)] <- 0; os$time <- pmin(os$time, 120)
xena <- file.path(root, "01_raw/TCGA_EBpp_RNASeqV2.geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0)); smp <- hdr[-1]
si <- data.frame(sample = smp, patient = substr(smp, 1, 12), code = substr(smp, 14, 15)); si$type <- ptype[si$patient]; si <- si[!is.na(si$type), ]
prim <- si[si$code == "01" | (si$type == "LAML" & si$code == "03"), ]; prim <- prim[!duplicated(prim$patient), ]
ok_ct <- names(table(prim$type))[table(prim$type) >= 30]; prim <- prim[prim$type %in% ok_ct, ]
E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(prim$sample, hdr)), data.table = FALSE)
g0 <- E[[1]]; E <- E[!duplicated(g0) & !is.na(g0), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; al <- c(PPPDE2 = "DESI1", C11orf10 = "TMEM258", GALNTL2 = "GALNT15", GALNTL1 = "GALNT16", WBSCR17 = "GALNT17", GALNTL4 = "GALNT18", KIAA0776 = "UFL1"); rownames(E)[rownames(E) %in% names(al)] <- al[rownames(E)[rownames(E) %in% names(al)]]; E <- as.matrix(E); cat("matrix", dim(E), "\n")
gk <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"); t2g <- unique(gk[, c("gs_name", "gene_symbol")])
for (sys in args) {
  genes <- sets[[sys]]; res <- file.path(root, "04_results/2026-10-12_Figure3_multisystem_run01", sys); dir.create(res, recursive = TRUE, showWarnings = FALSE)
  cat("== ", sys, " (", length(genes), " genes; absent from matrix: ", paste(setdiff(genes, rownames(E)), collapse = ","), ")\n", sep = "")
  summ <- list(); enz <- list(); ora <- list(); cl_all <- list()
  for (ct in sort(ok_ct)) {
    p <- prim[prim$type == ct, ]; X <- E[, p$sample]; X <- X[rowSums(!is.finite(X)) == 0, ]; X <- X[apply(X, 1, var) > 0, ]
    gu <- intersect(genes, rownames(X)); gu <- gu[apply(X[gu, , drop = FALSE], 1, median) >= 1]
    if (length(gu) < 6) { summ[[ct]] <- data.table(type = ct, n = ncol(X), n_genes = length(gu), n_C1 = NA_integer_, n_C2 = NA_integer_, logrank_p = NA_real_, HR_C2_vs_C1 = NA_real_, n_deg = NA_integer_, minority_frac = NA_real_, genes_used = paste(gu, collapse = ";")); next }
    Z <- t(scale(t(X[gu, ]))); sink("/dev/null"); cc <- ConsensusClusterPlus(Z, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "km", distance = "euclidean", seed = 666, plot = NULL); sink()
    cls <- cc[[2]]$consensusClass; names(cls) <- p$sample
    d <- data.frame(sample = p$sample, patient = p$patient, cl = cls); d <- merge(d, os, by = "patient"); d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$event), ]
    sf <- survdiff(Surv(time, event) ~ cl, d); lp <- 1 - pchisq(sf$chisq, 1); hr2 <- tryCatch(exp(coef(coxph(Surv(time, event) ~ I(cl == 2), d))), warning = function(w) NA_real_)
    good <- if (is.na(hr2)) 1 else if (hr2 > 1) 1 else 2; d$C <- ifelse(d$cl == good, "C1", "C2"); hr <- if (is.na(hr2)) NA_real_ else if (hr2 > 1) hr2 else 1 / hr2
    cl_all[[ct]] <- data.table(type = ct, d[, c("sample", "patient", "time", "event", "C")])
    grp <- factor(cls[colnames(X)], levels = c(good, 3 - good), labels = c("C1", "C2")); fit <- eBayes(lmFit(X, model.matrix(~ grp))); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none"); tt$gene <- rownames(tt)
    e <- tt[gu, ]; enz[[ct]] <- data.table(type = ct, gene = gu, logFC = e$logFC, P = e$P.Value, BH = e$adj.P.Val)
    deg <- tt$gene[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5]
    eo <- tryCatch(as.data.frame(enricher(deg, TERM2GENE = t2g, universe = tt$gene, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)), error = function(e) NULL)
    if (!is.null(eo) && nrow(eo) > 0) ora[[ct]] <- data.table(type = ct, pathway = eo$ID, count = eo$Count, GeneRatio = eo$GeneRatio, p = eo$pvalue, FDR = eo$p.adjust)
    mf <- min(table(d$C)) / nrow(d)
    summ[[ct]] <- data.table(type = ct, n = nrow(d), n_genes = length(gu), n_C1 = sum(d$C == "C1"), n_C2 = sum(d$C == "C2"), logrank_p = lp, HR_C2_vs_C1 = hr, n_deg = length(deg), minority_frac = mf, genes_used = paste(gu, collapse = ";"))
  }
  summ <- rbindlist(summ, fill = TRUE); summ$system <- sys; summ$balanced <- summ$minority_frac >= 0.10; summ$BH_within <- p.adjust(summ$logrank_p, "BH")
  fwrite(summ, file.path(res, "cluster_survival_summary.csv")); fwrite(rbindlist(enz), file.path(res, "enzyme_limma_C2_vs_C1.csv")); fwrite(rbindlist(ora), file.path(res, "kegg_ora_all.csv")); fwrite(rbindlist(cl_all), file.path(res, "cluster_assignments.csv.gz"))
  print(summ[order(logrank_p)][1:8, .(type, n, n_genes, n_C1, n_C2, p = signif(logrank_p, 3), BHw = signif(BH_within, 3), HR = round(HR_C2_vs_C1, 2), balanced)])
}
