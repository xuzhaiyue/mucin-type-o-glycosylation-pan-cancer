# Figure 2 data for GALNT (freeze v5): enzyme-pathway Spearman (raw and proliferation-partial) per cancer, reusing ssGSEA scores saved by 07_figure2_pathway_scores.R (pathway scores do not depend on the gene set). Spec: 00_meta/figure2_spec_v1.md
suppressPackageStartupMessages({library(data.table); library(GSVA); library(msigdbr); library(readxl)})
set.seed(666); root <- "${PROJECT_ROOT}"
res <- file.path(root, "04_results/2026-10-12_Figure2_GALNT_run01"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
inter <- file.path(root, "02_intermediate/figure2"); dir.create(inter, recursive = TRUE, showWarnings = FALSE)
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1"); al <- c(GALNTL2 = "GALNT15", GALNTL1 = "GALNT16", WBSCR17 = "GALNT17", GALNTL4 = "GALNT18")
prolif <- c("MKI67","TOP2A","PCNA","MCM2","MCM6","CCNB1","CDK1","AURKA")
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
ptype <- setNames(cdr$type, cdr$bcr_patient_barcode)
xena <- file.path(root, "01_raw/TCGA_EBpp_RNASeqV2.geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0)); smp <- hdr[-1]
si <- data.frame(sample = smp, patient = substr(smp, 1, 12), code = substr(smp, 14, 15)); si$type <- ptype[si$patient]; si <- si[!is.na(si$type), ]
prim <- si[si$code == "01" | (si$type == "LAML" & si$code == "03"), ]; prim <- prim[!duplicated(prim$patient), ]
ok_ct <- names(table(prim$type))[table(prim$type) >= 30]; prim <- prim[prim$type %in% ok_ct, ]; cat("cancers:", length(ok_ct), "tumors:", nrow(prim), "\n")
E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(prim$sample, hdr)), data.table = FALSE)
g0 <- E[[1]]; E <- E[!duplicated(g0) & !is.na(g0), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; E <- as.matrix(E); rownames(E)[rownames(E) %in% names(al)] <- al[rownames(E)[rownames(E) %in% names(al)]]; cat("matrix", dim(E), "\n")
gs_h <- msigdbr(species = "Homo sapiens", collection = "H"); gs_k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
gsl <- c(split(gs_h$gene_symbol, gs_h$gs_name), split(gs_k$gene_symbol, gs_k$gs_name)); gsl <- lapply(gsl, unique)
src <- setNames(c(rep("Hallmark", length(unique(gs_h$gs_name))), rep("KEGG", length(unique(gs_k$gs_name)))), c(unique(gs_h$gs_name), unique(gs_k$gs_name)))
member <- sapply(gsl, function(s) genes %in% s); rownames(member) <- genes
cor_t <- function(r, n) { t <- r * sqrt((n - 2) / (1 - r^2)); 2 * pt(-abs(t), n - 2) }
out <- list(); prl <- list()
for (ct in sort(ok_ct)) {
  t0 <- Sys.time(); s0 <- readRDS(file.path(inter, paste0("scores_", ct, ".rds"))); sc <- s0$scores; X <- E[, colnames(sc)]; X <- rbind(X, s0$prolif[, colnames(sc)]); X <- X[rowSums(!is.finite(X)) == 0, ]
  gu <- genes[genes %in% rownames(X)]; gu <- gu[apply(X[gu, , drop = FALSE], 1, median) >= 1]; if (length(gu) < 6) next
  n <- ncol(X); pz <- colMeans(t(scale(t(X[prolif, ])))); 
  R <- cor(t(X[gu, ]), t(sc), method = "spearman"); rz_p <- cor(t(X[gu, ]), pz, method = "spearman")[, 1]; ry_p <- cor(t(sc), pz, method = "spearman")[, 1]
  Rp <- (R - outer(rz_p, ry_p)) / sqrt(outer(1 - rz_p^2, 1 - ry_p^2))
  P <- cor_t(R, n); Pp <- 2 * pt(-abs(Rp * sqrt((n - 3) / (1 - Rp^2))), n - 3)
  m <- member[gu, colnames(R)]; R[m] <- NA; Rp[m] <- NA; P[m] <- NA; Pp[m] <- NA
  d <- data.table(type = ct, n = n, enzyme = rep(gu, ncol(R)), pathway = rep(colnames(R), each = length(gu)), rho = as.vector(R), p = as.vector(P), rho_partial = as.vector(Rp), p_partial = as.vector(Pp))
  d$BH <- p.adjust(d$p, "BH"); d$BH_partial <- p.adjust(d$p_partial, "BH"); out[[ct]] <- d
  prl[[ct]] <- data.table(type = ct, enzyme = gu, rho_prolif = rz_p)
  cat(ct, n, "done", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
out <- rbindlist(out); out$source <- src[out$pathway]; fwrite(out, file.path(res, "enzyme_pathway_correlations.csv.gz")); fwrite(rbindlist(prl), file.path(res, "enzyme_proliferation_rho.csv"))
cat("rows", nrow(out), " BH<0.05:", mean(out$BH < 0.05, na.rm = TRUE), " partial BH<0.05:", mean(out$BH_partial < 0.05, na.rm = TRUE), "\n")
