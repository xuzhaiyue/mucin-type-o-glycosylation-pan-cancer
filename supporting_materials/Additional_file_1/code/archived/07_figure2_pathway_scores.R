# Figure 2 data: ssGSEA pathway scores + enzyme-pathway Spearman (raw and proliferation-partial) per cancer. Spec: 00_meta/figure2_spec_v1.md
suppressPackageStartupMessages({library(data.table); library(GSVA); library(msigdbr); library(readxl)})
set.seed(666); root <- "${PROJECT_ROOT}"
res <- file.path(root, "04_results/2026-10-12_Figure2_pathway_run01"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
inter <- file.path(root, "02_intermediate/figure2"); dir.create(inter, recursive = TRUE, showWarnings = FALSE)
genes <- c("NEDD8","NAE1","UBA3","UBE2M","UBE2F","RBX1","RNF7","DCUN1D1","DCUN1D2","DCUN1D3","DCUN1D4","DCUN1D5","COPS5","COPS6","SENP8","UCHL3")
prolif <- c("MKI67","TOP2A","PCNA","MCM2","MCM6","CCNB1","CDK1","AURKA")
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
ptype <- setNames(cdr$type, cdr$bcr_patient_barcode)
xena <- file.path(root, "01_raw/TCGA_EBpp_RNASeqV2.geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", xena), nrows = 0)); smp <- hdr[-1]
si <- data.frame(sample = smp, patient = substr(smp, 1, 12), code = substr(smp, 14, 15)); si$type <- ptype[si$patient]; si <- si[!is.na(si$type), ]
prim <- si[si$code == "01" | (si$type == "LAML" & si$code == "03"), ]; prim <- prim[!duplicated(prim$patient), ]
ok_ct <- names(table(prim$type))[table(prim$type) >= 30]; prim <- prim[prim$type %in% ok_ct, ]; cat("cancers:", length(ok_ct), "tumors:", nrow(prim), "\n")
E <- fread(cmd = paste("gzip -dc", xena), select = c(1, match(prim$sample, hdr)), data.table = FALSE)
g0 <- E[[1]]; E <- E[!duplicated(g0) & !is.na(g0), ]; rownames(E) <- E[[1]]; E[[1]] <- NULL; E <- as.matrix(E); cat("matrix", dim(E), "\n")
gs_h <- msigdbr(species = "Homo sapiens", collection = "H"); gs_k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
gsl <- c(split(gs_h$gene_symbol, gs_h$gs_name), split(gs_k$gene_symbol, gs_k$gs_name)); gsl <- lapply(gsl, unique)
src <- setNames(c(rep("Hallmark", length(unique(gs_h$gs_name))), rep("KEGG", length(unique(gs_k$gs_name)))), c(unique(gs_h$gs_name), unique(gs_k$gs_name)))
member <- sapply(gsl, function(s) genes %in% s); rownames(member) <- genes
cor_t <- function(r, n) { t <- r * sqrt((n - 2) / (1 - r^2)); 2 * pt(-abs(t), n - 2) }
out <- list(); prl <- list()
for (ct in sort(ok_ct)) {
  t0 <- Sys.time(); s <- prim$sample[prim$type == ct]; X <- E[, s]; X <- X[rowSums(!is.finite(X)) == 0, ]; X <- X[apply(X, 1, var) > 0, ]
  sc <- gsva(ssgseaParam(X, gsl, minSize = 10, maxSize = 500, normalize = TRUE), verbose = FALSE)
  saveRDS(list(scores = sc, enz = X[genes, ], prolif = X[prolif, ]), file.path(inter, paste0("scores_", ct, ".rds")))
  n <- ncol(X); pz <- colMeans(t(scale(t(X[prolif, ])))); 
  R <- cor(t(X[genes, ]), t(sc), method = "spearman"); rz_p <- cor(t(X[genes, ]), pz, method = "spearman")[, 1]; ry_p <- cor(t(sc), pz, method = "spearman")[, 1]
  Rp <- (R - outer(rz_p, ry_p)) / sqrt(outer(1 - rz_p^2, 1 - ry_p^2))
  P <- cor_t(R, n); Pp <- 2 * pt(-abs(Rp * sqrt((n - 3) / (1 - Rp^2))), n - 3)
  m <- member[genes, colnames(R)]; R[m] <- NA; Rp[m] <- NA; P[m] <- NA; Pp[m] <- NA
  d <- data.table(type = ct, n = n, enzyme = rep(genes, ncol(R)), pathway = rep(colnames(R), each = length(genes)), rho = as.vector(R), p = as.vector(P), rho_partial = as.vector(Rp), p_partial = as.vector(Pp))
  d$BH <- p.adjust(d$p, "BH"); d$BH_partial <- p.adjust(d$p_partial, "BH"); out[[ct]] <- d
  prl[[ct]] <- data.table(type = ct, enzyme = genes, rho_prolif = rz_p)
  cat(ct, n, "done", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
out <- rbindlist(out); out$source <- src[out$pathway]; fwrite(out, file.path(res, "enzyme_pathway_correlations.csv.gz")); fwrite(rbindlist(prl), file.path(res, "enzyme_proliferation_rho.csv"))
cat("rows", nrow(out), " BH<0.05:", mean(out$BH < 0.05, na.rm = TRUE), " partial BH<0.05:", mean(out$BH_partial < 0.05, na.rm = TRUE), "\n")
