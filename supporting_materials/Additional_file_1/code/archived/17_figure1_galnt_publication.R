# Publication-style Figure 1 for GALNT (freeze v5). Layout follows the template Fig. 1; tables from 16_figure1_galnt_data.R.
# Reads tables from run01 where possible; recomputes E and F (F now includes CPTAC-OV).
suppressPackageStartupMessages({library(data.table); library(limma); library(readxl); library(ggplot2); library(patchwork); library(cowplot); library(scales)})
set.seed(666)
root <- "${PROJECT_ROOT}"; raw <- file.path(root, "01_raw")
r1 <- file.path(root, "04_results/2026-10-12_Figure1_GALNT_run01")
out <- file.path(root, "05_figures/2026-10-12_Figure1_GALNT_run01"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
res2 <- file.path(root, "04_results/2026-10-12_Figure1_GALNT_run01"); dir.create(res2, recursive = TRUE, showWarnings = FALSE)
fam <- "Arial"
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1")
cls <- setNames(c(rep("init", 20), rep("ext", 2)), genes)
col_w <- "#D9942B"; col_e <- "#7B52A6"; lab_col <- ifelse(cls[genes] == "init", col_w, col_e); names(lab_col) <- genes
theme_pub <- function(bs = 6) theme_classic(base_size = bs, base_family = fam) +
  theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25), plot.title = element_text(size = bs + 0.5, face = "bold"),
        legend.key.size = unit(2.6, "mm"), legend.title = element_text(size = bs), legend.text = element_text(size = bs - 0.5), plot.margin = margin(2, 2, 2, 2))
mm <- function(pt) pt / .pt   # font size in pt -> geom_text size

## ---------------------------------------------------------------- A: overview with localization
hpa <- fread(file.path(raw, "uniprot/hpa_subcellular_galnt22.tsv"), sep = "\t")
loc_map <- function(txt) { x <- tolower(txt); r <- c()
  if (grepl("golgi", x)) r <- c(r, "Golgi apparatus"); if (grepl("endoplasmic", x)) r <- c(r, "ER")
  if (grepl("cytosol|vesicles", x)) r <- c(r, "Cytoplasm"); if (grepl("nucleoplasm|nuclear bodies|nucleoli", x)) r <- c(r, "Nucleus"); if (grepl("plasma membrane", x)) r <- c(r, "Plasma membrane"); r }
locs <- rbindlist(lapply(genes, function(g) { h <- hpa[Gene == g]; l <- loc_map(paste(h$main, h$additional)); if (!length(l)) return(NULL); data.frame(gene = g, loc = l) }))
fwrite(locs, file.path(res2, "A_subcellular_localization_HPA.csv"))
loc_cols <- c("Golgi apparatus" = "#F4A460", ER = "#A6D854", Cytoplasm = "#E5535A", Nucleus = "#6BAED6", "Plasma membrane" = "#BDBDBD")
rr <- function(x0, x1, y0, y1, r = 1.6, id = 1, ax = 0.39) { th <- seq(0, pi / 2, length.out = 8); rx <- r * ax; ry <- r
  data.frame(x = c(x1 - rx + rx * cos(th), x0 + rx - rx * sin(th), x0 + rx - rx * cos(th), x1 - rx + rx * sin(th)),
             y = c(y1 - ry + ry * sin(th), y1 - ry + ry * cos(th), y0 + ry - ry * sin(th), y0 + ry - ry * cos(th)), id = id) }
boxes <- list(
  list(t = "Initiation (GalNAc-Ts)", g = paste0("GALNT", 1:7), x0 = 2, x1 = 21, row = 1), list(t = "Initiation (cont.)", g = paste0("GALNT", 8:14), x0 = 23, x1 = 42, row = 1),
  list(t = "Initiation (cont.)", g = c(paste0("GALNT", 15:18), "GALNTL5", "GALNTL6"), x0 = 44, x1 = 63, row = 1),
  list(t = "Extension (core 1)", g = c("C1GALT1", "C1GALT1C1"), x0 = 2, x1 = 21, row = 2))
y_top <- c(`1` = 53, `2` = 21.5); line_h <- 2.9
polys <- list(); txts <- list(); dots <- list(); k <- 0
for (b in boxes) { k <- k + 1; yt <- y_top[[as.character(b$row)]]; yb <- yt - 3.6 - length(b$g) * line_h - 1
  polys[[k]] <- rr(b$x0, b$x1, yb, yt, id = k); txts[[length(txts) + 1]] <- data.frame(x = (b$x0 + b$x1) / 2, y = yt - 1.8, label = b$t, face = "bold", hj = 0.5)
  for (i in seq_along(b$g)) { yy <- yt - 3.6 - (i - 0.5) * line_h; txts[[length(txts) + 1]] <- data.frame(x = b$x0 + 1.4, y = yy, label = b$g[i], face = "plain", hj = 0)
    l <- locs$loc[locs$gene == b$g[i]]; for (j in seq_along(l)) dots[[length(dots) + 1]] <- data.frame(x = b$x1 - 1.2 - (length(l) - j) * 2.0, y = yy, loc = l[j]) } }
polys <- rbindlist(polys); txts <- rbindlist(txts); dots <- rbindlist(dots)
bg <- rr(0, 100, 0, 60, r = 3.5); leg <- data.frame(loc = names(loc_cols), y = 51 - 3.4 * (0:4))
pA <- ggplot() + geom_polygon(data = bg, aes(x, y), fill = "#E6F1F7", colour = NA) +
  geom_polygon(data = polys, aes(x, y, group = id), fill = "#ECECF3", colour = "grey25", linewidth = 0.3) +
  geom_text(data = txts, aes(x, y, label = label, fontface = face, hjust = hj), size = mm(5.2), family = fam) +
  geom_point(data = dots, aes(x, y, colour = loc), size = 1.3) + scale_colour_manual(values = loc_cols) +
  annotate("text", x = 2, y = 57.3, label = "Mucin-type O-glycosylation initiation (GalNAc-transferases)", fontface = "bold", size = mm(6.5), hjust = 0, family = fam) +
  annotate("text", x = 2, y = 24.6, label = "Core 1 extension", fontface = "bold", size = mm(6.5), hjust = 0, family = fam) +
  annotate("text", x = 26, y = 14, label = "No dedicated eraser enzyme is known\nfor secretory-pathway O-GalNAc glycans", size = mm(5.4), hjust = 0, vjust = 0.5, lineheight = 1.05, colour = "grey25", family = fam) +
  annotate("text", x = 70, y = 54.5, label = "Subcellular localization", fontface = "plain", size = mm(5.5), hjust = 0, family = fam) +
  geom_point(data = leg, aes(x = 71, y = y - 2.8, colour = loc), size = 1.8) + geom_text(data = leg, aes(x = 73.5, y = y - 2.8, label = loc), hjust = 0, size = mm(5.2), family = fam) +
  annotate("text", x = 70, y = 29, label = "Initiation genes", colour = col_w, fontface = "bold", size = mm(5.5), hjust = 0, family = fam) +
  annotate("text", x = 70, y = 25, label = "Extension genes", colour = col_e, fontface = "bold", size = mm(5.5), hjust = 0, family = fam) +
  annotate("text", x = 70, y = 4, label = "Localization: Human Protein Atlas\n(no dot = no HPA annotation)", size = mm(4.5), hjust = 0, vjust = 0.5, lineheight = 1, colour = "grey30", family = fam) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 60), expand = FALSE) + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

## ---------------------------------------------------------------- B: mutation classes
cdr <- as.data.frame(read_excel("${TCGA_CDR_ROOT}/01_raw/pancan_endpoints/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1, .name_repair = "minimal"))
ptype <- setNames(cdr$type, cdr$bcr_patient_barcode)
hm <- names(fread(cmd = paste("gzip -dc", file.path(raw, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz")), nrows = 0))[-1]
ms <- data.frame(sample = hm, patient = substr(hm, 1, 12)); ms$type <- ptype[ms$patient]; ms <- ms[!is.na(ms$type) & !duplicated(ms$patient), ]
maf <- fread(cmd = paste("gzip -dc", file.path(raw, "mc3.v0.2.8.PUBLIC.xena.gz")), select = c("sample", "gene", "effect"))
maf$gene[maf$gene == "WBSCR17"] <- "GALNT17"; maf <- maf[gene %in% genes & sample %in% ms$sample]
ns <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","Translation_Start_Site")
maf <- maf[effect %in% ns]; N <- nrow(ms)
cnt <- maf[, .N, by = .(gene, effect)]; fwrite(cnt, file.path(res2, "B_mutation_counts_by_class.csv"))
tot <- maf[, .(n_mut = .N, n_samples = uniqueN(sample)), by = gene][order(-n_mut)]; tot$pct <- 100 * tot$n_samples / N; fwrite(tot, file.path(res2, "B_mutation_totals.csv"))
top <- head(tot, 10); cnt <- cnt[gene %in% top$gene]; cnt$gene <- factor(cnt$gene, levels = rev(top$gene)); cnt$effect <- factor(cnt$effect, levels = ns)
mut_cols <- c(Missense_Mutation = "#A48BC2", Nonsense_Mutation = "#E3342F", Frame_Shift_Del = "#3A7DC1", Splice_Site = "#F28E2B", Frame_Shift_Ins = "#5B2A86",
              In_Frame_Del = "#F3E88A", Nonstop_Mutation = "#8FD0C4", In_Frame_Ins = "#E86A8F", Translation_Start_Site = "#D95F02")
top$gene <- factor(top$gene, levels = rev(top$gene))
pB <- ggplot(cnt, aes(N, gene, fill = effect)) + geom_col(width = 0.8) +
  geom_text(data = top, aes(x = n_mut, y = gene, label = sprintf("%.1f%%", pct)), inherit.aes = FALSE, hjust = -0.15, size = mm(5.2), family = fam) +
  scale_fill_manual(values = mut_cols, labels = function(x) gsub("_", " ", x), name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.62))) +
  labs(x = "No. of mutations", y = NULL) + guides(fill = guide_legend(ncol = 1)) + theme_pub() +
  theme(axis.text.y = element_text(colour = lab_col[levels(top$gene)], face = "plain", size = 6), legend.position = c(0.80, 0.36), legend.key.size = unit(1.6, "mm"), legend.key.spacing.y = unit(0.1, "mm"), legend.text = element_text(size = 4.2), legend.margin = margin(0, 0, 0, 0), legend.background = element_blank())

## ---------------------------------------------------------------- C: frequency heatmap
fc <- fread(file.path(r1, "C_mutation_frequency_by_cancer.csv")); ld <- fread(file.path(r1, "C_mutation_load_proxy_by_cancer.csv"))
top_ct <- sort(head(ld[order(-load)]$type, 18)); fc <- fc[type %in% top_ct]; fc$type <- factor(fc$type, levels = top_ct); fc$gene <- factor(fc$gene, levels = genes)
pC <- ggplot(fc, aes(type, gene, fill = pct)) + geom_tile(colour = "grey70", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", pct)), size = mm(4.4), family = fam) +
  scale_fill_gradient(low = "white", high = "#7E5AA6", limits = c(0, max(fc$pct)), name = "Frequency%") +
  scale_y_discrete(limits = rev(genes), position = "right") + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(colour = rev(lab_col[genes]), size = 5.2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.8), plot.margin = margin(13, 2, 2, 7), legend.position = "right", legend.key.height = unit(5, "mm"), legend.key.width = unit(2.5, "mm"), legend.title = element_text(size = 5))

## ---------------------------------------------------------------- D: expression + prognosis
tn <- fread(file.path(r1, "D_tumor_vs_normal_combined.csv")); sv <- fread(file.path(r1, "D_os_logrank_median_split.csv"))
all_ct <- sort(union(tn$type, sv$type)); cat("D cancers:", length(all_ct), "\n")
ratio <- tn[, .(up = 100 * mean(sig == "up"), down = 100 * mean(sig == "down")), by = gene]; ratio$gene <- factor(ratio$gene, levels = rev(genes))
ratio$star <- ifelse(ratio$up > 25, "up", ifelse(ratio$down > 25, "down", NA))
rb <- melt(ratio, id.vars = c("gene", "star"), measure.vars = c("up", "down"), variable.name = "dir", value.name = "pct"); rb$dir <- factor(rb$dir, levels = c("down", "up"))
dir_cols <- c(up = "#EC7E67", down = "#3FB08E")
pDb <- ggplot(rb, aes(pct, gene, fill = dir)) + geom_col(width = 0.8) + scale_fill_manual(values = dir_cols, labels = c(down = "DOWN", up = "UP"), name = "Regulation") +
  geom_text(data = ratio[!is.na(star)], aes(x = up + down + 6, y = gene, label = "*", colour = star), inherit.aes = FALSE, size = 4, family = fam, show.legend = FALSE) +
  scale_colour_manual(values = dir_cols) + scale_x_reverse(expand = expansion(mult = c(0.05, 0)), limits = c(110, 0)) + labs(x = "% cancers", y = NULL) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.title.x = element_text(size = 5.5))
td <- CJ(type = all_ct, gene = genes); td <- merge(td, tn[, .(type, gene, log2FC)], by = c("type", "gene"), all.x = TRUE)
sp <- sv[p_logrank < 0.05]; sp$cat <- ifelse(sp$HR_high_vs_low < 1, "Favorable", "Unfavorable"); sp$lp <- pmin(-log10(sp$p_logrank), 6)
td$type <- factor(td$type, levels = all_ct); td$gene <- factor(td$gene, levels = rev(genes)); sp$type <- factor(sp$type, levels = all_ct); sp$gene <- factor(sp$gene, levels = rev(genes))
pDm <- ggplot(td, aes(type, gene)) + geom_tile(aes(fill = pmax(pmin(log2FC, 2.5), -2.5)), colour = "grey85", linewidth = 0.2) +
  scale_fill_gradient2(low = "#3FB08E", mid = "white", high = "#EC7E67", midpoint = 0, limits = c(-2.5, 2.5), na.value = "white", name = "logFC") +
  geom_point(data = sp, aes(size = lp, colour = cat), alpha = 0.85) + scale_colour_manual(values = c(Favorable = "#C45BC8", Unfavorable = "#7C8DC9"), name = "Prognostic impact") +
  scale_size_continuous(range = c(0.6, 3.4), limits = c(1.3, 6), breaks = -log10(c(0.05, 0.01, 0.001, 1e-6)), labels = c("0.05", "0.01", "0.001", "1e-06"), name = "Log-rank p", trans = "identity") +
  scale_y_discrete(position = "right") + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(colour = rev(lab_col[genes]), size = 6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.2),
        panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3))
pDb <- pDb + theme(legend.position = "top")
pD <- wrap_elements(full = (pDb | pDm) + plot_layout(widths = c(0.2, 0.8), guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal", legend.box = "horizontal", legend.box.just = "top", legend.justification = "right", legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.4), legend.title = element_text(size = 4.8, face = "plain"), legend.title.position = "top",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 14), legend.box.spacing = unit(1, "mm"), legend.spacing.x = unit(1, "mm")))

## ---------------------------------------------------------------- E, F: CPTAC
rd <- function(f) { d <- fread(f, data.table = FALSE, quote = "\""); g <- as.character(d[[1]]); d[[1]] <- NULL; d <- as.matrix(d); k <- !duplicated(g) & !is.na(g) & g != ""; d <- d[k, ]; rownames(d) <- g[k]; storage.mode(d) <- "double"; d }
pr <- rd(file.path(raw, "linkedomics/BRCA_prot.cct")); rn <- rd(file.path(raw, "linkedomics/BRCA_rna.cct")); cs <- intersect(colnames(pr), colnames(rn)); cg <- intersect(rownames(pr), rownames(rn))
sp_f <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 20) return(NA_real_); suppressWarnings(cor(a[ok], b[ok], method = "spearman")) }
enz <- data.frame(gene = genes, rho = sapply(genes, function(g) if (g %in% rownames(rn) && g %in% rownames(pr)) sp_f(rn[g, cs], pr[g, cs]) else NA_real_)); fwrite(enz, file.path(res2, "E_enzyme_mrna_protein_rho.csv"))
pairs <- data.frame(a = sample(cg, 6000, TRUE), b = sample(cg, 6000, TRUE)); pairs <- pairs[pairs$a != pairs$b, ]
rand <- mapply(function(a, b) sp_f(rn[a, cs], pr[b, cs]), pairs$a, pairs$b); rand <- rand[is.finite(rand)]
wp <- wilcox.test(enz$rho, rand)$p.value; fwrite(data.frame(stat = c("n_random", "enzyme_median", "random_median", "wilcoxon_p"), value = c(length(rand), median(enz$rho, na.rm = TRUE), median(rand), wp)), file.path(res2, "E_summary.csv"))
de <- rbind(data.frame(R = rand, grp = "Random"), data.frame(R = enz$rho[is.finite(enz$rho)], grp = "GALNT-family enzymes"))
cat("E: enzymes with CPTAC mRNA-protein pair:", sum(is.finite(enz$rho)), "of", length(genes), "\n")
pE <- ggplot(de, aes(R, fill = grp, colour = grp, linetype = grp)) + geom_density(alpha = 0.55, linewidth = 0.35, adjust = 1) +
  scale_fill_manual(values = c(Random = "#F4A38F", "GALNT-family enzymes" = "#8FD3BC"), name = "From") + scale_colour_manual(values = c(Random = "grey35", "GALNT-family enzymes" = "grey35"), name = "From") +
  scale_linetype_manual(values = c(Random = "solid", "GALNT-family enzymes" = "dashed"), name = "From") +
  annotate("text", x = -0.6, y = 2.65, label = sprintf("Wilcoxon\np = %.1e", wp), size = mm(5), hjust = 0, vjust = 1, lineheight = 0.9, family = fam) +
  labs(x = "R", y = "Density") + coord_cartesian(xlim = c(-0.6, 0.8)) + theme_pub() + theme(legend.position = c(0.8, 0.8), legend.background = element_blank(), legend.key.size = unit(2.4, "mm"))

cohorts <- list(KIRC = "CCRCC", COAD = "COAD", HNSC = "HNSCC", LUSC = "LSCC", LUAD = "LUAD", PAAD = "PDAC", UCEC = "UCEC", OV = "OV")
pf <- list()
for (tc in names(cohorts)) { nm <- cohorts[[tc]]
  Tm <- rd(file.path(raw, sprintf("linkedomics/%s_prot_T.cct", nm))); Nm <- rd(file.path(raw, sprintf("linkedomics/%s_prot_N.cct", nm)))
  cg2 <- intersect(rownames(Tm), rownames(Nm)); X <- cbind(Nm[cg2, ], Tm[cg2, ]); g <- factor(c(rep("N", ncol(Nm)), rep("T", ncol(Tm))), levels = c("N", "T"))
  kg <- rowSums(is.finite(X[, g == "T"])) >= 10 & rowSums(is.finite(X[, g == "N"])) >= 5; X <- X[kg, ]
  tt <- topTable(eBayes(lmFit(X, model.matrix(~ g))), coef = 2, number = Inf, sort.by = "none"); tt$BH <- p.adjust(tt$P.Value, "BH")
  i <- match(genes, rownames(tt)); pf[[tc]] <- data.frame(type = tc, cohort = nm, gene = genes, log2FC = tt$logFC[i], p = tt$P.Value[i], BH = tt$BH[i], nT = ncol(Tm), nN = ncol(Nm)) }
pf <- rbindlist(pf); mr <- tn[, .(type, gene, mrna_sig = sig)]; pf <- merge(pf, mr, by = c("type", "gene"), all.x = TRUE)
pf$dir <- ifelse(pf$log2FC > 0, "up", "down"); pf$prot_sig <- !is.na(pf$BH) & pf$BH < 0.05
pf$aligned <- pf$prot_sig & !is.na(pf$mrna_sig) & pf$mrna_sig == pf$dir; fwrite(pf, file.path(res2, "F_cptac_protein_tumor_vs_normal.csv"))
pf$lp <- pmin(-log10(pf$p), 6); pf$gene <- factor(pf$gene, levels = rev(genes)); pf$type <- factor(pf$type, levels = names(cohorts))
pF <- ggplot(pf[!is.na(pf$p), ], aes(type, gene)) + geom_point(aes(size = lp, colour = dir), alpha = 0.9) +
  geom_text(data = pf[pf$aligned %in% TRUE, ], aes(label = "*"), size = mm(11), vjust = 0.35, fontface = "bold", family = fam) +
  scale_colour_manual(values = dir_cols, labels = c(down = "DOWN", up = "UP"), name = "Regulation") +
  scale_size_continuous(range = c(0.2, 2.2), limits = c(0, 6), breaks = -log10(c(0.05, 0.01, 0.001, 1e-6)), labels = c("0.05", "0.01", "0.001", "1e-06"), name = "p value") +
  scale_y_discrete(limits = rev(genes), position = "right") + scale_x_discrete(drop = FALSE) + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(colour = rev(lab_col[genes]), size = 6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.2), panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
        legend.position = "right", legend.box.spacing = unit(0, "mm"), legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.8), legend.title = element_text(size = 5.2)) +
  guides(colour = guide_legend(override.aes = list(size = 2.2)))

## ---------------------------------------------------------------- compose
design <- "AAAAAAAAAAAAAAAAAAAA
BBBBBBBBCCCCCCCCCCCC
DDDDDDDDDDDDDEEEEEEE
DDDDDDDDDDDDDFFFFFFF"
comp <- pA + pB + pC + pD + pE + pF + plot_layout(design = design, heights = c(1.45, 2.4, 1.9, 3.9)) &
  theme(text = element_text(family = fam), plot.tag = element_text(face = "bold", size = 10, family = fam), plot.tag.position = c(0, 1))
comp <- comp + plot_annotation(tag_levels = "A")
W <- 7.2; H <- 9.5
ggsave(file.path(out, "Figure1.svg"), comp, width = W, height = H, device = svglite::svglite)
ggsave(file.path(out, "Figure1.png"), comp, width = W, height = H, dpi = 400, device = ragg::agg_png)
ggsave(file.path(out, "Figure1.tiff"), comp, width = W, height = H, dpi = 600, device = ragg::agg_tiff, compression = "lzw")
for (nm in c("A", "B", "C", "D", "E", "F")) { p <- get(paste0("p", nm)); ggsave(file.path(out, paste0("panel_", nm, ".png")), p, width = 4, height = 3, dpi = 300, device = ragg::agg_png) }
cat("done\n")
