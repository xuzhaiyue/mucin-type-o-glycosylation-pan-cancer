# Figure 3 (GALNT, freeze v5). Spec: 00_meta/freeze_v5_galnt.md. Uses tables from 13_figure3_multisystem.R and 14_galnt_confounders.R.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(survival); library(ggnewscale); library(ggVennDiagram)})
root <- "${PROJECT_ROOT}"; r3 <- file.path(root, "04_results/2026-10-12_Figure3_multisystem_run01/GALNT")
out <- file.path(root, "05_figures/2026-10-12_Figure3_GALNT_run01"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
fam <- "Arial"; mm <- function(pt) pt / .pt
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1")
grp <- setNames(c(rep("init", 20), rep("ext", 2)), genes); col_w <- "#D9942B"; col_e <- "#7B52A6"; lab_col <- ifelse(grp == "init", col_w, col_e)
cC1 <- "#E0A526"; cC2 <- "#4A7FB5"; FOCUS <- "PAAD"
theme_pub <- function(bs = 6) theme_classic(base_size = bs, base_family = fam) + theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
  plot.title = element_text(size = bs + 0.5, face = "bold"), legend.key.size = unit(2.6, "mm"), legend.title = element_text(size = bs), legend.text = element_text(size = bs - 0.5), plot.margin = margin(2, 2, 2, 2))
nice <- function(x) { y <- tools::toTitleCase(tolower(gsub("_", " ", sub("^KEGG_", "", x))))
  rep <- c("Atp" = "ATP", "Ecm" = "ECM", "Dna" = "DNA", "Hcm" = "HCM", "Arvc" = "ARVC", "Cams" = "CAMs", "Nod" = "NOD", "Rna" = "RNA", "Tgf" = "TGF", "Mapk" = "MAPK", "Pi3k" = "PI3K", "Akt" = "AKT", "Jak" = "JAK", "Stat" = "STAT")
  for (k in names(rep)) y <- gsub(paste0("\\b", k, "\\b"), rep[[k]], y); y <- gsub("Cytochrome P450", "CYP450", y); y <- gsub("^[Oo] Glycan", "O-Glycan", y); y }
sm <- fread(file.path(r3, "cluster_survival_summary.csv")); cl <- fread(file.path(r3, "cluster_assignments.csv.gz")); en <- fread(file.path(r3, "enzyme_limma_C2_vs_C1.csv")); ora <- fread(file.path(r3, "kegg_ora_all.csv"))
sig <- sm[BH_within < 0.05][order(logrank_p)]$type; show <- unique(c(sig, "COAD", "STAD")); show <- show[order(match(show, sm[order(logrank_p)]$type))]
conf <- fread(file.path(r3, "confounder_cox.csv"))
cat("shown:", show, "\n"); tum <- show

## A: Kaplan-Meier
km <- function(ct) {
  d <- cl[type == ct]; sf <- survfit(Surv(time, event) ~ C, d); s <- summary(sf, censored = TRUE); st <- sub("C=", "", as.character(s$strata))
  df <- rbind(data.frame(time = 0, surv = 1, C = c("C1", "C2"), nc = 0), data.frame(time = s$time, surv = s$surv, C = st, nc = s$n.censor))
  r <- sm[type == ct]; ptxt <- ifelse(r$logrank_p < 0.001, "p < 0.001", paste0("p = ", formatC(r$logrank_p, format = "f", digits = 3))); cj <- conf[type == ct]; atxt <- paste0("adj. p = ", ifelse(cj$p_M2 < 0.001, "<0.001", formatC(cj$p_M2, format = "f", digits = 3)))
  lab <- c(C1 = paste0("C1 (n=", r$n_C1, ")"), C2 = paste0("C2 (n=", r$n_C2, ")")); is_f <- r$BH_within >= 0.05
  ggplot(df, aes(time, surv, colour = C)) + geom_step(linewidth = 0.45) + geom_point(data = df[df$nc > 0 & df$time <= 120, ], shape = 3, size = 0.7, stroke = 0.25) +
    annotate("text", x = 2, y = 0.22, label = ct, hjust = 0, size = mm(7), fontface = "bold", family = fam) + annotate("text", x = 2, y = 0.12, label = ptxt, hjust = 0, size = mm(6), family = fam, colour = ifelse(is_f, "grey40", "black")) + annotate("text", x = 2, y = 0.03, label = atxt, hjust = 0, size = mm(5.4), family = fam, colour = "grey35") +
    scale_colour_manual(values = c(C1 = cC1, C2 = cC2), labels = lab, name = NULL) + scale_x_continuous(breaks = seq(0, 120, 24), limits = c(0, 120), expand = c(0.01, 0)) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
    labs(x = "Time in months", y = "Overall survival") + theme_pub(6) + theme(legend.position = c(0.97, 0.99), legend.justification = c(1, 1), legend.key.height = unit(2.4, "mm"), legend.key.width = unit(4, "mm"), legend.background = element_blank(), legend.text = element_text(size = 5.5),
    axis.title = element_text(size = 5.8), axis.text = element_text(size = 5.3))
}
kms <- lapply(tum, km)
nA <- length(kms)
pA <- wrap_plots(kms, ncol = 3) + plot_annotation(title = NULL)

## B: enzyme logFC (C2 vs C1)
eb <- en[type %in% tum]; eb$type <- factor(eb$type, levels = tum); eb <- merge(CJ(type = factor(tum, levels = tum), gene = genes), eb, by = c("type", "gene"), all.x = TRUE); eb$gene <- factor(eb$gene, levels = rev(genes)); eb$star <- ifelse(!is.na(eb$BH) & eb$BH < 0.05, "*", "")
dirn <- eb[, .(med = median(logFC, na.rm = TRUE)), by = type]; dirn$dir <- ifelse(dirn$med > 0, "Up", "Down"); lim <- ceiling(quantile(abs(eb$logFC), 0.98, na.rm = TRUE) * 10) / 10; print(lim)
dirn$dir <- factor(dirn$dir)
pB1 <- ggplot(dirn, aes(type, 1, fill = dir)) + geom_tile(colour = "white", linewidth = 0.4, height = 0.9) + scale_fill_manual(values = c(Up = "#E8706B", Down = "#7FCDB9"), name = "Direction") +
  scale_y_continuous(breaks = 1, labels = "Direction") + labs(x = NULL, y = NULL, title = "GALNT-family expression, C2 vs C1") + theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5.5), legend.position = "right", legend.title = element_text(size = 5.5), legend.text = element_text(size = 5.3))
pB2 <- ggplot(eb, aes(type, gene)) + geom_tile(aes(fill = logFC), colour = "white", linewidth = 0.3) + geom_text(aes(label = star), size = mm(7), vjust = 0.75, family = fam) +
  scale_fill_gradient2(low = "#3B7BB5", mid = "white", high = "#D6604D", limits = c(-lim, lim), oob = squish, name = "log2FC", na.value = "grey88") + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.4, colour = rev(lab_col[genes])), axis.text.x = element_text(size = 5.8, angle = 90, vjust = 0.5, hjust = 1, face = ifelse(tum == FOCUS, "bold", "plain")), legend.position = "right", legend.title = element_text(size = 5.5), legend.text = element_text(size = 5.3), legend.key.height = unit(5, "mm"))
pB <- (pB1 / pB2) + plot_layout(heights = c(1, 14), guides = "collect")
fwrite(eb[, .(type, gene, logFC, P, BH)], file.path(r3, "B_enzyme_logFC_shown.csv"))

## C: KEGG ORA heatmap
cat_map <- c(ECM_RECEPTOR_INTERACTION = "Signaling molecules and interaction", NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION = "Signaling molecules and interaction", CALCIUM_SIGNALING_PATHWAY = "Signal transduction", CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION = "Signaling molecules and interaction",
  FOCAL_ADHESION = "Cellular community", O_GLYCAN_BIOSYNTHESIS = "Glycan biosynthesis", CELL_ADHESION_MOLECULES_CAMS = "Signaling molecules and interaction", COMPLEMENT_AND_COAGULATION_CASCADES = "Immune system", DILATED_CARDIOMYOPATHY = "Cardiovascular disease",
  ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC = "Cardiovascular disease", DRUG_METABOLISM_CYTOCHROME_P450 = "Xenobiotics metabolism", STEROID_HORMONE_BIOSYNTHESIS = "Lipid metabolism", HEMATOPOIETIC_CELL_LINEAGE = "Immune system",
  HYPERTROPHIC_CARDIOMYOPATHY_HCM = "Cardiovascular disease", METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450 = "Xenobiotics metabolism", ARACHIDONIC_ACID_METABOLISM = "Lipid metabolism", AXON_GUIDANCE = "Development",
  SYSTEMIC_LUPUS_ERYTHEMATOSUS = "Immune disease", TYPE_I_DIABETES_MELLITUS = "Endocrine disease", INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION = "Immune system")
o <- ora[type %in% tum]; o$nl <- -log10(pmax(o$FDR, 1e-10)); rk <- o[, .(n = sum(FDR < 0.05), m = mean(nl)), by = pathway][order(-n, -m)]; top <- head(rk$pathway, 20)
oc <- CJ(pathway = top, type = tum); oc <- merge(oc, o[, .(pathway, type, nl, FDR)], by = c("pathway", "type"), all.x = TRUE); oc$nl[is.na(oc$nl)] <- 0
cats <- unname(cat_map[sub("^KEGG_", "", top)]); cats[is.na(cats)] <- "Other"; cdf <- data.frame(pathway = top, cat = cats)
dag <- function(x) ifelse(grepl("^O-Glycan", x), paste0(x, "\u2020"), x); oc$pathway_l <- factor(dag(nice(oc$pathway)), levels = rev(dag(nice(top)))); oc$type <- factor(oc$type, levels = tum); cdf$pathway_l <- factor(dag(nice(cdf$pathway)), levels = rev(dag(nice(top)))); cdf$type <- "Cat"
catcol <- c("Signaling molecules and interaction" = "#B3DE69", "Signal transduction" = "#7FD6D6", "Cellular community" = "#E5C494", "Glycan biosynthesis" = "#8DA0CB", "Immune system" = "#8FD14F", "Cardiovascular disease" = "#F4A6A6", "Xenobiotics metabolism" = "#B3B3B3", "Lipid metabolism" = "#FFD92F", "Development" = "#FC8D62", "Immune disease" = "#E78AC3", "Endocrine disease" = "#66C2A5", "Other" = "grey85")
fwrite(merge(oc, cdf[, c("pathway", "cat")], by = "pathway"), file.path(r3, "C_kegg_ora_shown.csv"))
lv <- c("Cat", tum); oc$x <- as.numeric(factor(as.character(oc$type), levels = lv)); cdf$x <- 1
pC <- ggplot() + geom_tile(data = cdf, aes(x = x, y = pathway_l, fill = cat), colour = "white", linewidth = 0.3, width = 0.9) + scale_fill_manual(values = catcol, name = "Pathway category", guide = guide_legend(order = 2, ncol = 3, keyheight = unit(2.4, "mm"), title.position = "top")) +
  new_scale_fill() + geom_tile(data = oc, aes(x = x, y = pathway_l, fill = nl), colour = "grey85", linewidth = 0.2) + geom_text(data = oc[oc$FDR < 0.05 & !is.na(oc$FDR), ], aes(x = x, y = pathway_l, label = "*"), size = mm(7), vjust = 0.75, family = fam) +
  scale_fill_gradientn(colours = c("white", "#FFF3B0", "#F4A25A", "#E0402B", "#C0151B"), limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10), oob = squish, name = expression(-log[10](FDR)), guide = guide_colourbar(order = 1, barwidth = unit(22, "mm"), barheight = unit(2.2, "mm"), title.position = "top")) +
  scale_x_continuous(breaks = seq_along(lv), labels = c("", tum), expand = c(0.01, 0)) + labs(x = NULL, y = NULL, title = "KEGG enrichment (C2 vs C1 DEGs)", subtitle = "\u2020 pathway contains GALNT genes (partly circular)") + theme_pub() +
  theme(plot.subtitle = element_text(size = 5, colour = "grey30"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.5), axis.text.x = element_text(size = 5.8, angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 5.5), legend.text = element_text(size = 5.2), legend.box.just = "top", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))

## D: Venn (focus + 3 most significant other cancers)
vt <- unique(c(FOCUS, head(setdiff(sm[order(logrank_p)]$type, FOCUS), 3))); sets <- lapply(vt, function(ct) ora[type == ct & FDR < 0.05]$pathway); names(sets) <- vt; print(lengths(sets))
pD <- ggVennDiagram(sets, label = "count", label_alpha = 0, label_size = 2.3, set_size = 2.6, edge_lty = "dashed", edge_size = 0.3) + scale_fill_gradient(low = "white", high = "#B9CDE5", guide = "none") + scale_x_continuous(expand = expansion(mult = 0.12)) + labs(title = "Shared significant KEGG pathways") +
  theme(plot.title = element_text(size = 6.5, face = "bold", family = fam), plot.margin = margin(2, 2, 2, 2), text = element_text(family = fam))
inter_all <- Reduce(intersect, sets); writeLines(c(paste("sets:", paste(names(sets), lengths(sets), collapse = "; ")), paste("shared by all:", paste(inter_all, collapse = ", "))), file.path(r3, "D_venn_summary.txt"))

design <- "AAAAAAAAAAABBBBBBB\nAAAAAAAAAAABBBBBBB\nAAAAAAAAAAABBBBBBB\nCCCCCCCCCCDDDDDDDD"
comp <- wrap_elements(full = pA) + wrap_elements(full = pB) + wrap_elements(full = pC) + wrap_elements(full = pD) + plot_layout(design = design, heights = c(2.1, 2.1, 2.1, 3.7)) & theme(text = element_text(family = fam), plot.tag = element_text(face = "bold", size = 10, family = fam), plot.tag.position = c(0, 1))
comp <- comp + plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(4, 2, 2, 8)))
W <- 7.2; Hh <- 9.5
ggsave(file.path(out, "Figure3.svg"), comp, width = W, height = Hh, device = svglite::svglite)
ggsave(file.path(out, "Figure3.png"), comp, width = W, height = Hh, dpi = 400, device = ragg::agg_png)
ggsave(file.path(out, "Figure3.tiff"), comp, width = W, height = Hh, dpi = 600, device = ragg::agg_tiff, compression = "lzw")
cat("done\n")
