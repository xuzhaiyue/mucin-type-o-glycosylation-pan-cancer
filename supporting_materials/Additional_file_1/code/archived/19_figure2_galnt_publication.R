# Figure 2 for GALNT (freeze v5; enzyme-only analogue of the template's pathway analysis). Spec: 00_meta/figure2_spec_v1.md. Uses tables from 18_figure2_galnt_scores.R.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
root <- "${PROJECT_ROOT}"; r2 <- file.path(root, "04_results/2026-10-12_Figure2_GALNT_run01")
out <- file.path(root, "05_figures/2026-10-12_Figure2_GALNT_run01"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
fam <- "Arial"
genes <- c(paste0("GALNT", 1:18), "GALNTL5", "GALNTL6", "C1GALT1", "C1GALT1C1")
grp <- setNames(c(rep("init", 20), rep("ext", 2)), genes); col_w <- "#D9942B"; col_e <- "#7B52A6"; lab_col <- ifelse(grp == "init", col_w, col_e)
mm <- function(pt) pt / .pt
theme_pub <- function(bs = 6) theme_classic(base_size = bs, base_family = fam) + theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
  plot.title = element_text(size = bs + 0.5, face = "bold"), legend.key.size = unit(2.6, "mm"), legend.title = element_text(size = bs), legend.text = element_text(size = bs - 0.5), plot.margin = margin(2, 2, 2, 2))
d <- fread(file.path(r2, "enzyme_pathway_correlations.csv.gz")); pr <- fread(file.path(r2, "enzyme_proliferation_rho.csv"))
d$strong <- d$BH_partial < 0.05 & abs(d$rho_partial) >= 0.3; nice <- function(x) { y <- tools::toTitleCase(tolower(gsub("_", " ", sub("^(HALLMARK|KEGG)_", "", x))))
  rep <- c("Myc" = "MYC", "Dna" = "DNA", "Uv" = "UV", "Tgf" = "TGF", "Pi3k" = "PI3K", "Akt" = "AKT", "Mtorc1" = "mTORC1", "Mtor" = "mTOR", "Il2" = "IL2", "Il6" = "IL6", "Stat5" = "STAT5", "Stat3" = "STAT3", "Jak" = "JAK",
           "Tnfa" = "TNF-\u03b1", "Nfkb" = "NF-\u03baB", "Kras" = "KRAS", "E2f" = "E2F", "G2m" = "G2M", "Dn" = "DN", "Wnt" = "WNT", "Beta" = "\u03b2", "Nfkb" = "NF-\u03baB", "Interferon Gamma" = "IFN-\u03b3", "Interferon Alpha" = "IFN-\u03b1")
  for (k in names(rep)) y <- gsub(paste0("\\b", k, "\\b"), rep[[k]], y); y <- gsub("\u03b2 Catenin", "\u03b2-Catenin", y); y }
FOCUS <- "PAAD"   # focus cancer: digestive-system restriction (post hoc, disclosed); was KIRC under the pan-cancer ranking
H <- d[source == "Hallmark"]; n_ct <- uniqueN(d$type)

## B: pan-cancer Hallmark landscape (median partial rho, dot = #cancers with strong significant)
M <- dcast(H[, .(rho = median(rho_partial, na.rm = TRUE)), by = .(pathway, enzyme)], pathway ~ enzyme, value.var = "rho"); mat <- as.matrix(M[, -1]); rownames(mat) <- M$pathway
mat[is.na(mat)] <- 0; ord <- hclust(dist(mat), "average")$order; pw_levels <- rev(nice(rownames(mat)[ord]))
bd <- H[, .(rho = median(rho_partial, na.rm = TRUE), n_strong = sum(strong, na.rm = TRUE)), by = .(pathway, enzyme)]; bd$pathway <- factor(nice(bd$pathway), levels = pw_levels); bd$enzyme <- factor(bd$enzyme, levels = genes)
fwrite(bd, file.path(r2, "B_hallmark_landscape.csv"))
pB <- ggplot(bd, aes(enzyme, pathway)) + geom_tile(aes(fill = rho), colour = "grey92", linewidth = 0.15) + geom_point(aes(size = n_strong), shape = 21, fill = NA, colour = "grey35", stroke = 0.3) +
  scale_fill_gradient2(low = "#3B7BB5", mid = "white", high = "#D6604D", limits = c(-0.45, 0.45), oob = squish, name = "Median\npartial rho") +
  scale_size_continuous(range = c(0, 2.6), breaks = c(5, 10, 20), name = "No. of cancers\n(BH<0.05,|rho|>=0.3)") + labs(x = NULL, y = NULL, title = "Pan-cancer enzyme-pathway correlation (33 cancers)", subtitle = "Grey tile: enzyme is a pathway member (excluded)") + theme_pub() +
  theme(plot.subtitle = element_text(size = 5, colour = "grey30"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = lab_col[genes], size = 6), axis.text.y = element_text(size = 5.3),
        legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 5.2), legend.text = element_text(size = 5), legend.key.size = unit(2.2, "mm"), legend.key.width = unit(6, "mm"))

## E: KIRC
K <- H[type == FOCUS]; K$pathway <- factor(nice(K$pathway), levels = pw_levels); K$enzyme <- factor(K$enzyme, levels = genes)
pE <- ggplot(K, aes(enzyme, pathway, fill = rho_partial)) + geom_tile(colour = "grey92", linewidth = 0.15) +
  geom_text(data = K[BH_partial < 0.05 & abs(rho_partial) >= 0.3], label = "*", size = mm(7), vjust = 0.7, family = fam) +
  scale_fill_gradient2(low = "#3B7BB5", mid = "white", high = "#D6604D", limits = c(-0.6, 0.6), oob = squish, name = "Partial rho") +
  labs(x = NULL, y = NULL, title = sprintf("%s (n = %d tumors)", FOCUS, K$n[1]), subtitle = "* BH < 0.05 and |partial rho| >= 0.3") + theme_pub() +
  theme(plot.title = element_text(size = 6.5, face = "bold", margin = margin(l = 10, b = 2)), plot.subtitle = element_text(size = 5, colour = "grey30", margin = margin(l = 10)), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = lab_col[genes], size = 6),
        legend.position = "bottom", legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm"), legend.title = element_text(size = 5.2, vjust = 1), legend.text = element_text(size = 5))

## C: top pathways per group (enzyme x cancer pairs with strong partial correlation)
d$group <- ifelse(grp[d$enzyme] == "init", "Initiation (GALNT)", "Extension (C1GALT1/C1)"); H2 <- d[source == "Hallmark"]
tot <- H2[, .(pos = sum(strong & rho_partial > 0, na.rm = TRUE), neg = sum(strong & rho_partial < 0, na.rm = TRUE), np = sum(!is.na(rho_partial))), by = .(group, pathway)]; tot$pos <- 100 * tot$pos / tot$np; tot$neg <- -100 * tot$neg / tot$np
sel <- unique(c(tot[order(-(pos - neg))][, head(.SD, 1), by = group]$pathway, tot[, .(s = sum(pos - neg)), by = pathway][order(-s)][1:6]$pathway, tot[, .(s = sum(pos - neg)), by = pathway][order(s)][1:6]$pathway))
cd <- melt(tot[pathway %in% sel], id.vars = c("group", "pathway", "np"), measure.vars = c("pos", "neg"), variable.name = "dir", value.name = "pct"); cd$pathway <- nice(cd$pathway)
lv <- unique(nice(tot[pathway %in% sel][, .(s = sum(pos)), by = pathway][order(s)]$pathway)); cd$pathway <- factor(cd$pathway, levels = lv)
fwrite(cd, file.path(r2, "C_top_pathways_by_group.csv"))
pC <- ggplot(cd, aes(pct, pathway, fill = dir)) + geom_col(width = 0.75) + facet_wrap(~ group, nrow = 1) + geom_vline(xintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = c(pos = "#D6604D", neg = "#3B7BB5"), labels = c(pos = "Positive", neg = "Negative"), name = NULL) + labs(x = "% of enzyme-cancer pairs\n(BH<0.05, |partial rho|>=0.3)", y = NULL, title = "Most frequent pathway associations") +
  theme_pub() + theme(strip.background = element_blank(), strip.text = element_text(size = 5.6, face = "bold"), axis.text.y = element_text(size = 5.5), legend.position = "bottom", legend.key.size = unit(2.2, "mm"), panel.spacing = unit(2, "mm"))

## D: proliferation dependence
pm <- pr[, .(rho = median(rho_prolif), lo = quantile(rho_prolif, 0.25), hi = quantile(rho_prolif, 0.75)), by = enzyme]; pm$enzyme <- factor(pm$enzyme, levels = pm$enzyme[order(pm$rho)])
pm$grp <- grp[as.character(pm$enzyme)]
pD <- ggplot(pm, aes(rho, enzyme, fill = grp)) + geom_col(width = 0.7) + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.35, linewidth = 0.25) +
  scale_fill_manual(values = c(init = col_w, ext = col_e), labels = c(init = "Initiation", ext = "Extension"), name = NULL) +
  labs(x = "Spearman rho with proliferation score\n(median and IQR across cancers)", y = NULL, title = "Enzyme-proliferation association") + theme_pub() +
  theme(legend.position = c(0.78, 0.18), axis.text.y = element_text(colour = lab_col[levels(pm$enzyme)], size = 6))

## F: strong associations per cancer
cn <- d[source == "Hallmark", .(n_strong = sum(strong, na.rm = TRUE), n = sum(!is.na(rho_partial)), N = max(n)), by = type]; cn$pct <- 100 * cn$n_strong / cn$n; cn$type <- factor(cn$type, levels = cn$type[order(cn$pct)])
cn$hl <- cn$type == FOCUS; fwrite(cn, file.path(r2, "F_strong_pairs_per_cancer.csv"))
pF <- ggplot(cn, aes(pct, type, fill = hl)) + geom_col(width = 0.75) + scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"), guide = "none") +
  labs(x = "% pairs with strong association\n(BH<0.05, |partial rho|>=0.3)", y = NULL, title = "Strength by cancer") + theme_pub() + theme(plot.margin = margin(2, 4, 2, 9), axis.text.y = element_text(size = 5.3, face = ifelse(levels(cn$type) == FOCUS, "bold", "plain"), colour = ifelse(levels(cn$type) == FOCUS, "#B2182B", "black")))

## A: workflow (designed schematic; 0.1-in units, 72 x 12 canvas)
rr <- function(x0, y0, x1, y1, r = 0.8, n = 8) { th <- seq(0, pi / 2, length.out = n); cx <- c(x1 - r, x0 + r, x0 + r, x1 - r); cy <- c(y1 - r, y1 - r, y0 + r, y0 + r)
  data.frame(x = unlist(lapply(1:4, function(i) cx[i] + r * cos(th + (i - 1) * pi / 2))), y = unlist(lapply(1:4, function(i) cy[i] + r * sin(th + (i - 1) * pi / 2)))) }
circ <- function(cx, cy, r, n = 40) { th <- seq(0, 2 * pi, length.out = n); data.frame(x = cx + r * cos(th), y = cy + r * sin(th)) }
hdr <- "#7B5EA7"; body <- "#F4F0F9"; edge <- "#CFC2E0"; xs <- c(1, 19.5, 38, 56.5); cw <- 16.5
PL <- list(); add <- function(g) PL[[length(PL) + 1]] <<- g
for (i in 1:4) {
  x0 <- xs[i]
  add(geom_polygon(data = rr(x0, 0.4, x0 + cw, 9.6), aes(x, y), fill = body, colour = edge, linewidth = 0.25))
  add(geom_polygon(data = rbind(data.frame(x = c(x0, x0 + cw, x0 + cw), y = c(8.0, 8.0, 8.8)), rr(x0, 8.0, x0 + cw, 9.6)[rr(x0, 8.0, x0 + cw, 9.6)$y > 8.0, ]), aes(x, y), fill = hdr, colour = NA))
  add(geom_polygon(data = circ(x0 + 1.5, 8.8, 1.0), aes(x, y), fill = "white", colour = NA))
  add(annotate("text", x = x0 + 1.5, y = 8.8, label = as.character(i), size = mm(6), fontface = "bold", colour = hdr, family = fam))
  add(annotate("text", x = x0 + 3.0, y = 8.8, label = c("Activity scoring", "Association", "Confounder control", "Summary")[i], size = mm(5.8), fontface = "bold", colour = "white", hjust = 0, family = fam))
  add(annotate("text", x = x0 + cw / 2, y = 2.0, label = c("ssGSEA activity per cancer\n50 Hallmark + 186 KEGG sets", "Spearman correlation: 22 enzymes\nx 33 cancers (BH-adjusted)", "Rank partial correlation\ncontrolling for proliferation", "Strong: BH<0.05, |rho|>=0.3\nPAAD as focus cancer")[i], size = mm(5), lineheight = 0.95, family = fam, colour = "grey15"))
  if (i < 4) add(geom_polygon(data = data.frame(x = c(x0 + cw + 0.45, x0 + cw + 1.55, x0 + cw + 0.45), y = c(4.1, 4.9, 5.7) - 0.0 + c(0.0, 0.0, 0.0)), aes(x, y), fill = hdr, colour = NA))
}
cx <- xs + cw / 2
## icon 1: score heatmap (genes x pathways)
set.seed(3); v <- c(0.9, -0.3, 0.5, -0.8, 0.2, 0.7, -0.5, 0.1, 0.8, -0.2, 0.4, -0.6, 0.6, 0.3, -0.7, 0.0, 0.9, -0.4)
hm <- expand.grid(i = 0:5, j = 0:2); hm$v <- v; add(geom_tile(data = hm, aes(x = cx[1] - 3.7 + i * 1.48, y = 7.1 - j * 1.3, fill = v), width = 1.42, height = 1.24))
## icon 2: scatter + fit
sc <- data.frame(x = c(0.5, 1.1, 1.6, 2.2, 2.6, 3.3, 3.8, 4.4, 5.0, 5.6), y = c(4.2, 4.9, 4.5, 5.6, 5.4, 6.2, 6.0, 6.9, 6.6, 7.2)); sc$x <- cx[2] - 3.1 + sc$x
add(annotate("segment", x = cx[2] - 3.4, xend = cx[2] + 3.0, y = 3.9, yend = 3.9, linewidth = 0.25)); add(annotate("segment", x = cx[2] - 3.4, xend = cx[2] - 3.4, y = 3.9, yend = 7.6, linewidth = 0.25))
add(geom_point(data = sc, aes(x, y), colour = col_w, size = 0.9)); add(annotate("segment", x = cx[2] - 3.0, xend = cx[2] + 2.8, y = 4.2, yend = 7.2, colour = "grey25", linewidth = 0.35))
## icon 3: partial correlation (two circles + proliferation)
add(geom_polygon(data = circ(cx[3] - 1.2, 6.5, 1.6), aes(x, y), fill = col_w, alpha = 0.55, colour = NA)); add(geom_polygon(data = circ(cx[3] + 1.2, 6.5, 1.6), aes(x, y), fill = "#3B7BB5", alpha = 0.45, colour = NA))
add(geom_polygon(data = circ(cx[3], 5.2, 1.6), aes(x, y), fill = "grey55", alpha = 0.55, colour = NA))
add(annotate("text", x = cx[3] - 3.0, y = 6.9, label = "enzyme", size = mm(4), colour = col_w, family = fam, fontface = "bold", hjust = 1)); add(annotate("text", x = cx[3] + 3.0, y = 6.9, label = "pathway", size = mm(4), colour = "#3B7BB5", family = fam, fontface = "bold", hjust = 0))
add(annotate("text", x = cx[3] + 2.0, y = 4.5, label = "proliferation", size = mm(4), colour = "grey35", family = fam, fontface = "bold", hjust = 0))
## icon 4: mini bars with PAAD highlight
bb <- data.frame(y = 7.3 - 0:3 * 1.05, w = c(5.8, 4.6, 3.6, 2.4), f = c("hl", "n", "n", "n")); add(geom_tile(data = bb, aes(x = cx[4] - 3.4 + w / 2, y = y, width = w, height = 0.8, colour = NULL), fill = c("#B2182B", rep("grey65", 3))))
add(annotate("text", x = cx[4] - 3.4 + 5.8 + 0.4, y = 7.3, label = "PAAD", size = mm(4.2), colour = "#B2182B", hjust = 0, fontface = "bold", family = fam))
pA <- ggplot() + PL + scale_fill_gradient2(low = "#3B7BB5", mid = "white", high = "#D6604D", limits = c(-1, 1), guide = "none") +
  annotate("text", x = 4.5, y = 11.3, label = "Pathway analysis workflow", size = mm(6.5), fontface = "bold", hjust = 0, family = fam) +
  annotate("rect", xmin = 45.5, xmax = 46.7, ymin = 10.9, ymax = 11.7, fill = col_w) + annotate("text", x = 47.1, y = 11.3, label = "Initiation (GALNT)", size = mm(5), hjust = 0, family = fam) +
  annotate("rect", xmin = 59.5, xmax = 60.7, ymin = 10.9, ymax = 11.7, fill = col_e) + annotate("text", x = 61.1, y = 11.3, label = "Extension (C1GALT1/C1)", size = mm(5), hjust = 0, family = fam) +
  coord_cartesian(xlim = c(0, 73), ylim = c(0, 12), expand = FALSE) + theme_void() + theme(plot.margin = margin(1, 1, 1, 1))

design <- "AAAAAAAAAAAAAAAAAAAA
BBBBBBBBBBBBCCCCCCCC
DDDDDDDDDEEEEEEFFFFF"
comp <- wrap_elements(full = pA) + pB + pE + pC + pD + pF + plot_layout(design = design, heights = c(1.15, 5.15, 2.6)) & theme(text = element_text(family = fam), plot.tag = element_text(face = "bold", size = 10, family = fam), plot.tag.position = c(0, 1))
comp <- comp + plot_annotation(tag_levels = "A")
W <- 7.2; Hh <- 9.4
ggsave(file.path(out, "Figure2.svg"), comp, width = W, height = Hh, device = svglite::svglite)
ggsave(file.path(out, "Figure2.png"), comp, width = W, height = Hh, dpi = 400, device = ragg::agg_png)
ggsave(file.path(out, "Figure2.tiff"), comp, width = W, height = Hh, dpi = 600, device = ragg::agg_tiff, compression = "lzw")
cat("done\n")
