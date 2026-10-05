#!/usr/bin/env Rscript
# =============================================================================
# Isoform-selective PI3K inhibition in NCI-H69 small-cell lung cancer cells
# Re-analysis of GEO GSE40564 (Wojtalla et al., Clin Cancer Res 2013; PMID 23172887)
# Platform: Affymetrix HG-U133 Plus 2.0 (GPL570); n = 3 per group
#   Ctr  = DMSO vehicle, 24 h
#   PIK  = PIK-75 100 nM (p110alpha-selective inhibitor), 24 h
#   TGX  = TGX-221 (p110beta-selective inhibitor), 24 h
#
# Outputs (in ./results):
#   Table_S1_DE_<contrast>.csv       full gene-level limma tables
#   Table_S2_GSEA_Hallmark.csv       fgsea results, all contrasts
#   Table_S3_GSEA_Reactome.csv
#   Table_S4_PROGENy.csv, Table_S5_TF_activity.csv
#   Table_1_DE_summary.csv           DEG counts per contrast
#   Figure1_heatmap.(pdf|png)        row z-score heatmap, top DEGs, 9 samples
#   Figure2_PCA_volcano.(pdf|png)    PCA + volcano plots (same axes)
#   Figure3_GSEA_hallmark.(pdf|png)  Hallmark NES dot plot
#   Figure4_PROGENy_TF.(pdf|png)     pathway and TF activity
#   sessionInfo.txt
#
# Requires: R >= 4.3; Bioconductor: GEOquery, limma, fgsea, progeny, dorothea,
#           decoupleR; CRAN: msigdbr (>= 10), pheatmap, ggplot2, ggrepel, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(fgsea); library(msigdbr)
  library(pheatmap); library(ggplot2); library(ggrepel); library(patchwork)
  library(progeny); library(dorothea); library(decoupleR)
})
set.seed(20260930)
dir.create("results", showWarnings = FALSE)
out <- function(f) file.path("results", f)

# ---- 0. Analysis parameters (report these in Methods) -----------------------
FDR_CUT   <- 0.05
LFC_CUT   <- 1          # |log2FC| >= 1 (2-fold) for "DEG" counts and heatmap
N_HEAT    <- 60         # top genes (moderated F) shown in Figure 1
N_PCA     <- 1000       # most variable genes used for PCA

# ---- 1. Download series matrix + GPL570 annotation ---------------------------
gse  <- getGEO("GSE40564", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = ".")[[1]]
expr <- exprs(gse)                 # submitter-processed RMA, log2 scale
pd   <- pData(gse)
fd   <- fData(gse)
stopifnot(max(expr, na.rm = TRUE) < 20)   # confirm log2 scale

grp <- factor(ifelse(grepl("Ctr", pd$title), "DMSO",
               ifelse(grepl("PIK75", pd$title), "PIK75", "TGX221")),
              levels = c("DMSO", "PIK75", "TGX221"))
names(grp) <- colnames(expr)
short <- sub("^H69_", "", pd$title)             # e.g. Ctr_1, PIK75_3
print(data.frame(GSM = colnames(expr), title = pd$title, group = grp))

# ---- 2. Probe filtering and collapsing -----------------------------------------
sym <- as.character(fd[rownames(expr), "Gene symbol"])
keep <- !grepl("^AFFX", rownames(expr)) &            # control probes
        !is.na(sym) & sym != "" &                     # unannotated
        !grepl("///", sym, fixed = TRUE)              # multi-gene (cross-hybridising)
expr <- expr[keep, ]; sym <- sym[keep]

# Expression filter: probe above the 25th percentile of all intensities
# in at least 3 samples (= size of the smallest group)
thr  <- quantile(expr, 0.25)
expr <- expr[rowSums(expr > thr) >= 3, ]
sym  <- as.character(fd[rownames(expr), "Gene symbol"])

# Collapse to one probe per gene: the probe with the highest mean intensity
# (WGCNA::collapseRows "MaxMean"). Chosen BEFORE and independently of any
# contrast, so the same probe represents a gene in every comparison.
o     <- order(sym, -rowMeans(expr))
first <- !duplicated(sym[o])
gexpr <- expr[o[first], ]
rownames(gexpr) <- sym[o[first]]
probe_used <- setNames(rownames(expr)[o[first]], sym[o[first]])
cat(sprintf("Probes: %d on array -> %d after filtering -> %d genes\n",
            nrow(exprs(gse)), nrow(expr), nrow(gexpr)))

# ---- 3. limma: all three pairwise contrasts ----------------------------------
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
fit  <- lmFit(gexpr, design)
cm   <- makeContrasts(PIK75_vs_DMSO   = PIK75  - DMSO,
                      TGX221_vs_DMSO  = TGX221 - DMSO,
                      TGX221_vs_PIK75 = TGX221 - PIK75, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)

tabs <- lapply(colnames(cm), function(k) {
  tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt$probe <- probe_used[tt$gene]
  tt[order(tt$P.Value), c("gene", "probe", "logFC", "AveExpr", "t",
                          "P.Value", "adj.P.Val")]
})
names(tabs) <- colnames(cm)
for (k in names(tabs))
  write.csv(tabs[[k]], out(paste0("Table_S1_DE_", k, ".csv")), row.names = FALSE)

de_summary <- do.call(rbind, lapply(names(tabs), function(k) {
  t <- tabs[[k]]; s <- t$adj.P.Val < FDR_CUT
  data.frame(contrast = k, genes_tested = nrow(t),
             FDR05 = sum(s),
             FDR05_up_lfc1 = sum(s & t$logFC >=  LFC_CUT),
             FDR05_down_lfc1 = sum(s & t$logFC <= -LFC_CUT))
}))
write.csv(de_summary, out("Table_1_DE_summary.csv"), row.names = FALSE)
print(de_summary)

# ---- 4. Figure 1: row z-score heatmap, top genes by moderated F ------------
ftab <- topTable(fit2, coef = 1:2, number = Inf)     # F-test: any change vs DMSO
ftab <- ftab[ftab$adj.P.Val < FDR_CUT &
             apply(abs(ftab[, c("PIK75_vs_DMSO", "TGX221_vs_DMSO")]), 1, max) >= LFC_CUT, ]
top  <- head(rownames(ftab), N_HEAT)
z    <- t(scale(t(gexpr[top, ])))
colnames(z) <- short
ann  <- data.frame(Treatment = grp, row.names = short)
ann_col <- list(Treatment = c(DMSO = "#8c8c8c", PIK75 = "#c0392b", TGX221 = "#2e86c1"))
pal  <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
hm <- function(file, ...) pheatmap(z, color = pal, breaks = seq(-2.5, 2.5, length.out = 102),
         annotation_col = ann, annotation_colors = ann_col,
         cluster_cols = TRUE, clustering_distance_rows = "correlation",
         show_rownames = TRUE, fontsize_row = 6.5, border_color = NA,
         main = sprintf("Top %d DEGs vs DMSO (row z-score)", length(top)), filename = file, ...)
hm(out("Figure1_heatmap.pdf"), width = 6.5, height = 9)
hm(out("Figure1_heatmap.png"), width = 6.5, height = 9)

# ---- 5. Figure 2: PCA + volcanoes on identical axes ------------------------
v   <- apply(gexpr, 1, var)
pcx <- prcomp(t(gexpr[order(-v)[1:N_PCA], ]), center = TRUE, scale. = FALSE)
ve  <- 100 * pcx$sdev^2 / sum(pcx$sdev^2)
pcd <- data.frame(pcx$x[, 1:2], group = grp, label = short)
cols <- ann_col$Treatment
p_pca <- ggplot(pcd, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3.2) + geom_text_repel(aes(label = label), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = sprintf("A  PCA, top %d variable genes", N_PCA)) +
  theme_classic(base_size = 11) + theme(legend.position = "top")

volc <- function(k, ttl, xl, yl, lab_genes) {
  d <- tabs[[k]]; d$y <- -log10(d$adj.P.Val)
  d$cls <- ifelse(d$adj.P.Val < FDR_CUT & d$logFC >=  LFC_CUT, "Up",
           ifelse(d$adj.P.Val < FDR_CUT & d$logFC <= -LFC_CUT, "Down", "NS"))
  n_up <- sum(d$cls == "Up"); n_dn <- sum(d$cls == "Down")
  ggplot(d, aes(logFC, y, colour = cls)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = 2, linewidth = 0.3) +
    geom_text_repel(data = subset(d, gene %in% lab_genes), aes(label = gene),
                    colour = "black", size = 2.8, max.overlaps = 30, min.segment.length = 0) +
    scale_colour_manual(values = c(Up = "#b2182b", Down = "#2166ac", NS = "grey75"), guide = "none") +
    coord_cartesian(xlim = xl, ylim = yl) +
    labs(title = sprintf("%s  (up %d / down %d)", ttl, n_up, n_dn),
         x = "log2 fold change", y = "-log10 FDR") + theme_classic(base_size = 11)
}
xl <- range(c(tabs$PIK75_vs_DMSO$logFC, tabs$TGX221_vs_DMSO$logFC))
yl <- c(0, max(-log10(tabs$PIK75_vs_DMSO$adj.P.Val)))
lab <- c("BRCA1", "MCL1", "E2F1", "MYC", "CCNE2", "MCM10", "RAD51", "DDIT4", "TXNIP",
         "HSPA5", "CTNNB1", "BCL2", "PIK3IP1", "IRS2", "CDKN1B")
p_v1 <- volc("PIK75_vs_DMSO",  "B  PIK-75 vs DMSO",  xl, yl, lab)
p_v2 <- volc("TGX221_vs_DMSO", "C  TGX-221 vs DMSO", xl, yl, lab)
fig2 <- p_pca + p_v1 + p_v2 + plot_layout(ncol = 3, widths = c(1, 1, 1))
ggsave(out("Figure2_PCA_volcano.pdf"), fig2, width = 15, height = 5)
ggsave(out("Figure2_PCA_volcano.png"), fig2, width = 15, height = 5, dpi = 300)

# ---- 6. GSEA: Hallmark + Reactome, ranked by moderated t -------------------
H  <- msigdbr(species = "Homo sapiens", collection = "H")
RE <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
toList <- function(df) split(df$gene_symbol, df$gs_name)
run_gsea <- function(sets) do.call(rbind, lapply(names(tabs), function(k) {
  r <- setNames(tabs[[k]]$t, tabs[[k]]$gene)
  g <- fgsea(pathways = sets, stats = sort(r, decreasing = TRUE),
             minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000)
  g$contrast <- k
  g$leadingEdge <- vapply(g$leadingEdge, function(x) paste(head(x, 15), collapse = ";"), "")
  as.data.frame(g[order(g$padj), ])
}))
gH <- run_gsea(toList(H));  write.csv(gH, out("Table_S2_GSEA_Hallmark.csv"), row.names = FALSE)
gR <- run_gsea(toList(RE)); write.csv(gR, out("Table_S3_GSEA_Reactome.csv"), row.names = FALSE)

# Figure 3: Hallmark dot plot (sets significant in >= 1 contrast)
sig <- unique(gH$pathway[gH$padj < FDR_CUT])
dd  <- subset(gH, pathway %in% sig)
dd$pathway  <- gsub("_", " ", sub("^HALLMARK_", "", dd$pathway))
ord <- with(subset(dd, contrast == "PIK75_vs_DMSO"), pathway[order(NES)])
dd$pathway  <- factor(dd$pathway, levels = ord)
dd$contrast <- factor(dd$contrast, levels = names(tabs),
                      labels = c("PIK-75 vs DMSO", "TGX-221 vs DMSO", "TGX-221 vs PIK-75"))
dd$sig <- dd$padj < FDR_CUT
p3 <- ggplot(dd, aes(contrast, pathway)) +
  geom_point(aes(size = -log10(padj), fill = NES, alpha = sig), shape = 21, colour = "grey30") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), name = "FDR < 0.05") +
  scale_size_continuous(range = c(1, 6), name = "-log10 FDR") +
  labs(x = NULL, y = NULL, title = "MSigDB Hallmark GSEA (fgsea, ranked by moderated t)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(out("Figure3_GSEA_hallmark.pdf"), p3, width = 7, height = 0.22 * length(sig) + 2)
ggsave(out("Figure3_GSEA_hallmark.png"), p3, width = 7, height = 0.22 * length(sig) + 2, dpi = 300)

# ---- 7. Pathway (PROGENy) and TF (DoRothEA A-C) activity -------------------
tmat <- sapply(tabs, function(t) setNames(t$t, t$gene)[rownames(gexpr)])
# PROGENy footprint weights shipped with the progeny package (no web access needed)
pm    <- as.matrix(progeny::getModel("Human", top = 500))
net_p <- data.frame(source = rep(colnames(pm), each = nrow(pm)),
                    target = rep(rownames(pm), times = ncol(pm)),
                    weight = as.vector(pm))
net_p <- subset(net_p, weight != 0)
pa <- run_mlm(mat = tmat, net = net_p, .source = "source", .target = "target",
              .mor = "weight", minsize = 5)
write.csv(pa, out("Table_S4_PROGENy.csv"), row.names = FALSE)

data(dorothea_hs, package = "dorothea")
net_t <- subset(dorothea_hs, confidence %in% c("A", "B", "C"))
ta <- run_ulm(mat = tmat, net = net_t, .source = "tf", .target = "target",
              .mor = "mor", minsize = 10)
ta$p_adj <- ave(ta$p_value, ta$condition, FUN = function(p) p.adjust(p, "BH"))
write.csv(ta, out("Table_S5_TF_activity.csv"), row.names = FALSE)

lbl <- c(PIK75_vs_DMSO = "PIK-75 vs DMSO", TGX221_vs_DMSO = "TGX-221 vs DMSO",
         TGX221_vs_PIK75 = "TGX-221 vs PIK-75")
pa$condition <- factor(lbl[pa$condition], levels = lbl)
pa$source <- factor(pa$source, levels = with(subset(pa, condition == lbl[1]), source[order(score)]))
p4a <- ggplot(subset(pa, condition != lbl[3]), aes(score, source, fill = score)) +
  geom_col() + facet_wrap(~ condition) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", guide = "none") +
  labs(x = "PROGENy activity score (MLM t)", y = NULL, title = "A  Pathway activity") +
  theme_bw(base_size = 10)
tt1 <- subset(ta, condition == "PIK75_vs_DMSO")
top_tf <- c(head(tt1$source[order(tt1$score)], 12), head(tt1$source[order(-tt1$score)], 12))
tb <- subset(ta, source %in% top_tf & condition != "TGX221_vs_PIK75")
tb$condition <- factor(lbl[tb$condition], levels = lbl[1:2])
tb$source <- factor(tb$source, levels = tt1$source[order(tt1$score)][tt1$source[order(tt1$score)] %in% top_tf])
p4b <- ggplot(tb, aes(score, source, fill = score)) + geom_col() + facet_wrap(~ condition) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", guide = "none") +
  labs(x = "TF activity (ULM t, DoRothEA A-C)", y = NULL, title = "B  Transcription-factor activity") +
  theme_bw(base_size = 10)
fig4 <- p4a + p4b + plot_layout(ncol = 2)
ggsave(out("Figure4_PROGENy_TF.pdf"), fig4, width = 12, height = 6)
ggsave(out("Figure4_PROGENy_TF.png"), fig4, width = 12, height = 6, dpi = 300)

# ---- 8. Key-gene table for the text -----------------------------------------
key <- c("BRCA1", "BRCA2", "RAD51", "FANCD2", "E2F1", "MYBL2", "MKI67", "CCNE2", "MCM10",
         "MCL1", "BCL2", "BCL2L11", "PMAIP1", "MYC", "DDIT4", "TXNIP", "PIK3IP1", "IRS2",
         "CDKN1B", "HSPA5", "CTNNB1", "TP53", "PIK3CA", "PIK3CB", "ASCL1", "NEUROD1")
kt <- do.call(cbind, lapply(names(tabs), function(k) {
  t <- tabs[[k]]; rownames(t) <- t$gene; t <- t[intersect(key, t$gene), ]
  setNames(data.frame(round(t$logFC, 2), signif(t$adj.P.Val, 2)),
           paste0(k, c("_log2FC", "_FDR")))
}))
kt <- cbind(gene = intersect(key, tabs[[1]]$gene), probe = probe_used[intersect(key, tabs[[1]]$gene)], kt)
write.csv(kt, out("Table_2_key_genes.csv"), row.names = FALSE)


# ---- 9. Sensitivity: probe-collapsing rule (MaxMean vs limma::avereps) -------
# Several GPL570 genes have probes that disagree (3'-UTR isoforms, cross-
# hybridisation). Re-fit with probe averaging and report concordance.
aexpr <- avereps(expr, ID = sym)
afit  <- eBayes(contrasts.fit(lmFit(aexpr, design), cm), trend = TRUE, robust = TRUE)
sens <- do.call(rbind, lapply(colnames(cm), function(k) {
  a <- topTable(afit, coef = k, number = Inf, sort.by = "none")
  m <- tabs[[k]][match(rownames(a), tabs[[k]]$gene), ]
  data.frame(contrast = k,
             r_logFC = round(cor(a$logFC, m$logFC, use = "complete.obs"), 3),
             FDR05_avereps = sum(a$adj.P.Val < FDR_CUT),
             FDR05_maxmean = sum(m$adj.P.Val < FDR_CUT))
}))
write.csv(sens, out("Table_S6_probe_collapsing_sensitivity.csv"), row.names = FALSE)
print(sens)

# Probe-level view of genes whose direction depends on the probe chosen
pfit <- eBayes(contrasts.fit(lmFit(expr, design), cm), trend = TRUE, robust = TRUE)
pt   <- topTable(pfit, coef = "PIK75_vs_DMSO", number = Inf, sort.by = "none")
pt$gene <- sym; pt$probe <- rownames(pt); pt$AveExpr <- round(pt$AveExpr, 2)
chk <- c("MCL1", "BCL2", "HSPA5", "EGFR", "MDM2", "CTNNB1", "BRCA1", "TP53")
pp  <- pt[pt$gene %in% chk, c("gene", "probe", "AveExpr", "logFC", "adj.P.Val")]
pp$used_in_main <- pp$probe %in% probe_used
write.csv(pp[order(pp$gene, -pp$AveExpr), ], out("Table_S7_probe_level_key_genes.csv"), row.names = FALSE)

writeLines(c(sprintf("PCA variance explained: PC1 %.1f%%, PC2 %.1f%%", ve[1], ve[2]),
             capture.output(sessionInfo())), out("sessionInfo.txt"))
cat("Done.\n")
