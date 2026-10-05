#!/usr/bin/env Rscript
# =============================================================================
# Part 3: robustness audit of the GSE40564 analysis
#   A. Array QC and re-processing from CEL files (own RMA; spike-in-anchored
#      normalisation without the quantile step)
#   B. Statistical robustness (treat, arrayWeights, leave-one-out, probe rule)
#   C. Is the response "global transcriptional repression"?  mRNA half-life test
#   D. Are the "induced" sets real?  ISR genes, ribosomal genes, self-contained tests
#   E. PIK-75 vs TGX-221 effect sizes; SREBP, FOXO/feedback gene panels
#   F. Footprint nulls and regulon overlap
# Inputs: ./cel/*.CEL.gz (GSE40564_RAW.tar), ./ext/human_halflife_PC1.csv
#         (Agarwal & Kelley, Genome Biol 2022;23:245, Additional file 3, sheet 'human', columns
#          'Gene name' and 'half-life (PC1)', saved with header gene,halflife_PC1;
#          https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02811-x/MediaObjects/13059_2022_2811_MOESM3_ESM.xlsx)
# Run with:  R_THREADS=1 Rscript GSE40564_cel_preprocess.R && Rscript GSE40564_audit.R
# =============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(fgsea); library(msigdbr)
  library(decoupleR); library(dorothea); library(progeny); library(data.table)
})
set.seed(20260930)
dir.create("results", showWarnings = FALSE); out <- function(f) file.path("results", f)
FDR_CUT <- 0.05; LFC_CUT <- 1
say <- function(...) cat(sprintf(...), "\n")

# ---- shared helpers ---------------------------------------------------------
gse <- getGEO(filename = "GSE40564_series_matrix.txt.gz", getGPL = FALSE)
ann <- fread("GPL570.annot.gz", skip = "ID\tGene title", fill = TRUE, sep = "\t", quote = "")
symOf <- setNames(ann[["Gene symbol"]], ann$ID)
grp <- factor(rep(c("DMSO", "PIK75", "TGX221"), each = 3), levels = c("DMSO", "PIK75", "TGX221"))
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
cm <- makeContrasts(PIK75_vs_DMSO = PIK75 - DMSO, TGX221_vs_DMSO = TGX221 - DMSO, levels = design)

collapse <- function(e, thr_q = 0.25) {          # same rule as the main pipeline
  s <- symOf[rownames(e)]
  k <- !grepl("^AFFX", rownames(e)) & !is.na(s) & s != "" & !grepl("///", s, fixed = TRUE)
  e <- e[k, ]; e <- e[rowSums(e > quantile(e, thr_q)) >= 3, ]; s <- symOf[rownames(e)]
  o <- order(s, -rowMeans(e)); f <- !duplicated(s[o]); g <- e[o[f], ]; rownames(g) <- s[o[f]]
  list(g = g, probes = e, sym = s)
}
fitDE <- function(g, w = NULL) eBayes(contrasts.fit(lmFit(g, design, weights = w), cm), trend = TRUE, robust = TRUE)
tt <- function(fit, k) { x <- topTable(fit, coef = k, number = Inf, sort.by = "none"); x$gene <- rownames(x); x }
cnt <- function(x) c(FDR05 = sum(x$adj.P.Val < FDR_CUT),
                     up2 = sum(x$adj.P.Val < FDR_CUT & x$logFC >= LFC_CUT),
                     down2 = sum(x$adj.P.Val < FDR_CUT & x$logFC <= -LFC_CUT),
                     frac_negative = round(mean(x$logFC < 0), 3))

sub_e <- exprs(gse)                               # submitter RMA
main  <- collapse(sub_e); g <- main$g
fit   <- fitDE(g); C1 <- tt(fit, "PIK75_vs_DMSO"); C2 <- tt(fit, "TGX221_vs_DMSO")

# ---- A. Re-processing from CEL files ----------------------------------------
cp   <- readRDS("cel/cel_processed.rds")          # from GSE40564_cel_preprocess.R
own  <- cp$own[, colnames(sub_e)]; nonq <- cp$nonq[, colnames(sub_e)]; qc <- as.data.table(cp$qc)
spk  <- grep("^AFFX-(r2-)?(Ec-)?(Bio[BCD]|bio[BCD]|P1-cre|Cre)|^AFFX-(r2-Bs-)?(lys|phe|thr|dap|Lys|Phe|Thr|Dap)", rownames(nonq), value = TRUE)
hyb  <- grep("[Bb]io[BCD]|[Cc]re", spk, value = TRUE); polyA <- setdiff(spk, hyb)
qc[, hyb_ctrl := round(colMeans(nonq[hyb, ]), 2)][, polyA_ctrl := round(colMeans(nonq[polyA, ]), 2)]
qc[, ACTB_5p_minus_3p := round(nonq["AFFX-HSAC07/X00351_5_at", ] - nonq["AFFX-HSAC07/X00351_3_at", ], 2)]
anch <- sweep(nonq, 2, colMeans(nonq[spk, ]) - mean(nonq[spk, ]))   # anchor arrays on the spike-in controls
fwrite(qc, out("Table_S11_array_QC.csv")); print(qc)

reproc <- rbindlist(lapply(list(submitter_RMA = sub_e, own_RMA = own, spike_anchored_no_quantile = anch), function(e) {
  f <- fitDE(collapse(e)$g)
  rbindlist(lapply(colnames(cm), function(k) { x <- tt(f, k); m <- if (k == "PIK75_vs_DMSO") C1 else C2
    i <- intersect(x$gene, m$gene)
    data.table(contrast = k, genes = nrow(x), t(cnt(x)),
               r_logFC_vs_main = round(cor(x[i, "logFC"], m[i, "logFC"]), 3)) }))
}), idcol = "processing")
fwrite(reproc, out("Table_S12_reprocessing_sensitivity.csv")); print(reproc)
A <- tt(fitDE(collapse(anch)$g), "PIK75_vs_DMSO")           # spike-anchored C1, used below

# ---- B. Statistical robustness -----------------------------------------------
tr <- function(k, lfc) { x <- topTreat(treat(contrasts.fit(lmFit(g, design), cm), lfc = lfc, trend = TRUE, robust = TRUE),
                                       coef = k, number = Inf); c(up = sum(x$adj.P.Val < FDR_CUT & x$logFC > 0), down = sum(x$adj.P.Val < FDR_CUT & x$logFC < 0)) }
treat_tab <- rbindlist(lapply(colnames(cm), function(k) rbindlist(lapply(c(log2(1.5), 1), function(l)
  data.table(contrast = k, treat_lfc = round(l, 2), t(tr(k, l)))))))
aw <- arrayWeights(g, design); names(aw) <- qc$sample
fw <- fitDE(g, w = matrix(aw, nrow(g), 9, byrow = TRUE))
loo <- rbindlist(lapply(7:9, function(i) { d <- design[-i, ]; f <- eBayes(contrasts.fit(lmFit(g[, -i], d), cm), trend = TRUE, robust = TRUE)
  x <- tt(f, "TGX221_vs_DMSO"); s <- x$gene[x$adj.P.Val < FDR_CUT]; full <- C2$gene[C2$adj.P.Val < FDR_CUT]
  data.table(dropped = qc$sample[i], FDR05 = length(s), overlap_with_full_358 = length(intersect(s, full)),
             r_logFC = round(cor(x$logFC, C2[x$gene, "logFC"]), 3)) }))
av <- avereps(main$probes, ID = main$sym); fa <- fitDE(av)
conc <- rbindlist(lapply(colnames(cm), function(k) { a <- tt(fa, k); m <- tt(fit, k); i <- intersect(a$gene, m$gene)
  sa <- a$gene[a$adj.P.Val < FDR_CUT & abs(a$logFC) >= LFC_CUT]; sm <- m$gene[m$adj.P.Val < FDR_CUT & abs(m$logFC) >= LFC_CUT]
  both <- i[a[i, "adj.P.Val"] < FDR_CUT & m[i, "adj.P.Val"] < FDR_CUT]
  data.table(contrast = k, pearson = round(cor(a[i, "logFC"], m[i, "logFC"]), 3),
             spearman = round(cor(a[i, "logFC"], m[i, "logFC"], method = "spearman"), 3),
             sign_concordance_sig_both = round(mean(sign(a[both, "logFC"]) == sign(m[both, "logFC"])), 3),
             DEG_maxmean = length(sm), DEG_avereps = length(sa),
             jaccard_DEG = round(length(intersect(sa, sm)) / max(1, length(union(sa, sm))), 3)) }))
fwrite(treat_tab, out("Table_S13_treat.csv")); fwrite(loo, out("Table_S14_TGX_leave_one_out.csv"))
fwrite(conc, out("Table_S6b_collapse_concordance.csv"))
print(treat_tab); say("arrayWeights: %s", paste(names(aw), round(aw, 2), collapse = "  "))
say("C2 FDR<0.05 with arrayWeights: %d", sum(tt(fw, "TGX221_vs_DMSO")$adj.P.Val < FDR_CUT)); print(loo); print(conc)

# ---- C. mRNA half-life: is the response what transcription shut-off predicts? ----
hl <- fread("ext/human_halflife_PC1.csv")
hlf <- function(x, lab) { d <- merge(data.table(gene = x$gene, logFC = x$logFC, t = x$t, A = x$AveExpr), hl, by = "gene")
  ct <- cor.test(d$logFC, d$halflife_PC1, method = "spearman", exact = FALSE)
  q <- cut(d$halflife_PC1, quantile(d$halflife_PC1, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))
  data.table(analysis = lab, n = nrow(d), spearman_rho = round(ct$estimate, 3), p = signif(ct$p.value, 2),
             t(round(tapply(d$logFC, q, mean), 3))) }
half <- rbind(hlf(C1, "PIK-75 vs DMSO (main)"), hlf(C2, "TGX-221 vs DMSO (main)"),
              hlf(A, "PIK-75 vs DMSO (spike-anchored)"))
fwrite(half, out("Table_S15_halflife.csv")); print(half)
# half-life-residualised statistics for GSEA
d1 <- merge(data.table(gene = C1$gene, t = C1$t), hl, by = "gene")
d1[, t_resid := resid(lm(t ~ splines::ns(halflife_PC1, 4)))]
say("Variance in C1 moderated t explained by half-life: %.1f%%", 100 * summary(lm(t ~ splines::ns(halflife_PC1, 4), d1))$r.squared)

# ---- D. Are the 'induced' gene sets real? ------------------------------------
H  <- msigdbr(species = "Homo sapiens", collection = "H");  HL <- split(H$gene_symbol, H$gs_name)
RE <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"); RL <- split(RE$gene_symbol, RE$gs_name)
focus <- c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
           "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
           "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1",
           "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION", "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
           "REACTOME_CHOLESTEROL_BIOSYNTHESIS", "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS")
SETS <- c(HL, RL)[focus]; isRP <- function(x) grepl("^RP[LS][0-9]|^RPLP|^RPSA$", x)
idx <- ids2indices(SETS, rownames(g))
ro  <- lapply(colnames(cm), function(k) { r <- mroast(g, idx, design, contrast = cm[, k], nrot = 4999); r$set <- rownames(r); r$contrast <- k; r })
gs <- function(stat, sets) as.data.table(fgsea(sets, sort(stat, decreasing = TRUE), minSize = 10, maxSize = 600, eps = 0, nPermSimple = 10000))[, .(pathway, NES, padj)]
g_main  <- gs(setNames(C1$t, C1$gene), SETS); g_resid <- gs(setNames(d1$t_resid, d1$gene), SETS)
g_anch  <- gs(setNames(A$t, A$gene), SETS);   g_tgx   <- gs(setNames(C2$t, C2$gene), SETS)
eff <- rbindlist(lapply(names(SETS), function(s) { gg <- intersect(SETS[[s]], rownames(g)); r1 <- ro[[1]][s, ]; r2 <- ro[[2]][s, ]
  data.table(set = sub("^HALLMARK_|^REACTOME_", "", s), n = length(gg), frac_ribosomal = round(mean(isRP(gg)), 2),
    PIK_mean_log2FC = round(mean(C1[gg, "logFC"]), 3), TGX_mean_log2FC = round(mean(C2[gg, "logFC"]), 3),
    PIK_pct_down = round(100 * mean(C1[gg, "logFC"] < 0)), PIK_anchored_mean_log2FC = round(mean(A[intersect(gg, A$gene), "logFC"]), 3),
    PIK_NES = round(g_main[pathway == s, NES], 2), PIK_NES_halflife_resid = round(g_resid[pathway == s, NES], 2),
    PIK_NES_anchored = round(g_anch[pathway == s, NES], 2), TGX_NES = round(g_tgx[pathway == s, NES], 2),
    PIK_roast = paste0(r1$Direction, " p=", signif(r1$FDR, 2)), TGX_roast = paste0(r2$Direction, " p=", signif(r2$FDR, 2))) }))
fwrite(eff, out("Table_S16_set_effect_sizes.csv")); print(eff)
say("Background: mean log2FC of all genes, PIK-75 %.3f (%.0f%% negative); TGX-221 %.3f", mean(C1$logFC), 100 * mean(C1$logFC < 0), mean(C2$logFC))
say("Ribosomal protein genes (n=%d): mean log2FC PIK-75 %.3f, spike-anchored %.3f", sum(isRP(C1$gene)), mean(C1$logFC[isRP(C1$gene)]), mean(A$logFC[isRP(A$gene)]))

# ---- E. Gene panels ------------------------------------------------------------
panels <- list(
  ISR_ATF4 = c("ATF4", "DDIT3", "ASNS", "TRIB3", "CHAC1", "SLC7A11", "PPP1R15A", "ATF3", "SESN2", "PSAT1", "PHGDH", "SLC7A5", "VEGFA", "DDIT4", "HSPA5", "XBP1"),
  SREBP = c("SREBF1", "SREBF2", "HMGCR", "HMGCS1", "LDLR", "INSIG1", "SQLE", "FDFT1", "LSS", "DHCR7", "MVD", "IDI1", "FASN", "SCD", "ACACA", "ACLY"),
  FOXO_feedback = c("CDKN1B", "MYC", "IRS2", "INSR", "IGF1R", "ERBB3", "PIK3IP1", "SESN3", "GADD45A", "BCL6", "PDK4", "KLF2", "CITED2", "TXNIP", "FOXO1", "FOXO3"),
  Proliferation_HR = c("E2F1", "MKI67", "CCNE2", "CCNE1", "MCM2", "PCNA", "TOP2A", "BRCA1", "BRCA2", "RAD51", "FANCD2", "SLFN11", "MCL1", "BCL2"))
pan <- rbindlist(lapply(names(panels), function(p) { gg <- intersect(panels[[p]], C1$gene)
  data.table(panel = p, gene = gg, PIK_log2FC = round(C1[gg, "logFC"], 2), PIK_FDR = signif(C1[gg, "adj.P.Val"], 2),
             PIK_anchored_log2FC = round(A[gg, "logFC"], 2), TGX_log2FC = round(C2[gg, "logFC"], 2), TGX_FDR = signif(C2[gg, "adj.P.Val"], 2),
             AveExpr = round(C1[gg, "AveExpr"], 1)) }))
fwrite(pan, out("Table_S17_gene_panels.csv")); print(pan, nrows = 100)

# ---- F. Footprint nulls and regulon overlap -------------------------------------
tmat <- cbind(PIK75 = setNames(C1$t, C1$gene), TGX221 = setNames(C2$t, C2$gene)[C1$gene])
pm <- as.matrix(progeny::getModel("Human", top = 500))
net_p <- subset(data.frame(source = rep(colnames(pm), each = nrow(pm)), target = rep(rownames(pm), ncol(pm)), weight = as.vector(pm)), weight != 0)
obs <- as.data.table(run_mlm(tmat, net_p, .source = "source", .target = "target", .mor = "weight", minsize = 5))
nul <- rbindlist(lapply(1:200, function(i) { m <- tmat; rownames(m) <- sample(rownames(m))
  as.data.table(run_mlm(m, net_p, .source = "source", .target = "target", .mor = "weight", minsize = 5))[, .(source, condition, score)] }))
foot <- merge(obs[, .(source, condition, score)], nul[, .(null_mean = mean(score), null_sd = sd(score)), by = .(source, condition)], by = c("source", "condition"))
foot <- merge(foot, merge(obs, nul, by = c("source", "condition"), allow.cartesian = TRUE)[, .(emp_p = (1 + sum(abs(score.y) >= abs(score.x))) / (1 + .N)), by = .(source, condition)], by = c("source", "condition"))
foot[, `:=`(score = round(score, 2), null_mean = round(null_mean, 2), null_sd = round(null_sd, 2))]
fwrite(foot[order(condition, score)], out("Table_S18_PROGENy_permutation_null.csv")); print(foot[source %in% c("PI3K", "Hypoxia", "EGFR", "Androgen", "MAPK")][order(condition, score)])
# which genes drive the EGFR footprint score?
eg <- data.table(net_p)[source == "EGFR" & target %in% C1$gene]; eg[, t := C1[target, "t"]][, contrib := weight * t]
say("Top EGFR-footprint contributors (PIK-75): %s", paste(head(eg[order(-contrib), sprintf("%s(%+.1f)", target, contrib)], 12), collapse = ", "))
data(dorothea_hs, package = "dorothea"); net_t <- subset(dorothea_hs, confidence %in% c("A", "B", "C"))
ta <- as.data.table(run_ulm(tmat, net_t, .source = "tf", .target = "target", .mor = "mor", minsize = 10))
say("DoRothEA, PIK-75: %d TFs scored; %.0f%% have negative scores; median score %.2f", ta[condition == "PIK75", .N], 100 * ta[condition == "PIK75", mean(score < 0)], ta[condition == "PIK75", median(score)])
reg <- lapply(c("SREBF1", "SREBF2", "HIF1A", "ARNT", "FOXO3", "FOXO1", "E2F4", "ETS1"), function(x) intersect(net_t$target[net_t$tf == x & net_t$mor > 0], C1$gene))
names(reg) <- c("SREBF1", "SREBF2", "HIF1A", "ARNT", "FOXO3", "FOXO1", "E2F4", "ETS1")
jac <- function(a, b) round(length(intersect(a, b)) / length(union(a, b)), 2)
say("Regulon sizes: %s", paste(names(reg), lengths(reg), collapse = "  "))
say("Jaccard: SREBF1-SREBF2 %.2f | HIF1A-ARNT %.2f | SREBF1-HIF1A %.2f | SREBF2-ARNT %.2f", jac(reg$SREBF1, reg$SREBF2), jac(reg$HIF1A, reg$ARNT), jac(reg$SREBF1, reg$HIF1A), jac(reg$SREBF2, reg$ARNT))
say("ETS1 regulon contains BRCA1: %s, BRCA2: %s", "BRCA1" %in% net_t$target[net_t$tf == "ETS1"], "BRCA2" %in% net_t$target[net_t$tf == "ETS1"])
print(ta[source %in% c("FOXO1", "FOXO3", "FOXO4", "ATF4", "SREBF1", "SREBF2", "HIF1A", "E2F4", "ETS1", "MYC")][order(condition, source), .(source, condition, score = round(score, 2), p_value = signif(p_value, 2))])
cat("Done.\n")
