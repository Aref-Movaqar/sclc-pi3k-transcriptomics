#!/usr/bin/env Rscript
# =============================================================================
# Part 2: orthogonal in-silico context for the GSE40564 findings
#   (a) DepMap CRISPR (Chronos) dependency on PI3K class IA isoforms in SCLC lines
#   (b) UCologne SCLC primary tumours (George et al., Nature 2015; cBioPortal
#       study sclc_ucologne_2015): isoform expression, NE subtype markers, and a
#       score for the PIK-75-suppressed programme
# Run AFTER GSE40564_pipeline.R (uses results/Table_S1_DE_PIK75_vs_DMSO.csv)
#
# DepMap input: DepMap Public 24Q4 (figshare article 27993248)
#   CRISPRGeneEffect.csv  https://ndownloader.figshare.com/files/51064667
#   Model.csv             https://ndownloader.figshare.com/files/51065297
#   -> place in ./depmap/  (update to the newest public release before submission)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(jsonlite); library(httr); library(survival)
})
set.seed(20260930)
out <- function(f) file.path("results", f)

# ---- (a) DepMap -------------------------------------------------------------
genes <- c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "AKT1", "MTOR", "RPTOR",
           "CDK9", "PRKDC", "BRCA1")
hdr   <- names(fread("depmap/CRISPRGeneEffect.csv", nrows = 0))
cols  <- c(hdr[1], hdr[sub(" \\(.*", "", hdr) %in% genes])
ge    <- fread("depmap/CRISPRGeneEffect.csv", select = cols)
setnames(ge, c("ModelID", sub(" \\(.*", "", cols[-1])))
mdl   <- fread("depmap/Model.csv")[, .(ModelID, CellLineName, OncotreeCode)]
ge    <- merge(ge, mdl, by = "ModelID")
ge[, group := fifelse(OncotreeCode == "SCLC", "SCLC",
              fifelse(OncotreeCode %in% c("LUAD", "LUSC", "NSCLC"), "NSCLC", "Other"))]
long  <- melt(ge, id.vars = c("ModelID", "CellLineName", "OncotreeCode", "group"),
              variable.name = "gene", value.name = "chronos")
long  <- long[!is.na(chronos)]

dep_sum <- long[, .(n = .N, median_chronos = round(median(chronos), 3),
                    frac_dependent = round(mean(chronos < -0.5), 3)),
                by = .(gene, group)][order(gene, group)]
fwrite(dep_sum, out("Table_S8_DepMap_summary.csv")); print(dep_sum[gene %in% c("PIK3CA", "PIK3CB", "CDK9", "PRKDC")])

sclc <- dcast(long[group == "SCLC"], ModelID + CellLineName ~ gene, value.var = "chronos")
wt   <- wilcox.test(sclc$PIK3CA, sclc$PIK3CB, paired = TRUE)
cat(sprintf("SCLC lines with CRISPR data: %d; paired Wilcoxon PIK3CA vs PIK3CB p = %.2g\n",
            nrow(sclc), wt$p.value))
cat(sprintf("Genes present in CRISPR matrix: %s\n", paste(intersect(genes, names(ge)), collapse = ", ")))
cat(sprintf("NCI-H69 screened: %s\n", "NCI-H69" %in% sclc$CellLineName))
print(sclc[CellLineName == "NCI-H69"])
# SCLC vs NSCLC for each isoform
cmp <- rbindlist(lapply(c("PIK3CA", "PIK3CB"), function(g) {
  x <- long[gene == g]
  data.table(gene = g, p_SCLC_vs_NSCLC =
               signif(wilcox.test(chronos ~ group, data = x[group != "Other"])$p.value, 2))
}))
print(cmp)
fwrite(data.table(test = c("paired Wilcoxon PIK3CA vs PIK3CB (SCLC)", cmp$gene),
                  p = c(signif(wt$p.value, 2), cmp$p_SCLC_vs_NSCLC)),
       out("Table_S8b_DepMap_tests.csv"))

cols3 <- c(SCLC = "#c0392b", NSCLC = "#2e86c1", Other = "#bfbfbf")
p5a <- ggplot(long[gene %in% c("PIK3CA", "PIK3CB", "PIK3CD")],
              aes(group, chronos, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = -0.5, linetype = 2, linewidth = 0.3) +
  geom_violin(scale = "width", colour = NA, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  facet_wrap(~ gene) + scale_fill_manual(values = cols3, guide = "none") +
  labs(x = NULL, y = "CRISPR gene effect (Chronos)",
       title = "A  DepMap 24Q4: PI3K class IA isoform dependency",
       subtitle = "dashed line = -0.5 (commonly used dependency threshold); NCI-H69 was not screened") +
  coord_cartesian(ylim = c(-1.6, 0.4)) +
  theme_bw(base_size = 10)

# ---- (b) UCologne SCLC (cBioPortal) ------------------------------------------
api  <- "https://www.cbioportal.org/api"
post <- function(path, body) fromJSON(content(POST(paste0(api, path), body = toJSON(body, auto_unbox = TRUE),
                                                   content_type_json(), accept_json()), "text", encoding = "UTF-8"))
de   <- fread(out("Table_S1_DE_PIK75_vs_DMSO.csv"))
sig_dn <- de[adj.P.Val < 0.05 & logFC <= -1, gene]     # PIK-75-suppressed genes
sig_up <- de[adj.P.Val < 0.05 & logFC >=  1, gene]
qgenes <- unique(c("PIK3CA", "PIK3CB", "PIK3CD", "ASCL1", "NEUROD1", "POU2F3", "YAP1", sig_dn, sig_up))
map  <- as.data.table(post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", qgenes))[, .(entrezGeneId, hugoGeneSymbol)]
expr <- as.data.table(post("/molecular-profiles/sclc_ucologne_2015_rna_seq_mrna/molecular-data/fetch",
                           list(sampleListId = "sclc_ucologne_2015_all", entrezGeneIds = map$entrezGeneId)))
expr <- merge(expr[, .(sampleId, patientId, entrezGeneId, value)], map, by = "entrezGeneId")
X    <- dcast(expr, sampleId + patientId ~ hugoGeneSymbol, value.var = "value")
num  <- setdiff(names(X), c("sampleId", "patientId"))
X[, (num) := lapply(.SD, function(v) log2(as.numeric(v) + 1)), .SDcols = num]
cat(sprintf("UCologne samples with RNA-seq: %d\n", nrow(X)))

zs  <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
g_dn <- intersect(sig_dn, num); g_up <- intersect(sig_up, num)
Z   <- X[, lapply(.SD, zs), .SDcols = num]
X[, p110a_programme := rowMeans(Z[, g_dn, with = FALSE], na.rm = TRUE) - rowMeans(Z[, g_up, with = FALSE], na.rm = TRUE)]
cat(sprintf("Signature genes found: %d down / %d up\n", length(g_dn), length(g_up)))

# NE subtype by dominant marker (Rudin et al., Nat Rev Cancer 2019 nomenclature)
mk <- c(A = "ASCL1", N = "NEUROD1", P = "POU2F3", Y = "YAP1")
X[, subtype := paste0("SCLC-", names(mk)[apply(as.matrix(Z[, mk, with = FALSE]), 1, which.max)])]

# Clinical
clin <- as.data.table(fromJSON(content(GET(paste0(api,
         "/studies/sclc_ucologne_2015/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY")),
         "text", encoding = "UTF-8")))
clin <- dcast(clin[clinicalAttributeId %in% c("OS_MONTHS", "OS_STATUS", "UICC_TUMOR_STAGE", "AGE", "SEX")],
              patientId ~ clinicalAttributeId, value.var = "value")
D <- merge(X, clin, by = "patientId", all.x = TRUE)
D[, OS_MONTHS := as.numeric(OS_MONTHS)]
D[, event := as.integer(grepl("^1|DECEASED", OS_STATUS))]

cors <- D[, .(r_PIK3CA_prog = cor(PIK3CA, p110a_programme, method = "spearman"),
              r_PIK3CB_prog = cor(PIK3CB, p110a_programme, method = "spearman"),
              r_PIK3CA_PIK3CB = cor(PIK3CA, PIK3CB, method = "spearman"),
              median_PIK3CA = median(PIK3CA), median_PIK3CB = median(PIK3CB))]
print(round(cors, 3))
iso_p <- wilcox.test(D$PIK3CA, D$PIK3CB, paired = TRUE)$p.value
kw    <- kruskal.test(p110a_programme ~ subtype, data = D)$p.value
cox   <- coxph(Surv(OS_MONTHS, event) ~ scale(p110a_programme), data = D[!is.na(OS_MONTHS)])
cs    <- summary(cox)
# Stage-adjusted sensitivity model (UICC stage collapsed to I-II vs III-IV)
D[, stage2 := fifelse(grepl("^(I|II)[AaBb]?$", UICC_TUMOR_STAGE), "I-II",
              fifelse(grepl("^(III|IV)", UICC_TUMOR_STAGE), "III-IV", NA_character_))]
print(D[, .N, by = .(UICC_TUMOR_STAGE, stage2)])
cox_adj <- coxph(Surv(OS_MONTHS, event) ~ scale(p110a_programme) + stage2,
                 data = D[!is.na(OS_MONTHS) & !is.na(stage2)])
ca <- summary(cox_adj)
fwrite(data.table(term = rownames(ca$conf.int), HR = round(ca$conf.int[, 1], 2),
                  lo = round(ca$conf.int[, 3], 2), hi = round(ca$conf.int[, 4], 2),
                  p = signif(ca$coefficients[, 5], 2), n = ca$n, events = ca$nevent),
       out("Table_S9b_UCologne_cox_stage_adjusted.csv"))
print(D[, .(median_score = round(median(p110a_programme), 3), .N), by = subtype][order(-median_score)])
res_b <- data.table(n_samples = nrow(D), n_OS = cs$n, n_events = cs$nevent,
                    HR_per_SD = round(cs$conf.int[1, 1], 2),
                    HR_lo = round(cs$conf.int[1, 3], 2), HR_hi = round(cs$conf.int[1, 4], 2),
                    cox_p = signif(cs$coefficients[1, 5], 2),
                    kruskal_subtype_p = signif(kw, 2), paired_PIK3CA_vs_PIK3CB_p = signif(iso_p, 2))
res_b <- cbind(res_b, round(cors, 3))
fwrite(res_b, out("Table_S9_UCologne_summary.csv")); print(res_b)
print(D[, .N, by = subtype])
fwrite(D[, .(sampleId, subtype, PIK3CA, PIK3CB, p110a_programme, OS_MONTHS, event)],
       out("Table_S10_UCologne_per_sample.csv"))

p5b <- ggplot(melt(D[, .(sampleId, PIK3CA, PIK3CB, PIK3CD)], id.vars = "sampleId"),
              aes(variable, value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  labs(x = NULL, y = "log2(expression + 1)", title = "B  Isoform mRNA, UCologne SCLC tumours") +
  theme_bw(base_size = 10)
p5c <- ggplot(D, aes(subtype, p110a_programme)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  labs(x = NULL, y = "PIK-75-suppressed programme score",
       title = "C  Programme score by NE subtype", subtitle = sprintf("Kruskal-Wallis p = %.2g", kw)) +
  theme_bw(base_size = 10)
fig5 <- p5a / (p5b | p5c) + plot_layout(heights = c(1, 1))
ggsave(out("Figure5_DepMap_UCologne.pdf"), fig5, width = 10, height = 8)
ggsave(out("Figure5_DepMap_UCologne.png"), fig5, width = 10, height = 8, dpi = 300)
cat("Done.\n")
