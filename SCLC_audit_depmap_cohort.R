#!/usr/bin/env Rscript
# =============================================================================
# Part 4: audit of the DepMap and UCologne analyses
#   (a) DepMap: effect size with CI (SCLC vs NSCLC), PIK3CA / PTEN genotype of SCLC
#       lines, genotype of NCI-H69
#   (b) UCologne: neuroendocrine (NE) confounding of the PIK-75 programme score,
#       random-signature nulls, NE-decoupled score, PI3K-pathway mutation status
# Extra DepMap 24Q4 files needed in ./depmap/:
#   OmicsSomaticMutationsMatrixHotspot.csv   https://ndownloader.figshare.com/files/51065750
#   OmicsSomaticMutationsMatrixDamaging.csv  https://ndownloader.figshare.com/files/51065747
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(httr); library(survival)})
set.seed(20260930); out <- function(f) file.path("results", f); say <- function(...) cat(sprintf(...), "\n")

# ---- (a) DepMap ----------------------------------------------------------------
rd <- function(f, genes) { h <- names(fread(f, nrows = 0)); cols <- c(h[1], h[sub(" \\(.*", "", h) %in% genes])
  x <- fread(f, select = cols); setnames(x, c("ModelID", sub(" \\(.*", "", cols[-1]))); x }
ge  <- rd("depmap/CRISPRGeneEffect.csv", c("PIK3CA", "PIK3CB"))
hot <- rd("depmap/OmicsSomaticMutationsMatrixHotspot.csv", c("PIK3CA", "PTEN", "AKT1")); setnames(hot, -1, paste0(names(hot)[-1], "_hotspot"))
dam <- rd("depmap/OmicsSomaticMutationsMatrixDamaging.csv", c("PTEN", "PIK3CA", "PIK3R1")); setnames(dam, -1, paste0(names(dam)[-1], "_damaging"))
mdl <- fread("depmap/Model.csv")[, .(ModelID, CellLineName, OncotreeCode)]
D <- Reduce(function(a, b) merge(a, b, by = "ModelID", all.x = TRUE), list(mdl, ge, hot, dam))
D[, group := fifelse(OncotreeCode == "SCLC", "SCLC", fifelse(OncotreeCode %in% c("LUAD", "LUSC", "NSCLC"), "NSCLC", "Other"))]
cat("NCI-H69 genotype in DepMap:\n"); print(D[CellLineName == "NCI-H69"])
S <- D[group == "SCLC" & !is.na(PIK3CA)]
w <- wilcox.test(D[group == "SCLC", PIK3CA], D[group == "NSCLC", PIK3CA], conf.int = TRUE)
say("PIK3CA Chronos, SCLC minus NSCLC: Hodges-Lehmann shift %.3f (95%% CI %.3f to %.3f), p = %.2f", w$estimate, w$conf.int[1], w$conf.int[2], w$p.value)
w2 <- wilcox.test(S$PIK3CA, S$PIK3CB, paired = TRUE, conf.int = TRUE)
say("Within SCLC, PIK3CA minus PIK3CB: paired shift %.3f (95%% CI %.3f to %.3f)", w2$estimate, w2$conf.int[1], w2$conf.int[2])
S[, PI3K_altered := (PIK3CA_hotspot > 0) | (PTEN_damaging > 0)]
gt <- S[, .(n = .N, median_PIK3CA = round(median(PIK3CA), 2), median_PIK3CB = round(median(PIK3CB), 2),
            pct_PIK3CA_dependent = round(100 * mean(PIK3CA < -0.5))), by = .(PIK3CA_hotspot = PIK3CA_hotspot > 0, PTEN_damaging = PTEN_damaging > 0)]
print(gt); fwrite(gt, out("Table_S19_DepMap_SCLC_genotype.csv"))
say("SCLC lines screened: %d; PIK3CA hotspot %d; PTEN damaging %d", nrow(S), sum(S$PIK3CA_hotspot > 0, na.rm = TRUE), sum(S$PTEN_damaging > 0, na.rm = TRUE))
wt <- S[PI3K_altered == FALSE]
say("Excluding PIK3CA-hotspot/PTEN-damaged lines (n = %d): median PIK3CA %.2f, PIK3CB %.2f, paired p = %.2g", nrow(wt), median(wt$PIK3CA), median(wt$PIK3CB), wilcox.test(wt$PIK3CA, wt$PIK3CB, paired = TRUE)$p.value)

# ---- (b) UCologne ---------------------------------------------------------------
api  <- "https://www.cbioportal.org/api"
post <- function(path, body) fromJSON(content(POST(paste0(api, path), body = toJSON(body, auto_unbox = TRUE), content_type_json(), accept_json()), "text", encoding = "UTF-8"))
de   <- fread(out("Table_S1_DE_PIK75_vs_DMSO.csv"))
cache <- "ext/ucologne_expr.rds"
if (!file.exists(cache)) {
  ne_genes <- c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3","SYT4","RTN1","MYT1L","SYP","KIF1A","TMSB15A","SYN1","SYT11","RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6","TMSB15B","CELF3",
                "RAB27B","TGFBR2","SLC16A5","S100A10","ITGB4","YAP1","LGALS3","EPHA2","S100A16","PLAU","ABCC3","ARHGDIB","CYR61","PTGES","CCND1","IFITM2","IFITM3","AHNAK","CAV2","TACSTD2","TGFBI","EMP1","CAV1","ANXA1","MYOF","NEUROD1","POU2F3","PIK3CB","PIK3CD")
  q <- unique(c(de$gene, ne_genes)); chunks <- split(q, ceiling(seq_along(q) / 400))
  map <- rbindlist(lapply(chunks, function(x) as.data.table(post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", x))[, .(entrezGeneId, hugoGeneSymbol)]))
  ex <- rbindlist(lapply(split(map$entrezGeneId, ceiling(seq_along(map$entrezGeneId) / 400)), function(ids)
    as.data.table(post("/molecular-profiles/sclc_ucologne_2015_rna_seq_mrna/molecular-data/fetch",
                       list(sampleListId = "sclc_ucologne_2015_all", entrezGeneIds = ids)))[, .(sampleId, patientId, entrezGeneId, value)]))
  ex <- merge(ex, map, by = "entrezGeneId")
  M <- dcast(ex, hugoGeneSymbol ~ sampleId, value.var = "value", fun.aggregate = mean)
  X <- log2(as.matrix(M[, -1]) + 1); rownames(X) <- M$hugoGeneSymbol
  saveRDS(list(X = X, pat = unique(ex[, .(sampleId, patientId)])), cache)
}
cc <- readRDS(cache); X <- cc$X; X <- X[apply(X, 1, sd, na.rm = TRUE) > 0 & rowMeans(is.na(X)) == 0, ]
Z <- t(scale(t(X))); say("UCologne matrix: %d genes x %d tumours", nrow(Z), ncol(Z))
dn <- intersect(de[adj.P.Val < 0.05 & logFC <= -1, gene], rownames(Z)); up <- intersect(de[adj.P.Val < 0.05 & logFC >= 1, gene], rownames(Z))
score <- function(d, u) colMeans(Z[d, , drop = FALSE]) - if (length(u)) colMeans(Z[u, , drop = FALSE]) else 0
NEpos <- c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3","SYT4","RTN1","MYT1L","SYP","KIF1A","TMSB15A","SYN1","SYT11","RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6","TMSB15B","CELF3")
NEneg <- c("RAB27B","TGFBR2","SLC16A5","S100A10","ITGB4","YAP1","LGALS3","EPHA2","S100A16","PLAU","ABCC3","ARHGDIB","CYR61","PTGES","CCND1","IFITM2","IFITM3","AHNAK","CAV2","TACSTD2","TGFBI","EMP1","CAV1","ANXA1","MYOF")
say("NE genes found: %d of 25 NE, %d of 25 non-NE", sum(NEpos %in% rownames(Z)), sum(NEneg %in% rownames(Z)))
NE <- colMeans(Z[intersect(NEpos, rownames(Z)), ]) - colMeans(Z[intersect(NEneg, rownames(Z)), ])
clin <- as.data.table(fromJSON(content(GET(paste0(api, "/studies/sclc_ucologne_2015/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY")), "text", encoding = "UTF-8")))
clin <- dcast(clin[clinicalAttributeId %in% c("OS_MONTHS", "OS_STATUS", "UICC_TUMOR_STAGE", "AGE", "SEX")], patientId ~ clinicalAttributeId, value.var = "value")
P <- merge(data.table(sampleId = colnames(Z), prog = score(dn, up), NE = NE, PIK3CA = X["PIK3CA", ], PIK3CB = X["PIK3CB", ]), cc$pat, by = "sampleId")
P <- merge(P, clin, by = "patientId", all.x = TRUE)
P[, OS := as.numeric(OS_MONTHS)][, ev := as.integer(grepl("^1|DECEASED", OS_STATUS))][, age := as.numeric(AGE)]
P[, stage := fifelse(grepl("^(I|II)[AaBb]?$", UICC_TUMOR_STAGE), "I-II", fifelse(grepl("^(III|IV)", UICC_TUMOR_STAGE), "III-IV", NA_character_))]
S2 <- P[!is.na(OS) & !is.na(stage)]
ct <- cor.test(P$prog, P$NE, method = "spearman", exact = FALSE)
say("Programme score vs NE score: Spearman rho = %.2f (p = %.1g)", ct$estimate, ct$p.value)
hr <- function(f, lab) { s <- summary(coxph(as.formula(f), S2)); i <- 1
  data.table(model = lab, n = s$n, events = s$nevent, HR = round(s$conf.int[i, 1], 2), lo = round(s$conf.int[i, 3], 2), hi = round(s$conf.int[i, 4], 2), p = signif(s$coefficients[i, 5], 2)) }
cox <- rbind(hr("Surv(OS, ev) ~ scale(prog) + stage", "programme + stage"),
             hr("Surv(OS, ev) ~ scale(prog) + stage + scale(NE)", "programme + stage + NE score"),
             hr("Surv(OS, ev) ~ scale(prog) + stage + scale(NE) + age + SEX", "programme + stage + NE + age + sex"),
             hr("Surv(OS, ev) ~ scale(NE) + stage", "NE score + stage (HR is for NE)"))
# NE-decoupled programme: drop signature genes correlated with NE score
rNE <- apply(Z, 1, function(v) cor(v, NE, method = "spearman"))
dn2 <- dn[abs(rNE[dn]) < 0.3]; up2 <- up[abs(rNE[up]) < 0.3]
S2[, prog_noNE := score(dn2, up2)[sampleId]]
cox <- rbind(cox, hr("Surv(OS, ev) ~ scale(prog_noNE) + stage", sprintf("NE-decoupled programme (%d down / %d up genes) + stage", length(dn2), length(up2))))
say("Signature genes with |rho to NE| >= 0.3: %d of %d", sum(abs(rNE[c(dn, up)]) >= 0.3), length(c(dn, up)))
# Random-signature nulls (stage-adjusted Cox z)
zcox <- function(v) summary(coxph(Surv(OS, ev) ~ scale(v[S2$sampleId]) + stage, S2))$coefficients[1, 4]
z_obs <- zcox(score(dn, up)); pool <- intersect(de$gene, rownames(Z))
z_rand <- replicate(1000, zcox(score(sample(pool, length(dn)), sample(pool, length(up)))))
dec <- cut(rNE[pool], quantile(rNE[pool], 0:10 / 10), include.lowest = TRUE, labels = FALSE); names(dec) <- pool
matched <- function(g) unlist(lapply(split(g, dec[g]), function(x) sample(pool[dec[pool] == dec[x[1]]], length(x))))
z_ne <- replicate(1000, zcox(score(matched(dn), matched(up))))
nulls <- data.table(null = c("random genes, size-matched", "random genes, matched on correlation with NE score"),
                    observed_z = round(z_obs, 2), null_mean_z = round(c(mean(z_rand), mean(z_ne)), 2), null_sd_z = round(c(sd(z_rand), sd(z_ne)), 2),
                    frac_null_nominal_p05 = round(c(mean(abs(z_rand) > 1.96), mean(abs(z_ne) > 1.96)), 3),
                    empirical_p = round(c(mean(c(abs(z_rand) >= abs(z_obs), TRUE)), mean(c(abs(z_ne) >= abs(z_obs), TRUE))), 3))
# Subtype association against the same random-signature null
mk <- c(A = "ASCL1", N = "NEUROD1", P = "POU2F3", Y = "YAP1")
sub <- factor(names(mk)[apply(Z[mk, ], 2, which.max)])
kw <- function(v) kruskal.test(v ~ sub)$statistic
kw_obs <- kw(score(dn, up)); kw_rand <- replicate(1000, kw(score(sample(pool, length(dn)), sample(pool, length(up)))))
nulls <- rbind(nulls, data.table(null = "subtype (Kruskal-Wallis chi2), random genes size-matched", observed_z = round(kw_obs, 2),
                                 null_mean_z = round(mean(kw_rand), 2), null_sd_z = round(sd(kw_rand), 2),
                                 frac_null_nominal_p05 = round(mean(kw_rand > qchisq(0.95, 3)), 3),
                                 empirical_p = round(mean(c(kw_rand >= kw_obs, TRUE)), 3)))
# PIK3CA vs PIK3CB correlation with the score: bootstrap difference
bd <- replicate(5000, { i <- sample(nrow(P), replace = TRUE); cor(P$PIK3CA[i], P$prog[i], method = "spearman") - cor(P$PIK3CB[i], P$prog[i], method = "spearman") })
say("rho(PIK3CA, score) - rho(PIK3CB, score) = %.2f (bootstrap 95%% CI %.2f to %.2f)", cor(P$PIK3CA, P$prog, method = "spearman") - cor(P$PIK3CB, P$prog, method = "spearman"), quantile(bd, 0.025), quantile(bd, 0.975))
# PI3K-pathway mutations
pg <- c("PIK3CA", "PTEN", "AKT1", "AKT2", "AKT3", "PIK3R1", "PIK3CB", "MTOR", "RICTOR", "TSC1", "TSC2")
pm <- as.data.table(post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", pg))
mu <- tryCatch(as.data.table(post("/molecular-profiles/sclc_ucologne_2015_mutations/mutations/fetch?projection=DETAILED",
               list(sampleListId = "sclc_ucologne_2015_all", entrezGeneIds = pm$entrezGeneId))), error = function(e) NULL)
if (!is.null(mu) && nrow(mu)) {
  mu <- merge(mu[, .(sampleId, entrezGeneId, proteinChange, mutationType)], pm[, .(entrezGeneId, hugoGeneSymbol)], by = "entrezGeneId")
  mu <- mu[!mutationType %in% c("Silent")]
  print(mu[, .N, by = hugoGeneSymbol][order(-N)]); print(mu[hugoGeneSymbol %in% c("PIK3CA", "PTEN"), .(sampleId, hugoGeneSymbol, proteinChange, mutationType)])
  P[, core_alt := sampleId %in% mu[hugoGeneSymbol %in% c("PIK3CA", "PTEN", "AKT1", "AKT2", "AKT3", "PIK3R1"), sampleId]]
  say("PIK3CA/PTEN/AKT/PIK3R1-mutated tumours with RNA-seq: %d of %d; median score %.2f vs %.2f; Wilcoxon p = %.2f",
      sum(P$core_alt), nrow(P), median(P[core_alt == TRUE, prog]), median(P[core_alt == FALSE, prog]), wilcox.test(prog ~ core_alt, P)$p.value)
}
fwrite(cox, out("Table_S20_UCologne_cox_NE_adjusted.csv")); fwrite(nulls, out("Table_S21_UCologne_random_signature_null.csv"))
print(cox); print(nulls); cat("Done.\n")
