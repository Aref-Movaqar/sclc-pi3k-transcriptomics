#!/usr/bin/env Rscript
# Reviewer checks (Computational Biology and Chemistry, first revision). Reviewer 2 points 1-6, 9, 10.
suppressPackageStartupMessages({ library(GEOquery); library(limma); library(msigdbr); library(data.table)
  library(survival); library(httr); library(jsonlite) })
set.seed(20261028)
say <- function(...) cat(sprintf(...), "\n")
out <- function(f) file.path("results_rev", f); dir.create("results_rev", showWarnings = FALSE)

gse <- getGEO(filename = "GSE40564_series_matrix.txt.gz", getGPL = FALSE)
ann <- fread("GPL570.annot.gz", skip = "ID\tGene title", fill = TRUE, sep = "\t", quote = "")
symOf <- setNames(ann[["Gene symbol"]], ann$ID)
grp <- factor(rep(c("DMSO", "PIK75", "TGX221"), each = 3), levels = c("DMSO", "PIK75", "TGX221"))
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
cm <- makeContrasts(C1 = PIK75 - DMSO, C2 = TGX221 - DMSO, C3 = TGX221 - PIK75, levels = design)
collapse <- function(e) {
  s <- symOf[rownames(e)]; k <- !grepl("^AFFX", rownames(e)) & !is.na(s) & s != "" & !grepl("///", s, fixed = TRUE)
  e <- e[k, ]; e <- e[rowSums(e > quantile(e, 0.25)) >= 3, ]; s <- symOf[rownames(e)]
  o <- order(s, -rowMeans(e)); f <- !duplicated(s[o]); g <- e[o[f], ]; rownames(g) <- s[o[f]]; g }
fitDE <- function(g, des = design, con = cm) eBayes(contrasts.fit(lmFit(g, des), con), trend = TRUE, robust = TRUE)
tt <- function(fit, k) { x <- topTable(fit, coef = k, number = Inf, sort.by = "none", confint = 0.90); x$gene <- rownames(x); as.data.table(x) }

sub_e <- exprs(gse); g <- collapse(sub_e); stopifnot(nrow(g) == 17458)
fit <- fitDE(g); C1 <- tt(fit, "C1"); C2 <- tt(fit, "C2"); C3 <- tt(fit, "C3")
cp <- readRDS("cel/cel_processed.rds"); nonq <- cp$nonq[, colnames(sub_e)]
spk <- grep("^AFFX-(r2-)?(Ec-)?(Bio[BCD]|bio[BCD]|P1-cre|Cre)|^AFFX-(r2-Bs-)?(lys|phe|thr|dap|Lys|Phe|Thr|Dap)", rownames(nonq), value = TRUE)
anch <- sweep(nonq, 2, colMeans(nonq[spk, ]) - mean(nonq[spk, ])); ga <- collapse(anch)
fitA <- fitDE(ga); A1 <- tt(fitA, "C1"); A2 <- tt(fitA, "C2")

# ---- gene sets --------------------------------------------------------------
H <- msigdbr(species = "Homo sapiens", collection = "H"); Hs <- split(H$gene_symbol, H$gs_name)
gm <- list(); for (l in readLines("ext/Reactome_2022.gmt")) { f <- strsplit(l, "\t")[[1]]; gm[[f[1]]] <- f[-(1:2)][f[-(1:2)] != ""] }
rx <- function(p) gm[[grep(p, names(gm), fixed = TRUE)[1]]]
S <- list(chol_biosynthesis_R = rx("Cholesterol Biosynthesis R-HSA-191273"), SREBF_R = rx("R-HSA-2426168"),
          chol_homeostasis = Hs$HALLMARK_CHOLESTEROL_HOMEOSTASIS, mTORC1 = Hs$HALLMARK_MTORC1_SIGNALING,
          hypoxia = Hs$HALLMARK_HYPOXIA, glycolysis = Hs$HALLMARK_GLYCOLYSIS, E2F = Hs$HALLMARK_E2F_TARGETS,
          G2M = Hs$HALLMARK_G2M_CHECKPOINT, HDR_R = rx("R-HSA-5693538"), DNA_repair = Hs$HALLMARK_DNA_REPAIR,
          UPR = Hs$HALLMARK_UNFOLDED_PROTEIN_RESPONSE, PI3K_AKT_MTOR = Hs$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
S <- lapply(S, function(x) intersect(unique(x), rownames(g)))

# ---- R2.1 / R2.2: direct C3 test, correlation-aware camera, overlap-free sets --
core <- c("chol_homeostasis", "mTORC1", "hypoxia", "glycolysis")
shared <- unique(unlist(lapply(core, function(a) unlist(lapply(setdiff(core, a), function(b) intersect(S[[a]], S[[b]]))))))
srebp_all <- unique(c(S$SREBF_R, S$chol_biosynthesis_R, S$chol_homeostasis))
S_u <- lapply(S[core], function(x) setdiff(x, unique(c(shared, srebp_all)))); names(S_u) <- paste0(core, "_unique")
ov <- outer(core, core, Vectorize(function(a, b) length(intersect(S[[a]], S[[b]])))); dimnames(ov) <- list(core, core)
print(ov); say("genes shared among the four core Hallmark sets: %d", length(shared))
say("mTORC1 genes that are also SREBP/cholesterol genes: %d of %d", length(intersect(S$mTORC1, srebp_all)), length(S$mTORC1))
allS <- c(S, S_u)
cam <- function(k) { r <- camera(g, allS, design, contrast = cm[, k], inter.gene.cor = NA, trend.var = TRUE, sort = FALSE)
  data.table(set = rownames(r), contrast = k, n = r$NGenes, corr = round(r$Correlation, 3), direction = r$Direction, p = signif(r$PValue, 2), FDR = signif(r$FDR, 2)) }
camr <- rbindlist(lapply(c("C1", "C2", "C3"), cam))
mean_fc <- function(D) sapply(allS, function(x) round(mean(D[gene %in% x, logFC]), 3))
eff <- data.table(set = names(allS), n = lengths(allS), C1_mean = mean_fc(C1), C2_mean = mean_fc(C2), C3_mean = mean_fc(C3),
                  ratio_C1_over_C2 = round(mean_fc(C1) / mean_fc(C2), 1))
camw <- dcast(camr, set ~ contrast, value.var = c("p", "FDR", "corr"))
res1 <- merge(eff, camw, by = "set", sort = FALSE); print(res1); fwrite(res1, out("R2_1_2_camera_C3_overlap.csv"))

# ---- R2.3: equivalence (90% CI within +/-0.32 log2, i.e. 1.25-fold) for genes called unchanged; spike-in values for modest changes
G3 <- c("SREBF2", "SREBF1", "E2F1", "MKI67", "PCNA", "TOP2A", "DUSP6", "ETV4", "BRCA1", "BRCA2", "RAD51", "CDKN1B", "SESN3", "GADD45A", "CITED2", "TRIB3", "DDIT4", "INSIG1", "FASN", "ACACA", "MCL1")
eq <- merge(C1[gene %in% G3, .(gene, logFC, CI.L, CI.R, FDR = signif(adj.P.Val, 2))], A1[gene %in% G3, .(gene, anch_logFC = logFC, anch_CI.L = CI.L, anch_CI.R = CI.R, anch_FDR = signif(adj.P.Val, 2))], by = "gene")
eq[, equivalent_main := CI.L > -0.32 & CI.R < 0.32][, equivalent_anch := anch_CI.L > -0.32 & anch_CI.R < 0.32]
tr <- topTreat(treat(fit, lfc = 0.32, trend = TRUE, robust = TRUE), coef = "C1", number = Inf, sort.by = "none"); tr$gene <- rownames(tr)
eq <- merge(eq, as.data.table(tr)[, .(gene, treat_1.25fold_FDR = signif(adj.P.Val, 2))], by = "gene")
num <- names(eq)[sapply(eq, is.numeric)]; eq[, (num) := lapply(.SD, round, 3), .SDcols = num]
print(eq[match(G3, gene)]); fwrite(eq[match(G3, gene)], out("R2_3_equivalence_spikein.csv"))

# ---- R2.4: rank of SREBP targets; composition of Figure 1 genes ----------------
C1[, rank_logFC := frank(logFC)][, rank_t := frank(t)]
t2 <- c("SQLE", "HMGCS1", "HMGCR", "LDLR", "SCD", "LSS", "FDFT1", "ACACA", "INSIG1", "FASN", "STARD4", "MVD", "IDI1", "FDPS", "DHCR7", "SC5D", "CYP51A1", "MSMO1", "NSDHL", "HSD17B7", "TMEM97", "LDLRAP1")
rk <- C1[gene %in% t2, .(gene, logFC = round(logFC, 2), rank_by_logFC = rank_logFC, rank_by_t = rank_t, of = nrow(C1))][order(rank_by_logFC)]
print(rk); fwrite(rk, out("R2_4_SREBP_ranks.csv"))
trd <- topTreat(treat(fit, lfc = 1, trend = TRUE, robust = TRUE), coef = "C1", number = Inf, sort.by = "none"); trd$gene <- rownames(trd)
down316 <- trd$gene[trd$adj.P.Val < 0.05 & trd$logFC < 0]; say("treat 2-fold down genes: %d", length(down316))
for (nm in c("SREBF_R", "chol_biosynthesis_R", "chol_homeostasis", "mTORC1", "hypoxia", "glycolysis")) {
  x <- S[[nm]]; k <- sum(down316 %in% x); e <- length(down316) * length(x) / nrow(g)
  say("%s: %d of %d set genes among the %d two-fold-down genes (expected %.1f; hypergeometric p = %.2g)", nm, k, length(x), length(down316), e,
      phyper(k - 1, length(x), nrow(g) - length(x), length(down316), lower.tail = FALSE)) }
# Figure 1 genes
ft <- topTable(fit, coef = 1:2, number = Inf)
ft <- ft[ft$adj.P.Val < 0.05 & apply(abs(ft[, c("C1", "C2")]), 1, max) >= 1, ]; top60 <- head(rownames(ft), 60)
say("Figure 1 top-60: SREBP/cholesterol genes: %s", paste(intersect(top60, srebp_all), collapse = ", "))
GOn <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
neur <- unique(GOn$gene_symbol[GOn$gs_name %in% c("GOBP_NEUROTRANSMITTER_SECRETION", "GOBP_SYNAPTIC_SIGNALING", "GOBP_NEURON_PROJECTION_DEVELOPMENT", "GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS")])
say("Figure 1 top-60 in neuronal GO sets (synaptic signalling / neurotransmitter / neuron projection): %d: %s", sum(top60 %in% neur), paste(intersect(top60, neur), collapse = ", "))
k <- sum(down316 %in% neur); say("two-fold-down genes in neuronal GO sets: %d of %d (expected %.1f; p = %.2g)", k, length(down316),
  length(down316) * mean(rownames(g) %in% neur), phyper(k - 1, sum(rownames(g) %in% neur), nrow(g) - sum(rownames(g) %in% neur), length(down316), lower.tail = FALSE))
cam_neur <- camera(g, list(neuronal = intersect(neur, rownames(g))), design, contrast = cm[, "C1"], inter.gene.cor = NA, trend.var = TRUE)
say("neuronal GO union (n = %d): mean log2FC PIK-75 %.3f, TGX-221 %.3f; camera C1 p = %.2g (%s)", sum(rownames(g) %in% neur),
    mean(C1[gene %in% neur, logFC]), mean(C2[gene %in% neur, logFC]), cam_neur$PValue, cam_neur$Direction)

# ---- R2.5: leave out PIK75_5 (extra 5' signal loss) ---------------------------
keep <- colnames(g)[-6]; stopifnot(grepl("GSM996648", colnames(g)[6]))
d5 <- model.matrix(~ 0 + droplevels(grp[-6])); colnames(d5) <- levels(grp)
f5 <- eBayes(contrasts.fit(lmFit(g[, keep], d5), makeContrasts(C1 = PIK75 - DMSO, levels = d5)), trend = TRUE, robust = TRUE)
L5 <- as.data.table(topTable(f5, coef = 1, number = Inf, sort.by = "none"), keep.rownames = "gene")
say("Without PIK75_5: FDR<0.05 %d; 2-fold down %d / up %d; r(log2FC with all arrays) = %.3f", sum(L5$adj.P.Val < 0.05),
    sum(L5$adj.P.Val < 0.05 & L5$logFC <= -1), sum(L5$adj.P.Val < 0.05 & L5$logFC >= 1), cor(L5$logFC, C1$logFC[match(L5$gene, C1$gene)]))
print(L5[gene %in% t2[1:10], .(gene, logFC = round(logFC, 2), FDR = signif(adj.P.Val, 2))])
say("Set means without PIK75_5: %s", paste(sapply(c("chol_biosynthesis_R", "chol_homeostasis", "mTORC1", "hypoxia", "glycolysis"), function(s) sprintf("%s %.2f", s, mean(L5[gene %in% S[[s]], logFC]))), collapse = "; "))
# per-array 3'/5' ratios for all AFFX housekeeping probe pairs
hk <- grep("^AFFX-HUMGAPDH/M33197_[53M]_at$|^AFFX-HSAC07/X00351_[53M]_at$", rownames(nonq), value = TRUE)
print(round(nonq[hk, ], 2))
qc <- as.data.table(cp$qc); print(qc)

# ---- R2.6: half-life test under spike-in anchoring; MCL1 probes ---------------
hl <- fread("ext/human_halflife_PC1.csv")[!duplicated(gene)]
hfun <- function(D, lab) { d <- merge(D[, .(gene, logFC, adj.P.Val)], hl, by = "gene"); d[, pct := frank(halflife_PC1) / .N * 100]
  dn <- d$logFC <= -1 & d$adj.P.Val < 0.05; s5 <- d$pct <= 5
  say("%s: n=%d rho=%.3f; shortest 5%%: %%down2fold %.1f vs %.1f; mean log2FC %.3f vs %.3f; %%up2fold %.2f vs %.2f; MW p=%.2g", lab, nrow(d),
      cor(d$logFC, d$halflife_PC1, method = "spearman"), 100 * mean(dn[s5]), 100 * mean(dn[!s5]), mean(d$logFC[s5]), mean(d$logFC[!s5]),
      100 * mean((d$logFC >= 1 & d$adj.P.Val < 0.05)[s5]), 100 * mean((d$logFC >= 1 & d$adj.P.Val < 0.05)[!s5]), wilcox.test(d$logFC[s5], d$logFC[!s5])$p.value)
  d }
dm <- hfun(C1, "main"); da <- hfun(A1, "spike-anchored")
s5 <- dm[pct <= 5][order(-logFC)]; say("shortest-5%% genes most increased: %s", paste(head(sprintf("%s %+.2f", s5$gene, s5$logFC), 15), collapse = ", "))
ieg <- c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "ATF3", "DUSP1", "DUSP6", "MYC", "KLF6", "NR4A1", "ZFP36")
print(dm[gene %in% ieg, .(gene, logFC = round(logFC, 2), pct = round(pct, 1))][order(pct)])
mp <- grep("^MCL1$", symOf[rownames(sub_e)]); mpr <- rownames(sub_e)[mp]
pf <- eBayes(contrasts.fit(lmFit(sub_e[mpr, , drop = FALSE], design), cm[, 1:2]), trend = FALSE)
pa <- eBayes(contrasts.fit(lmFit(anch[mpr, , drop = FALSE], design), cm[, 1:2]), trend = FALSE)
print(data.table(probe = mpr, mean_intensity = round(rowMeans(sub_e[mpr, ]), 2), PIK75_log2FC = round(pf$coefficients[, 1], 2), PIK75_anch = round(pa$coefficients[, 1], 2), TGX221_log2FC = round(pf$coefficients[, 2], 2)))
fwrite(rbind(dm[, .(gene, logFC, pct, norm = "main")], da[, .(gene, logFC, pct, norm = "spike")]), out("R2_6_halflife_both.csv"))

# ---- R2.9 / R2.10: cohort nulls matched on expression level and proliferation; NE genes ----
cc <- readRDS("ext/ucologne_expr.rds"); X <- cc$X; X <- X[apply(X, 1, sd, na.rm = TRUE) > 0 & rowMeans(is.na(X)) == 0, ]
Z <- t(scale(t(X)))
dn <- intersect(C1[adj.P.Val < 0.05 & logFC <= -1, gene], rownames(Z)); up <- intersect(C1[adj.P.Val < 0.05 & logFC >= 1, gene], rownames(Z))
score <- function(d, u) colMeans(Z[d, , drop = FALSE]) - if (length(u)) colMeans(Z[u, , drop = FALSE]) else 0
NEpos <- c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3","SYT4","RTN1","MYT1L","SYP","KIF1A","TMSB15A","SYN1","SYT11","RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6","TMSB15B","CELF3")
NEneg <- c("RAB27B","TGFBR2","SLC16A5","S100A10","ITGB4","YAP1","LGALS3","EPHA2","S100A16","PLAU","ABCC3","ARHGDIB","CYR61","PTGES","CCND1","IFITM2","IFITM3","AHNAK","CAV2","TACSTD2","TGFBI","EMP1","CAV1","ANXA1","MYOF")
NE <- colMeans(Z[intersect(NEpos, rownames(Z)), ]) - colMeans(Z[intersect(NEneg, rownames(Z)), ])
prolif_genes <- intersect(unique(c(Hs$HALLMARK_E2F_TARGETS, Hs$HALLMARK_G2M_CHECKPOINT)), rownames(Z))
PROL <- colMeans(Z[prolif_genes, ]); say("proliferation meta-gene from %d E2F/G2M genes", length(prolif_genes))
api <- "https://www.cbioportal.org/api"
clin <- as.data.table(fromJSON(content(GET(paste0(api, "/studies/sclc_ucologne_2015/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY")), "text", encoding = "UTF-8")))
clin <- dcast(clin[clinicalAttributeId %in% c("OS_MONTHS", "OS_STATUS", "UICC_TUMOR_STAGE")], patientId ~ clinicalAttributeId, value.var = "value")
P <- merge(data.table(sampleId = colnames(Z)), cc$pat, by = "sampleId"); P <- merge(P, clin, by = "patientId", all.x = TRUE)
P[, OS := as.numeric(OS_MONTHS)][, ev := as.integer(grepl("^1|DECEASED", OS_STATUS))]
P[, stage := fifelse(grepl("^(I|II)[AaBb]?$", UICC_TUMOR_STAGE), "I-II", fifelse(grepl("^(III|IV)", UICC_TUMOR_STAGE), "III-IV", NA_character_))]
S2 <- P[!is.na(OS) & !is.na(stage)]; say("cohort for Cox: %d patients, %d deaths", nrow(S2), sum(S2$ev))
zcox <- function(v) summary(coxph(Surv(OS, ev) ~ scale(v[S2$sampleId]) + stage, S2))$coefficients[1, 4]
pool <- intersect(C1$gene, rownames(Z)); obs <- score(dn, up); z_obs <- zcox(obs)
say("observed: stage-adjusted Cox z = %.2f; rho(score, proliferation) = %.2f; rho(score, NE) = %.2f", z_obs, cor(obs, PROL, method = "spearman"), cor(obs, NE, method = "spearman"))
binby <- function(v) { b <- cut(v, quantile(v, 0:10 / 10), include.lowest = TRUE, labels = FALSE); names(b) <- names(v); b }
mexp <- rowMeans(X[pool, ]); rP <- apply(Z[pool, ], 1, cor, y = PROL, method = "spearman"); rN <- apply(Z[pool, ], 1, cor, y = NE, method = "spearman")
bE <- binby(mexp); bP <- binby(rP); bN <- binby(rN)
bEP <- paste(binby(mexp), cut(rP, quantile(rP, 0:5 / 5), include.lowest = TRUE, labels = FALSE)); names(bEP) <- pool
matcher <- function(b) function(gs) unlist(lapply(split(gs, b[gs]), function(x) { cand <- pool[b[pool] == b[x[1]]]; sample(cand, length(x), replace = length(cand) < length(x)) }))
nullrun <- function(lab, fn, B = 1000) { z <- replicate(B, zcox(score(fn(dn), fn(up))))
  data.table(null = lab, observed_z = round(z_obs, 2), null_mean_z = round(mean(z), 2), null_sd_z = round(sd(z), 2),
             frac_null_nominal_p05 = round(mean(abs(z) > 1.96), 3), empirical_p = round(mean(c(abs(z) >= abs(z_obs), TRUE)), 3)) }
nulls <- rbind(nullrun("size-matched", function(gs) sample(pool, length(gs))),
               nullrun("matched on mean expression (deciles)", matcher(bE)),
               nullrun("matched on correlation with proliferation (deciles)", matcher(bP)),
               nullrun("matched on correlation with NE score (deciles)", matcher(bN)),
               nullrun("matched on expression x proliferation (10 x 5 bins)", matcher(bEP)))
print(nulls); fwrite(nulls, out("R2_9_matched_nulls.csv"))
cx <- summary(coxph(Surv(OS, ev) ~ scale(obs[S2$sampleId]) + stage + scale(PROL[S2$sampleId]), S2))
say("Cox score + stage + proliferation: HR %.2f (%.2f-%.2f), p = %.2g", cx$conf.int[1, 1], cx$conf.int[1, 3], cx$conf.int[1, 4], cx$coefficients[1, 5])
# NE content of the PIK-75 signature and within-subtype correlations
say("PIK-75 down genes in cohort: %d; NE-panel genes among them: %s; non-NE panel among up: %s", length(dn), paste(intersect(dn, NEpos), collapse = ","), paste(intersect(up, NEneg), collapse = ","))
say("down genes positively correlated with NE score (rho >= 0.3): %d of %d; negatively (<= -0.3): %d", sum(rN[dn] >= 0.3), length(dn), sum(rN[dn] <= -0.3))
say("up genes positively correlated with NE score (rho >= 0.3): %d of %d; negatively: %d", sum(rN[up] >= 0.3), length(up), sum(rN[up] <= -0.3))
mk <- c(A = "ASCL1", N = "NEUROD1", P = "POU2F3", Y = "YAP1"); sub <- factor(names(mk)[apply(Z[mk, ], 2, which.max)])
dn2 <- dn[abs(rN[dn]) < 0.3]; up2 <- up[abs(rN[up]) < 0.3]; s2 <- score(dn2, up2)
for (s in levels(sub)) say("SCLC-%s (n=%d): rho(score, NE) = %.2f; NE-decoupled score (%d/%d genes) rho = %.2f", s, sum(sub == s),
  cor(obs[sub == s], NE[sub == s], method = "spearman"), length(dn2), length(up2), cor(s2[sub == s], NE[sub == s], method = "spearman"))
say("rho(score, NE) overall %.2f; decoupled %.2f", cor(obs, NE, method = "spearman"), cor(s2, NE, method = "spearman"))
cat("Done.\n")
