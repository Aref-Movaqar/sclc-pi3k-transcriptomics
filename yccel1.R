suppressPackageStartupMessages({library(limma); library(msigdbr); library(data.table)})
x <- as.matrix(read.csv("GSE271761_counts.csv", row.names = 1, check.names = FALSE))
grp <- factor(c("DMSO","DMSO","ALP","ALP","sgC","sgC","sgP","sgP"), levels = c("DMSO","ALP","sgC","sgP"))
cpm <- t(t(x) / colSums(x)) * 1e6; d <- x[rowSums(cpm >= 1) >= 2, ]
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
v <- voom(d, des); fit <- lmFit(v, des); cm <- makeContrasts(ALP = ALP - DMSO, KO = sgP - sgC, levels = des)
f2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE)
cat("genes:", nrow(v), "\n"); for (k in c("ALP","KO")) { t <- topTable(f2, coef = k, number = Inf); cat(k, ": FDR<0.05 =", sum(t$adj.P.Val < 0.05), "\n") }
g <- c("PIK3CA","SQLE","HMGCS1","HMGCR","LDLR","SCD","LSS","FDFT1","INSIG1","ACACA","FASN","STARD4","SESN3","DDIT4","TXNIP","MCL1")
tab <- data.table(gene = g)[, c("ALP_log2FC","ALP_FDR","KO_log2FC","KO_FDR") := {
  a <- topTable(f2, coef = "ALP", number = Inf)[g, ]; b <- topTable(f2, coef = "KO", number = Inf)[g, ]
  list(round(a$logFC, 2), signif(a$adj.P.Val, 2), round(b$logFC, 2), signif(b$adj.P.Val, 2)) }]
print(tab)
gm <- list(); for (l in readLines("ext/Reactome_2022.gmt", warn = FALSE)) { f <- strsplit(l, "\t")[[1]]; gm[[f[1]]] <- f[-(1:2)][f[-(1:2)] != ""] }
H <- msigdbr(species = "Homo sapiens", collection = "H"); Hs <- split(H$gene_symbol, H$gs_name)
S <- list(SREBF_R = gm[[grep("R-HSA-2426168", names(gm))]], chol_biosynthesis_R = gm[[grep("Cholesterol Biosynthesis R-HSA-191273", names(gm), fixed = TRUE)]],
          chol_homeostasis = Hs$HALLMARK_CHOLESTEROL_HOMEOSTASIS, mTORC1 = Hs$HALLMARK_MTORC1_SIGNALING, hypoxia = Hs$HALLMARK_HYPOXIA, glycolysis = Hs$HALLMARK_GLYCOLYSIS, PI3K_AKT_MTOR = Hs$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
idx <- ids2indices(S, rownames(v))
for (k in c("ALP","KO")) { r <- camera(v, idx, des, contrast = cm[, k], inter.gene.cor = NA, sort = FALSE); tt <- topTable(f2, coef = k, number = Inf)
  print(data.table(contrast = k, set = rownames(r), n = r$NGenes, mean_log2FC = round(sapply(S, function(s) mean(tt[intersect(s, rownames(tt)), "logFC"])), 3), dir = r$Direction, p = signif(r$PValue, 2), FDR = signif(r$FDR, 2))) }
# H69 PIK-75 two-fold-down SREBP genes: same direction?
out <- rbindlist(lapply(c("ALP","KO"), function(k) { r <- camera(v, idx, des, contrast = cm[, k], inter.gene.cor = NA, sort = FALSE); tt <- topTable(f2, coef = k, number = Inf)
  data.table(contrast = ifelse(k == "ALP", "alpelisib 5 uM 6 h vs DMSO", "sgPIK3CA vs sgControl"), set = rownames(r), n = r$NGenes, mean_log2FC = round(sapply(S, function(s) mean(tt[intersect(s, rownames(tt)), "logFC"])), 3), direction = r$Direction, camera_p = signif(r$PValue, 2)) }))
fwrite(out, "results_rev/Table_S29b_GSE271761_camera.csv"); fwrite(tab, "results_rev/Table_S29a_GSE271761_genes.csv")
