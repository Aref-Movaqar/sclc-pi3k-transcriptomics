suppressPackageStartupMessages({ library(GEOquery); library(limma); library(msigdbr); library(data.table) })
gse <- getGEO(filename = "GSE40564_series_matrix.txt.gz", getGPL = FALSE)
ann <- fread("GPL570.annot.gz", skip = "ID\tGene title", fill = TRUE, sep = "\t", quote = ""); symOf <- setNames(ann[["Gene symbol"]], ann$ID)
grp <- factor(rep(c("DMSO", "PIK75", "TGX221"), each = 3), levels = c("DMSO", "PIK75", "TGX221"))
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
cm <- makeContrasts(C1 = PIK75 - DMSO, C2 = TGX221 - DMSO, C3 = TGX221 - PIK75, levels = design)
e <- exprs(gse); s <- symOf[rownames(e)]; k <- !grepl("^AFFX", rownames(e)) & !is.na(s) & s != "" & !grepl("///", s, fixed = TRUE)
e <- e[k, ]; e <- e[rowSums(e > quantile(e, 0.25)) >= 3, ]; s <- symOf[rownames(e)]; o <- order(s, -rowMeans(e)); f <- !duplicated(s[o]); g <- e[o[f], ]; rownames(g) <- s[o[f]]
fit <- eBayes(contrasts.fit(lmFit(g, design), cm), trend = TRUE, robust = TRUE)
H <- msigdbr(species = "Homo sapiens", collection = "H"); RE <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
L <- c(split(H$gene_symbol, H$gs_name), split(RE$gene_symbol, RE$gs_name))
focus <- c("REACTOME_CHOLESTEROL_BIOSYNTHESIS", "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS")
stopifnot(all(focus %in% names(L)))
S <- lapply(L[focus], function(x) intersect(unique(x), rownames(g)))
core <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS")
srebp <- unique(c(S$REACTOME_CHOLESTEROL_BIOSYNTHESIS, S$REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP, S$HALLMARK_CHOLESTEROL_HOMEOSTASIS))
shared <- unique(unlist(lapply(core, function(a) unlist(lapply(setdiff(core, a), function(b) intersect(S[[a]], S[[b]]))))))
for (a in core) S[[paste0(a, "_UNIQUE")]] <- setdiff(S[[a]], c(shared, srebp))
cat("shared among mTORC1/hypoxia/glycolysis:", length(shared), "; SREBP/cholesterol genes in mTORC1:", length(intersect(S$HALLMARK_MTORC1_SIGNALING, srebp)), "\n")
cam <- function(kk) { r <- camera(g, S, design, contrast = cm[, kk], inter.gene.cor = NA, trend.var = TRUE, sort = FALSE); r }
r1 <- cam("C1"); r2 <- cam("C2"); r3 <- cam("C3")
tt <- function(kk) { x <- topTable(fit, coef = kk, number = Inf, sort.by = "none"); setNames(x$logFC, rownames(x)) }
l1 <- tt("C1"); l2 <- tt("C2"); l3 <- tt("C3")
res <- data.table(set = names(S), n = lengths(S), PIK75_mean = sapply(S, function(x) round(mean(l1[x]), 3)), PIK75_pct_down = sapply(S, function(x) round(100 * mean(l1[x] < 0))),
  TGX221_mean = sapply(S, function(x) round(mean(l2[x]), 3)), C3_mean = sapply(S, function(x) round(mean(l3[x]), 3)),
  inter_gene_cor = round(r1$Correlation, 3),
  camera_C1_dir = r1$Direction, camera_C1_p = signif(r1$PValue, 2), camera_C1_FDR = signif(p.adjust(r1$PValue[1:16], "BH")[seq_along(S)], 2),
  camera_C2_dir = r2$Direction, camera_C2_p = signif(r2$PValue, 2),
  camera_C3_dir = r3$Direction, camera_C3_p = signif(r3$PValue, 2), camera_C3_FDR = signif(p.adjust(r3$PValue[1:16], "BH")[seq_along(S)], 2))
res[set %like% "_UNIQUE", c("camera_C1_FDR", "camera_C3_FDR") := NA]
print(res[, .(set = substr(set, 1, 45), n, PIK75_mean, PIK75_pct_down, TGX221_mean, camera_C1_p, camera_C1_FDR, camera_C2_p, camera_C3_dir, camera_C3_p, camera_C3_FDR)])
dir.create("results_rev", showWarnings = FALSE); fwrite(res, "results_rev/Table_S24_camera_C1_C2_C3.csv")
