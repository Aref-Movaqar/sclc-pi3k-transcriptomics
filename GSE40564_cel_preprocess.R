#!/usr/bin/env Rscript
# Part 3a: re-process GSE40564 from CEL files (run before GSE40564_audit.R)
# Input: unpack https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40564/suppl/GSE40564_RAW.tar into ./cel/
# RMA is implemented step by step (affy::bg.adjust convolution background ->
# optional limma::normalizeQuantiles -> Tukey median polish) so that the quantile
# step can be switched off. On most systems affy::rma(ab) / rma(ab, normalize = FALSE)
# gives the same result in one call.
suppressPackageStartupMessages({library(affy); library(limma)})
cels <- sort(list.files("cel", pattern = "CEL.gz$", full.names = TRUE))
ab <- ReadAffy(filenames = cels)
gsm <- sub("_.*", "", basename(cels))
qc <- data.frame(sample = sub(".*H69_", "", sub(".CEL.gz", "", basename(cels))), gsm = gsm,
                 scan = ab@protocolData@data$ScanDate,
                 raw_PM_median = round(apply(log2(pm(ab)), 2, median), 2),
                 raw_PM_q95 = round(apply(log2(pm(ab)), 2, quantile, 0.95), 2),
                 RNAdeg_slope = round(AffyRNAdeg(ab)$slope, 2))
P  <- pm(ab); pn <- probeNames(ab)
bg <- log2(apply(P, 2, bg.adjust))                       # RMA convolution background, per array
summarise <- function(m) {                                # median polish per probe set
  rows <- split(seq_len(nrow(m)), pn)
  r <- vapply(rows, function(i) { mp <- medpolish(m[i, , drop = FALSE], trace.iter = FALSE); mp$overall + mp$col }, numeric(ncol(m)))
  r <- t(r); colnames(r) <- gsm; r
}
own  <- summarise(normalizeQuantiles(bg))                 # RMA with quantile normalisation
nonq <- summarise(bg)                                     # same, without the quantile step
rle <- sweep(own, 1, apply(own, 1, median))               # relative log expression
qc$RLE_median <- round(apply(rle, 2, median), 3); qc$RLE_IQR <- round(apply(rle, 2, IQR), 3)
saveRDS(list(own = own, nonq = nonq, qc = qc), "cel/cel_processed.rds")
print(qc)
