################################################################################
# A39 — Does the uncapped ratio follow translation or the decay rate? (Larry,   #
# 2026-11-20: "Is it possible that the less-translated genes are decayed        #
# faster, so there are fewer PARE and GMUCT signals?"). Not part of the         #
# manuscript.                                                                   #
# Across the single-TuORF genes and the genes without TuORFs together (the gene #
# table of module A31), rank-based partial Spearman correlations of the         #
# uncapped ratio (GMUCT; PARE) with mORF TE given half-life, and with half-life #
# given mORF TE, each also with exon length and RNA TPM as covariates; and the  #
# correlation of half-life with mORF TE.                                        #
# If faster decay made more decay intermediates per mRNA, the uncapped ratio    #
# would fall with half-life; if the intermediates mainly come from decay during #
# translation, it would rise with TE when half-life is held constant.           #
# Output: output/A39_uncapped_TE_halflife_summary.txt                           #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A39_uncapped_TE_halflife_partial.R
################################################################################
suppressPackageStartupMessages(library(data.table))
options(warn = 1); set.seed(20261120)
OUT_DIR <- "output"
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A39 uncapped ratio, mORF TE and half-life — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
a11 <- fread(file.path(OUT_DIR, "A11_Table_genes.csv")); a12 <- fread(file.path(OUT_DIR, "A12_Table_genes.csv")); a12 <- a12[match(a11$gene_id, gene_id)]
len <- fread(file.path(OUT_DIR, "A12_gene_counts.csv"), select = c("gene_id", "exon_length"))
G <- data.table(gene_id = a11$gene_id, TE = a11$mORF_TE, RNA_TPM = a11$RNA_TPM_this_study, HL = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_))
G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]
pr <- function(x, y, Z = NULL) {   # rank-based partial correlation with a bootstrap 95% interval
  ok <- is.finite(x) & is.finite(y); if (!is.null(Z)) ok <- ok & Reduce(`&`, lapply(Z, is.finite))
  rx <- rank(x[ok]); ry <- rank(y[ok])
  if (!is.null(Z)) { RZ <- as.data.frame(lapply(Z, function(z) rank(z[ok]))); rx <- resid(lm(rx ~ ., data = RZ)); ry <- resid(lm(ry ~ ., data = RZ)) }
  r <- cor(rx, ry); b <- replicate(500, { s <- sample.int(length(rx), replace = TRUE); cor(rx[s], ry[s]) })
  sprintf("rho = %6.3f (%.3f to %.3f), n = %s", r, quantile(b, 0.025), quantile(b, 0.975), format(sum(ok), big.mark = ",")) }
d <- G[is.finite(TE) & TE > 0 & is.finite(HL) & HL > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]
.add("Across single-TuORF and No-TuORF genes with TE > 0, half-life and RNA TPM:")
.add("  half-life vs mORF TE:                                   ", pr(d$HL, d$TE))
for (v in c(U_gmuct = "GMUCT", U_pare = "PARE")) {
  u <- d[[names(which(c(U_gmuct = "GMUCT", U_pare = "PARE") == v))]]
  .add(sprintf("  uncapped ratio (%s) vs mORF TE:                     %s", v, pr(u, d$TE)))
  .add(sprintf("  uncapped ratio (%s) vs half-life:                   %s", v, pr(u, d$HL)))
  .add(sprintf("  uncapped ratio (%s) vs mORF TE, given half-life:    %s", v, pr(u, d$TE, list(d$HL))))
  .add(sprintf("  uncapped ratio (%s) vs half-life, given mORF TE:    %s", v, pr(u, d$HL, list(d$TE))))
  .add(sprintf("  uncapped ratio (%s) vs mORF TE, given half-life, exon length and RNA TPM: %s", v, pr(u, d$TE, list(d$HL, d$exon_length, d$RNA_TPM))))
  .add(sprintf("  uncapped ratio (%s) vs half-life, given mORF TE, exon length and RNA TPM: %s", v, pr(u, d$HL, list(d$TE, d$exon_length, d$RNA_TPM))))
}
writeLines(.log, file.path(OUT_DIR, "A39_uncapped_TE_halflife_summary.txt")); .add("written: output/A39_uncapped_TE_halflife_summary.txt")
