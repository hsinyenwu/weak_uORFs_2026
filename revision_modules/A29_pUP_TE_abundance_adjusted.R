################################################################################
# A29 (2026-11-07, V22; logic review) — Is the %uORF+ versus mORF TE relation of the single-TuORF cohort a by-product of mRNA
#   abundance? TE rises with RNA abundance and %uORF+ falls with it, so the relation is re-examined (i) as a partial Spearman
#   correlation given log2 RNA TPM (and given RNA TPM and leader length), computed from the residuals of rank regressions, and
#   (ii) within thirds of the cohort by RNA abundance (Spearman rho and the median TE of the <30% and 80-100% uORF+ groups).
# Input: weak_uORF_revision/revision_tables/TableS3_single_TuORF_genes_revised.xlsx (genes with TE > 0 and RNA TPM > 0).
# Outputs: output/A29_Table_pUP_TE_abundance.csv, output/A29_summary.txt.   Run from Weak_uORF_submission/:  Rscript revision_modules/A29_pUP_TE_abundance_adjusted.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
S3 <- as.data.table(read.xlsx("/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS3_single_TuORF_genes_revised.xlsx", sheet = 1))
d <- S3[is.finite(mORF_TE) & mORF_TE > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]
# V23 (2026-11-08): leader length = annotated 5' leader of the TuORF's isoform (Araport11+CTRL); before, the column TL_region_nt was used, which is the uORF length
suppressPackageStartupMessages({ library(GenomicFeatures) }); txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF("Files/Araport11+CTRL_20181206.gtf", format = "gtf") else GenomicFeatures::makeTxDbFromGFF("Files/Araport11+CTRL_20181206.gtf", format = "gtf")
ll <- sum(width(fiveUTRsByTranscript(txdb, use.names = TRUE))); d[, leader_len := as.numeric(ll[tx_id])]; stopifnot(mean(is.finite(d$leader_len)) > 0.99)
psp <- function(x, y, Z) { Z <- as.data.frame(lapply(as.data.frame(Z), rank)); rx <- resid(lm(rank(x) ~ ., data = Z)); ry <- resid(lm(rank(y) ~ ., data = Z)); cor(rx, ry) }
sp <- function(x, y) cor(x, y, method = "spearman")
out <- list(); addt <- function(test, stat, value, n, note = "") out[[length(out) + 1]] <<- data.table(test = test, statistic = stat, value = round(value, 4), n = n, note = note)
.add("A29: single-TuORF cohort genes with TE > 0 and RNA TPM > 0: n = %d", nrow(d))
r0 <- sp(d$pct_uORF_plus, d$mORF_TE); rA <- sp(d$pct_uORF_plus, d$RNA_TPM); rT <- sp(d$mORF_TE, d$RNA_TPM)
addt("unadjusted", "Spearman rho(%uORF+, TE)", r0, nrow(d)); addt("covariates", "Spearman rho(%uORF+, RNA TPM)", rA, nrow(d)); addt("covariates", "Spearman rho(TE, RNA TPM)", rT, nrow(d))
.add("  Spearman rho(%%uORF+, TE) = %.3f | rho(%%uORF+, RNA TPM) = %.3f | rho(TE, RNA TPM) = %.3f", r0, rA, rT)
p1 <- psp(d$pct_uORF_plus, d$mORF_TE, data.frame(a = log2(d$RNA_TPM))); addt("partial", "partial Spearman rho(%uORF+, TE | log2 RNA TPM)", p1, nrow(d))
p2 <- psp(d$pct_uORF_plus[is.finite(d$leader_len)], d$mORF_TE[is.finite(d$leader_len)], data.frame(a = log2(d$RNA_TPM[is.finite(d$leader_len)]), b = log2(d$leader_len[is.finite(d$leader_len)]))); addt("partial", "partial Spearman rho(%uORF+, TE | log2 RNA TPM, log2 leader length)", p2, nrow(d))
.add("  partial rho given log2 RNA TPM = %.3f; given RNA TPM and leader length = %.3f", p1, p2)
d[, abundance_third := cut(rank(RNA_TPM, ties.method = "first"), 3, labels = c("lowest third", "middle third", "highest third"))]
tb <- d[, .(n = .N, RNA_TPM_range = sprintf("%.1f-%.1f", min(RNA_TPM), max(RNA_TPM)), rho_pUP_TE = sp(pct_uORF_plus, mORF_TE), median_TE_low_pUP = median(mORF_TE[pct_uORF_plus < 30]), n_low = sum(pct_uORF_plus < 30),
            median_TE_high_pUP = median(mORF_TE[pct_uORF_plus >= 80]), n_high = sum(pct_uORF_plus >= 80), p_wilcoxon = wilcox.test(mORF_TE[pct_uORF_plus < 30], mORF_TE[pct_uORF_plus >= 80])$p.value), by = abundance_third][order(abundance_third)]
for (i in seq_len(nrow(tb))) { .add("  %s (n = %d, RNA TPM %s): rho(%%uORF+, TE) = %.3f; median TE <30%% uORF+ %.2f (n = %d) vs 80-100%% %.2f (n = %d), ratio %.2f, Wilcoxon p = %.1e",
    tb$abundance_third[i], tb$n[i], tb$RNA_TPM_range[i], tb$rho_pUP_TE[i], tb$median_TE_low_pUP[i], tb$n_low[i], tb$median_TE_high_pUP[i], tb$n_high[i], tb$median_TE_low_pUP[i] / tb$median_TE_high_pUP[i], tb$p_wilcoxon[i])
  addt(paste("within", tb$abundance_third[i]), "Spearman rho(%uORF+, TE)", tb$rho_pUP_TE[i], tb$n[i]); addt(paste("within", tb$abundance_third[i]), "median TE <30% uORF+", tb$median_TE_low_pUP[i], tb$n_low[i]); addt(paste("within", tb$abundance_third[i]), "median TE 80-100% uORF+", tb$median_TE_high_pUP[i], tb$n_high[i]) }
fwrite(rbindlist(out), "output/A29_Table_pUP_TE_abundance.csv"); fwrite(tb, "output/A29_Table_pUP_TE_by_abundance_third.csv")
writeLines(c(sprintf("A29 %%uORF+ vs TE adjusted for mRNA abundance — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), "output/A29_summary.txt"); cat("A29 done\n")
