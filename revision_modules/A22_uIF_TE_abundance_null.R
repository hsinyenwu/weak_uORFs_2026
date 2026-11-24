################################################################################
# A22 — Is the part of the uIF–TE correlation that exceeds the shared-term null  #
#       (Supplemental Fig. S3D) also expected by construction?                   #
#                                                                               #
# uIF = uORF P-site density / mORF P-site density; TE = Ribo TPM / RNA TPM.     #
# Null A (A3, Fig. S3D): the uORF density is shuffled across genes; each gene   #
#   keeps its mORF density and TE.  This breaks the link between the uORF       #
#   density and the gene's mRNA abundance, although the uORF density of a gene  #
#   must scale with its mRNA abundance whatever the uORF does.                  #
# Null B (abundance-preserving): the uORF density PER RNA (uORF_TL / RNA_TPM)   #
#   is shuffled across genes and multiplied by each gene's own RNA_TPM.  Under   #
#   this null the ribosome occupancy of the uORF per transcript is unrelated to #
#   TE, but the uORF density still rises with mRNA abundance.                   #
# For both nulls: r(log2 uIF, log2 TE), the partial correlation at fixed mORF   #
#   density, and the median TE of the high / medium / low uIF classes (Fig 1C). #
# Input: output/Table_S2_single_TuORF_genes.xlsx (2,995 single-TuORF genes).   #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A22_uIF_TE_abundance_null.R
# Not part of the manuscript (check for Larry's question of 2026-11-04).        #
################################################################################
suppressPackageStartupMessages({ library(openxlsx); library(dplyr) })
set.seed(20261104); N_PERM <- 2000; NO_TUORF_MEDIAN_TE <- 1.22   # No-TuORF genes, n = 12,724 (A3 summary, section 7)
S2 <- read.xlsx("output/Table_S2_single_TuORF_genes.xlsx", sheet = 1)
C <- S2 %>% filter(is.finite(uORF_TL), is.finite(mORF_TL), is.finite(RNA_TPM), is.finite(mORF_TE), uORF_TL > 0, mORF_TL > 0, RNA_TPM > 0, mORF_TE > 0)
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A22 uIF–TE nulls — %s; single-TuORF genes with TE and positive densities: n = %d", format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(C)))
lu <- log2(C$uORF_TL); lm_ <- log2(C$mORF_TL); lte <- log2(C$mORF_TE); la <- log2(C$RNA_TPM); luif <- lu - lm_
q <- C$uORF_TL / C$RNA_TPM            # uORF density per RNA
pc <- function(x, y, z) cor(resid(lm(x ~ z)), resid(lm(y ~ z)))
cls <- function(uif) factor(ifelse(uif > 1, "High", ifelse(uif > 0.2, "Medium", "Low")), levels = c("High", "Medium", "Low"))
med_by <- function(uif) tapply(C$mORF_TE, cls(uif), median)
.add(sprintf("observed: r(log2 uIF, log2 TE) = %.3f; partial r at fixed mORF density = %.3f", cor(luif, lte), pc(luif, lte, lm_)))
.add(sprintf("observed pieces: r(log2 uORF density, log2 RNA TPM) = %.3f; r(log2 mORF density, log2 RNA TPM) = %.3f; r(log2 uORF density, log2 mORF density) = %.3f", cor(lu, la), cor(lm_, la), cor(lu, lm_)))
.add(sprintf("observed: r(log2 uORF density per RNA, log2 TE) = %.3f (shared RNA denominator biases this towards +); r(log2 RNA TPM, log2 TE) = %.3f", cor(log2(q), lte), cor(la, lte)))
m_obs <- med_by(C$uIF); .add(sprintf("observed median TE by uIF class: High %.2f (n=%d) | Medium %.2f (n=%d) | Low %.2f (n=%d) | No-TuORF %.2f", m_obs[1], sum(cls(C$uIF) == "High"), m_obs[2], sum(cls(C$uIF) == "Medium"), m_obs[3], sum(cls(C$uIF) == "Low"), NO_TUORF_MEDIAN_TE))
run_null <- function(make_uif, label) {
  res <- t(replicate(N_PERM, { u <- make_uif(); lu2 <- log2(u); m <- med_by(u); c(r = cor(lu2, lte), partial = pc(lu2, lte, lm_), High = m[["High"]], Medium = m[["Medium"]], Low = m[["Low"]], nHigh = sum(cls(u) == "High"), nLow = sum(cls(u) == "Low")) }))
  res <- as.data.frame(res)
  .add(sprintf("%s: r mean %.3f [%.3f, %.3f]; partial r mean %.3f [%.3f, %.3f]", label, mean(res$r), quantile(res$r, .025), quantile(res$r, .975), mean(res$partial), quantile(res$partial, .025), quantile(res$partial, .975)))
  .add(sprintf("%s: median TE by class, mean over permutations: High %.2f | Medium %.2f | Low %.2f (class sizes High %.0f, Low %.0f); Low >= No-TuORF (%.2f) in %.1f%% of permutations",
               label, mean(res$High), mean(res$Medium), mean(res$Low), mean(res$nHigh), mean(res$nLow), NO_TUORF_MEDIAN_TE, 100 * mean(res$Low >= NO_TUORF_MEDIAN_TE)))
  .add(sprintf("%s: observed r more negative than the null in %.1f%% of permutations; observed partial r more negative in %.1f%%", label, 100 * mean(res$r > cor(luif, lte)), 100 * mean(res$partial > pc(luif, lte, lm_))))
  res }
resA <- run_null(function() sample(C$uORF_TL) / C$mORF_TL, "null A (uORF density shuffled; A3/S3D)")
resB <- run_null(function() sample(q) * C$RNA_TPM / C$mORF_TL, "null B (uORF density per RNA shuffled, abundance kept)")
.add(sprintf("link of the uORF density to mRNA abundance, r(log2 uORF density, log2 RNA TPM): observed %.3f; null B %.3f (null A 0 by construction); r(log2 uORF density per RNA, log2 RNA TPM) observed = %.3f", cor(lu, la), mean(replicate(200, cor(log2(sample(q) * C$RNA_TPM), la))), cor(log2(q), la)))
# How much of the observed correlation each null reproduces
.add(sprintf("share of the observed r reproduced: null A %.0f%%, null B %.0f%%", 100 * mean(resA$r) / cor(luif, lte), 100 * mean(resB$r) / cor(luif, lte)))
# Within-abundance check: uIF vs TE inside RNA-abundance quintiles (the abundance term cannot act within a narrow abundance band)
C$abin <- cut(rank(C$RNA_TPM), 5, labels = paste0("Q", 1:5))
w <- C %>% group_by(abin) %>% summarise(n = n(), r = cor(log2(uORF_TL / mORF_TL), log2(mORF_TE)), r_q_TE = cor(log2(uORF_TL / RNA_TPM), log2(mORF_TE)), .groups = "drop")
.add("uIF vs TE within RNA-abundance quintiles (r; and r of uORF density per RNA vs TE):"); for (i in seq_len(nrow(w))) .add(sprintf("  %s n=%d r=%.3f r(q,TE)=%.3f", w$abin[i], w$n[i], w$r[i], w$r_q_TE[i]))
write.csv(data.frame(perm = seq_len(N_PERM), r_nullA = resA$r, partial_nullA = resA$partial, r_nullB = resB$r, partial_nullB = resB$partial, lowTE_nullA = resA$Low, lowTE_nullB = resB$Low), "output/A22_uIF_TE_nulls.csv", row.names = FALSE)
writeLines(.log, "output/A22_summary.txt"); .add("written: output/A22_summary.txt, output/A22_uIF_TE_nulls.csv")
