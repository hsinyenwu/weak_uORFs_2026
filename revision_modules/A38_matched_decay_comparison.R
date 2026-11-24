################################################################################
# A38 — Matched comparison of mRNA half-life and uncapped ratios (Larry,        #
# 2026-11-20: "You use other similar length but no TuORF genes to do this?",    #
# then "Yes, run it"). A check of the linear-model adjustment of module A31.    #
# Not part of the manuscript.                                                   #
# For each single-TuORF gene of a %uORF+ group, the 5 genes without TuORFs that #
# are closest in transcript length (log10 exon length) and RNA abundance        #
# (log10 RNA TPM of this study) are taken as its controls (distance on the two  #
# standardized log values; controls may be used for several genes). The group  #
# is then compared with its matched controls:                                   #
#   ratio of medians (group / matched controls), 95% bootstrap interval over    #
#     the group genes (1,000 resamples, each gene with its own controls);       #
#   geometric mean of the per-gene ratios (gene / geometric mean of its 5       #
#     controls), 95% interval from the t distribution on the log values, and    #
#     the Wilcoxon signed-rank p of the log ratios.                             #
# A second matching uses transcript length only.                                #
# Balance: geometric means of length and TPM of the group / matched controls.   #
# Inputs as module A31: output/A11_Table_genes.csv, output/A12_Table_genes.csv, #
# output/A12_gene_counts.csv. Output: output/A38_matched_decay_summary.txt,     #
# output/A38_matched_decay.csv                                                  #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A38_matched_decay_comparison.R
################################################################################
suppressPackageStartupMessages(library(data.table))
options(warn = 1); set.seed(20261120)
OUT_DIR <- "output"; K <- 5L; NBOOT <- 1000L
PUP <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A38 matched comparison of the decay measures — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
a11 <- fread(file.path(OUT_DIR, "A11_Table_genes.csv")); a12 <- fread(file.path(OUT_DIR, "A12_Table_genes.csv"))
a12 <- a12[match(a11$gene_id, gene_id)]; stopifnot(identical(a11$gene_id, a12$gene_id))
len <- fread(file.path(OUT_DIR, "A12_gene_counts.csv"), select = c("gene_id", "exon_length"))
G <- data.table(gene_id = a11$gene_id, g = a11$pct_uORF_plus_group, RNA_TPM = a11$RNA_TPM_this_study, HL = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_))
G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]; stopifnot(!anyNA(G$exon_length), all(G$g %in% PUP))
MEAS <- c(HL = "mRNA half-life", U_pare = "uncapped ratio (PARE)", U_gmuct = "uncapped ratio (GMUCT)")
gm <- function(x) exp(mean(log(x)))
match_ctrl <- function(Xg, Xc, k) {   # indices of the k nearest controls for every group gene (Euclidean distance on standardized columns)
  t(vapply(seq_len(nrow(Xg)), function(i) { dd <- colSums((t(Xc) - Xg[i, ])^2); order(dd)[seq_len(k)] }, integer(k))) }
res <- list()
for (v in names(MEAS)) for (mode in c("length and abundance", "length only")) {
  d <- G[is.finite(get(v)) & get(v) > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]
  ctl <- d[g == "No TuORF"]; cols <- if (mode == "length only") c("exon_length") else c("exon_length", "RNA_TPM")
  Xc0 <- as.matrix(log10(ctl[, ..cols])); mu <- colMeans(Xc0); sdv <- apply(Xc0, 2, sd); Xc <- sweep(sweep(Xc0, 2, mu), 2, sdv, "/")
  for (gr in PUP[1:3]) {
    dg <- d[g == gr]; Xg <- sweep(sweep(as.matrix(log10(dg[, ..cols])), 2, mu), 2, sdv, "/")
    idx <- match_ctrl(Xg, Xc, K)
    yg <- dg[[v]]; yc <- matrix(ctl[[v]][idx], ncol = K)
    dist1 <- sqrt(rowSums((Xg - Xc[idx[, K], , drop = FALSE])^2))   # distance to the 5th control
    ratio_med <- median(yg) / median(as.vector(yc))
    boot <- replicate(NBOOT, { s <- sample.int(length(yg), replace = TRUE); median(yg[s]) / median(as.vector(yc[s, ])) })
    lr <- log(yg) - rowMeans(log(yc)); tt <- t.test(lr); wx <- wilcox.test(lr)
    res[[length(res) + 1]] <- data.table(measure = MEAS[[v]], matching = mode, group = gr, n_group = length(yg), n_controls_available = nrow(ctl), n_distinct_controls_used = length(unique(as.vector(idx))),
      length_ratio_before = gm(dg$exon_length) / gm(ctl$exon_length), length_ratio_after = gm(dg$exon_length) / gm(ctl$exon_length[as.vector(idx)]),
      TPM_ratio_before = gm(dg$RNA_TPM) / gm(ctl$RNA_TPM), TPM_ratio_after = gm(dg$RNA_TPM) / gm(ctl$RNA_TPM[as.vector(idx)]),
      median_distance_to_5th_control_SD = median(dist1),
      ratio_of_medians_unmatched = median(yg) / median(ctl[[v]]), ratio_of_medians_matched = ratio_med, rm_lo = unname(quantile(boot, 0.025)), rm_hi = unname(quantile(boot, 0.975)),
      geo_mean_ratio_unmatched = gm(yg) / gm(ctl[[v]]), geo_mean_ratio_matched = exp(mean(lr)), gmr_lo = exp(tt$conf.int[1]), gmr_hi = exp(tt$conf.int[2]), wilcoxon_p = wx$p.value)
  }
}
R <- rbindlist(res)
for (i in seq_len(nrow(R))) with(R[i], .add(sprintf("%-24s %-21s %-18s n = %4d (controls used %4d of %5d) | balance length %.2f -> %.2f, TPM %.2f -> %.2f | median ratio %.2f -> %.2f (%.2f-%.2f) | geometric mean ratio %.2f -> %.2f (%.2f-%.2f), Wilcoxon p = %.2g",
  measure, matching, group, n_group, n_distinct_controls_used, n_controls_available, length_ratio_before, length_ratio_after, TPM_ratio_before, TPM_ratio_after,
  ratio_of_medians_unmatched, ratio_of_medians_matched, rm_lo, rm_hi, geo_mean_ratio_unmatched, geo_mean_ratio_matched, gmr_lo, gmr_hi, wilcoxon_p)))
# (V53, 2026-11-23) the values of this note are read from output/A31_csv_adj.csv (run module A31 first); until V52 they were typed here and went out of date
adj31 <- fread(file.path(OUT_DIR, "A31_csv_adj.csv")); adj31 <- adj31[adj31$group == "80 to 100% uORF+", ]
.v31 <- function(ms, md) { r <- adj31[adj31$measure == ms & adj31$model == md, ]; stopifnot(nrow(r) == 1); r }
.p31 <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
ms31 <- c("mRNA half-life (h)", "uncapped ratio (PARE)", "uncapped ratio (GMUCT)")
b31 <- lapply(ms31, .v31, md = "adjusted for transcript length and RNA abundance"); l31 <- lapply(ms31, .v31, md = "adjusted for transcript length")
.add(sprintf("For comparison, module A31 (linear models of the log10 value, adjusted for log10 exon length and log10 RNA TPM), 80-100%% uORF+ vs genes without TuORFs: half-life %.2f (%s), PARE %.2f (%s), GMUCT %.2f (%s); length only: %.2f, %.2f, %.2f.",
             b31[[1]]$fold_vs_No_TuORF, .p31(b31[[1]]$p), b31[[2]]$fold_vs_No_TuORF, .p31(b31[[2]]$p), b31[[3]]$fold_vs_No_TuORF, .p31(b31[[3]]$p), l31[[1]]$fold_vs_No_TuORF, l31[[2]]$fold_vs_No_TuORF, l31[[3]]$fold_vs_No_TuORF))
fwrite(R, file.path(OUT_DIR, "A38_matched_decay.csv")); writeLines(.log, file.path(OUT_DIR, "A38_matched_decay_summary.txt"))
.add("written: output/A38_matched_decay_summary.txt, output/A38_matched_decay.csv")
