################################################################################
# A41 — The decay results with a second half-life dataset (Larry, 2026-11-20:   #
# "if you use Sorenson et al. (2018) ... do you see the same result").         #
# Not part of the manuscript.                                                   #
# Half-lives: Sorenson et al. 2018 (PNAS; cordycepin transcription shut-off,    #
# seedlings), wild-type decay rate alpha (per minute) from                      #
# ~/Desktop/uORFs_miRNA/Sorenson_et_al_pnas2017_DR.xlsx (lab copy, 18,674 genes; #
# genes with every alpha > 0 kept, as in the lab script Sorenson_2017_half_lives.R); #
# half-life = ln(2) / alpha in minutes, capped at 1,440 min (24 h; very stable  #
# genes have alpha near 0), as in the half-life prediction project.             #
# Repeats with these half-lives: (1) the %uORF+ groups against genes without    #
# TuORFs (medians, Kolmogorov-Smirnov, geometric means unadjusted and adjusted  #
# for log10 exon length and log10 RNA TPM, length only, and the matched         #
# comparison of module A38); (2) the uncapped ratio against half-life and mORF  #
# TE (module A40: correlations, partial correlations, medians in fifths of      #
# half-life within thirds of TE); (3) agreement with the 5-EU half-lives.       #
# Outputs: output/A41_sorenson_summary.txt, output/A41_Figure_sorenson_halflife.pdf
# Run from Weak_uORF_submission/:  Rscript revision_modules/A41_sorenson_halflife_check.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(openxlsx); library(ggplot2); library(patchwork) })
options(warn = 1); set.seed(20261120)
OUT_DIR <- "output"; CAP <- 1440; K <- 5L; NBOOT <- 1000L
PUP <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A41 Sorenson et al. 2018 half-lives — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
DR <- as.data.table(read.xlsx("/Users/wu/Desktop/uORFs_miRNA/Sorenson_et_al_pnas2017_DR.xlsx")); stopifnot(nrow(DR) == 18674L, "WT" %in% names(DR))
DR <- DR[WT > 0 & sov > 0 & vcs > 0 & vcs.sov > 0]; stopifnot(nrow(DR) == 18576L)
DR[, HLs := pmin(log(2) / WT, CAP)]
.add(sprintf("Sorenson genes: %s; half-life median %.0f min; capped at %d min: %s genes", format(nrow(DR), big.mark = ","), median(DR$HLs), CAP, format(sum(DR$HLs >= CAP), big.mark = ",")))
a11 <- fread(file.path(OUT_DIR, "A11_Table_genes.csv")); a12 <- fread(file.path(OUT_DIR, "A12_Table_genes.csv")); a12 <- a12[match(a11$gene_id, gene_id)]
len <- fread(file.path(OUT_DIR, "A12_gene_counts.csv"), select = c("gene_id", "exon_length"))
G <- data.table(gene_id = a11$gene_id, g = factor(a11$pct_uORF_plus_group, PUP), TE = a11$mORF_TE, RNA_TPM = a11$RNA_TPM_this_study, HL5 = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_))
G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]; G[, HLs := DR$HLs[match(gene_id, DR$gene_id)]]
stopifnot(!anyNA(G$g))
.add(sprintf("genes of the study (single-TuORF + No TuORF) with a Sorenson half-life: %s of %s; with both half-lives: %s; Spearman rho (Sorenson, 5-EU) = %.2f",
             format(sum(is.finite(G$HLs)), big.mark = ","), format(nrow(G), big.mark = ","), format(sum(is.finite(G$HLs) & is.finite(G$HL5)), big.mark = ","),
             cor(G$HLs, G$HL5, method = "spearman", use = "complete.obs")))
# (1) groups
d <- G[is.finite(HLs)]; ref <- d[g == "No TuORF", HLs]
.add(""); .add("(1) half-life (Sorenson) by %uORF+ group: median, ratio of medians vs No TuORF, KS p")
for (gr in PUP) { x <- d[g == gr, HLs]; .add(sprintf("  %-18s n = %5d, median %5.0f min, %.2fx, KS p = %.2g", gr, length(x), median(x), median(x) / median(ref), if (gr == "No TuORF") NA else suppressWarnings(ks.test(x, ref)$p.value))) }
dm <- d[is.finite(RNA_TPM) & RNA_TPM > 0]; dm[, g := relevel(g, ref = "No TuORF")]
for (mm in list(c("unadjusted"), c("length", "log10(exon_length)"), c("abundance", "log10(RNA_TPM)"), c("length and abundance", "log10(exon_length)", "log10(RNA_TPM)"))) {
  f <- as.formula(paste("log10(HLs) ~ g", if (length(mm) > 1) paste("+", paste(mm[-1], collapse = " + ")) else "")); m <- lm(f, dm); s <- summary(m)$coefficients
  .add(sprintf("  geometric mean ratio, %-21s: %s", mm[1], paste(sprintf("%s %.2fx (p = %.2g)", PUP[1:3], 10^s[paste0("g", PUP[1:3]), 1], s[paste0("g", PUP[1:3]), 4]), collapse = "; ")))
}
ctl <- dm[g == "No TuORF"]; Xc0 <- as.matrix(log10(ctl[, .(exon_length, RNA_TPM)])); mu <- colMeans(Xc0); sdv <- apply(Xc0, 2, sd); Xc <- sweep(sweep(Xc0, 2, mu), 2, sdv, "/")
for (gr in PUP[1:3]) {
  dg <- dm[g == gr]; Xg <- sweep(sweep(as.matrix(log10(dg[, .(exon_length, RNA_TPM)])), 2, mu), 2, sdv, "/")
  idx <- t(vapply(seq_len(nrow(Xg)), function(i) order(colSums((t(Xc) - Xg[i, ])^2))[seq_len(K)], integer(K)))
  yc <- matrix(ctl$HLs[idx], ncol = K); rm <- median(dg$HLs) / median(as.vector(yc))
  b <- replicate(NBOOT, { s <- sample.int(nrow(dg), replace = TRUE); median(dg$HLs[s]) / median(as.vector(yc[s, ])) })
  .add(sprintf("  matched (5 genes without TuORFs of similar length and abundance): %-18s n = %4d, ratio of medians %.2fx (%.2f-%.2f)", gr, nrow(dg), rm, quantile(b, 0.025), quantile(b, 0.975)))
}
# (2) uncapped ratio vs Sorenson half-life and TE
prho <- function(x, y, Z) { rx <- rank(x); ry <- rank(y); RZ <- as.data.frame(lapply(Z, rank)); cor(resid(lm(rx ~ ., data = RZ)), resid(lm(ry ~ ., data = RZ))) }
cut_q <- function(x, k) cut(x, unique(quantile(x, seq(0, 1, length.out = k + 1))), include.lowest = TRUE, labels = FALSE)
.add(""); .add("(2) uncapped ratio vs half-life (Sorenson) and mORF TE, single-TuORF + No-TuORF genes with TE > 0, RNA TPM > 0")
plots <- list(); dg <- position_dodge(width = 0.3)
for (v in c(U_gmuct = "GMUCT", U_pare = "PARE")) {
  key <- names(which(c(U_gmuct = "GMUCT", U_pare = "PARE") == v))
  e <- G[is.finite(get(key)) & get(key) > 0 & is.finite(HLs) & is.finite(TE) & TE > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]; e[, U := get(key)]
  .add(sprintf("  %s: n = %s; rho with half-life %.2f (partial, TE, length and TPM constant: %.2f); rho with TE %.2f (partial, half-life, length and TPM constant: %.2f)", v, format(nrow(e), big.mark = ","),
               cor(e$U, e$HLs, method = "spearman"), prho(e$U, e$HLs, list(e$TE, e$exon_length, e$RNA_TPM)), cor(e$U, e$TE, method = "spearman"), prho(e$U, e$TE, list(e$HLs, e$exon_length, e$RNA_TPM))))
  e[, hl5 := cut_q(HLs, 5)]; e[, te3 := factor(cut_q(TE, 3), 1:3, c("low TE", "middle TE", "high TE"))]
  A <- e[, .(med = median(U), lo = quantile(replicate(NBOOT, median(sample(U, replace = TRUE))), 0.025), hi = quantile(replicate(NBOOT, median(sample(U, replace = TRUE))), 0.975)), by = .(te3, hl5)]
  for (t3 in levels(e$te3)) .add(sprintf("    %-9s medians by half-life fifth (short to long): %s", t3, paste(sprintf("%.2f", A[te3 == t3][order(hl5), med]), collapse = ", ")))
  xl <- sprintf("%.0f", e[, median(HLs), by = hl5][order(hl5), V1])
  plots[[v]] <- ggplot(A, aes(factor(hl5), med, colour = te3, group = te3)) + geom_line(linewidth = 0.5, position = dg) + geom_point(size = 1.4, position = dg) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.35, position = dg) + scale_x_discrete(labels = xl) +
    scale_colour_manual(values = c("#9ecae1", "#4292c6", "#08306b"), name = NULL) +
    labs(x = "half-life, Sorenson (min; median of fifth)", y = sprintf("uncapped ratio (%s), median", v), subtitle = sprintf("%s versus Sorenson half-life within thirds of TE", v))
}
# figure: half-life by group (cumulative) + the two uncapped panels
cdf <- ggplot(d, aes(HLs, colour = g)) + stat_ecdf(linewidth = 0.5) + scale_x_log10() + scale_colour_manual(values = c("#D55E00", "#999999", "#0072B2", "black"), name = NULL) +
  labs(x = "half-life, Sorenson (min; capped at 1,440)", y = "cumulative fraction of genes", subtitle = "half-life by %uORF+ group")
th <- theme_bw(base_size = 8) + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = element_blank(), legend.key.size = unit(0.3, "cm"),
  legend.text = element_text(size = 8), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.minor = element_blank())
P <- (cdf + plots$GMUCT + plots$PARE) + plot_layout(ncol = 3) + plot_annotation(tag_levels = "A") & th & theme(plot.tag = element_text(face = "bold", size = 10))
P[[1]] <- P[[1]] + theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0))
ggsave(file.path(OUT_DIR, "A41_Figure_sorenson_halflife.pdf"), P, width = 9, height = 3.4, device = cairo_pdf)
writeLines(.log, file.path(OUT_DIR, "A41_sorenson_summary.txt")); .add("written: output/A41_sorenson_summary.txt, output/A41_Figure_sorenson_halflife.pdf")
