################################################################################
# A40 — Supplemental Fig. S22 (draft): the uncapped ratio against mRNA half-life #
# and mORF TE (Larry, 2026-11-20: "This part is very confusing and need some    #
# data to backup. Do not just say supplemental Method, put the data in a        #
# figure"). Genes: single-TuORF genes and genes without TuORFs together, with   #
# TE > 0, a half-life (Szabo et al. 2020) and RNA TPM > 0, and an uncapped ratio #
# (GMUCT, Yu et al. 2021; PARE, Nagarajan et al. 2019; the sets of modules A11   #
# and A12).                                                                      #
# (A, B) median uncapped ratio in fifths of half-life (short to long), drawn    #
#   separately for thirds of mORF TE: if the ratio measured the decay rate, it  #
#   would be highest at the shortest half-lives.                                #
# (C, D) median uncapped ratio in fifths of mORF TE, drawn separately for       #
#   thirds of half-life.                                                         #
# Error bars: 95% bootstrap intervals of the median (1,000 resamples).          #
# Panel labels: Spearman rho, and the partial rho with the other measure,       #
#   exon length and RNA TPM held constant (rank-based, as module A39).          #
# Output: output/A40_Figure_uncapped_halflife_TE.pdf, output/A40_summary.txt,   #
#         output/A40_Table_bins.csv                                             #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A40_uncapped_halflife_TE_figure.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
options(warn = 1); set.seed(20261120)
OUT_DIR <- "output"; NBOOT <- 1000L
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A40 uncapped ratio vs half-life and mORF TE — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
a11 <- fread(file.path(OUT_DIR, "A11_Table_genes.csv")); a12 <- fread(file.path(OUT_DIR, "A12_Table_genes.csv")); a12 <- a12[match(a11$gene_id, gene_id)]
len <- fread(file.path(OUT_DIR, "A12_gene_counts.csv"), select = c("gene_id", "exon_length"))
G <- data.table(gene_id = a11$gene_id, TE = a11$mORF_TE, RNA_TPM = a11$RNA_TPM_this_study, HL = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_))
G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]
G <- G[is.finite(TE) & TE > 0 & is.finite(HL) & HL > 0 & is.finite(RNA_TPM) & RNA_TPM > 0 & is.finite(exon_length)]
prho <- function(x, y, Z) { rx <- rank(x); ry <- rank(y); RZ <- as.data.frame(lapply(Z, rank)); cor(resid(lm(rx ~ ., data = RZ)), resid(lm(ry ~ ., data = RZ))) }
medci <- function(x) { b <- replicate(NBOOT, median(sample(x, replace = TRUE))); c(n = length(x), med = median(x), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975))) }
cut_q <- function(x, k) cut(x, unique(quantile(x, seq(0, 1, length.out = k + 1))), include.lowest = TRUE, labels = FALSE)
SETS <- c(U_gmuct = "GMUCT", U_pare = "PARE")
bins <- list(); plots <- list(); labs <- list(); dg <- position_dodge(width = 0.3)
for (v in names(SETS)) {
  d <- G[is.finite(get(v)) & get(v) > 0]; d[, U := get(v)]
  d[, hl5 := cut_q(HL, 5)]; d[, te5 := cut_q(TE, 5)]
  d[, te3 := factor(cut_q(TE, 3), 1:3, c("low TE", "middle TE", "high TE"))]; d[, hl3 := factor(cut_q(HL, 3), 1:3, c("short half-life", "middle half-life", "long half-life"))]
  r_hl <- cor(d$U, d$HL, method = "spearman"); r_te <- cor(d$U, d$TE, method = "spearman")
  p_hl <- prho(d$U, d$HL, list(d$TE, d$exon_length, d$RNA_TPM)); p_te <- prho(d$U, d$TE, list(d$HL, d$exon_length, d$RNA_TPM))
  .add(sprintf("  rho(half-life, mORF TE) in the %s set: %.2f", SETS[[v]], cor(d$HL, d$TE, method = "spearman")))
  .add(sprintf("%s: n = %s genes; rho with half-life %.2f (partial, TE, exon length and RNA TPM constant: %.2f); rho with mORF TE %.2f (partial, half-life, exon length and RNA TPM constant: %.2f)",
               SETS[[v]], format(nrow(d), big.mark = ","), r_hl, p_hl, r_te, p_te))
  A <- d[, as.list(medci(U)), by = .(stratum = te3, bin = hl5)][, x_med := d[, median(HL), by = hl5][match(bin, hl5), V1]]
  C <- d[, as.list(medci(U)), by = .(stratum = hl3, bin = te5)][, x_med := d[, median(TE), by = te5][match(bin, te5), V1]]
  bins[[v]] <- rbind(cbind(set = SETS[[v]], panel = "by half-life fifth within TE third", A), cbind(set = SETS[[v]], panel = "by TE fifth within half-life third", C))
  hl_lab <- sprintf("%.1f", sort(unique(A$x_med))); te_lab <- sprintf("%.2f", sort(unique(C$x_med)))
  pa <- ggplot(A, aes(factor(bin), med, colour = stratum, group = stratum)) + geom_line(linewidth = 0.5, position = dg) + geom_point(size = 1.4, position = dg) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.35, position = dg) + scale_x_discrete(labels = hl_lab) +
    scale_colour_manual(values = c("#9ecae1", "#4292c6", "#08306b"), name = NULL) +
    labs(x = "mRNA half-life (h; median of each fifth)", y = sprintf("uncapped ratio (%s), median", SETS[[v]]),
         subtitle = sprintf("%s: rho with half-life %.2f\nTE, length, abundance held constant: %.2f", SETS[[v]], r_hl, p_hl))
  pc <- ggplot(C, aes(factor(bin), med, colour = stratum, group = stratum)) + geom_line(linewidth = 0.5, position = dg) + geom_point(size = 1.4, position = dg) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.35, position = dg) + scale_x_discrete(labels = te_lab) +
    scale_colour_manual(values = c("#fdae6b", "#e6550d", "#7f2704"), name = NULL) +
    labs(x = "mORF TE (median of each fifth)", y = sprintf("uncapped ratio (%s), median", SETS[[v]]),
         subtitle = sprintf("%s: rho with mORF TE %.2f\nhalf-life, length, abundance held constant: %.2f", SETS[[v]], r_te, p_te))
  plots[[v]] <- list(pa, pc)
}
th <- theme_bw(base_size = 8) + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = element_blank(),
  legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.minor = element_blank())
P <- (plots$U_gmuct[[1]] + plots$U_pare[[1]]) / (plots$U_gmuct[[2]] + plots$U_pare[[2]]) + plot_annotation(tag_levels = "A") & th & theme(plot.tag = element_text(face = "bold", size = 10))
P[[1]][[1]] <- P[[1]][[1]] + theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0))   # panel A: key at the bottom right, where it covers no data
# panel order: A = GMUCT by half-life, B = PARE by half-life, C = GMUCT by TE, D = PARE by TE
ggsave(file.path(OUT_DIR, "A40_Figure_uncapped_halflife_TE.pdf"), P, width = 6.5, height = 6.2, device = cairo_pdf)
fwrite(rbindlist(bins), file.path(OUT_DIR, "A40_Table_bins.csv"))
writeLines(.log, file.path(OUT_DIR, "A40_summary.txt")); .add("written: output/A40_Figure_uncapped_halflife_TE.pdf, output/A40_summary.txt, output/A40_Table_bins.csv")
