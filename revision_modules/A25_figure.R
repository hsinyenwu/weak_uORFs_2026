################################################################################
# A25_figure.R — Supplemental Fig. S3 (translated overlapping uORFs versus TuORFs): sourced by A25 at its end, or run alone:
#   A25_FIGURE_ONLY=1 Rscript revision_modules/A25_figure.R   (loads output/A2_Table_TuORF_leader_splicing.csv and output/A25_Table_ouORF_leader_splicing.csv)
# V22 (2026-11-07): panel subtitles wrapped so nothing is clipped at the page edge; tick labels of D on two rows; shorter strip labels.
################################################################################
if (!exists("O")) {   # stand-alone mode
  suppressPackageStartupMessages(library(data.table)); OUT_DIR <- "output"; t0 <- Sys.time(); .log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
  A2 <- fread(file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.csv")); O <- fread(file.path(OUT_DIR, "A25_Table_ouORF_leader_splicing.csv")); O[, pct_uORF_plus := pct_uORF_plus_CAGE] }
suppressPackageStartupMessages(library(ggplot2)); source("revision_modules/fig_compact.R")
tu <- A2[, .(pct_uORF_plus, S_TIS, cov_TIS, annot_class, kind = "TuORFs (n = 7,423)")]
ou <- as.data.table(O)[, .(pct_uORF_plus, S_TIS, cov_TIS, annot_class, kind = sprintf("Overlapping uORFs (n = %s)", format(nrow(O), big.mark = ",")))]
both <- rbind(tu, ou); both[, kind := factor(kind, levels = c("TuORFs (n = 7,423)", unique(ou$kind)))]
colk <- setNames(c("#2166AC", "#E66101"), levels(both$kind))
pA <- ggplot(both, aes(pct_uORF_plus, colour = kind)) + stat_ecdf(linewidth = 0.6) + scale_colour_manual(values = colk, name = NULL) +
  labs(x = "%uORF+ (CAGE-seq)", y = "Cumulative fraction of ORFs", title = "A  Inclusion by TSS selection",
       subtitle = sprintf("median %%uORF+: TuORFs %.0f%%, overlapping uORFs %.0f%%; <= 50%% uORF+: %.0f%% and %.0f%%", median(tu$pct_uORF_plus, na.rm = TRUE), median(ou$pct_uORF_plus, na.rm = TRUE), 100 * mean(tu$pct_uORF_plus <= 50, na.rm = TRUE), 100 * mean(ou$pct_uORF_plus <= 50, na.rm = TRUE))) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
cls <- both[, .(n = .N), by = .(kind, cls = substr(annot_class, 1, 2))][, pct := 100 * n / sum(n), by = kind]
cls[, cls := factor(sub("\\.$", "", cls), levels = c("A", "B", "C1", "C2", "D"))]
pB <- ggplot(cls, aes(cls, pct, fill = kind)) + geom_col(position = position_dodge(0.7), width = 0.65) + scale_fill_manual(values = colk, name = NULL) +
  geom_text(aes(label = sprintf("%.0f", pct)), position = position_dodge(0.7), vjust = -0.3, size = 7 / .pt, family = FIG_FONT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = "Annotation class", y = "% of ORFs", title = "B  Annotated leader introns", subtitle = "A none; B elsewhere; C1/C2 in the ORF body; D start codon in an intron of another isoform") +
  theme_bw(base_size = 9) + theme(legend.position = "none", panel.grid.minor = element_blank())
cov <- both[!is.na(S_TIS)]
pC <- ggplot(cov, aes(S_TIS, colour = kind)) + stat_ecdf(linewidth = 0.6) + scale_colour_manual(values = colk, name = NULL) + geom_vline(xintercept = 0.9, linetype = 3, linewidth = 0.3) +
  labs(x = "Fraction of RNA-seq reads keeping the start codon (S)", y = "Cumulative fraction of ORFs", title = "C  Leader splicing at the start codon",
       subtitle = sprintf("ORFs with >= 10 reads: TuORFs %s (%.0f%%), overlapping uORFs %s (%.0f%%); S < 0.90: %.0f%% and %.0f%%", format(sum(!is.na(tu$S_TIS)), big.mark = ","), 100 * mean(!is.na(tu$S_TIS)), sum(!is.na(ou$S_TIS)), 100 * mean(!is.na(ou$S_TIS)),
                         100 * mean(tu$S_TIS[!is.na(tu$S_TIS)] < 0.9), 100 * mean(ou$S_TIS[!is.na(ou$S_TIS)] < 0.9))) +
  theme_bw(base_size = 9) + theme(legend.position = "none", panel.grid.minor = element_blank())
cuts <- c(0, 10, 50, 80); lab <- c("> 0%", ">= 10%", ">= 50%", ">= 80%")
shr <- rbindlist(lapply(split(both, both$kind), function(d) {
  p <- d$pct_uORF_plus; s <- d$S_TIS; ok <- !is.na(s); lost <- (p[ok] / 100) * (1 - s[ok])
  rbind(data.table(kind = d$kind[1], mechanism = "TSS selection", cut = lab, pct = sapply(cuts, function(k) 100 * mean(if (k == 0) 100 - p > 1e-9 else 100 - p >= k - 1e-9, na.rm = TRUE))),
        data.table(kind = d$kind[1], mechanism = "Leader splicing", cut = lab, pct = sapply(cuts, function(k) 100 * mean(if (k == 0) 100 * lost > 1e-9 else 100 * lost >= k - 1e-9))))
}))
shr[, `:=`(cut = factor(cut, levels = lab), mechanism = factor(mechanism, levels = c("TSS selection", "Leader splicing")))]
pD <- ggplot(shr, aes(cut, pct, fill = kind)) + geom_col(position = position_dodge(0.7), width = 0.65) + facet_wrap(~ mechanism) + scale_fill_manual(values = colk, name = NULL) +
  geom_text(aes(label = sprintf("%.0f", pct)), position = position_dodge(0.7), vjust = -0.3, size = 7 / .pt, family = FIG_FONT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = "Percentage of the gene's transcripts that lack the ORF", y = "% of ORFs", title = "D  ORFs excluded from part of the gene's transcripts") + scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_bw(base_size = 9) + theme(legend.position = "none", panel.grid.minor = element_blank(), strip.background = element_blank())
fwrite(shr, file.path(OUT_DIR, "A25_Table_ouORF_vs_TuORF_cutoffs.csv"))
g <- fig_grid(list(legend_in(compact(pA), 0.03, 0.99), compact(pB), free_xlab(compact(pC)), compact(pD)), ncol = 2, h = 6.2,   # V30 (2026-11-09, Larry): the legend of A inside the plot; aligned plotting areas (fig_grid)   # V22: compact() wraps titles/subtitles for half-width panels (compact_wide() is for full-width panels and let the text run off the page)
  top = fig_title("Figure S5. Translated overlapping uORFs: TSS selection and leader splicing"))   # V21 numbering (2026-11-07; S16 until the schematic of the leader-splicing audit became S8)
fig_save(file.path(OUT_DIR, "A25_Figure_ouORF_vs_TuORF.pdf"), g, h = 6.2)
.add(sprintf("  figure: %s", file.path(OUT_DIR, "A25_Figure_ouORF_vs_TuORF.pdf")))
for (k in seq_len(nrow(shr))) .add(sprintf("  cutoffs: %-28s %-40s %-7s %5.1f%%", shr$kind[k], shr$mechanism[k], shr$cut[k], shr$pct[k]))
