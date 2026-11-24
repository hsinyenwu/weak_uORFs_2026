# Revised Figure 3C: three abundance-weighted "ladders" (gene count -> mRNA mass -> after TSS weighting) from the
# transcriptome-impact tables in output/ (CAGE), with the TSS-seq endpoint in the legend. Reviewer 2 §2.6.
# LADDER_SOURCE=TSSseq draws the same ladders from the TSS-seq table (Supplemental Fig. S13 of the revision).
# V20 (2026-11-06, Larry): the table gains two rows with the shares after TSS weighting and leader splicing (module A26: the TSS-weighted share minus the
# splicing loss, i.e. the mRNA whose uORF start codon is removed by leader splicing), for detected TuORFs and for TuORFs + overlapping uORFs.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(gridExtra); library(grid)})
OUT <- "output_revised_figs_20261028"
# V23 (2026-11-08): every value comes from module A26 (output/A26_Table_prevalence_ladder.csv), which computes the three ladders, the shares after
# TSS weighting and leader splicing and the TPM table from one per-transcript table per TSS source, with one leader definition for numerator and
# denominator and the same inclusion rule for every set (before: the canonical metric table, whose TuORF rows required %uORF+ > 0 and whose
# numerators included 20 carriers outside the gene-list denominator).
A26 <- read.csv("output/A26_Table_prevalence_ladder.csv", check.names = FALSE)
SRC <- Sys.getenv("LADDER_SOURCE", "CAGE")      # "CAGE" = Figure 3C; "TSSseq" = Supplemental Fig. S13 (the same ladders from the TSS-seq table)
if (SRC == "TSSseq") { OTHER_SRC <- "CAGE"; OTHER <- "CAGE-seq"; TTL <- "Figure S13. uORF prevalence in the steady-state transcriptome estimated with TSS-seq data"; FN <- "FigureS13_revised_prevalence_ladder_TSSseq.pdf"; VALS <- "FigureS13_revised_prevalence_ladder_TSSseq_values.csv"; SRC_LAB <- "TSS-seq"
} else { OTHER_SRC <- "TSSseq"; OTHER <- "TSS-seq"; TTL <- "C"; FN <- "Figure3_revised_panel_C_prevalence_ladder.pdf"; VALS <- "Figure3_revised_panel_C_values.csv"; SRC_LAB <- "CAGE" }
SETS <- c("Detected TuORFs", "Detected TuORFs + overlapping uORFs", "Any predicted AUG uORF")
STEPS <- c("% of genes", "% of nuclear mRNA (unweighted)", "% of nuclear mRNA after TSS weighting")
val <- function(src, set, step) { v <- A26$value[A26$source == src & A26$set == set & A26$step == step]; stopifnot(length(v) == 1); v }
tpm <- function(src, set, step) { v <- A26$TPM[A26$source == src & A26$set == set & A26$step == step]; stopifnot(length(v) == 1); v }
lad <- data.frame(set = factor(rep(SETS, each = 3), levels = rev(SETS)), step = factor(rep(c("% of 5'UTR genes", "% of nuclear mRNA (unweighted)", "% of nuclear mRNA after TSS weighting"), 3), levels = c("% of 5'UTR genes", "% of nuclear mRNA (unweighted)", "% of nuclear mRNA after TSS weighting")),
                  value = unlist(lapply(SETS, function(s) sapply(STEPS, function(st) val(SRC, s, st)))))
tss_end <- sapply(SETS, function(s) val(OTHER_SRC, s, "% of nuclear mRNA after TSS weighting"))
stopifnot(!any(is.na(lad$value)))
print(lad)
p <- ggplot(lad, aes(x = step, y = value, fill = set)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", value)), position = position_dodge(0.75), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("Any predicted AUG uORF" = "#BDBDBD", "Detected TuORFs + overlapping uORFs" = "#74A9CF", "Detected TuORFs" = "#2166AC"), name = NULL) +   # the other source's endpoints are given in the subtitle
  scale_x_discrete(labels = c("% of genes\nwith the ORF", "% of nuclear mRNA\n(unweighted)", "% of nuclear mRNA\nafter TSS weighting")) +
  coord_cartesian(ylim = c(0, 65)) + theme_classic(base_size = 11) + labs(title = TTL, subtitle = sprintf("the same three fractions after TSS weighting with %s: any predicted AUG uORF %.1f%%, TuORFs + overlapping uORFs %.1f%%, TuORFs %.1f%%", OTHER, tss_end[3], tss_end[2], tss_end[1]), x = NULL, y = "Percent") +
  theme(plot.title = element_text(face = "bold", size = 16), legend.position = c(0.68, 0.85), legend.text = element_text(size = 8.5))
spl_step <- "% of nuclear mRNA after TSS weighting and leader splicing"
tab <- data.frame(Category = c("Plastid encoded", "Mitochondria encoded", "Nuclear encoded, no 5′ leader", "Nuclear encoded, with 5′ leader", "  with detected TuORFs (unweighted)", "  with TuORFs, TSS-weighted", "  with TuORFs, TSS-weighted and splicing-adjusted", "  with TuORFs or overlapping uORFs, TSS-weighted", "  with TuORFs or overlapping uORFs, TSS-weighted and splicing-adjusted", "  with any predicted AUG uORF, TSS-weighted"),
                  TPM = c(tpm(SRC, "pool", "Plastid encoded"), tpm(SRC, "pool", "Mitochondria encoded"), tpm(SRC, "pool", "Nuclear encoded, no 5' leader"), tpm(SRC, "pool", "Nuclear encoded, with 5' leader"),
                          tpm(SRC, SETS[1], STEPS[2]), tpm(SRC, SETS[1], STEPS[3]), tpm(SRC, SETS[1], spl_step), tpm(SRC, SETS[2], STEPS[3]), tpm(SRC, SETS[2], spl_step), tpm(SRC, SETS[3], STEPS[3])),
                  `% of nuclear mRNA with 5′ leader` = c("", "", "", "100", sprintf("%.1f", val(SRC, SETS[1], STEPS[2])), sprintf("%.1f", val(SRC, SETS[1], STEPS[3])), sprintf("%.1f", val(SRC, SETS[1], spl_step)),
                                                         sprintf("%.1f", val(SRC, SETS[2], STEPS[3])), sprintf("%.1f", val(SRC, SETS[2], spl_step)), sprintf("%.1f", val(SRC, SETS[3], STEPS[3]))), check.names = FALSE)
tab$TPM <- format(round(as.numeric(tab$TPM)), big.mark = ",")
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28): plot above the table; plot and table also saved alone for the Figure 3 composite
tg <- tableGrob(tab, rows = NULL, theme = ttheme_minimal(base_size = 8, base_family = FIG_FONT, core = list(fg_params = list(hjust = 0, x = 0.02)), colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))
pc <- legend_in(compact_wide(p), 0.995, 0.99)   # V30 (2026-11-09, Larry): legend inside the plot, in the free room above the lowest bars
cairo_pdf(file.path(OUT, FN), width = FIG_W, height = 6.9, family = FIG_FONT)   # V20: two more table rows
grid.arrange(pc, tg, ncol = 1, heights = c(1.25, 1.1))
dev.off()
fig_save(file.path(OUT, sub(".pdf", "_plot.pdf", FN, fixed = TRUE)), pc, w = 6.5, h = 2.9)   # full-width row of the Figure 3 composite
cairo_pdf(file.path(OUT, sub(".pdf", "_table.pdf", FN, fixed = TRUE)), width = 6.5, height = 2.8, family = FIG_FONT); grid.draw(tg); dev.off()   # V20: 10 rows   # full-width row of the Figure 3 composite
write.csv(lad, file.path(OUT, VALS), row.names = FALSE)
cat("saved", FN, "\n")
