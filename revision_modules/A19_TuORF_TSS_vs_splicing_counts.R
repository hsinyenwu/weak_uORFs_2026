################################################################################
# A19 — How many translated uORFs are affected by TSS selection and how many   #
#        by leader splicing (Larry, 2026-11-01; Reviewer 2 §2.1 follow-up)      #
#                                                                              #
# Per TuORF (from the A2 audit table):                                         #
#   share of the gene's transcripts that start downstream of the uORF          #
#       = 1 - p,            p = %uORF+/100 (CAGE)                              #
#   share that start upstream but lose the start codon by splicing             #
#       = p x (1 - S_TIS),  S_TIS from the RNA-seq junction reads              #
#   Splicing can only be measured for TuORFs with >= 10 reads across the start  #
#   codon (n = 5,828 of 7,423), so the TSS share is given for all TuORFs and,   #
#   like for like, for the same 5,828.                                         #
# Outputs (output/, prefix A19_): the cutoff table (csv), a summary txt with   #
#   the per-uORF averages (the 73%/27% split of the Results), and the table    #
#   panel drawn for the main figure (PDF, 6.5 in wide, Arial 8 pt).            #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A19_TuORF_TSS_vs_splicing_counts.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(gridExtra); library(grid) })
source("revision_modules/fig_compact.R")
OUT <- "output"; MIN_READS <- 10
a2 <- fread(file.path(OUT, "A2_Table_TuORF_leader_splicing.csv"))
a2[, p := pct_uORF_plus / 100]
m  <- a2[!is.na(S_TIS) & cov_TIS >= MIN_READS]
m[, `:=`(loss_tss = 1 - p, loss_spl = p * (1 - S_TIS))]
a2[, loss_tss := 1 - p]

cuts <- c(0, 0.10, 0.50, 0.80)
lab  <- c("more than 0%", "at least 10%", "at least 50%", "at least 80%")
sh <- function(x, k) 100 * mean(if (k == 0) x > 1e-12 else x >= k - 1e-9)   # V23: tolerance for floating-point boundaries
tab <- data.table(cutoff = lab,
  TSS_all      = sapply(cuts, function(k) sh(a2$loss_tss, k)),
  TSS_measured = sapply(cuts, function(k) sh(m$loss_tss, k)),
  splicing_measured = sapply(cuts, function(k) sh(m$loss_spl, k)))
tab[, `:=`(n_all = nrow(a2), n_measured = nrow(m))]
fwrite(tab, file.path(OUT, "A19_Table_TuORF_affected_TSS_vs_splicing.csv"))

s <- c(
  sprintf("A19 — translated uORFs affected by TSS selection and by leader splicing (%s)", format(Sys.Date())),
  sprintf("TuORFs: %d; with >= %d junction reads across the start codon (splicing measurable): %d", nrow(a2), MIN_READS, nrow(m)),
  "", "Share of TuORFs missing from the given share of their gene's transcripts:",
  sprintf("  %-18s TSS (all %d): %5.1f%%   TSS (same %d): %5.1f%%   splicing (%d): %5.1f%%", tab$cutoff, nrow(a2), tab$TSS_all, nrow(m), tab$TSS_measured, nrow(m), tab$splicing_measured),
  "", "Per-uORF averages over the measured TuORFs (the split quoted in the Results):",
  sprintf("  transcripts lacking the uORF: %.1f%% = %.1f%% start downstream of it (TSS) + %.1f%% lose the start codon by splicing", 100 * mean(1 - m$p * m$S_TIS), 100 * mean(m$loss_tss), 100 * mean(m$loss_spl)),
  sprintf("  of the transcripts that lack the uORF, %.0f%% lack it because of the TSS and %.0f%% because of splicing", 100 * mean(m$loss_tss) / mean(1 - m$p * m$S_TIS), 100 * mean(m$loss_spl) / mean(1 - m$p * m$S_TIS)),
  sprintf("  TuORFs losing the start codon by splicing in >= 10%% of transcripts: %d of %d (%.1f%%)", sum(m$loss_spl >= 0.10), nrow(m), 100 * mean(m$loss_spl >= 0.10)),
  sprintf("  share of the summed splicing loss carried by class D (start codon inside an intron of another isoform): %.0f%% (%d TuORFs)", 100 * sum(m[grepl("^D", annot_class)]$loss_spl) / sum(m$loss_spl), sum(grepl("^D", m$annot_class))),
  sprintf("  TuORFs with both TSS loss >= 10%% and splicing loss >= 10%%: %.1f%%", 100 * mean(m$loss_tss >= 0.10 & m$loss_spl >= 0.10)),
  "", "Annotation (all TuORFs):",
  sprintf("  start codon inside an intron of another isoform (class D): %.1f%%; body or stop codon altered (C1 + C2): %.1f%%; either: %.1f%%", 100 * mean(grepl("^D", a2$annot_class)), 100 * mean(grepl("^C", a2$annot_class)), 100 * mean(grepl("^[CD]", a2$annot_class))))
writeLines(s, file.path(OUT, "A19_summary.txt")); cat(s, sep = "\n")

## table panel for the main figure (Figure 2G in the V11 draft): 6.5 in wide, Arial, 8 pt
pt <- data.frame(
  a = tab$cutoff,
  b = sprintf("%.0f%%", tab$TSS_all),
  c = sprintf("%.0f%%", tab$TSS_measured),
  d = sprintf("%.0f%%", tab$splicing_measured), check.names = FALSE)
names(pt) <- c("Share of the gene's\ntranscripts that lack\nthe uORF",                                     # <= 24 characters per line at 8 pt bold in a 1.45-1.65 in column. V54 (2026-11-23, Larry: "Figure 2G label should not have ..."; he chose this label): was "Translated uORF excluded from ... of its gene's transcripts"
               sprintf("Transcription starts\ndownstream of the uORF\n(TSS selection),\nall uORFs (n = %s)", format(nrow(a2), big.mark = ",")),
               sprintf("Transcription starts\ndownstream of the uORF,\nuORFs with ≥ 10 reads\n(n = %s)", format(nrow(m), big.mark = ",")),
               sprintf("Leader splicing removes\nthe start codon,\nuORFs with ≥ 10 reads\n(n = %s)", format(nrow(m), big.mark = ",")))
th <- ttheme_minimal(base_size = 8, base_family = FIG_FONT,
  core = list(fg_params = list(hjust = 0, x = 0.04), padding = unit(c(5, 3.2), "mm")),
  colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold", lineheight = 0.9), padding = unit(c(5, 2.5), "mm")))
tg <- tableGrob(pt, rows = NULL, theme = th)
tg$widths <- unit(c(1.45, 1.65, 1.65, 1.65), "in")
tg <- gtable::gtable_add_grob(tg, grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 0.8)), t = 1, b = 1, l = 1, r = ncol(tg))   # rule under the header
tg <- gtable::gtable_add_grob(tg, grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 0.8)), t = nrow(tg), b = nrow(tg), l = 1, r = ncol(tg))   # rule under the last row
# V18 (2026-11-05, Larry): a note under the table says what TSS selection means and what the 5,828 TuORFs with junction reads are
note <- paste(strwrap(sprintf(paste0("TSS selection: the transcript starts downstream of the uORF and therefore lacks it. Leader splicing can be measured only where at least 10 ",
  "RNA-seq reads cover the uORF start codon (%s of the %s TuORFs). The middle column shows TSS selection for these same uORFs. Read coverage ",
  "favors uORFs that are present in many transcripts, so TSS exclusion is less common among them."), format(nrow(m), big.mark = ","), format(nrow(a2), big.mark = ",")), width = 108), collapse = "\n")
NOTE_H <- 0.62; H <- 1.9 + NOTE_H
cairo_pdf(file.path(OUT, "A19_Figure2G_TSS_vs_splicing_table.pdf"), width = FIG_W, height = H, family = FIG_FONT)
grid.newpage()
grid.text("D", x = unit(0.05, "in"),   # V54 (2026-11-23, Larry: "Reletter"): the table is cited before the plots of %uORF+ versus uIF and TE, so it is Figure 2D (was 2G); the file keeps its name
 y = unit(H - 0.08, "in"), just = c("left", "top"), gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FIG_FONT))
pushViewport(viewport(x = unit(0.08, "in"), y = unit(H - 0.1, "in"), width = unit(6.4, "in"), height = unit(H - 0.2 - NOTE_H, "in"), just = c("left", "top")))
grid.draw(tg); popViewport()
grid.text(note, x = unit(0.12, "in"), y = unit(NOTE_H - 0.04, "in"), just = c("left", "top"), gp = gpar(fontsize = 8, lineheight = 1.05, fontfamily = FIG_FONT))
dev.off()
cat("saved", file.path(OUT, "A19_Figure2G_TSS_vs_splicing_table.pdf"), "\n")
