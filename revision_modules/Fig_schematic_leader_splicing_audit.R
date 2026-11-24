################################################################################
# Schematic of the leader-splicing audit (2026-11-07, Larry: 'a figure that explains how you analyze TuORF and translated ouORF
# retention by splicing'). Vector PDF drawn with grid; Arial, >= 7 pt; 6.5 in wide.
#   A  one gene, three transcripts: uORF retained; uORF lost by TSS selection (%uORF+); start codon lost by leader splicing (S)
#   B  the junction-read count at a start codon: R (aligned block covers the position), E (splice gap spans it), S = R / (R + E)
#   C  annotation classes A, B, C1, C2, D (introns of the gene's isoforms relative to the uORF)
# V26 (2026-11-09): panel D states the per-transcript splicing factor (s over the uORFs the transcript carries; the gene-level wording of V23-V25 was stale after the V25 change in splicing_rule.R).
#   D  what is audited for TuORFs (start and stop codons) and for translated overlapping uORFs (start codon only), and how S enters
#      the effective inclusion and the mRNA pool share
# V23 (2026-11-08): %uORF+ defined on the leader window; the junction-consistent rule for a start codon inside an intron of another expressed isoform.
# Output: output/Fig_schematic_leader_splicing_audit.pdf (= Supplemental Fig. S8 of V21, 2026-11-07; the former S8-S21 become S9-S22)      Run from Weak_uORF_submission/:  Rscript revision_modules/Fig_schematic_leader_splicing_audit.R
################################################################################
suppressPackageStartupMessages(library(grid))
FONT <- "Arial"; OUT <- "output/Fig_schematic_leader_splicing_audit.pdf"
gp_t <- function(size = 8, face = "plain", col = "black") gpar(fontfamily = FONT, fontsize = size, fontface = face, col = col, lineheight = 0.95)
COL_UORF <- "#F5D800"; COL_CDS <- "black"; COL_LEADER <- "grey55"; COL_RET <- "#2166AC"; COL_SPL <- "#D6604D"; COL_CAGE <- "#CC66CC"; COL_INTRON <- "grey35"
wrap <- function(s, n) paste(strwrap(s, n), collapse = "\n")
txt <- function(x, y, lab, size = 8, face = "plain", just = "left", col = "black", w = NA) { if (!is.na(w)) lab <- wrap(lab, w); grid.text(lab, x = unit(x, "native"), y = unit(y, "native"), just = just, gp = gp_t(size, face, col)) }
seg <- function(x0, y0, x1, y1, col = "black", lwd = 1, lty = 1, arrow = NULL) grid.segments(unit(x0, "native"), unit(y0, "native"), unit(x1, "native"), unit(y1, "native"), gp = gpar(col = col, lwd = lwd, lty = lty), arrow = arrow)
box <- function(x0, x1, y, h, fill, col = fill) grid.rect(x = unit(x0, "native"), y = unit(y, "native"), width = unit(x1 - x0, "native"), height = unit(h, "native"), just = c("left", "centre"), gp = gpar(fill = fill, col = col, lwd = 0.6))
intron <- function(x0, x1, y, h = 4, col = COL_INTRON) { xm <- (x0 + x1) / 2; seg(x0, y, xm, y + h, col = col, lwd = 0.9); seg(xm, y + h, x1, y, col = col, lwd = 0.9) }
tss <- function(x, y, col = "black") { seg(x, y, x, y + 6, col = col, lwd = 1); seg(x, y + 6, x + 3, y + 6, col = col, lwd = 1, arrow = arrow(length = unit(1.3, "mm"), type = "closed")) }
panel <- function(x, y, w, h, letter) { pushViewport(viewport(x = unit(x, "in"), y = unit(y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top"), xscale = c(0, 100), yscale = c(0, 100))); txt(0.5, 96, letter, 10, "bold") }
endp <- function() popViewport()

cairo_pdf(OUT, width = 6.5, height = 8.5, family = FONT)
grid.newpage()
grid.text("Figure S6. How uORF retention is measured: TSS selection and leader splicing", x = unit(0.02, "in"), y = unit(8.5 - 0.12, "in"), just = c("left", "top"), gp = gp_t(10, "bold"))   # Supplemental Fig. S8 of V21 (2026-11-07)
## ---------------- A: one gene, three transcripts (drawing in x 3-58, text in x 61-100) ------------------------------------
panel(0.0, 8.2, 6.5, 3.1, "A")
txt(4, 96, "One gene: dispersed transcription start sites (TSSs) and a leader intron", 9, "bold")
yg <- 76
seg(4, yg, 58, yg, col = "grey80", lwd = 0.8)
box(5, 22, yg, 3, COL_LEADER); box(9, 15, yg, 5, COL_UORF, "black"); txt(12, yg + 6.5, "uORF", 7, just = "centre")
intron(22, 30, yg, 4); box(30, 38, yg, 3, COL_LEADER); box(38, 56, yg, 7, COL_CDS); txt(47, yg, "mORF", 7, "bold", just = "centre", col = "white")
txt(1, yg, "gene", 7, col = "grey30")
for (k in 1:4) seg(5 + k, yg + 3, 5 + k, yg + 3 + c(3, 6, 4, 2)[k], col = COL_CAGE, lwd = 1.5)
for (k in 1:5) seg(30 + k, yg + 3, 30 + k, yg + 3 + c(2, 6, 9, 5, 2)[k], col = COL_CAGE, lwd = 1.5)
txt(37, yg + 10, "CAGE signal (TSSs)", 6.5, col = COL_CAGE)
seg(9, yg - 60, 9, yg + 5, col = "black", lwd = 0.6, lty = 2); txt(9.5, yg - 58, "uORF start codon", 6.5, col = "grey20")
# transcripts
y1 <- 58; y2 <- 43; y3 <- 28
tr <- function(y, x0, uorf, note, lab) { tss(x0, y + 2.5); box(x0, 38, y, 3, COL_LEADER); if (uorf) box(9, 15, y, 5, COL_UORF, "black"); box(38, 56, y, 7, COL_CDS)
  txt(61, y + 3, lab, 7.5, "bold"); txt(61, y - 4, note, 7, col = "grey20", w = 46) }
tr(y1, 5, TRUE, "starts upstream of the uORF; the leader intron is spliced as annotated: uORF retained", "transcript 1")
box(22, 30, y1, 1, "white", "white"); intron(22, 30, y1, 3.5)
tr(y2, 33, FALSE, "starts downstream of the uORF start codon: uORF absent (TSS selection, measured by %uORF+)", "transcript 2")
box(5, 38, y3, 1, "white", "white"); tss(5, y3 + 2.5); box(5, 7, y3, 3, COL_LEADER); intron(7, 30, y3, 4); box(30, 38, y3, 3, COL_LEADER); box(38, 56, y3, 7, COL_CDS)
txt(61, y3 + 3, "transcript 3", 7.5, "bold"); txt(61, y3 - 4, "starts upstream, but a splice site before the uORF removes its start codon: uORF lost by leader splicing (measured by S)", 7, col = "grey20", w = 46)
txt(4, 11.5, "%uORF+ = percentage of the CAGE signal in the leader window (annotated leader plus 100 nt upstream) at or upstream of the uORF start codon (transcripts 1 and 3 start upstream, transcript 2 downstream).", 7, w = 118)
txt(4, 3.5, "S = fraction of the RNA-seq reads across the start codon that keep it (transcript 1 versus transcript 3), counted per uORF from the junction reads (B).", 7, w = 118)
endp()
## ---------------- B: junction reads at a start codon -----------------------------------------------------------------------
panel(0.0, 5.0, 3.3, 2.35, "B")
txt(5, 96, "Junction reads at the start codon", 9, "bold")
yb <- 76; box(5, 95, yb, 4, COL_LEADER); box(38, 52, yb, 6, COL_UORF, "black"); txt(45, yb + 8, "ATG", 7, just = "centre")
seg(40, 14, 40, yb + 4, col = "black", lwd = 0.6, lty = 2)
rr <- c(62, 54, 46, 38); for (k in seq_along(rr)) { y <- rr[k]; x0 <- c(16, 24, 30, 12)[k]; seg(x0, y, x0 + 36, y, col = COL_RET, lwd = 2.2) }
txt(97, 50, "R = 4 reads: an aligned\nblock covers the codon", 7, col = COL_RET, just = "right")
sr <- c(29, 21); for (k in seq_along(sr)) { y <- sr[k]; seg(6, y, 26, y, col = COL_SPL, lwd = 2.2); seg(26, y, 52, y, col = COL_SPL, lwd = 0.6, lty = 3); seg(52, y, 64, y, col = COL_SPL, lwd = 2.2) }
txt(97, 25, "E = 2 reads: a splice\ngap spans the codon", 7, col = COL_SPL, just = "right")
txt(5, 9, "S = R / (R + E) = 0.67", 8, "bold"); txt(5, 2.5, "counted only where R + E >= 10; S < 0.90 = spliced over", 7)
endp()
## ---------------- C: annotation classes ---------------------------------------------------------------------------------------
panel(3.3, 5.0, 3.2, 2.35, "C")
txt(5, 96, "Annotation classes", 9, "bold"); txt(45, 96, "(introns of all isoforms of the gene)", 7.5, col = "grey20")
cls <- list(list("A", "no leader intron in any isoform", NULL, NULL), list("B", "leader intron elsewhere in the gene", c(8, 20), "own"), list("C1", "intron inside the uORF body, own isoform (start codon kept)", c(38, 46), "own"),
            list("C2", "intron of another isoform alters the body or stop (start kept)", c(44, 60), "other"), list("D", "start codon inside an intron of another isoform", c(22, 40), "other"))
yy <- c(80, 65, 50, 35, 20)
for (k in seq_along(cls)) { y <- yy[k]; c1 <- cls[[k]]
  box(8, 70, y, 2.5, COL_LEADER); box(30, 52, y, 4.5, COL_UORF, "black"); box(70, 95, y, 5.5, COL_CDS)
  if (!is.null(c1[[3]])) intron(c1[[3]][1], c1[[3]][2], y, 4.5, col = if (c1[[4]] == "own") COL_INTRON else COL_SPL)
  txt(1, y, c1[[1]], 7.5, "bold"); txt(8, y - 6, c1[[2]], 6.5, col = "grey20") }
txt(5, 5, "grey caret: intron of the uORF's own isoform; red: of another isoform. Only class D removes the start codon.", 6.5, col = "grey20", w = 70)
endp()
## ---------------- D: TuORFs versus overlapping uORFs, and how S is used ------------------------------------------------------
panel(0.0, 2.55, 6.5, 2.55, "D")
txt(4, 96, "What is audited, and how S enters the estimates", 9, "bold")
yd1 <- 78; yd2 <- 58
box(5, 36, yd1, 3, COL_LEADER); box(12, 22, yd1, 5, COL_UORF, "black"); box(36, 56, yd1, 7, COL_CDS)
seg(12, yd1 - 5, 12, yd1 + 5, lwd = 0.6, lty = 2); seg(22, yd1 - 5, 22, yd1 + 5, lwd = 0.6, lty = 2)
txt(12, yd1 + 8, "start", 6.5, just = "centre"); txt(22, yd1 + 8, "stop", 6.5, just = "centre")
txt(61, yd1 + 3, "TuORF", 7.5, "bold"); txt(61, yd1 - 4, "start and stop codons audited; spliced over if S < 0.90 at either codon", 7, col = "grey20", w = 46)
box(5, 36, yd2, 3, COL_LEADER); box(36, 56, yd2, 7, COL_CDS); box(28, 44, yd2, 5, COL_UORF, "black")
seg(28, yd2 - 5, 28, yd2 + 5, lwd = 0.6, lty = 2); txt(28, yd2 + 8, "start", 6.5, just = "centre"); txt(44, yd2 - 8, "stop inside the mORF", 6.5, just = "centre", col = "grey20")
txt(61, yd2 + 3, "translated overlapping uORF", 7.5, "bold"); txt(61, yd2 - 4, "start codon only, because its stop codon is located in the coding sequence", 7, col = "grey20", w = 46)
txt(4, 41, "Effective inclusion of a uORF = %uORF+ x S   (S = 1 where fewer than 10 reads cover the codon)", 8)
txt(4, 31, "Fraction of the nuclear leader-containing mRNA that contains a detected uORF after both filters", 7.5, "bold")
txt(4, 16, "= sum of TPM x p x s over the transcripts that contain a uORF, divided by the TPM of all leader-containing transcripts. p = %uORF+/100 of the most included uORF of the transcript. s = max(p x S) / max(p) over the uORFs that the transcript contains. When the start codon is located in an intron of another expressed isoform, S is divided by the ratio of the TPM of the isoforms containing the uORF to the summed TPM of both kinds of isoform (the result is limited to at most 1). TSS selection acts through p and leader splicing through s.", 7, w = 118)
endp()
dev.off(); cat("saved", OUT, "\n")
