# fig_S10_S12_compact.R — Supplemental Figs. S9 (genes with two, three or four TuORFs: %uORF+ and uIF by uORF position) and   [the file name keeps old numbers]
# S10 (undetected uORFs) on 6.5-in pages (Larry, 2026-11-01). Sourced by the figure driver (build_revised_figure_driver.py)
# after the canonical FIGURE 11 block and the undetected-uORF module (canonical FIGURE 4), whose objects it reuses:
#   S10_rows      list (N = 2, 3, 4) of n_uorfs, data, pUP_ecdf, pUP_box, uIF_ecdf, uIF_box, letters_pUP, letters_uIF
#                 (letters: pairwise Wilcoxon rank-sum tests, Benjamini-Hochberg, p < 0.05; the panels of the submitted Figure S3)
#   p_S5D_hist    %uORF+ histogram of the undetected-uORF genes (panel D)
#   TE_Uu_Final_4, letters_TE_4, letters_RNA_4, countsBC_4   mORF TE and RNA of the five gene groups (panel E; letters: pairwise
#                 Wilcoxon rank-sum tests, Benjamini-Yekutieli, p < 0.05)
# S10 is drawn entirely in R (vector PDF). S12 = the browser panels A-C cropped from the submitted figure file (raster, 300 dpi)
# above the regenerated D and E (vector). Both files are written to OUT_DIR.
# Standalone (layout work without rerunning the driver): Rscript revision_modules/fig_S10_S12_compact.R DRIVER_OUT_DIR
# loads DRIVER_OUT_DIR/S10_S12_objects.rds, which the driver writes before sourcing this file.
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra); library(png) })
if (!exists("S10_rows")) {
  a <- commandArgs(trailingOnly = TRUE); OUT_DIR <- if (length(a)) a[1] else "_to_delete/driver_compact_20261028"
  source("revision_modules/fig_compact.R"); o <- readRDS(file.path(OUT_DIR, "S10_S12_objects.rds")); for (nm in names(o)) assign(nm, o[[nm]]) }
SUBMITTED_FIGS <- "/Users/wu/Desktop/weak_uORF_revision/Submit_Genome Research/Figures 20260826"
gp_txt <- function(size, bold = FALSE) gpar(fontsize = size, fontface = if (bold) "bold" else "plain", fontfamily = FIG_FONT)
pal_order <- function(N) RColorBrewer::brewer.pal(max(3, N), "Dark2")[seq_len(N)]     # one palette for all four panels of a row

## ---- Supplemental Fig. S9 ------------------------------------------------------------------------------------------------
s10_header <- function(N, cols, letter) {                 # panel letter, row title, uORF-order schematic, column headers (0.65 in tall)
  x0 <- 2.55; y <- 0.33; bw <- 0.17; gap <- 0.05; xm <- x0 + N * (bw + gap)
  g <- gList(textGrob(letter, x = unit(0.06, "in"), y = unit(0.63, "in"), just = c("left", "top"), gp = gp_txt(10, TRUE)),
             textGrob(sprintf("genes with %d uORFs", N), x = unit(3.25, "in"), y = unit(0.63, "in"), just = c("centre", "top"), gp = gp_txt(9, TRUE)),
             linesGrob(x = unit(c(x0 - 0.2, xm + 0.85), "in"), y = unit(c(y, y), "in"), gp = gpar(lwd = 1.2)),
             textGrob("uORFs", x = unit(x0 - 0.27, "in"), y = unit(y, "in"), just = c("right", "centre"), gp = gp_txt(8)),
             rectGrob(x = unit(xm, "in"), y = unit(y, "in"), width = unit(0.7, "in"), height = unit(0.12, "in"), just = "left", gp = gpar(fill = "black", col = NA)),
             textGrob("mORF", x = unit(xm + 0.35, "in"), y = unit(y - 0.09, "in"), just = c("centre", "top"), gp = gp_txt(8)),
             textGrob("%uORF+", x = unit(1.63, "in"), y = unit(0.02, "in"), just = c("centre", "bottom"), gp = gp_txt(9)),
             textGrob("uIF", x = unit(4.88, "in"), y = unit(0.02, "in"), just = c("centre", "bottom"), gp = gp_txt(9)))
  for (i in seq_len(N)) { xi <- x0 + (i - 1) * (bw + gap)
    g <- gList(g, rectGrob(x = unit(xi, "in"), y = unit(y, "in"), width = unit(bw, "in"), height = unit(0.1, "in"), just = "left", gp = gpar(fill = cols[i], col = NA)),
               textGrob(c("1st", "2nd", "3rd", "4th")[i], x = unit(xi + bw / 2, "in"), y = unit(y - 0.09, "in"), just = c("centre", "top"), gp = gp_txt(8))) }
  gTree(children = g) }
s10_row <- function(r) {
  N <- r$n_uorfs; d <- r$data; cols <- pal_order(N)
  med1 <- function(v) stats::median(d[[v]][d$Rank_Label == "1st"], na.rm = TRUE)   # the red dashed lines of the legend
  notitle <- function(p) p + labs(title = NULL, subtitle = NULL)
  small_letters <- function(p) { for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomText")) { p$layers[[i]]$aes_params$size <- 3; p$layers[[i]]$aes_params$fontface <- "italic" }; p }
  legend_in <- function(p) p + theme(legend.position = "inside", legend.position.inside = c(0.04, 0.99), legend.justification.inside = c(0, 1),
                                     legend.background = element_blank(), legend.key.size = unit(0.25, "cm"), legend.margin = margin(0, 0, 0, 0))
  p1 <- suppressMessages(notitle(r$pUP_ecdf) + scale_color_manual(values = cols, name = "Order") + labs(x = "%uORF+ (CAGE)", y = "Fraction"))
  p2 <- suppressMessages(small_letters(notitle(r$pUP_box)) + scale_fill_manual(values = cols) + labs(x = NULL, y = "%uORF+ (CAGE)") +
          geom_hline(yintercept = med1("pUP"), linetype = "dashed", colour = "red", linewidth = 0.4))
  p3 <- suppressMessages(notitle(r$uIF_ecdf) + scale_color_manual(values = cols, name = "Order") + labs(x = "log2(uIF)", y = "Fraction"))
  p4 <- suppressMessages(small_letters(notitle(r$uIF_box)) + scale_fill_manual(values = cols) + labs(x = NULL, y = "log2(uIF)") +
          geom_hline(yintercept = med1("log2_uIF"), linetype = "dashed", colour = "red", linewidth = 0.4))
  pl <- list(legend_in(compact(p1, geom_scale = 0.4)), compact(p2, geom_scale = 0.4), legend_in(compact(p3, geom_scale = 0.4)), compact(p4, geom_scale = 0.4))
  panels <- ggpubr::ggarrange(plotlist = pl, ncol = 4, nrow = 1, align = "h")
  arrangeGrob(s10_header(N, cols, LETTERS[N - 1]), ggplotGrob(panels), ncol = 1, heights = unit(c(0.65, 1.9), "in")) }
rows10 <- lapply(S10_rows, s10_row)
S10 <- arrangeGrob(grobs = c(list(fig_title("Figure S9. Distal uORFs are excluded more often than proximal uORFs in multi-uORF genes, as their position predicts")), rows10),
                   ncol = 1, heights = unit(c(0.3, rep(2.55, length(rows10))), "in"))
fig_save(file.path(OUT_DIR, "FigureS11_revised_multi_uORF_genes.pdf"), S10, h = 0.3 + 2.55 * length(rows10))
for (r in S10_rows) message(sprintf("S10: %d-uORF genes n = %d | letters %%uORF+: %s | letters uIF: %s", r$n_uorfs, dplyr::n_distinct(r$data$gene_id),
                                    paste(r$letters_pUP$Rank_Label, r$letters_pUP$Letters, collapse = " "), paste(r$letters_uIF$Rank_Label, r$letters_uIF$Letters, collapse = " ")))

## ---- Supplemental Fig. S10 ------------------------------------------------------------------------------------------------
S12_SUB <- file.path(SUBMITTED_FIGS, "FigS5_UuORFs_v20260507.pdf"); wk <- file.path(OUT_DIR, "S11_work"); dir.create(wk, showWarnings = FALSE)
SC <- 300 / 72                                              # the submitted page is 468 x 544 pt; A-C occupy the top 414 pt (D and E start at 418)
system2("pdftoppm", c("-r", 300, "-x", 0, "-y", 0, "-W", round(468 * SC), "-H", round(414 * SC), "-png", "-singlefile", shQuote(S12_SUB), shQuote(file.path(wk, "S12_ABC"))))
abc <- readPNG(file.path(wk, "S12_ABC.png")); h_abc <- 6.5 * dim(abc)[1] / dim(abc)[2]
# panel E: the five gene groups in the order of the submitted figure, labelled as in the manuscript ("Other" = genes without a TuORF or
# a detected overlapping uORF that are not in the undetected-uORF set)
grp_name <- c("Other" = "No TuORF", "0~30" = "Undetected, <30%", "30~80" = "Undetected, 30–80%", "80~100" = "Undetected, 80–100%", "Detected uORFs" = "TuORF")   # V40 (2026-11-11): the abbreviation UuORF is written out. The percentage is the %uORF+ group; the longer "Undetected uORF, 80–100% uORF+" does not fit in two columns at 8 pt
ord <- names(grp_name); lv <- countsBC_4$Group_Label[match(ord, as.character(countsBC_4$Group))]
leg_lab <- function(x) { g <- sub("\n.*$", "", x); n <- as.integer(sub(".*\\(n=(\\d+)\\).*", "\\1", x)); paste0(grp_name[g], " (n = ", formatC(n, big.mark = ",", format = "d"), ")") }
pal5 <- c("grey45", "#9ECAE1", "#3182BD", "#08519C", "#E41A1C")
dE <- TE_Uu_Final_4; dE$Group_Label <- factor(as.character(dE$Group_Label), levels = lv)
relevel_letters <- function(l) { l$Group_Label <- factor(as.character(l$Group_Label), levels = lv); l }
box_E <- function(y, ylab, ylim, lets, ly) ggplot(dE, aes(x = Group_Label, y = .data[[y]], colour = Group_Label)) +
  geom_hline(yintercept = stats::median(dE[[y]][dE$Group == "Other"], na.rm = TRUE), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.65) + coord_cartesian(ylim = ylim) +
  scale_colour_manual(values = pal5, labels = leg_lab, name = NULL, drop = FALSE) + guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  geom_text(data = relevel_letters(lets), aes(x = Group_Label, y = ly, label = Letters), inherit.aes = FALSE, size = 3, fontface = "italic", family = FIG_FONT, nudge_x = 0.27) +   # V41 (2026-11-12): to the right of the box centre, so that the whiskers of the RNA panel do not cross the letters
  theme_classic() + labs(x = NULL, y = ylab)
no_x <- function(p) p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())   # after compact(), which sets the axis text again
pE1 <- no_x(compact(box_E("TE", "mORF TE", c(0, 3.2), letters_TE_4, 3.05), geom_scale = 0.5))
pE2 <- no_x(compact(box_E("RNA", "RNA (TPM)", c(0, 50), letters_RNA_4, 47), geom_scale = 0.5))
E <- ggpubr::ggarrange(pE1, pE2, ncol = 2, nrow = 1, align = "h", common.legend = TRUE, legend = "bottom")
pD <- p_S5D_hist + labs(title = NULL)
for (i in seq_along(pD$layers)) { l <- pD$layers[[i]]; d <- l$data
  if (inherits(l$geom, "GeomText")) {
    pD$layers[[i]]$aes_params$size <- 8 / ggplot2::.pt                                            # all labels 8 pt in the 2.6-in panel
    if (is.data.frame(d) && nrow(d) == 1 && identical(d$x, 15)) pD$layers[[i]]$data$x <- 17       # the 0-30% label, clear of the clipped-bar note
    if (is.data.frame(d) && nrow(d) == 1 && identical(d$y, 200)) pD$layers[[i]]$data$y <- 170 }   # the clipped-bar note below the 0-30% label
  if (inherits(l$geom, "GeomSegment") && is.data.frame(d) && nrow(d) == 1 && identical(d$y, 200)) pD$layers[[i]]$data$y <- 182 }   # its arrow starts at the note
pD <- compact(pD, geom_scale = 0.5)
TT <- 0.3                                                   # title row: the figure number on the page (V14, 2026-11-03)
H <- TT + h_abc + 0.1 + 2.45
cairo_pdf(file.path(OUT_DIR, "FigureS12_revised_undetected_uORFs.pdf"), width = 6.5, height = H, family = FIG_FONT); grid.newpage()
grid.raster(abc, x = unit(0, "in"), y = unit(H - TT, "in"), width = unit(6.5, "in"), height = unit(h_abc, "in"), just = c("left", "top"), interpolate = TRUE)
vp_draw <- function(g, x, y, w, h) { pushViewport(viewport(x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top"))); grid.draw(g); popViewport() }
y0 <- TT + h_abc + 0.1
vp_draw(fig_title("Figure S10. Undetected uORFs are generally excluded from the mRNA sequences"), 0, 0, 6.5, TT)   # the legend title of the docx
vp_draw(ggplotGrob(pD), 0.1, y0 + 0.05, 2.6, 2.35); vp_draw(ggplotGrob(E), 2.75, y0 + 0.05, 3.75, 2.35)
grid.text("D", x = unit(0.06, "in"), y = unit(H - y0, "in"), just = c("left", "top"), gp = gp_txt(10, TRUE))
grid.text("E", x = unit(2.72, "in"), y = unit(H - y0, "in"), just = c("left", "top"), gp = gp_txt(10, TRUE))
invisible(dev.off())
message(sprintf("S11: page %.2f in tall; groups in E: %s", H, paste(leg_lab(lv), collapse = " | ")))
message("S11 letters TE: ", paste(letters_TE_4$Group, letters_TE_4$Letters, collapse = " "), " | RNA: ", paste(letters_RNA_4$Group, letters_RNA_4$Letters, collapse = " "))
