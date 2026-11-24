# compose_main_figures.R — 6.5-in-wide composites of the main figures for the revised manuscript (2026-10-28, Larry).
# Regenerated panels (compact PDFs from the figure driver / ladder script, Arial >= 8 pt) are combined with the schematic and
# single-gene browser panels cropped from the submitted figure PDFs (Submit_Genome Research/Figures 20260826), which already fit
# a 6.5 x 9 in page. Everything is rasterised at 300 dpi and placed with grid at exact inch positions; output PDF (raster) + PNG.
# usage: Rscript revision_modules/compose_main_figures.R DRIVER_OUT_DIR     (run from Weak_uORF_submission/)
suppressPackageStartupMessages({ library(png); library(grid) })
args <- commandArgs(trailingOnly = TRUE); DRV <- if (length(args)) args[1] else "_to_delete/driver_compact_20261028"
OUT <- "output_revised_figs_20261028"; WK <- file.path(OUT, "composite_work"); dir.create(WK, showWarnings = FALSE)
SUB <- "/Users/wu/Desktop/weak_uORF_revision/Submit_Genome Research/Figures 20260826"
DPI <- 300; S <- DPI / 72                                   # submitted PDFs: 1 pt = 1 px at 72 dpi
render <- function(pdf, out, page = 1) { system2("pdftoppm", c("-r", DPI, "-f", page, "-l", page, "-png", "-singlefile", shQuote(pdf), shQuote(file.path(WK, out)))); readPNG(file.path(WK, paste0(out, ".png"))) }
crop <- function(pdf, out, x, y, w, h) {                     # x, y, w, h in 72-dpi points of the submitted page
  system2("pdftoppm", c("-r", DPI, "-x", round(x * S), "-y", round(y * S), "-W", round(w * S), "-H", round(h * S), "-png", "-singlefile", shQuote(pdf), shQuote(file.path(WK, out))))
  readPNG(file.path(WK, paste0(out, ".png"))) }
place <- function(img, x, y, w, H) { h <- w * dim(img)[1] / dim(img)[2]                  # keep the aspect ratio; y measured from the top
  grid.raster(img, x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top"), interpolate = TRUE); h }
page <- function(file, H) cairo_pdf(file, width = 6.5, height = H, family = "Arial")
done <- function(file) { dev.off(); cat("saved", file, "\n") }   # PDF only (Larry, 2026-11-01); the PNG renders for the docx are made by the run script

## Figure 1 (V16, 2026-11-04): A schematic (submitted), B = uIF distribution (compact; formerly Supplemental Fig. S2), C-F single-gene browser panels (submitted rows D-G, re-lettered by reletter_fig1_rows.py); the uIF-TE panels are now Supplemental Fig. S2
f1 <- file.path(WK, "Fig1_submitted_relabelled_v16.pdf")   # the submitted file with the uIF-class labels (relabel_submitted_figure_labels.py, 2026-11-01) and the row letters C-F (V16); same geometry
A  <- crop(f1, "f1_A", 0, 0, 141, 134); Dg <- crop(f1, "f1_D", 0, 134, 468, 129); Eg <- crop(f1, "f1_E", 0, 273, 468, 119); Fg <- crop(f1, "f1_F", 0, 404, 468, 117); Gg <- crop(f1, "f1_G", 0, 532, 468, 118)
Dg[seq_len(round(3 * S)), round(100 * S):round(420 * S), ] <- 1   # the row-D crop starts 2 pt under the submitted B/C x-axis titles: blank their descenders (the D letter, x < 12 pt, is kept)
B <- render(file.path(DRV, "Figure1_compact_B_dist.pdf"), "f1_B")   # V16: the distribution panel (3.4 x 2.75 in)
H <- 8.3; page(file.path(OUT, "Figure1_compact_composite.pdf"), H)   # V21 (2026-11-07): B is 1.9 in tall, like A
place(A, 0.00, 0.05, 1.96, H); place(B, 2.00, 0.00, 3.40, H)      # B is 1.9 in tall (V21; 2.75 in V16-V20)
y <- 2.05; for (im in list(Dg, Eg, Fg, Gg)) { h <- place(im, 0.3, y, 5.9, H); y <- y + h + 0.02 }   # 5.9 in wide so the four browser rows end above 9 in
done(file.path(OUT, "Figure1_compact_composite.pdf")); cat("  Figure 1 bottom at", round(y, 2), "in\n")

## Figure 2 (V30, 2026-11-09, Larry): B-F are one aligned 3 x 2 grid (fig_main_layout.R: legends of C and F at the right, plotting areas of a row aligned);
## its top-left cell is empty and takes the schematic A (submitted), with the letter A at the height of the letters B and C
f2 <- file.path(SUB, "Fig2_uORF+_v20260414.pdf"); A2 <- crop(f2, "f2_A", 0, 0, 155, 80)
BF <- render(file.path(DRV, "Figure2_compact_panels_B-F.pdf"), "f2_BF")
G2 <- render("output/A19_Figure2G_TSS_vs_splicing_table.pdf", "f2_G")   # panel G: translated uORFs affected by TSS selection and by leader splicing (A19; Larry, 2026-11-01)
F2 <- list(H_BF = 4.4, A = c(0.08, 0.05, 1.85))                          # height of the B-F grid (in); x, y and width of the schematic A (in)
## V54 (2026-11-23, Larry: "Reletter"): the table (A19) is cited before the plots of the second row, so it is panel D and sits between the two rows of
## the grid; the plots of the second row are E-G (fig_main_layout.R). The grid is cut at 2.2 in, in the white band between its rows (2.14-2.26 in).
SPLIT <- 2.2; k <- round(SPLIT / F2$H_BF * dim(BF)[1]); BF1 <- BF[seq_len(k), , , drop = FALSE]; BF2 <- BF[(k + 1):dim(BF)[1], , , drop = FALSE]
H <- F2$H_BF + 0.05 + 2.52 + 0.05; page(file.path(OUT, "Figure2_compact_composite.pdf"), H)   # the table D is 2.52 in tall with its note
h1 <- place(BF1, 0, 0, 6.5, H); place(A2, F2$A[1], F2$A[2], F2$A[3], H)
place(G2, 0.0, h1 + 0.05, 6.5, H)              # full-width table row with its note (panel D)
place(BF2, 0, h1 + 0.05 + 2.52 + 0.05, 6.5, H)  # panels E-G
done(file.path(OUT, "Figure2_compact_composite.pdf"))

## Figure 3: A-B compact (ATG frequency) above C (prevalence ladder + TPM table)
AB <- render(file.path(DRV, "Figure3_compact_panels_A-B.pdf"), "f3_AB"); Cp <- render(file.path(OUT, "Figure3_revised_panel_C_prevalence_ladder_plot.pdf"), "f3_C"); Ct <- render(file.path(OUT, "Figure3_revised_panel_C_prevalence_ladder_table.pdf"), "f3_Ct")
H <- 8.35; page(file.path(OUT, "Figure3_compact_composite.pdf"), H)   # V30 (2026-11-09): A-B 2.45 in, C plot 2.9 in, table 2.8 in   # V20 (2026-11-06): the TPM table has two more rows (2.8 in)
h <- place(AB, 0, 0, 6.5, H); h2 <- place(Cp, 0.0, h + 0.1, 6.5, H); place(Ct, 0, h + 0.1 + h2 + 0.05, 6.5, H)   # C plot and the TPM table as full-width rows
done(file.path(OUT, "Figure3_compact_composite.pdf"))

## Figure 4: A-C single-gene browser panels (submitted, scaled to 6.5 in), D-F compact
f4 <- file.path(SUB, "Fig4_CpuORF_v20260511pdf.pdf")
A4 <- crop(f4, "f4_A", 0, 0, 403, 105); B4 <- crop(f4, "f4_B", 0, 108, 403, 96); C4 <- crop(f4, "f4_C", 0, 207, 403, 121)   # stop above the D-F letters of the submitted layout
DF4 <- render(file.path(DRV, "Figure4_compact_panels_D-F.pdf"), "f4_DEF")   # V30 (2026-11-09): D-F as one aligned row with the legends inside the plots (fig_main_layout.R), 2.1 in tall
H <- 7.6; page(file.path(OUT, "Figure4_compact_composite.pdf"), H)
y <- 0; for (im in list(A4, B4, C4)) { h <- place(im, 0, y, 6.5, H); y <- y + h + 0.05 }
place(DF4, 0, y, 6.5, H)
done(file.path(OUT, "Figure4_compact_composite.pdf")); cat("  Figure 4 D-F row at", round(y, 2), "in; bottom at", round(y + 2.1, 2), "in\n")
