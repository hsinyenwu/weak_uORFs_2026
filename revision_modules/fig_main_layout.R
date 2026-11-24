# fig_main_layout.R — V30 (2026-11-09, Larry): layout of the data panels of Figures 2-5 and of Supplemental Fig. S2.
# Legends sit inside the plotting area where there is free room and at the right of the panel otherwise (Figure 2C and 2F, as Larry asked);
# the panels of one row share the top and the bottom of their plotting areas (fig_grid(), patchwork). Text stays at 8 pt or more.
# Sourced at the end of the revised figure driver. Stand-alone, for layout work without rerunning the driver (the driver saves the plot objects):
#   Rscript revision_modules/fig_main_layout.R [DRIVER_OUT_DIR]      (run from Weak_uORF_submission/)
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
if (!exists("OUT_DIR")) { args <- commandArgs(trailingOnly = TRUE); OUT_DIR <- if (length(args)) args[1] else "_to_delete/driver_compact_20261028" }
if (!exists("fig_grid")) source("revision_modules/fig_compact.R")
list2env(readRDS(file.path(OUT_DIR, "main_layout_objects.rds")), globalenv())   # the plot objects that the driver saved (p3B_c ... pE; f3A and f3B = Figure 3A-B)
m4 <- theme(plot.margin = margin(3, 4, 3, 4))
one_line <- function(x) gsub("\n", " ", x, fixed = TRUE)

## Figure 2 B-F: a 3 x 2 grid; the top-left cell stays empty for the schematic (panel A), which the composite scripts place there.
## The legends of C and F are at the right of their panels (third column), so that the plotting areas of B-C and of D-E-F line up.
f2 <- function(p) compact(p, geom_scale = 0.45) + m4
B2 <- text_pt(f2(p3B_c)); for (l in B2$layers) if (!is.null(l$aes_params$label) && !is.null(l$data$x)) l$data$x <- ifelse(l$data$x > 80, 96, ifelse(l$data$x > 40, 55, l$data$x))
B2 <- suppressMessages(B2 + scale_y_continuous(expand = expansion(mult = c(0.02, 0.34))) + scale_x_continuous(breaks = seq(0, 100, 25)) + expand_limits(x = 108))   # the three count labels in free room above the bars; the third to the right of the line at 80
## V54 (2026-11-23, Larry: "Reletter"): the second row is lettered E, F, G (was D, E, F); the table of A19 is panel D and the composite scripts place it between the two rows
p3D <- p3D + labs(title = "E"); p3E <- p3E + labs(title = "F"); p3F <- p3F + labs(title = "G")
P2 <- list(patchwork::plot_spacer(), B2, legend_side(f2(p3C), wrap_paren), cor_left(text_pt(f2(p3D))), cor_left(text_pt(f2(p3E)), 0.01), legend_side(f2(p3F), wrap_paren))
fig_save(file.path(OUT_DIR, "Figure2_compact_panels_B-F.pdf"), fig_grid(P2, ncol = 3, w = 6.5, h = 4.4), w = 6.5, h = 4.4)

## Figure 3 A-B: the legend of A inside the plot (top right, one line per group, no legend title, 7 pt)
A3 <- free_xlab(legend_in(compact(f3A), 0.99, 0.99, wrap = one_line, pt = 7))   # legend text 7 pt (Larry, 2026-11-09: 'Make it 7 pt'; every other text of the figures stays at 8 pt or more); the axis title next to its axis (B has two-line tick labels)
B3 <- compact(f3B + scale_x_discrete(labels = function(x) sub("^No TuORF$", "No\nTuORF", sub("\n.*$", "", x)))) + .legend_small(title = TRUE)   # the chi-square paragraph is not part of the figure
fig_save(file.path(OUT_DIR, "Figure3_compact_panels_A-B.pdf"), fig_grid(list(A3, B3), ncol = 2, widths = c(1.7, 1), w = 6.5, h = 2.45), w = 6.5, h = 2.45)

## Figure 4 D-F: one aligned row, legends inside
f4 <- function(p) compact(p, geom_scale = 0.45) + m4
P4 <- list(legend_in(f4(p5D), 0.02, 0.99), legend_in(f4(p5E), 0.99, 0.03), legend_in(f4(p5F), 0.99, 0.03, wrap = wrap_n))
fig_save(file.path(OUT_DIR, "Figure4_compact_panels_D-F.pdf"), fig_grid(P4, ncol = 3, w = 6.5, h = 2.1), w = 6.5, h = 2.1)

## Supplemental Fig. S2: uIF versus TE (A) and TE by uIF class (B), legend of B inside the plot (bottom right)
fS2 <- function(p, ttl) compact(p + labs(title = ttl), geom_scale = 0.5)
S2 <- fig_grid(list(text_pt(fS2(pD_c, "A"), 9), legend_in(fS2(pE, "B"), 0.99, 0.03)), ncol = 2, top = fig_title("Figure S2. uIF and mORF translation efficiency"), w = 6.5, h = 3.3)
fig_save(file.path(OUT_DIR, "FigureS2_revised_uIF_vs_TE.pdf"), S2, w = 6.5, h = 3.3)
## Figure 5 A-F on one page: the legend of C inside the plot (top right); rows aligned
f5 <- function(p) compact(p + theme(axis.text.x = element_text(angle = 30, hjust = 1)), geom_scale = 0.5)
P5 <- list(f5(p8A), f5(p8B), legend_in(f5(p8C), 0.995, 0.995), f5(p8D), free_xlab(f5(p8E)), f5(p8H + labs(title = "F: uORF length before vs after TSS filtering", y = "uORF length (aa, log10 scale)")))
fig_save(file.path(OUT_DIR, "Figure5_compact_panels_A-F.pdf"), fig_grid(P5, ncol = 2, w = 6.5, h = 8.8), h = 8.8)

message("fig_main_layout.R: Figures 2 B-F, 3 A-B, 4 D-F, 5 and S2 written to ", OUT_DIR)
