# fig_compact.R — shared styling for the revision figures (Larry, 2026-10-28): pages 6.5 x <= 9 in, Arial, no text below 8 pt.
# Usage:  source("revision_modules/fig_compact.R")
#         g <- gridExtra::arrangeGrob(grobs = lapply(list(pA, pB, pC), compact), ncol = 2, top = fig_title("Figure S7. ...")); fig_save(file, g, h = 9)
#         ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = lapply(list(pA, pB), compact), ncol = 2, nrow = 1), top = fig_title("..."))
# compact(): Arial everywhere, axis/legend/strip text 8 pt, axis titles 9 pt, panel titles 10 pt bold, text layers raised to >= 8 pt,
# long titles/subtitles wrapped, and lines/points thinned by `geom_scale` (the factor by which the panel shrank from the size the
# plot was drawn at; default 0.55, the main-figure panels use 0.45 so they match the submitted figures).  fig_save(): cairo_pdf
# (embeds Arial) at 6.5 in width and at most 9 in height.
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
FIG_FONT <- "Arial"; FIG_W <- 6.5; FIG_H_MAX <- 9; FIG_PT_MIN <- 8; GEOM_SCALE <- 0.55
.mm_min <- FIG_PT_MIN / ggplot2::.pt                      # geom_text size (mm) that gives 8 pt
fig_theme <- function(base = FIG_PT_MIN) theme(
  text = element_text(family = FIG_FONT, size = base), axis.text = element_text(family = FIG_FONT, size = base), axis.text.x = element_text(family = FIG_FONT, size = base), axis.text.y = element_text(family = FIG_FONT, size = base),
  axis.title = element_text(family = FIG_FONT, size = base + 1), strip.text.x = element_text(family = FIG_FONT, size = base), strip.text.y = element_text(family = FIG_FONT, size = base),
  plot.title = element_text(family = FIG_FONT, size = base + 2, face = "bold"), plot.subtitle = element_text(family = FIG_FONT, size = base),
  plot.caption = element_text(family = FIG_FONT, size = base), legend.text = element_text(family = FIG_FONT, size = base, margin = margin(l = 2, r = 6)), legend.title = element_text(family = FIG_FONT, size = base),
  strip.text = element_text(family = FIG_FONT, size = base), legend.key.size = unit(0.3, "cm"), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(4, 10, 6, 8),
  plot.title.position = "plot", plot.caption.position = "plot", legend.box.spacing = unit(4, "pt"), legend.key.spacing.x = unit(6, "pt"), legend.key.spacing.y = unit(2, "pt"))
thin <- function(p, f = GEOM_SCALE) {                           # line widths and point sizes scaled by the panel shrink factor
  if (!inherits(p, "ggplot")) return(p)
  p <- p + theme(geom = element_geom(linewidth = 0.5 * f, borderwidth = 0.5 * f, pointsize = 1.5 * f),   # ggplot2 >= 4.0 defaults for layers without explicit sizes
                 line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), rect = element_rect(linewidth = 0.3))
  for (l in p$layers) { g <- class(l$geom)[1]; if (grepl("Text|Label|Logo", g)) next
    ap <- l$aes_params
    if (is.numeric(ap$linewidth)) l$aes_params$linewidth <- ap$linewidth * f
    if (is.numeric(ap$size)) l$aes_params$size <- ap$size * f                         # point size, or the old 'size' of a line layer
    if (is.numeric(ap$stroke)) l$aes_params$stroke <- ap$stroke * f
    if (is.numeric(l$geom_params$outlier.size)) l$geom_params$outlier.size <- l$geom_params$outlier.size * f
    if (is.list(l$geom_params$outlier_gp) && is.numeric(l$geom_params$outlier_gp$size)) l$geom_params$outlier_gp$size <- l$geom_params$outlier_gp$size * f }   # ggplot2 >= 4.0 boxplots
  p }
legend_out <- function(p) {                                     # a legend drawn inside the panel can overlap the data: move it below the plot
  lp <- p$theme$legend.position
  if (is.numeric(lp) || identical(lp, "inside")) p <- p + theme(legend.position = "bottom", legend.direction = "vertical", legend.justification = c(0, 1), legend.box.just = "left", legend.margin = margin(2, 0, 0, 0), legend.location = "plot")   # aligned to the whole panel, not the plot area
  p }
compact <- function(p, wrap_title = 32, wrap_sub = 46, geom_scale = GEOM_SCALE, wrap_x = 40) {      # characters per line at 10 pt bold / 8 pt in a 3.2-in panel
  if (!inherits(p, "ggplot")) return(p)
  p <- legend_out(thin(p, geom_scale)) + fig_theme()
  for (l in p$layers) if (grepl("Text|Label", class(l$geom)[1])) {                      # geom_text / geom_label / annotate("text")
    s <- l$aes_params$size; if (is.null(s) || s < .mm_min) l$aes_params$size <- .mm_min; l$aes_params$family <- FIG_FONT }
  for (k in c("title", "subtitle", "caption")) if (is.character(p$labels[[k]])) p$labels[[k]] <- paste(strwrap(gsub("\n", if (k == "title") " " else "; ", p$labels[[k]]), if (k == "title") wrap_title else wrap_sub), collapse = "\n")   # line breaks of subtitles become '; '
  if (is.character(p$labels$x) && nchar(p$labels$x) > wrap_x && !grepl("\n", p$labels$x, fixed = TRUE)) p$labels$x <- paste(strwrap(p$labels$x, wrap_x), collapse = "\n")
  if (is.character(p$labels$y) && nchar(p$labels$y) > 26 && !grepl("\n", p$labels$y, fixed = TRUE)) p$labels$y <- paste(strwrap(p$labels$y, 26), collapse = "\n")           # a rotated title longer than the panel is clipped
  p }
compact_wide <- function(p, wrap_title = 72, wrap_sub = 104, wrap_x = 90, ...) compact(p, wrap_title = wrap_title, wrap_sub = wrap_sub, wrap_x = wrap_x, ...)   # a panel spanning the full 6.5-in width
fig_title <- function(txt, size = 10, wrap = 82) textGrob(paste(strwrap(txt, wrap), collapse = "\n"), gp = gpar(fontsize = size, fontface = "bold", fontfamily = FIG_FONT, lineheight = 0.95))
fig_save <- function(file, g, w = FIG_W, h = FIG_H_MAX) ggsave(file, g, width = w, height = min(h, FIG_H_MAX), device = cairo_pdf, family = FIG_FONT)

## ---- V30 (2026-11-09, Larry: "put the legends in the plot"; "panels do not align"): legends inside the plotting area where there is free
## room, otherwise at the side; the panels of one row share the top and the bottom of their plotting areas (patchwork).
## Usage (after compact(), which still moves an inside legend below the plot):
##   legend_in(compact(p), x = 0.98, y = 0.03)            legend inside, anchored at (x, y) in panel coordinates (0-1); the corner is taken from x and y
##   legend_side(compact(p), wrap = wrap_paren)           legend at the right of the panel, labels broken before every " ("
##   fig_grid(list(pA, pB, pC, pD), ncol = 2, top = fig_title("Figure S.."))    aligned grid; design = "AB\nCC" for panels that span the page
.legend_small <- function(title = FALSE, pt = FIG_PT_MIN) theme(legend.key.size = unit(0.28, "cm"), legend.key = element_blank(), legend.key.spacing.y = unit(0, "pt"),
  legend.text = element_text(family = FIG_FONT, size = pt, margin = margin(l = 2, r = 2), lineheight = 0.9),
  legend.title = if (title) element_text(family = FIG_FONT, size = pt) else element_blank())
wrap_paren <- function(x) gsub(" (", "\n(", x, fixed = TRUE)                       # "High uIF (>1) (n=1116) a" -> three short lines
wrap_n <- function(x) sub(" \\(n ?=", "\n(n=", x)                                  # one break, before the count
legend_labels <- function(p, f) {                                                 # change the legend labels of the colour / fill scales without touching their colours
  p$scales <- p$scales$clone(); hit <- FALSE
  for (sc in p$scales$scales) if (any(sc$aesthetics %in% c("colour", "fill"))) { old <- sc$labels; hit <- TRUE
    sc$labels <- if (is.function(old)) (function(o) function(x) f(o(x)))(old) else if (is.character(old)) stats::setNames(f(old), names(old)) else f }
  if (!hit) for (a in intersect(c("colour", "fill"), unlist(lapply(c(list(p$mapping), lapply(p$layers, function(l) l$mapping)), names))))
    p <- p + (if (a == "colour") scale_colour_discrete(labels = f) else scale_fill_discrete(labels = f))
  p }
legend_in <- function(p, x = 0.98, y = 0.03, just = NULL, title = FALSE, horizontal = FALSE, wrap = NULL, bg = 0.75, pt = FIG_PT_MIN) {   # pt: legend text size; 8 pt unless Larry asks for less (Figure 3A: 7 pt, 2026-11-09)
  if (is.null(just)) just <- c(if (x > 0.5) 1 else 0, if (y > 0.5) 1 else 0)
  if (!is.null(wrap)) p <- legend_labels(p, wrap)
  p + .legend_small(title, pt) + theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification.inside = just,
        legend.direction = if (horizontal) "horizontal" else "vertical", legend.box.just = "left", legend.margin = margin(1, 2, 1, 2),
        legend.background = if (bg > 0) element_rect(fill = scales::alpha("white", bg), colour = NA) else element_blank()) }
legend_side <- function(p, wrap = NULL, title = FALSE) {
  if (!is.null(wrap)) p <- legend_labels(p, wrap)
  p + .legend_small(title) + theme(legend.position = "right", legend.direction = "vertical", legend.justification = c(0, 0.5), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"), legend.key.spacing.y = unit(3, "pt"), legend.background = element_blank()) }
legend_below <- function(p, horizontal = TRUE, title = FALSE) p + .legend_small(title) + theme(legend.position = "bottom", legend.direction = if (horizontal) "horizontal" else "vertical",
        legend.justification = c(0, 1), legend.location = "plot", legend.margin = margin(2, 0, 0, 0), legend.key.spacing.x = unit(6, "pt"))
text_pt <- function(p, pt = FIG_PT_MIN, geoms = "Text|Label") {                   # every text layer at exactly `pt` points (compact() only raises text below 8 pt)
  for (l in p$layers) if (grepl(geoms, class(l$geom)[1])) l$aes_params$size <- pt / ggplot2::.pt
  p }
with_cairo <- function(expr, w = FIG_W, h = FIG_H_MAX) {                          # sizes are measured on a cairo device with Arial (the base pdf device fails on the font)
  prev <- grDevices::dev.cur(); f <- tempfile(fileext = ".pdf"); grDevices::cairo_pdf(f, width = w, height = h, family = FIG_FONT)
  on.exit({ grDevices::dev.off(); unlink(f); if (prev > 1) grDevices::dev.set(prev) }); force(expr) }
fig_grid <- function(plots, ncol = NULL, nrow = NULL, design = NULL, heights = NULL, widths = NULL, top = NULL, w = FIG_W, h = FIG_H_MAX) {
  # Panels are aligned with patchwork (same top and bottom of the plotting areas within a row, same left and right within a column).
  # A row of the design that is one full-width panel ("CC") is laid out on its own, so that its axis labels do not push the panels of the
  # other rows; `heights` are then shares of the page for the blocks (as in gridExtra) and shares of the plotting areas inside a block.
  grob1 <- function(pl, ...) patchwork::patchworkGrob(patchwork::wrap_plots(pl, ...) + patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))))   # no extra margin around a block (the panels carry their own)
  with_cairo({
    rows <- if (is.character(design)) strsplit(strsplit(design, "\n")[[1]], "") else NULL
    full <- if (is.null(rows)) FALSE else vapply(rows, function(r) length(unique(r)) == 1, TRUE)
    if (!any(full)) g <- grob1(plots, ncol = ncol, nrow = nrow, design = design, heights = heights, widths = widths)
    else { if (is.null(heights)) heights <- rep(1, length(rows))
      blk <- cumsum(c(TRUE, full[-1] | full[-length(full)]))                                   # a full-width row starts a block and ends it
      gl <- lapply(split(seq_along(rows), blk), function(ix) { lets <- unique(unlist(rows[ix])); lets <- lets[lets != "#"]
        if (length(lets) == 1) return(grob1(plots[match(lets, LETTERS)]))
        grob1(plots[match(lets, LETTERS)], design = paste(vapply(rows[ix], function(r) paste(ifelse(r == "#", "#", LETTERS[match(r, lets)]), collapse = ""), ""), collapse = "\n"), heights = heights[ix], widths = widths) })
      g <- gridExtra::arrangeGrob(grobs = gl, ncol = 1, heights = vapply(split(heights, blk), sum, 0)) }
    if (!is.null(top)) g <- gridExtra::arrangeGrob(g, top = top); g }, w, h) }
cor_left <- function(p, x = 0.03) { for (l in p$layers) if (grepl("Cor", class(l$stat)[1])) l$stat_params$label.x.npc <- x; p }   # the correlation label starts at the left of a narrow panel
wrap_w <- function(w) function(x) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"), "")                                # legend labels wrapped at w characters
free_xlab <- function(p) patchwork::free(p, type = "label", side = "b")   # the x-axis title stays next to its axis when a neighbouring panel has taller tick labels (fig_grid aligns the titles otherwise)
fig_row <- function(plots, widths = 1) with_cairo(cowplot::plot_grid(plotlist = plots, nrow = 1, align = "h", axis = "tb", rel_widths = widths))   # one aligned row for a panel with facet strips or a legend below next to a plain panel: the plotting areas (strips included) share top and bottom, and the space below them is shared (use the result as one full-width panel of fig_grid)
