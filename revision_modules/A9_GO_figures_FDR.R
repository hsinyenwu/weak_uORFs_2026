################################################################################
# A9 — GO enrichment figures (manuscript Figs S8-S10) with the FDR mark          #
#      (Reviewer 2, §2.5 minor: provide background, test and multiple-testing   #
#       details; the submitted plots ranked by raw p without any FDR indication) #
#                                                                               #
# Reads the current GO module output (output/GO_enrichment_results_by_pct_uORF_plus.csv:
# one-sided Fisher's exact tests of TAIR GO terms for the %uORF+ groups of the   #
# single-TuORF cohort and the reference group ">=1 high-uIF uORF or >=2 TuORFs",#
# background = all uORF-containing genes with GO annotation, BH within group)   #
# and redraws the three bubble plots with the module's term ordering, adding:   #
#   * colour = -log10 p, size = gene count, black outline + asterisk = FDR<0.05; #
#   * group labels with n and the number of FDR<0.05 terms per group;           #
#   * root terms (biological_process etc.) not shown.                           #
# Also writes one panel with every FDR<0.05 term across ontologies.             #
# Outputs: output/A9_GO_{MF,BP,CC}_by_pct_uORF_plus_FDR.pdf, output/A9_GO_FDR_significant_terms.pdf, output/A9_summary.txt
# Run from Weak_uORF_submission/:  Rscript revision_modules/A9_GO_figures_FDR.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(stringr) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GO <- file.path(OUT_DIR, "GO_enrichment_results_by_pct_uORF_plus.csv"); stopifnot(file.exists(PATH_GO))
TOP_N <- 15L; FDR_CUT <- 0.05
ROOT <- c("GO:0008150", "GO:0003674", "GO:0005575")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
x <- fread(PATH_GO)
x[, Group := sub(">=1 strong uORF or >=2 uORFs", ">=1 high-uIF uORF or >=2 TuORFs", Group, fixed = TRUE)]   # old group name (runs before the 2026-10-28 renaming in the canonical script)
group_order <- c("1-uORF 0~30% uORF+", "1-uORF 30~80% uORF+", "1-uORF 80~100% uORF+", ">=1 high-uIF uORF or >=2 TuORFs")
stopifnot(all(group_order %in% x$Group))
pretty <- c(`1-uORF 0~30% uORF+` = "0–30% uORF+", `1-uORF 30~80% uORF+` = "30–80% uORF+", `1-uORF 80~100% uORF+` = "80–100% uORF+", `>=1 high-uIF uORF or >=2 TuORFs` = "≥1 high-uIF uORF\nor ≥2 TuORFs")
sizes <- x[, .(n = Group_size[1]), by = Group]; nmap <- setNames(sizes$n, sizes$Group)
.add("A9: GO results read: ", nrow(x), " tested terms; groups ", paste(sprintf("%s n=%d", names(nmap), nmap), collapse = "; "))
cnt <- x[!(GO_ID %in% ROOT), .(tested = .N, fdr05 = sum(p_adj_BH < FDR_CUT)), by = .(Group, Ontology)]   # root terms (biological_process etc.) are not counted or drawn (2026-11-01)
N_BG <- if (file.exists(file.path(OUT_DIR, "GO_universe_size.txt"))) as.integer(readLines(file.path(OUT_DIR, "GO_universe_size.txt"))[1]) else 4886L   # written by the main pipeline (run log: 4,886 genes with a detected TuORF and a GO term)
BG_TXT <- sprintf("all genes with a detected TuORF and Gene Ontology annotation (n = %s)", format(N_BG, big.mark = ","))
for (g in group_order) .add(sprintf("  %-32s FDR<%.2f terms: BP %d, MF %d, CC %d (tested %d)", g, FDR_CUT, cnt[Group == g & Ontology == "BP", fdr05], cnt[Group == g & Ontology == "MF", fdr05], cnt[Group == g & Ontology == "CC", fdr05], cnt[Group == g, sum(tested)]))
sig <- x[p_adj_BH < FDR_CUT][order(Group, Ontology, p_adj_BH)]
.add(sprintf("  FDR<%.2f terms in total: %d (%d in the low group, %d medium, %d high, %d reference; root terms among them: %d)", FDR_CUT, nrow(sig), sum(sig$Group == group_order[1]), sum(sig$Group == group_order[2]), sum(sig$Group == group_order[3]), sum(sig$Group == group_order[4]), sum(sig$GO_ID %in% ROOT)))

make_plot <- function(ont, ont_title, fig_label) {
  d <- x[Ontology == ont & !(GO_ID %in% ROOT)] %>% group_by(Group) %>% slice_min(order_by = p_value, n = TOP_N, with_ties = FALSE) %>% ungroup() %>% as.data.table()
  d[, Group := factor(Group, levels = group_order)]
  fdr_n <- cnt[Ontology == ont][match(group_order, Group), fdr05]
  labs_x <- sprintf("%s\nn = %s\n%d at\nFDR < %.2f", sapply(pretty[group_order], function(x) paste(strwrap(x, 12), collapse = "\n")), trimws(format(nmap[group_order], big.mark = ",")), fdr_n, FDR_CUT)   # wrapped so four columns fit 6.5 in
  d[, Group_Label := factor(as.character(Group), levels = group_order, labels = labs_x)]
  d[, Term_Label := paste0(str_trunc(Term, 52), " (", GO_ID, ")")]
  memb <- d[, .(n_groups = uniqueN(Group), first_group = min(as.integer(Group)), best_p = min(p_value)), by = Term_Label][, block := ifelse(n_groups == 1L, first_group, length(group_order) + 1L)][order(block, n_groups, first_group, best_p)]
  d[, Term_Label := factor(Term_Label, levels = rev(memb$Term_Label))]
  d[, sig := p_adj_BH < FDR_CUT]
  ggplot(d, aes(Group_Label, Term_Label)) +
    geom_point(data = d[sig == TRUE], aes(size = FG_Count), shape = 21, fill = NA, colour = "black", stroke = 1.1) +
    geom_point(aes(size = FG_Count, colour = -log10(p_value))) +
    geom_text(data = d[sig == TRUE], label = "*", nudge_x = 0.28, size = 4, fontface = "bold") +
    scale_colour_gradient(low = "#56B4E9", high = "#D55E00", name = expression(-log[10] * " p")) + scale_size_continuous(name = "genes in group\nwith the term", range = c(0.8, 3.6)) + scale_x_discrete(drop = FALSE) +
    labs(title = sprintf("%s. Gene Ontology enrichment by %%uORF+ group: %s", fig_label, ont_title),
         subtitle = sprintf("single-TuORF genes grouped by %%uORF+ (CAGE) plus the reference group\nbackground = %s; one-sided Fisher's exact test\n%d terms with the smallest p per group; black outline and * = Benjamini–Hochberg FDR < %.2f; root terms not shown", BG_TXT, TOP_N, FDR_CUT),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 0, face = "bold", size = 7, lineheight = 0.9), axis.text.y = element_text(size = 7.5), plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 6.8), legend.key.size = unit(0.4, "cm"), panel.grid.minor = element_blank())
}
plots <- list(MF = make_plot("MF", "Molecular Function", "Figure S16"), BP = make_plot("BP", "Biological Process", "Figure S17"), CC = make_plot("CC", "Cellular Component", "Figure S18"))   # +1 after the new S16 (2026-10-31); -1 after the mORF figure S11 was removed (V14, 2026-11-03))
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
go_compact <- function(p) compact(p, wrap_title = 80, wrap_sub = 105) + theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal", plot.title.position = "plot", plot.caption.position = "plot", legend.spacing.x = unit(0.15, "cm")) + scale_y_discrete(expand = expansion(add = 0.8)) +
  guides(colour = guide_colourbar(barwidth = unit(1.6, "cm"), barheight = unit(0.3, "cm"), order = 1, theme = theme(legend.text = element_text(margin = margin(t = 3)))), size = guide_legend(nrow = 1, order = 2))
for (ont in names(plots)) { fn <- file.path(OUT_DIR, sprintf("A9_GO_%s_by_pct_uORF_plus_FDR.pdf", ont)); fig_save(fn, go_compact(plots[[ont]]), h = 9); .add("  saved ", fn) }
# all FDR-significant terms in one panel
s2 <- copy(sig)[!(GO_ID %in% ROOT)]; s2[, Group := factor(Group, levels = group_order)]; s2[, Label := paste0(Term, " (", GO_ID, ")")]
s2[, Ontology := factor(Ontology, levels = c("MF", "BP", "CC"))]; setorder(s2, Group, Ontology, p_adj_BH); s2[, Label := factor(Label, levels = rev(unique(Label)))]
s2[, Group_Label := factor(as.character(Group), levels = group_order, labels = sprintf("%s (n = %s)", gsub("\n", " ", pretty[group_order]), format(nmap[group_order], big.mark = ",")))]
pS <- ggplot(s2, aes(-log10(p_adj_BH), Label, fill = Ontology)) + geom_col(width = 0.7) + geom_vline(xintercept = -log10(FDR_CUT), linetype = 3) + facet_grid(Group_Label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(MF = "#D55E00", BP = "#009E73", CC = "#56B4E9"), name = "Ontology") + geom_text(aes(label = sprintf("%d/%d", FG_Count, BG_Count)), hjust = -0.1, size = 2.4) + scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(-log[10] * " (Benjamini–Hochberg FDR)"), y = NULL, title = sprintf("All Gene Ontology terms at FDR < %.2f by %%uORF+ group", FDR_CUT), subtitle = "labels = genes in the group with the term / genes in the background with the term; dotted line = FDR 0.05; root terms not shown") +
  theme_bw(base_size = 9) + theme(strip.text.y = element_text(angle = 0, size = 7.5), axis.text.y = element_text(size = 7.5), plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
fig_save(file.path(OUT_DIR, "A9_GO_FDR_significant_terms.pdf"), go_compact(pS), h = 8); .add("  saved output/A9_GO_FDR_significant_terms.pdf")
fwrite(sig[, .(Group, Ontology, GO_ID, Term, FG_Count, BG_Count, Odds_Ratio = round(Odds_Ratio, 2), p_value = signif(p_value, 3), p_adj_BH = signif(p_adj_BH, 3))], file.path(OUT_DIR, "A9_GO_FDR_significant_terms.csv"))
writeLines(.log, file.path(OUT_DIR, "A9_summary.txt"))
