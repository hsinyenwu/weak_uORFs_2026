################################################################################
# A24 — Gene Ontology figures (Supplemental Figs. S15-S17) and Table S8 for              #
#       the gene-level grouping adopted on 2026-11-05 (Larry): all TuORF genes  #
#       grouped by the highest %uORF+ among their TuORFs                         #
#         only low            every TuORF below 30% uORF+                        #
#         at least one medium highest TuORF in 30-80%                            #
#         at least one high   highest TuORF at 80% or more                       #
# Input: output/A23_GO_all_TuORF_genes_by_max_pUP.csv (module A23: same         #
# annotation, universe, term filter, one-sided Fisher test and Benjamini-       #
# Hochberg correction as the canonical script). Drawing follows A9 (bubble      #
# plots, 15 smallest p per group, black outline + asterisk = FDR < 0.05, root   #
# terms not drawn) with three columns instead of four.                          #
# Outputs: output/A24_GO_{MF,BP,CC}_by_max_pUP_FDR.pdf, output/A24_GO_FDR_significant_terms.pdf,
#          weak_uORF_revision/revision_tables/TableS8_GO_enrichment_by_pct_uORF_plus.xlsx (rewritten; previous copy backed up),
#          output/A24_summary.txt
# Run from Weak_uORF_submission/:  Rscript revision_modules/A24_GO_by_max_pUP_figures_table.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(stringr); library(openxlsx) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GO <- file.path(OUT_DIR, "A23_GO_all_TuORF_genes_by_max_pUP.csv"); stopifnot(file.exists(PATH_GO))
TABS <- "/Users/wu/Desktop/weak_uORF_revision/revision_tables"; BAK <- file.path(TABS, "_replaced_by_V19_20261105"); dir.create(BAK, showWarnings = FALSE)
TOP_N <- 15L; FDR_CUT <- 0.05; ROOT <- c("GO:0008150", "GO:0003674", "GO:0005575")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
x <- fread(PATH_GO)
group_order <- c("all TuORF genes: only low %uORF+ (<30)", "all TuORF genes: at least one medium (30-80)", "all TuORF genes: at least one high (>=80)")
stopifnot(all(group_order %in% x$Group))
short <- c("only low %uORF+ (every TuORF <30%)", "at least one medium (highest TuORF 30-80%)", "at least one high (highest TuORF >=80%)")
names(short) <- group_order
pretty <- c("only low %uORF+\n(all TuORFs <30%)", "≥1 medium\n(highest 30–80%)", "≥1 high\n(highest ≥80%)"); names(pretty) <- group_order
sizes <- x[, .(n = Group_size[1]), by = Group]; nmap <- setNames(sizes$n, sizes$Group)
.add("A24: GO results read: ", nrow(x), " tested terms; groups ", paste(sprintf("%s n=%d", short[names(nmap)], nmap), collapse = "; "))
cnt <- x[!(GO_ID %in% ROOT), .(tested = .N, fdr05 = sum(p_adj_BH < FDR_CUT)), by = .(Group, Ontology)]
N_BG <- 4886L; BG_TXT <- sprintf("all genes with a detected TuORF and Gene Ontology annotation (n = %s)", format(N_BG, big.mark = ","))
for (g in group_order) .add(sprintf("  %-45s FDR<%.2f terms: BP %d, MF %d, CC %d (tested %d)", short[g], FDR_CUT, cnt[Group == g & Ontology == "BP", fdr05], cnt[Group == g & Ontology == "MF", fdr05], cnt[Group == g & Ontology == "CC", fdr05], cnt[Group == g, sum(tested)]))
sig <- x[p_adj_BH < FDR_CUT & !(GO_ID %in% ROOT)][order(Group, Ontology, p_adj_BH)]
.add(sprintf("  FDR<%.2f terms (root terms excluded): %d (%d only low, %d medium, %d high)", FDR_CUT, nrow(sig), sum(sig$Group == group_order[1]), sum(sig$Group == group_order[2]), sum(sig$Group == group_order[3])))

make_plot <- function(ont, ont_title, fig_label) {
  d <- x[Ontology == ont & !(GO_ID %in% ROOT)] %>% group_by(Group) %>% slice_min(order_by = p_value, n = TOP_N, with_ties = FALSE) %>% ungroup() %>% as.data.table()
  d[, Group := factor(Group, levels = group_order)]
  fdr_n <- cnt[Ontology == ont][match(group_order, Group), fdr05]; fdr_n[is.na(fdr_n)] <- 0L
  labs_x <- sprintf("%s\nn = %s\n%d at FDR < %.2f", pretty[group_order], trimws(format(nmap[group_order], big.mark = ",")), fdr_n, FDR_CUT)
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
    labs(title = sprintf("%s. Gene Ontology %s terms enriched in TuORF genes grouped by the highest %%uORF+ among their TuORFs", fig_label, ont_title),   # V23: = legend title
         subtitle = sprintf("all TuORF genes grouped by the highest %%uORF+ among their TuORFs (CAGE)\nbackground = %s; one-sided Fisher's exact test\n%d terms with the smallest p per group; black outline and * = Benjamini–Hochberg FDR < %.2f; root terms not shown", BG_TXT, TOP_N, FDR_CUT),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 0, face = "bold", size = 7, lineheight = 0.9), axis.text.y = element_text(size = 7.5), plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 6.8), legend.key.size = unit(0.4, "cm"), panel.grid.minor = element_blank())
}
plots <- list(MF = make_plot("MF", "Molecular Function", "Figure S15"), BP = make_plot("BP", "Biological Process", "Figure S16"), CC = make_plot("CC", "Cellular Component", "Figure S17"))
source("revision_modules/fig_compact.R")
go_compact <- function(p) compact(p, wrap_title = 80, wrap_sub = 105) + theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal", plot.title.position = "plot", plot.caption.position = "plot", legend.spacing.x = unit(0.15, "cm")) + scale_y_discrete(expand = expansion(add = 0.8)) +
  guides(colour = guide_colourbar(barwidth = unit(1.6, "cm"), barheight = unit(0.3, "cm"), order = 1, theme = theme(legend.text = element_text(margin = margin(t = 3)))), size = guide_legend(nrow = 1, order = 2))
for (ont in names(plots)) { fn <- file.path(OUT_DIR, sprintf("A24_GO_%s_by_max_pUP_FDR.pdf", ont)); fig_save(fn, go_compact(plots[[ont]]), h = 9); .add("  saved ", fn) }
# all FDR-significant terms in one panel
s2 <- copy(sig); s2[, Group := factor(Group, levels = group_order)]; s2[, Label := paste0(Term, " (", GO_ID, ")")]
s2[, Ontology := factor(Ontology, levels = c("MF", "BP", "CC"))]; setorder(s2, Group, Ontology, p_adj_BH); s2[, Label := factor(Label, levels = rev(unique(Label)))]
s2[, Group_Label := factor(as.character(Group), levels = group_order, labels = sprintf("%s (n = %s)", gsub("\n", " ", pretty[group_order]), format(nmap[group_order], big.mark = ",")))]
pS <- ggplot(s2, aes(-log10(p_adj_BH), Label, fill = Ontology)) + geom_col(width = 0.7) + geom_vline(xintercept = -log10(FDR_CUT), linetype = 3) + facet_grid(Group_Label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(MF = "#D55E00", BP = "#009E73", CC = "#56B4E9"), name = "Ontology") + geom_text(aes(label = sprintf("%d/%d", FG_Count, BG_Count)), hjust = -0.1, size = 2.4) + scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(-log[10] * " (Benjamini–Hochberg FDR)"), y = NULL, title = sprintf("All Gene Ontology terms at FDR < %.2f by %%uORF+ group", FDR_CUT), subtitle = "labels = genes in the group with the term / genes in the background with the term; dotted line = FDR 0.05; root terms not shown") +
  theme_bw(base_size = 9) + theme(strip.text.y = element_text(angle = 0, size = 7.5), axis.text.y = element_text(size = 7.5), plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
fig_save(file.path(OUT_DIR, "A24_GO_FDR_significant_terms.pdf"), go_compact(pS), h = 8); .add("  saved output/A24_GO_FDR_significant_terms.pdf")
# ---- Table S8 (same layout as before: data sheet + Legend sheet) --------------------------------------------------------
f11 <- file.path(TABS, "TableS8_GO_enrichment_by_pct_uORF_plus.xlsx")   # Supplemental Table S8
if (file.exists(f11)) { file.copy(f11, file.path(BAK, "TableS8_GO_enrichment_by_pct_uORF_plus_previous_copy.xlsx"), overwrite = TRUE); l11 <- read.xlsx(f11, sheet = "Legend") } else l11 <- data.frame(Column = character(0), Description = character(0))
d11 <- copy(x)[, Group := short[Group]][order(factor(Group, levels = short), Ontology, p_value)]
d11[, Odds_Ratio := ifelse(is.finite(Odds_Ratio), as.character(round(Odds_Ratio, 4)), "Inf")]
d11 <- d11[, .(Group, Group_size, Ontology, GO_ID, Term, FG_Count, BG_Count, Odds_Ratio, p_value, p_adj_BH)]
setdesc <- function(col, txt) { if (col %in% l11[[1]]) l11[[2]][l11[[1]] == col] <<- txt else l11 <<- rbind(l11, setNames(data.frame(col, txt), names(l11))) }
setdesc("Group", "all TuORF genes grouped by the highest %uORF+ among their TuORFs: only low (every TuORF <30% uORF+), at least one medium (highest TuORF 30-80%), at least one high (highest TuORF >=80%); CAGE-seq %uORF+")
setdesc("Group_size", "genes in the group (within the background)")
setdesc("BG_Count", sprintf("background genes with the term (background = all %s genes with a detected TuORF and Gene Ontology annotation)", format(N_BG, big.mark = ",")))
setdesc("p_adj_BH", "Benjamini-Hochberg FDR, computed within each group over the three ontologies; root terms (biological_process, molecular_function, cellular_component) are listed but not counted in the text or drawn in the figures")
setdesc("Odds_Ratio", "Fisher odds ratio ('Inf' when no background gene outside the group carries the term)")
setdesc("GO_ID", "Gene Ontology term identifier"); setdesc("Term", "Gene Ontology term"); setdesc("Ontology", "BP = biological process, MF = molecular function, CC = cellular component")   # V41 (2026-11-12): the abbreviation GO is written out in the legend and in the sheet name
wb <- createWorkbook(); addWorksheet(wb, "Gene Ontology by %uORF+ group"); writeData(wb, "Gene Ontology by %uORF+ group", d11); addWorksheet(wb, "Legend"); writeData(wb, "Legend", l11); saveWorkbook(wb, f11, overwrite = TRUE)
.add("  Table S8 rewritten: ", f11, " (", nrow(d11), " rows; previous copy in ", BAK, ")")
writeLines(.log, file.path(OUT_DIR, "A24_summary.txt"))
