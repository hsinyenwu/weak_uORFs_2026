#!/usr/bin/env python3
"""(2026-11-07, V21: the supplemental pages called S10, S11 and S12 below are S11, S12 and S13 since the schematic of the leader-splicing audit became S8.)
Build a driver that re-renders the manuscript's data panels (Figs 1B-C/S2, 2B-F, 4D-F, 5, 6 DLA + qPCR, S10, S11 D-E; supplemental numbers of V14, 2026-11-03: the mORF panels (S11 of V11-V13, submitted S4) are still drawn as former_S11_mORF_plus_panels_A-B.pdf but are no longer in the paper)
from the canonical script with (i) manuscript panel letters, (ii) the shared-term null annotated on Fig 1B,
(iii) outputs in a separate folder. Since 2026-10-28 the canonical script itself carries the high / medium / low
uIF class labels and the %uORF+ / No TuORF wording, so the label substitutions below are no-ops kept for older copies.
The canonical script itself is not modified.  usage: python3 build_revised_figure_driver.py SCRIPT DRIVER OUTDIR"""
import sys, re
src, dst, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
L = open(src, encoding='utf8').read().split('\n')
def block(a, b): return '\n'.join(L[a-1:b])          # 1-based inclusive
def sub1(txt, old, new, n=1):
    assert txt.count(old) == n, (old, txt.count(old)); return txt.replace(old, new)
def opt(txt, old, new):                                  # replace if present (label already renamed in the canonical script otherwise)
    return txt.replace(old, new)

base = block(1, 1223)
base = sub1(base, 'OUT_DIR_MAIN <- "output"', f'OUT_DIR_MAIN <- "{outdir}"')
base = sub1(base, 'OUT_DIR <- OUT_DIR_MAIN', 'OUT_DIR <- OUT_DIR_MAIN\nsource("revision_modules/fig_compact.R")   # 6.5 x <= 9 in pages, Arial, >= 8 pt (2026-10-28)')   # the base block ends inside a multi-line statement, so the helper is sourced here
fig1 = block(1224, 1303)
fig1 = sub1(fig1, 'ggplot2::labs(title = "C", x = "uIF", y = "Number of genes")', 'ggplot2::labs(title = "uIF distribution (Figure 1B)", x = "uIF", y = "Number of genes")')   # V40 (2026-11-11): this panel of a working page was titled "Figure S2"
fig1 = sub1(fig1, 'ggplot2::labs(title = "D", x = bquote(log[2]~"(uIF)"), y = bquote(log[2]~"(mORF TE)"))', 'ggplot2::labs(title = "B", x = bquote(log[2]~"(uIF)"), y = bquote(log[2]~"(mORF TE)"))')
fig1 = sub1(fig1, 'ggplot2::labs(title = "E", x = "mORF TE", y = "Fraction", color = "")', 'ggplot2::labs(title = "C", x = "mORF TE", y = "Fraction", color = "")')
# V15 (2026-11-04, Larry): the note on the built-in part of the uIF-TE correlation is no longer drawn on panel B (the null analysis was removed from the manuscript)
fig1 = sub1(fig1, 'save_pdf(Figure1_CDE, "Figure1_CDE.pdf", w = 10, h = 10)', 'save_pdf(Figure1_CDE, "Figure1_revised_panels_B_C_and_S2.pdf", w = 10, h = 10)\n# The uIF distribution as its own file (a working page; it is Figure 1B since V16) and the compact single panels B and C for the Figure 1 composite (2026-10-28; 6.5-in page, Arial, >= 8 pt)\nFigS2_alone <- ggpubr::annotate_figure(ggpubr::ggarrange(compact(pC + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) + ggplot2::labs(title = NULL, subtitle = sprintf("x axis truncated at uIF = 8; %d genes with uIF > 8 are not drawn", sum(OUT_df6$uIF > 8, na.rm = TRUE)))), ncol = 1, nrow = 1), top = fig_title("Not a manuscript figure: distribution of uIF among single-TuORF genes (Figure 1B)", wrap = 44))\nfig_save(file.path(OUT_DIR, "FigureS2_revised_uIF_distribution.pdf"), FigS2_alone, w = 4.2, h = 3.4)\npD_c <- pD; pD_c$layers <- Filter(function(l) !(inherits(l$geom, "GeomText") && is.character(l$aes_params$label) && grepl("shared mORF", l$aes_params$label)), pD_c$layers); pD_c <- pD_c + ggplot2::labs(subtitle = NULL)   # V15: no built-in-term note\nfig_save(file.path(OUT_DIR, "Figure1_compact_B.pdf"), compact(pD_c, geom_scale = 0.45, wrap_sub = 40), w = 2.25, h = 2.75)\nfig_save(file.path(OUT_DIR, "Figure1_compact_C.pdf"), compact(pE, geom_scale = 0.45), w = 2.25, h = 2.75)   # legend goes below the plot (compact)\n# V16 (2026-11-04, Larry): the uIF distribution becomes Figure 1B (compact single panel) and the uIF-TE scatter + class ECDF become Supplemental Fig. S2\nfig_save(file.path(OUT_DIR, "Figure1_compact_B_dist.pdf"), compact(pC + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) + ggplot2::labs(title = "B", subtitle = sprintf("x axis truncated at uIF = 8; %d genes with uIF > 8 are not drawn", sum(OUT_df6$uIF > 8, na.rm = TRUE))), geom_scale = 0.45, wrap_sub = 44), w = 3.4, h = 1.9)   # V21 (2026-11-07, Larry): same height as panel A\n# V30 (2026-11-09): Supplemental Fig. S2 (uIF versus TE) is drawn by revision_modules/fig_main_layout.R (legend inside the plot, aligned panels)')
fig2 = block(1304, 1532)
fig2 = fig2.replace('"%pUP (CAGE)"', '"%uORF+ (CAGE)"')
fig2 = fig2.replace('ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")', 'ggpubr::stat_cor(method = "pearson", label.x.npc = 0.3, label.y.npc = "top")', 1)    # D: not over the points at %uORF+ = 0
fig2 = fig2.replace('ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")', 'ggpubr::stat_cor(method = "pearson", label.x.npc = 0.12, label.y.npc = "top")', 1)   # E: the longer label has to end inside the panel
assert 'Figure3_BCDEFGH.pdf' in fig2
fig2 = fig2.replace('Figure3_BCDEFGH.pdf', 'Figure2_revised_panels_B-F_and_S11_G-H.pdf')
fig2 = sub1(fig2, 'save_pdf(Final_Fig3, "Figure2_revised_panels_B-F_and_S11_G-H.pdf", w = 10, h = 18)', 'save_pdf(Final_Fig3, "Figure2_revised_panels_B-F_and_S11_G-H.pdf", w = 10, h = 18)\n# The mORF figure that was removed from the manuscript in V14 as its own file (a working page) and compact single panels B-F for the Figure 2 composite (2026-10-28)\nFigS11_alone <- ggpubr::annotate_figure(ggpubr::ggarrange(compact(p3H + ggplot2::labs(title = "A")), compact(p3G + ggplot2::labs(title = "B")), ncol = 2, nrow = 1), top = fig_title("Not a manuscript figure (removed in V14): mORFs are less affected by TSS heterogeneity"))\nfig_save(file.path(OUT_DIR, "former_S11_mORF_plus_panels_A-B.pdf"), FigS11_alone, w = 6.5, h = 3.3)\np3B_c <- p3B; for (i in seq_along(p3B_c$layers)) if (!is.null(p3B_c$layers[[i]]$aes_params$label)) p3B_c$layers[[i]]$aes_params$label <- sub("(n=", "(", p3B_c$layers[[i]]$aes_params$label, fixed = TRUE)   # shorter count labels for the 2.15-in panel\nfor (i in seq_along(p3B_c$layers)) if (!is.null(p3B_c$layers[[i]]$aes_params$label) && !is.null(p3B_c$layers[[i]]$data$x)) p3B_c$layers[[i]]$data$x <- ifelse(p3B_c$layers[[i]]$data$x > 80, 86, ifelse(p3B_c$layers[[i]]$data$x > 40, 50, p3B_c$layers[[i]]$data$x))   # spread the three count labels across the 2.15-in panel\n# V30 (2026-11-09): the Figure 2 panels B-F are laid out as one aligned grid by revision_modules/fig_main_layout.R')
fig4 = block(2262, 2387)
fig4 = sub1(fig4, 'Plot_Data_uORF$ORF_Group <- factor(Plot_Data_uORF$ORF_Group, levels = c("CPuORF", "TuORF"))',
            'Plot_Data_uORF$ORF_Group <- factor(Plot_Data_uORF$ORF_Group, levels = c("CPuORF", "TuORF"))\n'
            'message("CPuORF module: share of genes with >= 90% uORF+ by group: ", paste(names(tapply(Plot_Data_uORF$pUP >= 90, Plot_Data_uORF$ORF_Group, mean)), sprintf("%d/%d = %.1f%%", tapply(Plot_Data_uORF$pUP >= 90, Plot_Data_uORF$ORF_Group, sum), tapply(Plot_Data_uORF$pUP >= 90, Plot_Data_uORF$ORF_Group, length), 100 * tapply(Plot_Data_uORF$pUP >= 90, Plot_Data_uORF$ORF_Group, mean)), collapse = " | "))')
assert fig4.count('Figure5_Robust_Matching_Stats.pdf') >= 1
fig4 = fig4.replace('Figure5_Robust_Matching_Stats.pdf', 'Figure4_revised_panels_D-F.pdf')
fig4 = sub1(fig4, 'save_pdf(Figure5, "Figure4_revised_panels_D-F.pdf", w=14, h=5)', 'save_pdf(Figure5, "Figure4_revised_panels_D-F.pdf", w=14, h=5)\n'
            '# V30 (2026-11-09): the Figure 4 panels D-F are laid out as one aligned row by revision_modules/fig_main_layout.R')
after5 = block(2388, 2422)
figATG = block(2424, 2818)   # ATG frequency across the TSS (manuscript Figure 3A-B) and the TSS sequence logos (Supplemental Figure S11 since V39)
figATG = sub1(figATG, '  save_pdf(fig6_ab_annotated, "Figure6_AB_50pct_CAGE_threshold.pdf", w = 14, h = 7)',
              '  save_pdf(fig6_ab_annotated, "Figure6_AB_50pct_CAGE_threshold.pdf", w = 14, h = 7)\n'
              '  # compact Figure 3A-B for the 6.5-in composite (2026-10-28)\n'
              '  # V30 (2026-11-09): Figure 3A-B is laid out by revision_modules/fig_main_layout.R (legend of A inside the plot, aligned panels)')
figATG = sub1(figATG, '  save_pdf(fig6_logos, "Figure6_All_SeqLogos.pdf", w = 12, h = 15)',
              '  save_pdf(fig6_logos, "Figure6_All_SeqLogos.pdf", w = 12, h = 15)\n'
              '  # Supplemental Figure S11 since V39 (S13 in V21-V38; the 50% CAGE site logos) on a 6.5-in page (2026-10-28)\n'
              '  logo_c <- function(p) suppressMessages(compact(p) + ggplot2::scale_x_continuous(breaks = seq(1, 51, 10), labels = c("-40", "-30", "-20", "-10", "+1", "+11")) + ggplot2::labs(x = "Position relative to the TSS (nt)"))   # the logos show window positions 21-71 = -40..+10 around the 50%-CAGE TSS (+1)   # the plots carry their own titles; fewer position ticks\n'
              '  S12_compact <- ggpubr::annotate_figure(ggpubr::ggarrange(logo_c(res_50$pC), logo_c(res_50$pD), ncol = 2), top = fig_title("Figure S11. Sequence context of transcription start sites in TuORF and No-TuORF genes", wrap = 58))\n'
              '  fig_save(file.path(OUT_DIR, "FigureS13_compact_logos_50pct.pdf"), S12_compact, h = 2.9)')
# V23 (2026-11-08): Figure 3A-B per gene. One transcript per gene (the one with the most CAGE signal, as module A7), so that the chi-square tests
# count genes; upstream versus downstream of the same genes is a paired comparison (Stuart-Maxwell test of marginal homogeneity; the V of the
# 2 x 3 table is kept as the effect size); axis labels name the TSS reference point and genes.
figATG = sub1(figATG, '  dt_valid <- dt_point[start >= 61 & (tx_len - start) >= 62]\n',
              '  dt_valid <- dt_point[start >= 61 & (tx_len - start) >= 62]\n'
              '  dt_valid[, gene_v23 := toupper(get_gene_id_from_tx(tx_name))]; data.table::setorder(dt_valid, gene_v23, -total_sum); dt_valid <- dt_valid[!duplicated(gene_v23)]   # V23: one transcript per gene (most CAGE signal)\n')
figATG = sub1(figATG, 'chisq_one <- function(tbl) {',
              'stuart_maxwell <- function(x, y, lev = c("0", "1", "2+")) {   # V23: paired 3-category comparison (marginal homogeneity)\n'
              '  tb <- table(factor(x, lev), factor(y, lev)); k <- nrow(tb); rs <- rowSums(tb); cs <- colSums(tb); d <- (rs - cs)[-k]\n'
              '  S <- matrix(0, k - 1, k - 1); for (i in 1:(k - 1)) for (j in 1:(k - 1)) S[i, j] <- if (i == j) rs[i] + cs[i] - 2 * tb[i, i] else -(tb[i, j] + tb[j, i])\n'
              '  st <- as.numeric(t(d) %*% solve(S) %*% d); list(chi2 = st, df = k - 1, p = stats::pchisq(st, k - 1, lower.tail = FALSE)) }\n'
              'chisq_one <- function(tbl) {')
figATG = sub1(figATG, '    out <- chisq_one(tbl)\n',
              '    out <- chisq_one(tbl)\n'
              '    sc1 <- sub("^.*_", "", gp[1]); sc2 <- sub("^.*_", "", gp[2])\n'
              '    if (sc1 == sc2) { m <- stats_cat_vec == sc1; cg <- function(v) as.character(cut(v, c(-Inf, 0, 1, Inf), labels = c("0", "1", "2+")))\n'
              '      sm <- stuart_maxwell(cg(up_counts[m]), cg(down_counts[m])); out$chi2 <- sm$chi2; out$df <- sm$df; out$p <- sm$p }   # V23: the same genes upstream and downstream: paired test\n')
figATG = sub1(figATG, 'x = "Position relative to threshold point", y = "% Transcripts")', 'x = "Position relative to the TSS reference point (nt)", y = "% of genes")')
figATG = sub1(figATG, 'ggplot2::labs(title = "B: ATG Count (Chi-square / Cramer\'s V)", y = "% Transcripts", x = "")', 'ggplot2::labs(title = "B: ATG count", y = "% of genes", x = "")')
# V23 (2026-11-08): a TuORF called from an AUG one codon upstream (one extra N-terminal Met) is the same CPuORF (SAMDC AT5G15950, ATHB1 AT3G01470)
fig4 = sub1(fig4, '    ORF_Group  = ifelse(key_strict %in% keys_strict | key_length %in% keys_length, "CPuORF", "TuORF"),',
            '    ORF_Group  = ifelse(key_strict %in% keys_strict | key_length %in% keys_length | (substr(pep_clean, 1, 2) == "MM" & paste0(trimws(gene_id), "_", substring(pep_clean, 2)) %in% keys_strict), "CPuORF", "TuORF"),   # V23: one extra N-terminal Met')
fig6 = block(2819, 3015)
fig6 = fig6.replace('"DLA_Analysis_Plot_WuORFs_UuORFs.pdf"', '"Figure6_revised_A_DLA_WT_vs_Mut.pdf"').replace('"DLA_Perfect_Kozak.pdf"', '"Figure6_revised_B_DLA_perfect_Kozak.pdf"')
fig6 = sub1(fig6, 'save_pdf(p7_4, "qPCR_perfect_Kozak.pdf", w = 10, h = 15)', 'save_pdf(p7_4, "qPCR_perfect_Kozak.pdf", w = 10, h = 15)\n'
            '# manuscript Figure 6 on one 6.5 x 9 in page (2026-10-28): the same boxplots with genes as columns and biological repeats as rows\n'
            'refacet <- function(p, ttl, yl) { nf <- nrow(unique(p$data[, c("Repeat", "Gene")])); compact(p + ggplot2::facet_wrap(~ Repeat + Gene, scales = "free_y", ncol = if (nf > 7) ceiling(nf / 2) else nf, labeller = ggplot2::labeller(.multi_line = FALSE)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.32))) + ggplot2::labs(title = ttl, subtitle = NULL, x = NULL, y = yl) + ggplot2::theme(legend.position = "none"), geom_scale = 0.4) }\n'
            'Fig6_compact <- ggpubr::ggarrange(refacet(p7_1, "A  Dual-luciferase assay: wild-type versus mutated uORF", "FLUC/REN activity"), refacet(p7_3, "B  Dual-luciferase assay: optimized Kozak context", "FLUC/REN activity"), refacet(p7_2, "C  Reporter mRNA (RT-qPCR): wild-type versus mutated uORF", "FLUC/REN mRNA"), refacet(p7_4, "D  Reporter mRNA (RT-qPCR): optimized Kozak context", "FLUC/REN mRNA"), ncol = 1, nrow = 4)\n'
            'fig_save(file.path(OUT_DIR, "Figure6_compact_panels_A-D.pdf"), Fig6_compact, h = 9)')
# V23 (2026-11-08): the points of Figure 6 are jittered horizontally only and with a fixed seed (before, every render moved them, also vertically);
# the RT-qPCR panels C-D carry no test (three technical replicates of one RNA sample per construct and experiment).
fig6 = fig6.replace('ggplot2::geom_jitter(width = 0.2, size = 1.5, alpha = 0.5)', 'ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.2, height = 0, seed = 1), size = 1.5, alpha = 0.5)')
assert fig6.count('position_jitter(width = 0.2, height = 0, seed = 1)') == 4
for nm, lst in (('p7_2', 'comparison_list_7'), ('p7_4', 'comparison_list_7_4')):
    i0 = fig6.index(nm + ' <- ggplot2::ggplot('); i1 = fig6.index('  ggpubr::stat_compare_means(', i0); i2 = fig6.index('  ) +\n', i1) + len('  ) +\n')
    assert ('comparisons = ' + lst + ',') in fig6[i1:i2]
    fig6 = fig6[:i1] + fig6[i2:]
fig5 = block(3016, 3695)
fig5 = sub1(fig5, 'save_pdf(Figure_8_Final, "Figure8_ABCDEFGH_Full_Violin.pdf", w = 12, h = 20)', 'save_pdf(Figure_8_Final, "Figure8_ABCDEFGH_Full_Violin.pdf", w = 12, h = 20)\n'
            '# manuscript Figure 5 = panels A-F on a 6.5-in page (2026-10-28)\n'
            'p8D <- p8D + ggplot2::scale_y_continuous(limits = c(0, 100)); p8D$layers <- lapply(p8D$layers, function(l) { if (inherits(l$geom, "GeomText")) l$mapping$label <- quote(sprintf("%.1f%%", Percent)); l })   # one decimal for every bar label\n# V30 (2026-11-09): Figure 5 is laid out by revision_modules/fig_main_layout.R (legend of C inside the plot, aligned panels)')
fig5 = opt(fig5, 'plot_title = "A: All uORFs vs No uORF"', 'plot_title = "A: All TuORFs vs No TuORF"')
fig5 = opt(fig5, 'plot_title = "B: Strong uORFs (pUP > 30) vs No uORF"', 'plot_title = "B: TuORFs with >30% uORF+ vs No TuORF"')
fig5 = fig5.replace('"No uORF"', '"No TuORF"')
fig5 = opt(fig5, '"D: % Removed by 30% pUP Cutoff"', '"D: % removed by the 30% uORF+ cutoff"')
fig5 = opt(fig5, '"E: pUP Distribution"', '"E: %uORF+ distribution"')
fig5 = opt(fig5, 'plot_title = "G: log2(uIF) vs Length (Strong uORFs >30% pUP)"', 'plot_title = "G: log2(uIF) vs Length (TuORFs >30% uORF+)"')
fig5 = fig5.replace('y = "% pUP"', 'y = "%uORF+ (CAGE)"')
assert 'Figure8_ABCDEFGH_Full_Violin.pdf' in fig5
fig5 = fig5.replace('Figure8_ABCDEFGH_Full_Violin.pdf', 'Figure5_revised_panels_A-H.pdf')
# Supplemental Figs S10 and S12 D-E regenerated from the canonical code (2026-11-01, Larry): the undetected-uORF module (canonical FIGURE 4,
# lines 1713-2263) and the uORF-order block (canonical FIGURE 11, lines 3936-4145) run after the main-figure blocks; the order function
# returns its panels, and revision_modules/fig_S10_S12_compact.R draws the 6.5-in pages (S10 vector; S12 = submitted browser panels A-C + new D-E).
uu = block(1713, 2263)
fig11 = block(3936, 4145)
fig11 = sub1(fig11, '  message("Saved: ", final_filename)\n}', '  message("Saved: ", final_filename)\n'
             '  invisible(list(n_uorfs = n_uorfs, data = df_subset, pUP_ecdf = pD, pUP_box = pE, uIF_ecdf = pA, uIF_box = pB, letters_pUP = stats_pup, letters_uIF = stats_uif))   # panels for Supplemental Fig. S10 (2026-11-01)\n}')
fig11 = sub1(fig11, 'generate_uorf_order_fig(df_ranked_final, 2, "2_uORFs")\ngenerate_uorf_order_fig(df_ranked_final, 3, "3_uORFs")\ngenerate_uorf_order_fig(df_ranked_final, 4, "4_uORFs")',
             'S10_rows <- list(generate_uorf_order_fig(df_ranked_final, 2, "2_uORFs"), generate_uorf_order_fig(df_ranked_final, 3, "3_uORFs"), generate_uorf_order_fig(df_ranked_final, 4, "4_uORFs"))')
supp = ('saveRDS(list(S10_rows = S10_rows, p_S5D_hist = p_S5D_hist, TE_Uu_Final_4 = TE_Uu_Final_4, letters_TE_4 = letters_TE_4, letters_RNA_4 = letters_RNA_4, countsBC_4 = countsBC_4),\n'
        '        file.path(OUT_DIR, "S10_S12_objects.rds"))   # for layout work on the two pages without rerunning the driver\n'
        'source("revision_modules/fig_S10_S12_compact.R")\n'
        '# V30 (2026-11-09, Larry): legends inside the plots or at the side and aligned panels for Figures 2-5 and Supplemental Fig. S2; the plot objects are saved so that the layout can be redrawn without rerunning the driver\n'
        'saveRDS(list(p3B_c = p3B_c, p3C = p3C, p3D = p3D, p3E = p3E, p3F = p3F, p5D = p5D, p5E = p5E, p5F = p5F, f3A = res_50$pA, f3B = res_50$pB, p8A = p8A, p8B = p8B, p8C = p8C, p8D = p8D, p8E = p8E, p8H = p8H, pD_c = pD_c, pE = pE, df_len_H = df_len_H),   # df_len_H: data that a text layer of p8H reads from the global environment\n'
        '        file.path(OUT_DIR, "main_layout_objects.rds"))\n'
        'source("revision_modules/fig_main_layout.R")')
driver = '\n'.join([base, fig1, fig2, fig4, after5, figATG, fig6, fig5, uu, fig11, supp, 'message("Revised figure driver: done.")'])
# global label renaming (uIF classes): no 'strong / weak' wording anywhere
for old, new in [('"Strong uORF (>1)"', '"High uIF (>1)"'), ('"Mid uORF (>0.2 - 1)"', '"Medium uIF (0.2-1)"'), ('"Weak uORF (0 - 0.2)"', '"Low uIF (<=0.2)"'),
                 ('"Strong (uIF > 1)"', '"High (uIF > 1)"'), ('"Weak (uIF <= 0.2)"', '"Low (uIF <= 0.2)"'),
                 ('"Strong (>1)"', '"High (>1)"'), ('"Mid (>0.2 - 1)"', '"Medium (0.2-1)"'), ('"Weak (0 - 0.2)"', '"Low (<=0.2)"')]:
    driver = driver.replace(old, new)
open(dst, 'w', encoding='utf8').write(driver)
print('driver written:', dst, len(driver.split('\n')), 'lines; remaining strong/weak label strings:', len(re.findall(r'"(Strong|Weak)[ (]', driver)))
