################################################################################
# A42 — Supplemental Figure S22 and its sheets in Supplemental Table S12        #
# (Larry, 2026-11-20: "put the data in a figure", then "Yes, add Sorenson").    #
# (A) mRNA half-life from Sorenson et al. 2018 (cordycepin shut-off; wild-type  #
#     decay rate alpha per minute; half-life = ln 2 / alpha, capped at 24 h) in #
#     single-TuORF genes by %uORF+ group and in genes without TuORFs (as        #
#     Supplemental Fig. S21A); letters from pairwise Kolmogorov-Smirnov tests.  #
# (B) Sorenson half-life versus the 5-EU half-life (Szabo et al. 2020).         #
# (C-F) median uncapped ratio in fifths of half-life (5-EU: C, D; Sorenson:     #
#     E, F), drawn separately for thirds of mORF TE.                             #
# (G, H) median uncapped ratio in fifths of mORF TE, for thirds of 5-EU         #
#     half-life. C, E, G = GMUCT (Yu et al. 2021); D, F, H = PARE (Nagarajan et  #
#     al. 2019). Genes: single-TuORF genes and genes without TuORFs with an      #
#     uncapped ratio, the half-life of the panel, mORF TE > 0 and RNA TPM > 0.   #
#     Error bars: 95% bootstrap intervals of the median (1,000 resamples).      #
# Partial correlations: rank-based, with the other measure, log exon length and #
# RNA TPM held constant (as modules A39-A41). Group ratios for A: medians;     #
# geometric means from linear models of log10 half-life, unadjusted and         #
# adjusted for log10 exon length and log10 RNA TPM (as module A31).             #
# Outputs: output_revised_figs_20261028/FigureS22_revised_halflife_Sorenson_uncapped.pdf, #
#          output_revised_figs_20261028/png_for_docx/S22.png (200 dpi),          #
#          new sheets S22A_Sorenson_groups, S22B_halflife_agreement, S22C-H_bins, S22_partial_correlations in #
#          weak_uORF_revision/revision_tables/TableS12_decay_inclusion_splicing.xlsx (previous copy backed up), #
#          output/A42_summary.txt                                                #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A42_FigureS22_Sorenson_uncapped.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(openxlsx); library(ggplot2); library(multcompView) })
options(warn = 1); set.seed(20261120)
OUT_DIR <- "output"; FIGS <- "output_revised_figs_20261028"; CAP_MIN <- 1440; NBOOT <- 1000L
T12 <- "/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS12_decay_inclusion_splicing.xlsx"
PUP <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF"); KEY <- setNames(c("0–30% uORF+", "30–80% uORF+", "80–100% uORF+", "No TuORF"), PUP)
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
big <- function(x) trimws(format(x, big.mark = ","))
.add(sprintf("A42 Supplemental Figure S22 — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
DR <- as.data.table(read.xlsx("/Users/wu/Desktop/uORFs_miRNA/Sorenson_et_al_pnas2017_DR.xlsx")); stopifnot(nrow(DR) == 18674L)
DR <- DR[WT > 0 & sov > 0 & vcs > 0 & vcs.sov > 0]; stopifnot(nrow(DR) == 18576L)
DR[, HLs_h := pmin(log(2) / WT, CAP_MIN) / 60]
a11 <- fread(file.path(OUT_DIR, "A11_Table_genes.csv")); a12 <- fread(file.path(OUT_DIR, "A12_Table_genes.csv")); a12 <- a12[match(a11$gene_id, gene_id)]
len <- fread(file.path(OUT_DIR, "A12_gene_counts.csv"), select = c("gene_id", "exon_length"))
G <- data.table(gene_id = a11$gene_id, g = factor(a11$pct_uORF_plus_group, PUP), TE = a11$mORF_TE, RNA_TPM = a11$RNA_TPM_this_study, HL5 = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_))
G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]; G[, HLs := DR$HLs_h[match(gene_id, DR$gene_id)]]; stopifnot(!anyNA(G$g), !anyNA(G$exon_length))
# ---- (A) Sorenson half-life by group ------------------------------------------------------------------------------------
dA <- G[is.finite(HLs)]
pw <- rbindlist(lapply(combn(PUP, 2, simplify = FALSE), function(p) data.table(group1 = p[1], group2 = p[2], KS_p = suppressWarnings(ks.test(dA[g == p[1], HLs], dA[g == p[2], HLs])$p.value))))
ids <- setNames(paste0("G", 1:4), PUP); pv <- setNames(pw$KS_p, paste(ids[pw$group1], ids[pw$group2], sep = "-"))
L <- multcompLetters(pv, compare = "<", threshold = 0.05, Letters = letters)$Letters[ids]; names(L) <- PUP
ch <- unique(unlist(strsplit(unname(L), ""))); mp <- setNames(letters[seq_along(ch)], ch); L[] <- vapply(unname(L), function(s) paste(mp[strsplit(s, "")[[1]]], collapse = ""), "")
nA <- table(dA$g)[PUP]; labA <- setNames(sprintf("%s\n(n = %s) %s", KEY[PUP], big(as.integer(nA)), L), PUP)
ref <- dA[g == "No TuORF", HLs]
dm <- dA[is.finite(RNA_TPM) & RNA_TPM > 0]; dm[, g := relevel(g, ref = "No TuORF")]
m0 <- summary(lm(log10(HLs) ~ g, dm))$coefficients; m2 <- summary(lm(log10(HLs) ~ g + log10(exon_length) + log10(RNA_TPM), dm))$coefficients
SA <- rbindlist(lapply(PUP, function(gr) { x <- dA[g == gr, HLs]; r <- paste0("g", gr)
  data.table(group = KEY[[gr]], n = length(x), median_half_life_h = median(x), ratio_of_medians_vs_No_TuORF = median(x) / median(ref), KS_letter = L[[gr]],
             geo_mean_ratio_unadjusted = if (gr == "No TuORF") 1 else 10^m0[r, 1], p_unadjusted = if (gr == "No TuORF") NA else m0[r, 4],
             geo_mean_ratio_adjusted_length_abundance = if (gr == "No TuORF") 1 else 10^m2[r, 1], p_adjusted = if (gr == "No TuORF") NA else m2[r, 4]) }))
.add("(A) Sorenson half-life by group:"); for (i in seq_len(nrow(SA))) with(SA[i], .add(sprintf("  %-14s n = %5d, median %.2f h, %.2fx, letter %s; geometric %.2fx (p = %.2g), adjusted %.2fx (p = %.2g)", group, n, median_half_life_h, ratio_of_medians_vs_No_TuORF, KS_letter, geo_mean_ratio_unadjusted, p_unadjusted, geo_mean_ratio_adjusted_length_abundance, p_adjusted)))
# ---- (B) agreement --------------------------------------------------------------------------------------------------------
dB <- G[is.finite(HLs) & is.finite(HL5) & HL5 > 0]; rB <- cor(dB$HLs, dB$HL5, method = "spearman")
.add(sprintf("(B) Sorenson vs 5-EU half-life: n = %s, rho = %.2f", big(nrow(dB)), rB))
# ---- (C-H) uncapped ratio vs half-life and TE -----------------------------------------------------------------------------
prho <- function(x, y, Z) { rx <- rank(x); ry <- rank(y); RZ <- as.data.frame(lapply(Z, rank)); cor(resid(lm(rx ~ ., data = RZ)), resid(lm(ry ~ ., data = RZ))) }
cut_q <- function(x, k) cut(x, unique(quantile(x, seq(0, 1, length.out = k + 1))), include.lowest = TRUE, labels = FALSE)
medci <- function(x) { b <- replicate(NBOOT, median(sample(x, replace = TRUE))); list(n = length(x), med = median(x), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975))) }
panels <- list(C = c("U_gmuct", "HL5", "hl"), D = c("U_pare", "HL5", "hl"), E = c("U_gmuct", "HLs", "hl"), F = c("U_pare", "HLs", "hl"), G = c("U_gmuct", "HL5", "te"), H = c("U_pare", "HL5", "te"))
SET <- c(U_gmuct = "GMUCT", U_pare = "PARE"); HLN <- c(HL5 = "5-EU", HLs = "Sorenson")
bins <- list(); pcor <- list(); plt <- list()
for (k in names(panels)) {
  u <- panels[[k]][1]; h <- panels[[k]][2]; mode <- panels[[k]][3]
  d <- G[is.finite(get(u)) & get(u) > 0 & is.finite(get(h)) & get(h) > 0 & is.finite(TE) & TE > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]
  d[, `:=`(U = get(u), H = get(h))]
  ph <- prho(d$U, d$H, list(d$TE, d$exon_length, d$RNA_TPM)); pt <- prho(d$U, d$TE, list(d$H, d$exon_length, d$RNA_TPM))
  pcor[[k]] <- data.table(panel = k, uncapped = SET[[u]], half_life = HLN[[h]], n = nrow(d), rho_half_life = cor(d$U, d$H, method = "spearman"), partial_rho_half_life = ph,
                          rho_TE = cor(d$U, d$TE, method = "spearman"), partial_rho_TE = pt, rho_half_life_vs_TE = cor(d$H, d$TE, method = "spearman"))
  if (mode == "hl") { d[, bin := cut_q(H, 5)]; d[, strat := factor(cut_q(TE, 3), 1:3, c("Low TE", "Middle TE", "High TE"))]; xmed <- d[, median(H), by = bin]; xl <- sprintf("%s half-life (h)", HLN[[h]]); pal <- c("#9ecae1", "#4292c6", "#08306b")
  } else { d[, bin := cut_q(TE, 5)]; d[, strat := factor(cut_q(H, 3), 1:3, c("Short", "Middle", "Long"))]; xmed <- d[, median(TE), by = bin]; xl <- "mORF TE"; pal <- c("#fdae6b", "#e6550d", "#7f2704") }
  B <- d[, medci(U), by = .(strat, bin)]; B[, x := xmed$V1[match(bin, xmed$bin)]]
  bins[[k]] <- cbind(panel = k, uncapped = SET[[u]], half_life = HLN[[h]], grouping = if (mode == "hl") "fifths of half-life within thirds of mORF TE" else "fifths of mORF TE within thirds of half-life", B)
  fm <- function(x) sub("^-", "\u2212", sprintf("%.2f", x))
  sub <- if (mode == "hl") sprintf("%s; %s half-lives; partial rho = %s", SET[[u]], HLN[[h]], fm(ph)) else sprintf("%s; thirds of 5-EU half-life; partial rho = %s", SET[[u]], fm(pt))
  xlab <- if (mode == "hl") sprintf("%.1f", sort(xmed$V1)) else sprintf("%.2f", sort(xmed$V1))
  dg <- position_dodge(width = 0.35)
  plt[[k]] <- ggplot(B, aes(factor(bin), med, colour = strat, group = strat)) + geom_line(position = dg, linewidth = 0.6) + geom_point(position = dg, size = 1.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.4, position = dg) + scale_x_discrete(labels = xlab) + scale_colour_manual(values = pal, name = if (mode == "hl") NULL else "Half-life:") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) + labs(title = k, subtitle = sub, x = paste0(if (mode == "hl") "Half-life (h)" else "mORF TE", ", median of each fifth"), y = "Uncapped ratio") + theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
}
PC <- rbindlist(pcor); .add("(C-H) correlations:"); for (i in seq_len(nrow(PC))) with(PC[i], .add(sprintf("  %s %-5s %-8s n = %6s, rho half-life %.2f (partial %.2f), rho TE %.2f (partial %.2f), rho half-life vs TE %.2f", panel, uncapped, half_life, big(n), rho_half_life, partial_rho_half_life, rho_TE, partial_rho_TE, rho_half_life_vs_TE)))
BB <- rbindlist(bins, fill = TRUE); .add(sprintf("  genes per point: GMUCT %d-%d, PARE %d-%d", min(BB[uncapped == "GMUCT", n]), max(BB[uncapped == "GMUCT", n]), min(BB[uncapped == "PARE", n]), max(BB[uncapped == "PARE", n])))
for (k in names(panels)) for (s in unique(BB[panel == k, strat])) .add(sprintf("    %s %-16s medians: %s", k, s, paste(sprintf("%.2f", BB[panel == k & strat == s][order(bin), med]), collapse = ", ")))
# ---- figure -----------------------------------------------------------------------------------------------------------------
source("revision_modules/fig_compact.R")
pal_pup <- c("#619CFF", "#00BA38", "#F8766D", "grey40")
dA[, lab := factor(labA[as.character(g)], labA[PUP])]
pA <- ggplot(dA, aes(HLs, colour = lab)) + stat_ecdf(geom = "step", linewidth = 0.8) + scale_colour_manual(values = setNames(pal_pup, labA[PUP]), name = NULL) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) + coord_cartesian(xlim = c(0.1, 24)) +
  labs(title = "A", subtitle = "Sorenson et al. 2018 half-lives of single-TuORF genes by %uORF+ group", x = "mRNA half-life (h)", y = "Fraction of genes") + theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
pB <- ggplot(dB, aes(HL5, HLs)) + geom_bin2d(bins = 45) + scale_fill_gradient(low = "grey85", high = "black", name = "Genes", trans = "log10") + scale_x_log10() + scale_y_log10() +
  labs(title = "B", subtitle = sprintf("Both half-lives (n = %s); Spearman's rho = %.2f", big(nrow(dB)), rB), x = "5-EU half-life (h)", y = "Sorenson half-life (h)") + theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
cc <- list(legend_in(compact(pA, wrap_sub = 60), 0.02, 0.97, wrap = function(x) x), legend_in(compact(pB, wrap_sub = 60), 0.99, 0.03, title = TRUE))
for (k in names(plt)) cc[[length(cc) + 1]] <- legend_in(compact(plt[[k]], wrap_sub = 60), 0.01, 0.99, horizontal = TRUE, title = panels[[k]][3] == "te")   # one-line key in the empty band at the top (the y axis has 30% headroom)
g <- fig_grid(cc, ncol = 2, top = fig_title("Figure S22. mRNA half-lives from a second dataset, and uncapped mRNA versus half-life and mORF TE", wrap = 64), h = 9)
fpdf <- file.path(FIGS, "FigureS22_revised_halflife_Sorenson_uncapped.pdf"); fig_save(fpdf, g, h = 9)
system2("pdftoppm", c("-r", "200", "-png", "-singlefile", fpdf, file.path(FIGS, "png_for_docx", "S22")))
.add("  saved ", fpdf, " and png_for_docx/S22.png")
# ---- Table S12 sheets -------------------------------------------------------------------------------------------------------
bak <- file.path(dirname(T12), "_replaced_by_V49_S22_20261120"); dir.create(bak, showWarnings = FALSE)
if (!file.exists(file.path(bak, basename(T12)))) file.copy(T12, file.path(bak, basename(T12)))
wb <- loadWorkbook(file.path(OUT_DIR, "A31_Table_decay_inclusion_splicing.xlsx"))   # V53 (2026-11-22): start from the table that module A31 wrote (until V52: from the copy of 2026-11-20 in the folder above, which a later run of A31 could not change)
add_sheet <- function(nm, df) { if (nm %in% names(wb)) removeWorksheet(wb, nm); addWorksheet(wb, nm); writeData(wb, nm, df) }
add_sheet("S22A_Sorenson_groups", SA)
add_sheet("S22B_halflife_agreement", data.table(n_genes = nrow(dB), spearman_rho = rB))
add_sheet("S22C-H_bins", BB[, .(panel, uncapped, half_life, grouping, stratum = strat, fifth = bin, x_median_of_fifth = x, n_genes = n, median_uncapped_ratio = med, ci95_lo = lo, ci95_hi = hi)])
add_sheet("S22_partial_correlations", PC)
rd <- read.xlsx(wb, sheet = "README", colNames = FALSE)[[1]]
rd <- rd[!grepl("^Supplemental Figure S22 sheets", rd)]
rd <- append(rd, paste0("Supplemental Figure S22 sheets: S22A_Sorenson_groups (half-lives of Sorenson et al. 2018: wild-type decay rate alpha per minute, half-life = ln 2 / alpha, capped at 24 h, in hours; ",
  "n, median, ratio of medians relative to genes without TuORFs, letters from pairwise Kolmogorov-Smirnov tests at p < 0.05, ratios of geometric means from linear models of log10 half-life, unadjusted and adjusted for log10 exon length and log10 RNA TPM); ",
  "S22B_halflife_agreement (Spearman rho of the two half-life datasets); S22C-H_bins (median uncapped ratio with 95% bootstrap interval, 1,000 resamples, in fifths of half-life within thirds of mORF TE, or in fifths of mORF TE within thirds of 5-EU half-life); ",
  "S22_partial_correlations (Spearman rho and rank-based partial correlations with the other measure, exon length and RNA TPM held constant). Code: revision_modules/A42_FigureS22_Sorenson_uncapped.R."), after = length(rd) - 1)
removeWorksheet(wb, "README"); addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = rd), colNames = FALSE)
sh <- names(wb); worksheetOrder(wb) <- c(which(sh == "README"), which(sh != "README"))
saveWorkbook(wb, T12, overwrite = TRUE); .add("  Table S12: sheets S22A_Sorenson_groups, S22B_halflife_agreement, S22C-H_bins, S22_partial_correlations; README line added (previous copy in ", bak, ")")
writeLines(.log, file.path(OUT_DIR, "A42_summary.txt"))
