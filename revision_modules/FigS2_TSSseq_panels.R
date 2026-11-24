################################################################################
# Figure S2 — the Figure 2B–F panels computed with TSS-seq instead of CAGE-seq   #
# (the "TSS_SOURCE = TSS" companion run of the canonical script, redrawn as a    #
# stand-alone module from the current tables: corrected %uORF+ window, latest    #
# TuORF catalogue, renamed uIF classes, manuscript panel letters A–E).           #
#                                                                                #
# Inputs: output/Table_S2_single_TuORF_genes.xlsx (cohort, uIF classes, TE),     #
#   output/Table_S1_all_TuORFs.xlsx sheet "Table S1B (TSS-seq)" (%uORF+ from     #
#   TSS-seq, Nielsen et al. 2019), RSEM gene tables + genes_w_5UTR + S1A/ouORF   #
#   tables for the No-TuORF genes of Figure 1C / 2F.                             #
# Outputs (output/): FigureS2_TSSseq_panels_A-E.pdf, FigS2_TSSseq_summary.txt,   #
#   FigS2_TSSseq_values.csv                                                      #
# Run from Weak_uORF_submission/:  Rscript revision_modules/FigS2_TSSseq_panels.R#
################################################################################
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(multcompView); library(ggpubr) })
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_S1 <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"); PATH_S2 <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_RSEM_RNA  <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_5UTR <- file.path(DATA_DIR, "genes_w_5UTR.tsv")
for (f in c(PATH_S1, PATH_S2, PATH_OUORF, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_5UTR)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
UIF_LV <- c("High uIF (>1)", "Medium uIF (0.2-1)", "Low uIF (<=0.2)"); PUP_LV <- c("80~100", "30~80", "0~30", "No TuORF")

# ---- data --------------------------------------------------------------------
S2 <- as.data.table(read.xlsx(PATH_S2, 1)); S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)"); S1B <- as.data.table(read.xlsx(PATH_S1, sheet = "Table S1B (TSS-seq)"))
ou <- read.xlsx(PATH_OUORF)
B <- S1B[gene_id %in% S2$gene_id, .(gene_id, pUP = pct_uORF_plus)]; stopifnot(!anyDuplicated(B$gene_id), nrow(B) == nrow(S2))
df <- merge(S2[, .(gene_id, uIF = round(ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, uIF), 3), uIF_raw = ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, uIF), TE = mORF_TE, pUP_cage = pct_uORF_plus)], B, by = "gene_id")   # V22: rounded to 3 decimals as in Table S3 and Figure 2C (classes 1,116 / 991 / 888)
df[, uIF_Group := factor(fcase(uIF > 1, UIF_LV[1], uIF > 0.2, UIF_LV[2], default = UIF_LV[3]), UIF_LV)]
df[, pUP_Group := factor(fcase(pUP >= 80, "80~100", pUP >= 30, "30~80", default = "0~30"), PUP_LV[1:3])]
rna <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); ribo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
ctrl <- data.table(gene_id = rna$gene_id, RNA_TPM = rna$TPM, Ribo_TPM = ribo$TPM[match(rna$gene_id, ribo$gene_id)])[, TE := Ribo_TPM / RNA_TPM]
ctrl <- ctrl[gene_id %in% read.delim(PATH_5UTR)$gene_id & RNA_TPM > 1 & is.finite(TE) & !gene_id %in% unique(c(S1A$gene_id, as.character(ou$gene_id)))]
dfF <- rbind(df[, .(gene_id, TE, pUP_Group = as.character(pUP_Group))], ctrl[, .(gene_id, TE, pUP_Group = "No TuORF")])[is.finite(TE)]
dfF[, pUP_Group := factor(pUP_Group, PUP_LV)]
.add(sprintf("Figure S8 (S10 in V46-V47, S9 in V45, S8 in V39-V44, S10 until V38, S2 as submitted; TSS-seq): %d single-TuORF genes (%d / %d / %d high / medium / low uIF); No-TuORF genes with TE: %d", nrow(df), sum(df$uIF_Group == UIF_LV[1]), sum(df$uIF_Group == UIF_LV[2]), sum(df$uIF_Group == UIF_LV[3]), nrow(ctrl)))

# ---- letters (pairwise KS, p < 0.05; letters re-mapped per character, shared letters kept) ----
ks_labels <- function(d, grp, val) {
  d <- d[is.finite(d[[val]]), ]; lv <- levels(d[[grp]]); cmb <- combn(lv, 2, simplify = FALSE)
  p <- sapply(cmb, function(g) suppressWarnings(ks.test(d[[val]][d[[grp]] == g[1]], d[[val]][d[[grp]] == g[2]])$p.value))
  ids <- setNames(paste0("G", seq_along(lv)), lv); names(p) <- sapply(cmb, function(g) paste(ids[g], collapse = "-"))
  L <- multcompLetters(p, compare = "<", threshold = 0.05, Letters = letters)$Letters[ids]; names(L) <- lv
  ch <- unique(unlist(strsplit(unname(L), ""))); mp <- setNames(letters[seq_along(ch)], ch); L[] <- vapply(unname(L), function(s) paste(mp[strsplit(s, "")[[1]]], collapse = ""), "")
  n <- table(d[[grp]])[lv]
  list(labels = setNames(sprintf("%s (n=%d) %s", lv, as.integer(n), L), lv), pairs = data.table(pair = sapply(cmb, paste, collapse = " vs "), KS_p = p), letters = L)
}
# ---- panels ------------------------------------------------------------------
th_title <- theme(plot.title = element_text(face = "bold", size = 16))
tot <- nrow(df); st <- df[, .(n = .N), by = pUP_Group][, pct := round(100 * n / tot, 1)]
lab <- function(g) sprintf("%s%%\n(n=%d)", st$pct[st$pUP_Group == g], st$n[st$pUP_Group == g])
pA <- ggplot(df, aes(x = pUP)) + geom_histogram(binwidth = 2, fill = "#3498db", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(30, 80), color = "black", linewidth = 0.5) +
  annotate("text", x = 15, y = Inf, label = lab("0~30"), vjust = 1.2, size = 3.5) + annotate("text", x = 55, y = Inf, label = lab("30~80"), vjust = 1.2, size = 3.5) +
  annotate("text", x = 91, y = Inf, label = lab("80~100"), vjust = 1.2, size = 3.5) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +   # V30: the three count labels in free room above the bars
  theme_classic() + labs(title = "A", x = "%uORF+ (TSS-seq)", y = "Number of genes") + th_title
kB <- ks_labels(as.data.frame(df), "uIF_Group", "pUP"); dB <- copy(df)[, lab := factor(kB$labels[as.character(uIF_Group)], kB$labels)]
pB <- ggplot(dB, aes(x = pUP, color = lab)) + stat_ecdf(geom = "step", linewidth = 1.2) + scale_color_manual(values = setNames(scales::hue_pal()(3), kB$labels)) +
  theme_classic() + labs(title = "B", x = "%uORF+ (TSS-seq)", y = "Fraction", color = "") + th_title + theme(legend.position = c(0.70, 0.4), legend.background = element_blank())
dC <- df[is.finite(uIF_raw) & uIF_raw > 0]   # V23: unrounded uIF for the log axis (four genes round to 0.000)
pC <- ggplot(dC, aes(x = pUP, y = log2(uIF_raw))) + geom_point(color = "#3498db", alpha = 0.5, size = 1) + geom_smooth(method = "lm", color = "red", linewidth = 1, se = FALSE, formula = y ~ x) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + theme_classic() + labs(title = "C", x = "%uORF+ (TSS-seq)", y = bquote(log[2]~"(uIF)")) + th_title
dD <- df[is.finite(TE) & TE > 0]
pD <- ggplot(dD, aes(x = pUP, y = log2(TE))) + geom_point(color = "#3498db", alpha = 0.5, size = 1) + geom_smooth(method = "lm", color = "red", linewidth = 1, se = FALSE, formula = y ~ x) +
  stat_cor(method = "pearson", label.x.npc = 0.35, label.y.npc = "top") + theme_classic() + labs(title = "D", x = "%uORF+ (TSS-seq)", y = bquote(log[2]~"(mORF TE)")) + th_title
kE <- ks_labels(as.data.frame(dfF), "pUP_Group", "TE"); dE <- copy(dfF)[, lab := factor(kE$labels[as.character(pUP_Group)], kE$labels)]
pE <- ggplot(dE[TE <= 3.5], aes(x = TE, color = lab)) + stat_ecdf(geom = "step", linewidth = 1.2) + scale_color_manual(values = setNames(c(scales::hue_pal()(3), "#999999"), kE$labels)) +
  coord_cartesian(xlim = c(0, 3.5)) + theme_classic() + labs(title = "E", x = "mORF TE", y = "Fraction", color = "") + th_title + theme(legend.position = c(0.70, 0.4), legend.background = element_blank())
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
# V30 (2026-11-09, Larry): two rows laid out like Figure 2 (A wide and B; C, D and E); the legends of B and E at the right of their panels; aligned plotting areas (fig_grid)
cS10 <- lapply(list(pA, pB, pC, pD, pE), compact); for (i in c(1, 3, 4)) cS10[[i]] <- cor_left(text_pt(cS10[[i]]))
cS10[[2]] <- legend_side(cS10[[2]], wrap_paren); cS10[[5]] <- legend_side(cS10[[5]], wrap_paren)
fig <- fig_grid(cS10, design = "AAB\nCDE", top = fig_title("Figure S8. uORF inclusion, uIF and mORF TE analyzed with TSS-seq data"), h = 5.4)
fig_save(file.path(OUT_DIR, "FigureS2_TSSseq_panels_A-E.pdf"), fig, h = 5.4)

# ---- numbers -----------------------------------------------------------------
rC <- cor.test(dC$pUP, log2(dC$uIF_raw)); rD <- cor.test(dD$pUP, log2(dD$TE))
.add(sprintf("  A: 100%% uORF+ (TSS-seq) in %d of %d genes (%.1f%%); below 90: %.1f%%; groups 0~30 / 30~80 / 80~100 = %s", sum(df$pUP >= 100), tot, 100 * mean(df$pUP >= 100), 100 * mean(df$pUP < 90), paste(sprintf("%.1f%% (n=%d)", st$pct[match(c("0~30", "30~80", "80~100"), st$pUP_Group)], st$n[match(c("0~30", "30~80", "80~100"), st$pUP_Group)]), collapse = " / ")))
.add("  B: medians of %uORF+ by uIF class: ", paste(sprintf("%s %.1f", UIF_LV, tapply(df$pUP, df$uIF_Group, median)), collapse = "; "), "; letters ", paste(kB$letters, collapse = " "))
for (i in seq_len(nrow(kB$pairs))) .add(sprintf("     %s: KS p = %.2g", kB$pairs$pair[i], kB$pairs$KS_p[i]))
.add(sprintf("  C: Pearson r(%%uORF+ TSS-seq, log2 uIF) = %.2f (p = %.1e), n = %d (genes with uIF > 0); CAGE in the same genes: %.2f", rC$estimate, rC$p.value, nrow(dC), cor(dC$pUP_cage, log2(dC$uIF_raw))))
.add(sprintf("  D: Pearson r(%%uORF+ TSS-seq, log2 TE) = %.2f (p = %.1e), n = %d; CAGE in the same genes: %.2f", rD$estimate, rD$p.value, nrow(dD), cor(dD$pUP_cage, log2(dD$TE))))
.add("  E: median TE by group: ", paste(sprintf("%s %.2f (n=%d)", PUP_LV, tapply(dfF$TE, dfF$pUP_Group, median), table(dfF$pUP_Group)), collapse = "; "), "; letters ", paste(kE$letters, collapse = " "))
for (i in seq_len(nrow(kE$pairs))) .add(sprintf("     %s: KS p = %.2g", kE$pairs$pair[i], kE$pairs$KS_p[i]))
.add(sprintf("  CAGE vs TSS-seq %%uORF+ per gene: Spearman %.2f; same group in %.0f%% of genes", cor(df$pUP_cage, df$pUP, method = "spearman"),
             100 * mean(fcase(df$pUP_cage >= 80, "80~100", df$pUP_cage >= 30, "30~80", default = "0~30") == as.character(df$pUP_Group))))
fwrite(df[, .(gene_id, uIF, uIF_class = uIF_Group, mORF_TE = TE, pct_uORF_plus_CAGE = pUP_cage, pct_uORF_plus_TSSseq = pUP, pUP_group_TSSseq = pUP_Group)], file.path(OUT_DIR, "FigS2_TSSseq_values.csv"))
writeLines(.log, file.path(OUT_DIR, "FigS2_TSSseq_summary.txt")); .add("done: output/FigureS2_TSSseq_panels_A-E.pdf")
