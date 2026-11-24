################################################################################
# A18 (2026-10-31, Larry; revised 2026-11-01) — mRNA-weighted loss of translated uORFs (Supplemental Fig. S14 of the revision)
# Of all leader-containing mRNA made from genes with a translated uORF (Figure 3C transcript set, same-gene mapping, RSEM TPM),
#   how much contains the TuORF, lacks it because transcription starts downstream of it (TSS), or lacks it because the start
#   codon is not in the mature mRNA although transcription starts upstream (leader splicing)?
# Transcripts that carry the TuORF in their exons: contains = TPM x p x s, lacks by TSS = TPM x (1 - p), lacks by splicing =
#   TPM x p x (1 - s), with p = %uORF+/100 of the transcript's best TuORF and s = max_k(p_k s_k) / max_k p_k over the TuORFs k that the
#   transcript carries (V25, 2026-11-09; over the gene's TuORFs in the audit table in V23/V24, which let TuORFs on no expressed isoform set the
#   factor; the best effective inclusion; s_k = fraction of RNA-seq reads across the start codon that keep it, s_k = 1 when
#   fewer than 10 reads; for class-D TuORFs with an expressed alternative isoform s_k = min(1, S_k / carrier share), V23; revision_modules/splicing_rule.R).
# Since 2026-11-01 the genes' other leader isoforms, which lack the TuORF by annotation, are included as well: an isoform whose
#   annotated start lies downstream of the TuORF start codon counts as lacking it by TSS choice; an isoform whose exons do not
#   contain the start codon (start codon in one of its introns, or the ORF not fully contained) counts as lacking it by splicing.
# Shown for all genes and after removing the 1% and 5% most highly expressed genes, with CAGE-seq and TSS-seq; panel C shows
#   the per-gene distribution (each gene counted once).
# Outputs: output/A18_Figure_mRNA_weighted_TSS_splicing.pdf, output/A18_Table_mRNA_weighted_TSS_splicing.csv,
#   output/A18_Table_per_gene_CAGE.csv, output/A18_Table_isoforms_without_TuORF.csv, output/A18_summary.txt
# Run from Weak_uORF_submission/:  Rscript revision_modules/A18_mRNA_weighted_TSS_splicing.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(gridExtra); library(grid); library(readxl); library(GenomicRanges); library(rtracklayer) })
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
PATH_GTF2 <- "Files/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf"
source("revision_modules/splicing_rule.R")   # V23 (2026-11-08): junction-consistent s_k; one leader definition per transcript (as A26)
ORFS <- orf_s_table(include_ouORF = FALSE)
CARR <- orf_carriers(ORFS)   # V25 (2026-11-09): per-transcript factor over the TuORFs the transcript carries
S_src <- list(CAGE = tx_s(CARR, ORFS, "p_CAGE")[, .(transcript_id, s_gene = s_tx, s_measured)], TSSseq = tx_s(CARR, ORFS, "p_TSSseq")[, .(transcript_id, s_gene = s_tx, s_measured)])   # V23: the TSS-seq panel weights s by TSS-seq inclusion; V25: per transcript (the column keeps its name s_gene)
LL <- leader_lengths(); leader_tx <- LL[leader_len > 0, transcript_id]

## V23 (2026-11-08): every leader isoform of a TuORF gene is classified relative to the gene's best TuORF, the most included TuORF that at
## least one expressed leader isoform carries (start and stop codon in its exons):
##   carriers of the best TuORF: contains = TPM x p x s, lacks by TSS = TPM x (1 - p), lacks by splicing = TPM x p x (1 - s)
##     (p = highest %uORF+/100 among the transcript's TuORFs, as in Figure 3C; s = s_t over the TuORFs the transcript carries, V25, junction-consistent;
##     the best TuORF is chosen by CAGE %uORF+ with ties broken by the order of Table S1A, and the TSS-seq panel keeps that choice, V25 review 2026-11-09);
##   isoforms whose annotated 5' end lies downstream of the best TuORF's start codon: lacks by TSS = TPM;
##   isoforms that start upstream but lack the best TuORF (start codon in an intron, or the ORF interrupted by an intron): lacks by splicing
##     = TPM x p_b, lacks by TSS = TPM x (1 - p_b), p_b = %uORF+/100 of the best TuORF.
## Before V23 these isoforms counted 100% as splicing, and isoforms that carried another, poorly included TuORF counted as carriers
## (which put 1.4 points of splicing loss under TSS).
tx0 <- fread("output/Table_Transcriptome_Impact_PerTranscript_CAGE.csv")[grepl("^Nuclear", Category) & transcript_id %in% leader_tx]
s1a <- as.data.table(read_excel("output/Table_S1_all_TuORFs.xlsx", sheet = "Table S1A (CAGE)")); s1b <- as.data.table(read_excel("output/Table_S1_all_TuORFs.xlsx", sheet = "Table S1B (TSS-seq)"))
s1a[, pT := s1b$pct_uORF_plus[match(paste(gene_id, uORF_genome_range), paste(s1b$gene_id, s1b$uORF_genome_range))]]
gtf <- import(PATH_GTF2); ex <- gtf[gtf$type == "exon"]; ex_by_tx <- split(ex, ex$transcript_id)
in_exons <- function(t, pos) { e <- ex_by_tx[[t]]; !is.null(e) && any(start(e) <= pos & end(e) >= pos) }
iso <- tx0[gene_id %in% s1a$gene_id, .(transcript_id, gene_id, TPM)]
car <- merge(iso[, .(transcript_id, gene_id)], s1a[, .(gene_id, uORF_genome_range, uORF_start_genomic, uORF_stop_genomic, pC = pct_uORF_plus, pT)], by = "gene_id", allow.cartesian = TRUE)
car[, carries := mapply(function(t, a, b) in_exons(t, a) && in_exons(t, b), transcript_id, uORF_start_genomic, uORF_stop_genomic)]
best <- car[carries == TRUE][order(-pC)][, .SD[1], by = gene_id][, .(gene_id, best_uORF = uORF_genome_range, best_start = uORF_start_genomic, pC_b = pC / 100, pT_b = fifelse(is.na(pT), 0, pT / 100))]
iso <- merge(iso, best, by = "gene_id")
best_carriers <- car[carries == TRUE & paste(gene_id, uORF_genome_range) %in% paste(best$gene_id, best$best_uORF), unique(transcript_id)]
iso[, carrier_best := transcript_id %in% best_carriers]
cls <- rep("carrier", nrow(iso))
for (i in which(!iso$carrier_best)) { t <- iso$transcript_id[i]; e <- ex_by_tx[[t]]; if (is.null(e)) { cls[i] <- "no model"; next }
  sp <- range(e); st <- as.character(strand(sp))[1]; pos <- iso$best_start[i]; tss <- if (st == "+") start(sp) else end(sp)
  cls[i] <- if (if (st == "+") pos < tss else pos > tss) "TSS (annotated 5' end downstream of the start codon)" else if (!in_exons(t, pos)) "splicing (start codon in an intron)" else "splicing (ORF interrupted by an intron)" }
iso[, reason := cls]
fwrite(iso[carrier_best == FALSE, .(transcript_id, gene_id, TPM, best_uORF, reason)], "output/A18_Table_isoforms_without_TuORF.csv")
genes <- unique(iso$gene_id)
.add("V23: %d TuORF genes with an expressed leader isoform that carries a TuORF; %d isoforms carry the gene's best TuORF, %d do not (%s)", length(genes), sum(iso$carrier_best), sum(!iso$carrier_best),
     paste(iso[carrier_best == FALSE, .(n = .N, TPM = round(sum(TPM))), by = reason][order(-TPM)][, sprintf("%s: %d isoforms, %.0f TPM", reason, n, TPM)], collapse = "; "))

run <- function(src) {
  tx <- fread(sprintf("output/Table_Transcriptome_Impact_PerTranscript_%s.csv", src))[, .(transcript_id, Max_pUP_TuORF)]
  d <- merge(iso, tx, by = "transcript_id", all.x = TRUE); d <- merge(d, S_src[[src]], by = "transcript_id", all.x = TRUE); d[is.na(s_gene), `:=`(s_gene = 1, s_measured = FALSE)]   # V25: s per transcript
  d[, pb := if (src == "CAGE") pC_b else pT_b]
  d[, p := fifelse(carrier_best, fifelse(is.na(Max_pUP_TuORF), 0, Max_pUP_TuORF / 100), NA_real_)]
  d[, `:=`(contains = fifelse(carrier_best, TPM * p * s_gene, 0),
           lack_tss = fifelse(carrier_best, TPM * (1 - p), fifelse(grepl("^splicing", reason), TPM * (1 - pb), TPM)),
           lack_spl = fifelse(carrier_best, TPM * p * (1 - s_gene), fifelse(grepl("^splicing", reason), TPM * pb, 0)))]
  g <- d[, .(TPM = sum(TPM), TPM_m = sum(TPM[carrier_best]), TPM_u = sum(TPM[!carrier_best]), contains = sum(contains), lack_tss = sum(lack_tss), lack_spl = sum(lack_spl),
             lack_tss_u = sum(lack_tss[!carrier_best]), lack_spl_u = sum(lack_spl[!carrier_best]), TPM_mm = sum(TPM[carrier_best & s_measured])), by = gene_id]   # V25: TPM_mm = carrier mRNA with a measured s (per transcript)
  stopifnot(all(abs(g$contains + g$lack_tss + g$lack_spl - g$TPM) < 1e-6 * pmax(1, g$TPM)))
  g <- g[order(-TPM)]                                   # most expressed genes first (the top 1% / 5% removal and the top-10 share rely on this order)
  n <- nrow(g); out <- list()
  for (drop in c(0, 0.01, 0.05)) {
    k <- ceiling(n * drop); dd <- if (k) g[-(1:k)] else g; T <- sum(dd$TPM)
    out[[length(out) + 1]] <- data.table(source = src, subset = c("All genes", "Without the top 1% of genes", "Without the top 5% of genes")[match(drop, c(0, 0.01, 0.05))],
      genes_removed = k, genes = nrow(dd), contains = 100 * sum(dd$contains) / T, lacks_TSS = 100 * sum(dd$lack_tss) / T, lacks_splicing = 100 * sum(dd$lack_spl) / T,
      lacks_TSS_annotated_isoforms = 100 * sum(dd$lack_tss_u) / T, lacks_splicing_annotated_isoforms = 100 * sum(dd$lack_spl_u) / T,
      isoforms_without_best_TuORF_TPM = 100 * sum(dd$TPM_u) / T, splicing_measured_TPM = 100 * sum(dd$TPM_mm) / sum(dd$TPM_m))
  }
  res <- rbindlist(out); res[, lacks_total := lacks_TSS + lacks_splicing]
  g[, `:=`(pct_lacking_TSS = 100 * lack_tss / TPM, pct_lacking_TSS_or_splicing = 100 * (lack_tss + lack_spl) / TPM)]
  list(res = res, genes = g, top10 = 100 * sum(g$TPM[1:10]) / sum(g$TPM))
}
R <- lapply(c(CAGE = "CAGE", TSSseq = "TSSseq"), run)
tab <- rbindlist(lapply(R, `[[`, "res"))
for (i in seq_len(nrow(tab))) .add("%-7s %-30s genes %5d (removed %3d): contains %.1f%% | lacks %.1f%% = TSS %.1f%% (of which isoforms without the best TuORF %.1f%%) + splicing %.1f%% (of which isoforms without the best TuORF %.1f%%); splicing measured for %.0f%% of the mRNA of the isoforms that carry the best TuORF; isoforms without it: %.1f%% of the mRNA",
  tab$source[i], tab$subset[i], tab$genes[i], tab$genes_removed[i], tab$contains[i], tab$lacks_total[i], tab$lacks_TSS[i], tab$lacks_TSS_annotated_isoforms[i], tab$lacks_splicing[i], tab$lacks_splicing_annotated_isoforms[i], tab$splicing_measured_TPM[i], tab$isoforms_without_best_TuORF_TPM[i])
gC <- R$CAGE$genes
.add("CAGE: the 10 most expressed genes make %.1f%% of this mRNA; per gene (each gene once, n = %d): median %% of transcripts lacking the TuORF %.1f (TSS only %.1f); mean %.1f (TSS only %.1f; splicing %.1f); genes with >= 80%% lacking %.1f%%; genes with >= 50%% lacking %.1f%%; genes with splicing loss >= 10 points %.1f%%",
  R$CAGE$top10, nrow(gC), median(gC$pct_lacking_TSS_or_splicing), median(gC$pct_lacking_TSS), mean(gC$pct_lacking_TSS_or_splicing), mean(gC$pct_lacking_TSS), mean(gC$pct_lacking_TSS_or_splicing - gC$pct_lacking_TSS),
  100 * mean(gC$pct_lacking_TSS_or_splicing >= 80), 100 * mean(gC$pct_lacking_TSS_or_splicing >= 50), 100 * mean(gC$pct_lacking_TSS_or_splicing - gC$pct_lacking_TSS >= 10))
fwrite(tab[, .(source, subset, genes_removed, genes, contains = round(contains, 2), lacks_total = round(lacks_total, 2), lacks_TSS = round(lacks_TSS, 2), lacks_splicing = round(lacks_splicing, 2),
               lacks_TSS_annotated_isoforms = round(lacks_TSS_annotated_isoforms, 2), lacks_splicing_annotated_isoforms = round(lacks_splicing_annotated_isoforms, 2), splicing_measured_TPM = round(splicing_measured_TPM, 1))], "output/A18_Table_mRNA_weighted_TSS_splicing.csv")
fwrite(gC[, .(gene_id, TPM = round(TPM, 2), TPM_isoforms_without_TuORF = round(TPM_u, 2), pct_lacking_TSS = round(pct_lacking_TSS, 2), pct_lacking_TSS_or_splicing = round(pct_lacking_TSS_or_splicing, 2))], "output/A18_Table_per_gene_CAGE.csv")
# ---- figure
source("revision_modules/fig_compact.R")
# 2026-11-19 (Larry: "Can you make S14 A and B the same style as D? move the legend on the right"): A and B as grouped bars (the three fractions on
# the x axis, one bar per gene subset), legends at the right of every panel
lv <- c("All genes", "Without the top 1% of genes", "Without the top 5% of genes")
parts <- c(contains = "contain\nthe uORF", lacks_TSS = "lack it: transcription\nstarts downstream (TSS)", lacks_splicing = "lack it: removed\nby leader splicing")
panel <- function(src, ttl) { d <- tab[source == src]
  long <- melt(d, id.vars = c("subset", "genes"), measure.vars = names(parts), variable.name = "part", value.name = "pct")
  long[, part := factor(parts[as.character(part)], levels = parts)]
  long[, grp := factor(sprintf("%s\n(n = %s)", subset, format(genes, big.mark = ",")), levels = sprintf("%s\n(n = %s)", lv, format(d$genes[match(lv, d$subset)], big.mark = ",")))]
  ggplot(long, aes(part, pct, fill = grp)) + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(0.8), vjust = -0.4, size = 8 / .pt, family = FIG_FONT) +
    scale_fill_manual(values = setNames(c("#08519C", "#4292C6", "#9ECAE1"), levels(long$grp)), name = NULL) + coord_cartesian(ylim = c(0, 72)) +
    labs(x = NULL, y = "% of mRNA", title = ttl) + theme_classic(base_size = 9) }
pA <- panel("CAGE", "A  mRNA-weighted fractions, CAGE-seq")
pB <- panel("TSSseq", "B  The same with TSS-seq")
dC <- rbind(gC[, .(v = pct_lacking_TSS, what = "TSS only")], gC[, .(v = pct_lacking_TSS_or_splicing, what = "TSS or splicing")])
pC <- ggplot(dC, aes(v, colour = what)) + stat_ecdf(linewidth = 0.6) + geom_vline(xintercept = median(gC$pct_lacking_TSS_or_splicing), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c("TSS only" = "grey45", "TSS or splicing" = "#E66101"), name = NULL) +
  labs(x = "% of the gene's transcripts that lack the uORF", y = "Cumulative fraction", title = "C  Per gene, each gene counted once (CAGE-seq)",
       subtitle = sprintf("n = %s genes; dotted line = median (%.1f%%)", format(nrow(gC), big.mark = ","), median(gC$pct_lacking_TSS_or_splicing))) +   # mean and the >= 80% share are in the summary (half-width panel, 2026-11-19)
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
# D (2026-11-19, Larry: "make a figure like Figure 3C for the sequence predicted uORF and ouORFs for the effects from TSS and splicing", then
#   "Panel of Fig. S13", "wait no. make it Figure S14", new panel of S14): the sequence-predicted sets of module A34 (CAGE-seq), in the style of Figure 3C
B34 <- fread("output/A34_Table_ladder.csv")[source == "CAGE"]
SNB <- c("Predicted uORFs ending in the leader", "Predicted overlapping uORFs", "All predicted AUG uORFs")
stepsB <- c("genes\nwith the uORF", "nuclear mRNA\n(unweighted)", "nuclear mRNA\nafter TSS\nweighting", "nuclear mRNA\nafter TSS weighting\nand splicing")
DB <- rbindlist(lapply(SNB, function(s_) { r <- B34[set == s_]; data.table(set = s_, step = stepsB, value = c(r$pct_genes, r$pct_mRNA, r$pct_TSS, r$pct_TSS_splicing)) }))
DB[, `:=`(set = factor(set, levels = SNB), step = factor(step, levels = stepsB))]
pD <- ggplot(DB, aes(step, value, fill = set)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", value)), position = position_dodge(0.8), vjust = -0.4, size = 8 / .pt, family = FIG_FONT) +
  scale_fill_manual(values = c(`Predicted uORFs ending in the leader` = "#BDBDBD", `Predicted overlapping uORFs` = "#9E9AC8", `All predicted AUG uORFs` = "#636363"), labels = c("Predicted uORFs\nending in the leader", "Predicted\noverlapping uORFs", "All predicted\nAUG uORFs"), name = NULL) +
  coord_cartesian(ylim = c(0, 72)) + labs(x = NULL, y = "%", title = "D  Sequence-predicted AUG uORFs: TSS usage and leader splicing") +
  theme_classic(base_size = 9)
rt <- function(p) compact_wide(p) + theme(legend.position = "right", legend.justification = c(0, 0.5), legend.key.spacing.y = unit(4, "pt"))   # 2026-11-19: legends at the right
rt_half <- function(p) compact(p) + theme(legend.position = "right", legend.justification = c(0, 0.5), legend.key.spacing.y = unit(4, "pt"))   # C at half width (Larry: "make C half size on x-axis")
pE <- ggplot() + theme_void()   # empty right half of the row of C (fig_grid lays out a row with one panel at full width)
g <- fig_grid(list(rt(pA), rt(pB), rt_half(pC), rt(pD), pE), design = "AA\nBB\nCE\nDD", heights = c(1, 1, 1, 1), h = 9,   # V30 (2026-11-09, Larry): the legend of C inside the plot (fig_grid); D added 2026-11-19
  top = fig_title("Figure S14. Translated and sequence-predicted uORFs: loss by TSS position and by leader splicing"))
fig_save("output/A18_Figure_mRNA_weighted_TSS_splicing.pdf", g, h = 9)
writeLines(.log, "output/A18_summary.txt")
