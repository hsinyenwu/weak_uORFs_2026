################################################################################
# A16 (2026-10-29, Larry): do the predicted-but-undetected uORFs that remain in the mRNA after TSS weighting
# (the ~1.6% of nuclear 5'-leader mRNA behind the 10.9% vs 10.0% of Figure 3C) go with translational repression?
# Gene level, nuclear genes with an annotated 5' leader (the Figure 3C transcript set). Results only.
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
tx <- fread("output/Table_Transcriptome_Impact_PerTranscript_CAGE.csv")[Category == "Nuclear_with_5UTR"]
rna <- fread("Files/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only", select = c("gene_id", "TPM")); setnames(rna, "TPM", "RNA_TPM")
ribo <- fread("Files/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only", select = c("gene_id", "TPM")); setnames(ribo, "TPM", "Ribo_TPM")
te <- merge(rna, ribo, by = "gene_id")[RNA_TPM > 0.1 & Ribo_TPM > 0, .(gene_id, TE = Ribo_TPM / RNA_TPM)]
tot <- sum(tx$TPM)
# mass components (per transcript, as in Figure 3C)
m <- tx[, .(pred_mass = sum(TPM * Max_pUP_pred / 100), tu_mass = sum(TPM * Max_pUP_TuORF / 100)), by = .(TuORF_mapped, has_pred = n_pred_uORFs > 0)]
cat(sprintf("TSS-weighted mass (%% of nuclear 5'-leader mRNA): predicted uORFs %.2f%% (on TuORF transcripts %.2f%%, on transcripts without a mapped TuORF %.2f%%); detected TuORFs %.2f%%\n",
    100 * sum(m$pred_mass) / tot, 100 * sum(m[TuORF_mapped == TRUE, pred_mass]) / tot, 100 * sum(m[TuORF_mapped == FALSE, pred_mass]) / tot, 100 * sum(m$tu_mass) / tot))
cat(sprintf("  TuORF transcripts without any predicted AUG uORF (no AUG uORF was predicted in the annotated leader of these transcripts): %d transcripts, TuORF-weighted mass %.2f%%\n",
    tx[TuORF_mapped == TRUE & n_pred_uORFs == 0, .N], 100 * tx[TuORF_mapped == TRUE & n_pred_uORFs == 0, sum(TPM * Max_pUP_TuORF / 100)] / tot))
# gene level
g <- tx[, .(TPM = sum(TPM), TuORF = any(TuORF_mapped), ouORF = any(ouORF_mapped), pred = any(n_pred_uORFs > 0),
            pUP_pred = max(Max_pUP_pred), pUP_tu = max(Max_pUP_TuORF), pred_mass = sum(TPM * Max_pUP_pred / 100), tu_mass = sum(TPM * Max_pUP_TuORF / 100)), by = gene_id]
g[, group := fifelse(TuORF, "TuORF gene", fifelse(ouORF, "ouORF gene (no TuORF)", fifelse(pred, "predicted-only uORF", "no uORF")))]
bin <- function(p) cut(p, c(-Inf, 30, 80, Inf), labels = c("<30", "30-80", "80-100"))
g[, pUP_bin := fifelse(group == "TuORF gene", as.character(bin(pUP_tu)), fifelse(group == "predicted-only uORF", as.character(bin(pUP_pred)), NA_character_))]
g <- merge(g, te, by = "gene_id", all.x = TRUE)
ref <- g[group == "no uORF" & !is.na(TE), TE]; med_ref <- median(ref)
cat(sprintf("\nGenes: no uORF %d, predicted-only uORF %d, TuORF %d, ouORF-only %d (nuclear genes with a 5' leader and expressed isoforms; TE for genes with RNA TPM > 0.1)\n",
    g[group == "no uORF", .N], g[group == "predicted-only uORF", .N], g[group == "TuORF gene", .N], g[group == "ouORF gene (no TuORF)", .N]))
cat(sprintf("Reference: no-uORF genes with TE, n = %d, median TE %.3f\n\n", length(ref), med_ref))
summ <- function(d, label) { v <- d[!is.na(TE), TE]
  cat(sprintf("%-42s n = %5d (with TE %5d) | median TE %.3f = %+.0f%% vs no-uORF | TE < no-uORF median: %.0f%% | KS p = %s | Wilcoxon p = %s | TSS-weighted mass %.2f%% of mRNA\n",
      label, nrow(d), length(v), median(v), 100 * (median(v) / med_ref - 1), 100 * mean(v < med_ref), format.pval(suppressWarnings(ks.test(v, ref))$p.value, 2), format.pval(wilcox.test(v, ref)$p.value, 2),
      100 * sum(if (grepl("TuORF", label)) d$tu_mass else d$pred_mass) / tot)) }
summ(g[group == "predicted-only uORF"], "predicted-only uORF genes, all")
for (b in c("80-100", "30-80", "<30")) summ(g[group == "predicted-only uORF" & pUP_bin == b], sprintf("  predicted-only, %%uORF+ %s", b))
summ(g[group == "predicted-only uORF" & pUP_pred >= 90], "  predicted-only, %uORF+ >= 90")
summ(g[group == "TuORF gene"], "TuORF genes, all (reference pattern)")
for (b in c("80-100", "30-80", "<30")) summ(g[group == "TuORF gene" & pUP_bin == b], sprintf("  TuORF genes, %%uORF+ %s", b))
# expression check: are the included predicted-only uORFs in lowly expressed genes (hence undetected)?
q <- g[group == "predicted-only uORF" & pUP_pred >= 80, .(median_RNA_TPM = median(TPM), share_TPM_below_1 = mean(TPM < 1))]
q2 <- g[group == "TuORF gene", .(median_RNA_TPM = median(TPM), share_TPM_below_1 = mean(TPM < 1))]
cat(sprintf("\nExpression: predicted-only genes with %%uORF+ 80-100: median gene TPM %.1f, %.0f%% below 1 TPM; TuORF genes: median %.1f, %.0f%% below 1 TPM\n", q$median_RNA_TPM, 100 * q$share_TPM_below_1, q2$median_RNA_TPM, 100 * q2$share_TPM_below_1))
# TE among the well-expressed included predicted-only genes (where Ribo-seq could have detected translation)
summ(g[group == "predicted-only uORF" & pUP_pred >= 80 & TPM >= 5], "  predicted-only, %uORF+ 80-100, gene TPM >= 5")
summ(g[group == "no uORF" & TPM >= 5], "  no uORF, gene TPM >= 5 (comparison)")
fwrite(g[order(group, -TPM)], "output/A16_Table_predicted_only_uORF_TE.csv")
d <- g[!is.na(TE) & (group == "no uORF" | (group == "predicted-only uORF" & !is.na(pUP_bin)) | (group == "TuORF gene" & !is.na(pUP_bin)))]
d[, set := fifelse(group == "no uORF", "no uORF", paste0(fifelse(group == "TuORF gene", "TuORF ", "predicted-only "), pUP_bin, "% uORF+"))]
p <- ggplot(d, aes(log10(TE), colour = set)) + stat_ecdf(linewidth = 0.5) + theme_bw(base_size = 9) + labs(x = "log10 mORF TE", y = "Cumulative fraction of genes", title = "A16: mORF TE by uORF status (nuclear genes with a 5' leader)", colour = NULL)
ggsave("output/A16_Figure_predicted_only_uORF_TE.pdf", p, width = 6.5, height = 4)
