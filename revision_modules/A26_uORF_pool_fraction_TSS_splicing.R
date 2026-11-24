################################################################################
# A26 (2026-11-06, Larry; rewritten for V23, 2026-11-08) — The prevalence ladder of Figure 3C and Supplemental Fig. S13: share of the
# nuclear leader-containing mRNA that carries a uORF, from one per-transcript table per TSS source, with one denominator and one
# inclusion rule for every set, and the share after TSS weighting and leader splicing for the detected sets.
#
# Pool (denominator): nuclear transcripts of the expressed isoform set (RSEM TPM, mean of the six libraries; at most two isoforms per gene,
#   isoform TPM >= 0.5) that have an annotated 5' leader (exonic bases upstream of the CDS start; splicing_rule.R::leader_lengths()).
# Carrier transcripts: transcripts whose exons contain the start and stop codon of an ORF of the set (same-gene mapping of the canonical
#   script: columns TuORF_mapped, ouORF_mapped, n_pred_uORFs of output/Table_Transcriptome_Impact_PerTranscript_{CAGE,TSSseq}.csv).
# For each set (detected TuORFs; detected TuORFs or overlapping uORFs; any predicted AUG uORF):
#   % of genes = genes with >= 1 carrier / genes with >= 1 leader-containing transcript;
#   % of mRNA (unweighted) = TPM of carriers / pool TPM (the same for both TSS sources);
#   % after TSS weighting = sum(TPM x p) / pool, p = highest %uORF+/100 among the transcript's ORFs of the set;
#   % after TSS weighting and leader splicing (detected sets) = sum(TPM x p x s_t) / pool, s_t = max_k(p_k s_k) / max_k(p_k) over the ORFs k
#   of the set that the transcript carries (splicing_rule.R: orf_carriers(), tx_s(); s_k from the junction-consistent rule).
# V23 changes: (1) leader defined per transcript for numerator and denominator (the canonical metric table summed carriers over all nuclear
#   transcripts but took the denominator from a gene list, so 20 TuORF carriers lay outside it); (2) the same inclusion rule (carrier, any p)
#   for every set (the canonical TuORF rows required p > 0 and the others did not, so the unweighted TuORF share depended on the TSS source);
#   (3) the junction-consistent s_k (the V21 and V22 losses are kept in the pool-fraction table for the record).
# V25 change (2026-11-09): the splicing factor is assigned per transcript over the ORFs it carries, not per gene over every ORF of the gene in
#   the audit tables (673 TuORFs lie on no expressed isoform and set the gene factor of 49 genes; see splicing_rule.R). 8.538% -> 8.546%.
# Outputs: output/A26_Table_prevalence_ladder.csv (all ladder values and the TPM table of Figure 3C / S13), output/A26_Table_uORF_pool_fractions.csv,
#   output/A26_Table_orf_s.csv, output/A26_Table_s_transcript.csv (V25; A26_Table_s_gene.csv until V24), output/A26_summary.txt
# Run from Weak_uORF_submission/ after A2 and A25:  Rscript revision_modules/A26_uORF_pool_fraction_TSS_splicing.R
################################################################################
suppressPackageStartupMessages({ library(data.table) })
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
source("revision_modules/splicing_rule.R")
ORFS <- orf_s_table(include_ouORF = TRUE); A2 <- ORFS[kind == "TuORF"]; A25 <- ORFS[kind == "ouORF"]
.add(sprintf("splicing rule (V23): %d ORFs; class D with an expressed alternative isoform and S measured: %d (carrier share computed for %d; median carrier share %.2f, median S %.3f); s_k < 1 for %d ORFs",
             nrow(ORFS), sum(ORFS$cls == "D" & ORFS$alt_ex & is.finite(ORFS$S)), sum(!is.na(ORFS$c_share)), median(ORFS$c_share, na.rm = TRUE), median(ORFS$S[!is.na(ORFS$c_share)], na.rm = TRUE), sum(ORFS$s < 1)))
fwrite(ORFS[, .(gene_id, kind, uid, tis, p_CAGE, p_TSSseq, S, cov, cls, alt_ex, c_share, s_V21, s_V22, s)], "output/A26_Table_orf_s.csv")
LL <- leader_lengths(); leader_tx <- LL[leader_len > 0, transcript_id]
CARR <- orf_carriers(ORFS)   # V25: the ORFs each expressed isoform carries (start and stop codon in its exons, same gene)
.add("V25: ORFs carried by at least one expressed isoform: %d of %d TuORFs and %d of %d overlapping uORFs; class-D ORFs with an expressed alternative isoform and S measured: carrier share > 0 (division applied) for %d, carrier share 0 (no expressed isoform has the start codon in an exon; never carried) for %d",
     uniqueN(CARR[uid %in% A2$uid, .(uid, gene_id)]), nrow(A2), uniqueN(CARR[uid %in% A25$uid, .(uid, gene_id)]), nrow(A25),
     sum(ORFS$cls == "D" & ORFS$alt_ex & is.finite(ORFS$S) & !is.na(ORFS$c_share) & ORFS$c_share > 0), sum(ORFS$cls == "D" & ORFS$alt_ex & is.finite(ORFS$S) & !is.na(ORFS$c_share) & ORFS$c_share == 0))

out <- list(); lad <- list(); tpmtab <- list(); extra <- list(); st_out <- NULL
for (src in c("CAGE", "TSSseq")) {
  all <- fread(sprintf("output/Table_Transcriptome_Impact_PerTranscript_%s.csv", src))
  all[, leader := transcript_id %in% leader_tx]
  nuc <- all[grepl("^Nuclear", Category)]; tx <- nuc[leader == TRUE]; pool <- sum(tx$TPM); G <- uniqueN(tx$gene_id)
  stopifnot(nuc[(TuORF_mapped | ouORF_mapped | n_pred_uORFs > 0) & !leader, .N] == 0)   # every carrier has a leader
  tx[, `:=`(pT = fifelse(is.na(Max_pUP_TuORF), 0, Max_pUP_TuORF) / 100, pO = fifelse(is.na(Max_pUP_ouORF), 0, Max_pUP_ouORF) / 100, pP = fifelse(is.na(Max_pUP_pred), 0, Max_pUP_pred) / 100)]
  pcol <- paste0("p_", src)
  # V25 check: the carrier mapping of splicing_rule.R against the per-transcript table (TuORF_mapped, ouORF_mapped) and p against max p_k over the carried ORFs
  carT <- tx_s(CARR, A2, pcol); carO <- tx_s(CARR, A25, pcol)
  .add("%-7s carrier check: TuORF carriers %d (per-transcript table) vs %d (orf_carriers; %d transcripts differ); overlapping-uORF carriers %d vs %d (%d differ); |p - max p_k| > 0.005 for %d of %d TuORF carriers in both",
       src, sum(tx$TuORF_mapped), sum(tx$transcript_id %in% carT$transcript_id), sum(xor(tx$TuORF_mapped, tx$transcript_id %in% carT$transcript_id)),
       sum(tx$ouORF_mapped), sum(tx$transcript_id %in% carO$transcript_id), sum(xor(tx$ouORF_mapped, tx$transcript_id %in% carO$transcript_id)),
       merge(tx[TuORF_mapped == TRUE, .(transcript_id, pT)], carT, by = "transcript_id")[abs(pT - p_max) > 0.005, .N], merge(tx[TuORF_mapped == TRUE, .(transcript_id, pT)], carT, by = "transcript_id")[, .N])
  sets <- list("Detected TuORFs" = list(carrier = tx$TuORF_mapped, p = tx$pT, orfs = A2),
               "Detected TuORFs + overlapping uORFs" = list(carrier = tx$TuORF_mapped | tx$ouORF_mapped, p = pmax(tx$pT, tx$pO), orfs = rbind(A2, A25)),
               "Any predicted AUG uORF" = list(carrier = tx$n_pred_uORFs > 0, p = tx$pP, orfs = NULL),
               "Overlapping uORFs only (transcripts without a TuORF)" = list(carrier = tx$ouORF_mapped & !tx$TuORF_mapped, p = tx$pO, orfs = A25))
  for (nm in names(sets)) {
    st <- sets[[nm]]; d <- tx[st$carrier]; d[, p := st$p[st$carrier]]; tss <- sum(d$TPM * d$p)
    row <- data.table(source = src, set = nm, pool_TPM = pool, pool_transcripts = nrow(tx), leader_genes = G, carrier_transcripts = nrow(d), carrier_genes = uniqueN(d$gene_id),
                      pct_genes = 100 * uniqueN(d$gene_id) / G, carrier_TPM = sum(d$TPM), pct_unweighted = 100 * sum(d$TPM) / pool, TSS_weighted_TPM = tss, pct_TSS = 100 * tss / pool)
    if (!is.null(st$orfs)) {
      st_tx <- tx_s(CARR, st$orfs, pcol); d <- merge(d, st_tx, by = "transcript_id", all.x = TRUE); d[is.na(s_tx), `:=`(s_tx = 1, s_tx_V21 = 1, s_tx_V22 = 1, s_measured = FALSE, s_quantified = FALSE)]   # V25: per transcript
      both <- sum(d$TPM * d$p * d$s_tx)
      row[, `:=`(TSS_splicing_TPM = both, pct_TSS_splicing = 100 * both / pool, splicing_loss_points = 100 * (tss - both) / pool, splicing_loss_pct_of_TSS_weighted = 100 * (tss - both) / tss,
                 splicing_loss_points_rule_V21 = 100 * (tss - sum(d$TPM * d$p * d$s_tx_V21)) / pool, splicing_loss_points_rule_V22 = 100 * (tss - sum(d$TPM * d$p * d$s_tx_V22)) / pool,
                 pct_TSS_weighted_mRNA_with_s_measured = 100 * sum(d$TPM[d$s_measured] * d$p[d$s_measured]) / tss)]
      if (src == "CAGE" && nm == "Detected TuORFs + overlapping uORFs") st_out <- merge(d[, .(transcript_id, gene_id, TPM, p)], st_tx, by = "transcript_id", all.x = TRUE) }
    out[[length(out) + 1]] <- row }
  # pieces of the ladder quoted in the Results
  extra[[src]] <- data.table(source = src,
    pct_predicted_only = 100 * tx[n_pred_uORFs > 0 & !TuORF_mapped & !ouORF_mapped, sum(TPM * pP)] / pool,
    pct_ouORF_without_TuORF = 100 * tx[ouORF_mapped & !TuORF_mapped, sum(TPM * pO)] / pool,
    pct_ouORF_raising_p_on_TuORF_transcripts = 100 * tx[TuORF_mapped == TRUE, sum(TPM * (pmax(pT, pO) - pT))] / pool)
  # TPM table of Figure 3C / S13 (CAGE TPM are the same in both tables)
  tpmtab[[src]] <- data.table(source = src, row = c("Plastid encoded", "Mitochondria encoded", "Nuclear encoded, no 5' leader", "Nuclear encoded, with 5' leader"),
    TPM = c(all[Category == "Chloroplast", sum(TPM)], all[Category == "Mitochondrial", sum(TPM)], nuc[leader == FALSE, sum(TPM)], pool))
}
res <- rbindlist(out, fill = TRUE); ext <- rbindlist(extra); tt <- rbindlist(tpmtab)
for (i in seq_len(nrow(res))) .add("%-7s %-54s genes %5d of %d (%.2f%%) | carriers %5d tx | unweighted %.2f%% | TSS-weighted %.3f%%%s", res$source[i], res$set[i], res$carrier_genes[i], res$leader_genes[i], res$pct_genes[i],
  res$carrier_transcripts[i], res$pct_unweighted[i], res$pct_TSS[i],
  if (is.na(res$pct_TSS_splicing[i])) "" else sprintf(" | TSS + splicing %.3f%% | loss %.3f points (%.1f%% of the TSS-weighted share; s measured for %.0f%% of it; V21 rule %.3f, V22 rule %.3f points)",
    res$pct_TSS_splicing[i], res$splicing_loss_points[i], res$splicing_loss_pct_of_TSS_weighted[i], res$pct_TSS_weighted_mRNA_with_s_measured[i], res$splicing_loss_points_rule_V21[i], res$splicing_loss_points_rule_V22[i]))
for (i in seq_len(nrow(ext))) .add("%-7s predicted uORFs on transcripts without a detected TuORF or overlapping uORF: %.3f%% | overlapping uORFs on transcripts without a TuORF: %.3f%% | overlapping uORFs more included than the TuORF on TuORF transcripts: %.3f%%",
  ext$source[i], ext$pct_predicted_only[i], ext$pct_ouORF_without_TuORF[i], ext$pct_ouORF_raising_p_on_TuORF_transcripts[i])
.add("pool: %s nuclear transcripts with a 5' leader, %s TPM (CAGE table); genes with a leader-containing transcript: %s", format(res[source == "CAGE"][1, pool_transcripts], big.mark = ","), format(round(res[source == "CAGE"][1, pool_TPM]), big.mark = ","), format(res[source == "CAGE"][1, leader_genes], big.mark = ","))
.add("Carrier transcripts of the combined set (CAGE): %d in %d genes; s_t < 1 for %d transcripts (%d genes); s measured for at least one carried ORF in %d transcripts", nrow(st_out), uniqueN(st_out$gene_id), sum(st_out$s_tx < 1, na.rm = TRUE), uniqueN(st_out[s_tx < 1, gene_id]), sum(st_out$s_measured, na.rm = TRUE))
# ladder table for the Figure 3C / S13 script: long format, one row per source x set x step, plus the TPM rows
L <- rbind(res[, .(source, set, step = "% of genes", value = pct_genes, TPM = NA_real_)], res[, .(source, set, step = "% of nuclear mRNA (unweighted)", value = pct_unweighted, TPM = carrier_TPM)],
           res[, .(source, set, step = "% of nuclear mRNA after TSS weighting", value = pct_TSS, TPM = TSS_weighted_TPM)],
           res[!is.na(pct_TSS_splicing), .(source, set, step = "% of nuclear mRNA after TSS weighting and leader splicing", value = pct_TSS_splicing, TPM = TSS_splicing_TPM)],
           tt[, .(source, set = "pool", step = row, value = NA_real_, TPM)])
fwrite(L, "output/A26_Table_prevalence_ladder.csv"); fwrite(ext, "output/A26_Table_ladder_text_values.csv")
fwrite(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], "output/A26_Table_uORF_pool_fractions.csv")
fwrite(st_out[, .(transcript_id, gene_id, TPM, p, s_tx, s_tx_V21, s_tx_V22, p_max, s_measured, s_quantified, n_orfs)], "output/A26_Table_s_transcript.csv")   # V25 (A26_Table_s_gene.csv until V24)
writeLines(c(sprintf("A26 prevalence ladder and pool shares after TSS weighting and leader splicing — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), "output/A26_summary.txt")
