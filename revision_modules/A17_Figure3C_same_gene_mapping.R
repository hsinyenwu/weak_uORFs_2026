################################################################################
# A17 (2026-10-29, Larry: "check again if 9.9% is correct"): Figure 3C maps each TuORF/ouORF to every expressed
# transcript whose exons contain both its ends, regardless of gene. For nested genes this counts the host gene's
# mRNA as uORF-containing (AT1G07590's 'TuORF' = the coding exons of AT1G07600/MT1A, 8,100 TPM). Here the mapping
# is restricted to transcripts of the uORF's own gene and all Figure 3C quantities are recomputed (CAGE and TSS-seq).
################################################################################
suppressPackageStartupMessages({library(GenomicFeatures); library(GenomicRanges); library(data.table); library(openxlsx)})
t0 <- Sys.time(); .log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF("Files/Araport11+CTRL_20181206.gtf", format = "gtf"))
ex_all <- exonsBy(txdb, by = "tx", use.names = TRUE); cds_all <- cdsBy(txdb, by = "tx", use.names = TRUE); utr5 <- names(fiveUTRsByTranscript(txdb, use.names = TRUE))
gene_of <- function(tx) sub("\\.\\d+$", "", tx)
S1 <- lapply(c(CAGE = "Table S1A (CAGE)", TSSseq = "Table S1B (TSS-seq)"), function(sh) { d <- as.data.table(read.xlsx("output/Table_S1_all_TuORFs.xlsx", sheet = sh))
  d[, .(gene_id, chr = as.character(chr), Strand, g1 = pmin(uORF_start_genomic, uORF_stop_genomic), g2 = pmax(uORF_start_genomic, uORF_stop_genomic), pUP = pct_uORF_plus)] })
oux <- as.data.table(read.xlsx("Files/CiPS_ouORF_all_Final.xlsx"))[, .(gene_id, chr = as.character(chr), Strand, g1 = pmin(ORF_gstart, ORF_gstop), g2 = pmax(ORF_gstart, ORF_gstop))]
oup <- fread("output/Table_ouORFs_pUP.csv")[, .(chr = as.character(chr), g1 = pmin(g_start, g_end), g2 = pmax(g_start, g_end), Strand, pUP_CAGE, pUP_TSS)]
OU <- merge(unique(oux), oup, by = c("chr", "Strand", "g1", "g2"))
.add("ouORFs: %d in the xlsx, %d with pUP", nrow(oux), nrow(OU))
# nested-gene check: TuORFs lying entirely inside the CDS of a transcript of another gene (same strand)
cds_u <- unlist(cds_all); mcols(cds_u)$tx <- rep(names(cds_all), elementNROWS(cds_all))
gr_u <- with(S1$CAGE, GRanges(chr, IRanges(g1, g2), Strand))
h <- findOverlaps(gr_u, cds_u, type = "within"); nest <- data.table(i = queryHits(h), tx = mcols(cds_u)$tx[subjectHits(h)])[gene_of(tx) != S1$CAGE$gene_id[i]]
nest <- unique(nest[, .(uORF_gene = S1$CAGE$gene_id[i], range = paste0(S1$CAGE$chr[i], "_", S1$CAGE$g1[i], "_", S1$CAGE$g2[i]), host_gene = gene_of(tx))])
.add("TuORFs lying inside an annotated CDS of ANOTHER gene (nested genes): %d TuORFs in %d genes", uniqueN(nest$range), uniqueN(nest$uORF_gene)); print(nest[order(uORF_gene)])
res <- list()
for (src in c("CAGE", "TSSseq")) {
  tx <- fread(sprintf("output/Table_Transcriptome_Impact_PerTranscript_%s.csv", src)); tx[, gene := gene_of(transcript_id)]
  ex <- ex_all[intersect(names(ex_all), tx$transcript_id)]; exu <- unlist(ex); mcols(exu)$tx <- rep(names(ex), elementNROWS(ex))
  map_same_gene <- function(d) {   # both ends inside exons of the same transcript, transcript of the uORF's own gene
    gs <- GRanges(d$chr, IRanges(d$g1, width = 1), d$Strand); ge <- GRanges(d$chr, IRanges(d$g2, width = 1), d$Strand)
    hs <- as.data.table(findOverlaps(gs, exu, type = "within"))[, .(i = queryHits, tx = mcols(exu)$tx[subjectHits])]
    he <- as.data.table(findOverlaps(ge, exu, type = "within"))[, .(i = queryHits, tx = mcols(exu)$tx[subjectHits])]
    m <- unique(merge(hs, he, by = c("i", "tx"))); m[, `:=`(pUP = d$pUP[i], same = gene_of(tx) == d$gene_id[i])]
    list(all = m[, .(Max_pUP_any = max(pUP)), by = tx], same = m[same == TRUE, .(Max_pUP_same = max(pUP)), by = tx], cross = m[same == FALSE]) }
  tu <- map_same_gene(S1[[src]]); ou <- map_same_gene(OU[, .(gene_id, chr, Strand, g1, g2, pUP = if (src == "CAGE") pUP_CAGE else pUP_TSS)])
  tx <- merge(tx, tu$all, by.x = "transcript_id", by.y = "tx", all.x = TRUE); tx <- merge(tx, tu$same, by.x = "transcript_id", by.y = "tx", all.x = TRUE)
  tx <- merge(tx, setnames(ou$all, "Max_pUP_any", "ou_any"), by.x = "transcript_id", by.y = "tx", all.x = TRUE); tx <- merge(tx, setnames(ou$same, "Max_pUP_same", "ou_same"), by.x = "transcript_id", by.y = "tx", all.x = TRUE)
  for (cc in c("Max_pUP_any", "Max_pUP_same", "ou_any", "ou_same")) tx[is.na(get(cc)), (cc) := 0]
  .add("[%s] reproduction of the canonical mapping: transcripts with a mapped TuORF %d (canonical %d); Max_pUP identical in %.1f%% of them", src, tx[Max_pUP_any > 0 | transcript_id %in% tu$all$tx, .N], tx[TuORF_mapped == TRUE, .N], 100 * tx[TuORF_mapped == TRUE, mean(abs(Max_pUP_any - Max_pUP_TuORF) < 0.01)])
  lost <- tx[TuORF_mapped == TRUE & !(transcript_id %in% tu$same$tx)][order(-TPM)]
  .add("[%s] transcripts whose only TuORF mapping is from ANOTHER gene: %d, TPM %.0f (%.2f%% of 5'-leader mRNA); TSS-weighted TuORF mass removed %.2f%%", src, nrow(lost), sum(lost$TPM), 100 * sum(lost$TPM) / tx[Category == "Nuclear_with_5UTR", sum(TPM)], 100 * lost[, sum(TPM * Max_pUP_TuORF / 100)] / tx[Category == "Nuclear_with_5UTR", sum(TPM)])
  if (src == "CAGE") print(lost[1:min(8, nrow(lost)), .(transcript_id, TPM = round(TPM, 1), Max_pUP_TuORF = round(Max_pUP_TuORF, 1))])
  lost_ou <- tx[ouORF_mapped == TRUE & !(transcript_id %in% ou$same$tx)]; .add("[%s] transcripts whose only ouORF mapping is from another gene: %d, TPM %.0f", src, nrow(lost_ou), sum(lost_ou$TPM))
  five <- tx[Category == "Nuclear_with_5UTR"]; den <- sum(five$TPM); nuc <- tx[Category %in% c("Nuclear_with_5UTR", "Nuclear_no_5UTR")]
  pc <- function(x) sprintf("%.1f%%", 100 * x); pc2 <- function(x) sprintf("%.2f%%", 100 * x)
  genes <- five[, .(tu_old = any(Max_pUP_TuORF > 0), tu_new = any(Max_pUP_same > 0), up_old = any(Max_pUP_TuORF > 0 | ouORF_mapped), up_new = any(Max_pUP_same > 0 | ou_same > 0), pred = any(n_pred_uORFs > 0)), by = gene_id]
  r <- data.table(source = src, quantity = c("genes with leader (n)", "% genes with TuORF", "% genes with TuORF or ouORF", "% genes with predicted uORF",
      "% mRNA from TuORF transcripts (unweighted)", "% mRNA, TuORF, TSS-weighted", "% mRNA from TuORF/ouORF transcripts (unweighted)", "% mRNA, TuORF/ouORF, TSS-weighted", "% mRNA from predicted-uORF transcripts (unweighted)", "% mRNA, predicted, TSS-weighted", "% mRNA, predicted-only (no TuORF/ouORF), TSS-weighted"),
    canonical_mapping = c(nrow(genes), pc(mean(genes$tu_old)), pc(mean(genes$up_old)), pc(mean(genes$pred)),
      pc(five[TuORF_mapped == TRUE, sum(TPM)] / den), pc2(nuc[, sum(TPM * Max_pUP_TuORF / 100)] / den), pc(five[TuORF_mapped | ouORF_mapped, sum(TPM)] / den), pc2(nuc[, sum(TPM * pmax(Max_pUP_TuORF, Max_pUP_ouORF) / 100)] / den),
      pc(five[n_pred_uORFs > 0, sum(TPM)] / den), pc2(nuc[, sum(TPM * Max_pUP_pred / 100)] / den), pc2(nuc[!TuORF_mapped & !ouORF_mapped, sum(TPM * Max_pUP_pred / 100)] / den)),
    same_gene_mapping = c(nrow(genes), pc(mean(genes$tu_new)), pc(mean(genes$up_new)), pc(mean(genes$pred)),
      pc(five[Max_pUP_same > 0 | transcript_id %in% tu$same$tx, sum(TPM)] / den), pc2(five[, sum(TPM * Max_pUP_same / 100)] / den), pc(five[transcript_id %in% c(tu$same$tx, ou$same$tx), sum(TPM)] / den), pc2(five[, sum(TPM * pmax(Max_pUP_same, ou_same) / 100)] / den),
      pc(five[n_pred_uORFs > 0, sum(TPM)] / den), pc2(five[, sum(TPM * Max_pUP_pred / 100)] / den), pc2(five[!(transcript_id %in% c(tu$same$tx, ou$same$tx)), sum(TPM * Max_pUP_pred / 100)] / den)))
  res[[src]] <- r; print(r)
  fwrite(tx[, .(transcript_id, gene_id, Category, TPM, TuORF_mapped, Max_pUP_TuORF, Max_pUP_TuORF_same_gene = Max_pUP_same, ouORF_mapped, Max_pUP_ouORF, Max_pUP_ouORF_same_gene = ou_same, n_pred_uORFs, Max_pUP_pred)], sprintf("output/A17_PerTranscript_same_gene_%s.csv", src))
}
fwrite(rbindlist(res), "output/A17_Figure3C_same_gene_vs_canonical.csv"); fwrite(nest, "output/A17_TuORFs_inside_other_gene_CDS.csv")
.add("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))); writeLines(.log, "output/A17_summary.txt")
