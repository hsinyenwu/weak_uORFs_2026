################################################################################
# A34 (2026-11-19, Larry: "Please also add leader splicing to the sequenced predicted set and make a figure like Figure 3C for the sequence
# predicted uORF and ouORFs for the effects from TSS and splicing")
# Leader splicing for the sequence-predicted AUG uORFs of module A33 (ending in the leader, and the predicted overlapping uORFs), with the rule
# that A26 uses for the detected sets (splicing_rule.R):
#   * S_k = fraction of RNA-seq reads across the start codon that keep it (the read counting of A2: merged control RNA-seq, unique mappers,
#     dUTP sense reads, splice gaps that are annotated introns of another gene only are not counted); S_k is used when >= 10 reads cover the
#     start codon, s_k = 1 otherwise.
#   * Class D (the start codon lies in an intron of another isoform of the gene, Araport11+CTRL) with the alternative isoform in the expressed
#     set: s_k = min(1, S_k / c_k), c_k = TPM of the expressed isoforms whose exons contain the start codon / TPM of those plus the expressed
#     isoforms whose intron contains it. Every other ORF: s_k = S_k.
#   * A predicted ORF is carried by the transcript in whose leader it was found (as in the canonical predicted-uORF count), so per transcript
#     s_t = max_k(p_k s_k) / max_k(p_k) over its ORFs of the set, and the share after TSS weighting and leader splicing = sum(TPM x p_t x s_t) / pool.
#   Check: the same read counting at the start codons of the 7,423 detected TuORFs must give the S_TIS of module A2.
# Outputs: output/A34_Table_predicted_ORF_s.csv, output/A34_Table_ladder.csv, output/A34_Figure_predicted_uORF_ladder.pdf, output/A34_summary.txt
# Run from Weak_uORF_submission/ after A33:  Rscript revision_modules/A34_predicted_uORFs_TSS_splicing.R
################################################################################
suppressPackageStartupMessages({ library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(GenomicAlignments); library(Rsamtools); library(rtracklayer); library(data.table); library(ggplot2); library(gridExtra); library(grid) })
options(warn = 1)
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
source("revision_modules/splicing_rule.R")   # leader_lengths()
DATA_DIR <- "Files"; OUT_DIR <- "output"; MIN_COV <- 10L
PATH_GTF  <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_GTF2 <- file.path(DATA_DIR, "Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf")
PATH_BAM  <- "/Users/wu/Desktop/CTRL_v1/RNA_CTRL_merged.bam"
ORF <- fread(file.path(OUT_DIR, "A33_Table_predicted_uORFs_with_overlapping.csv"))[type %in% c("ends in the leader", "overlapping the mORF")]
ORF[, `:=`(chr = as.character(chr), oid = .I)]
.add("predicted ORFs from A33: %s ending in the leader, %s overlapping the mORF", format(sum(ORF$type == "ends in the leader"), big.mark = ","), format(sum(ORF$type == "overlapping the mORF"), big.mark = ","))

# ---- annotation: introns of every isoform (full annotation), expressed isoforms ----
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
tx_all <- transcripts(txdb, columns = c("tx_name", "gene_id")); tx2gene <- setNames(as.character(unlist(mcols(tx_all)$gene_id)), mcols(tx_all)$tx_name)
cds_by <- cdsBy(txdb, by = "tx", use.names = TRUE); intr_by0 <- intronsByTranscript(txdb, use.names = TRUE); intr_by0 <- intr_by0[names(intr_by0) %in% names(cds_by)]   # coding isoforms only, as in A2
intr <- unlist(intr_by0, use.names = TRUE); intr$tx_id <- names(intr); names(intr) <- NULL; intr$gene_id <- unname(tx2gene[intr$tx_id])
seqlevels(intr) <- sub("^chr", "", seqlevels(intr), ignore.case = TRUE)
g2 <- import(PATH_GTF2); ex2 <- g2[g2$type == "exon"]; expr_tx <- unique(as.character(ex2$transcript_id))
seqlevels(ex2) <- sub("^chr", "", seqlevels(ex2), ignore.case = TRUE)
tis_gr <- GRanges(ORF$chr, IRanges(ORF$start_bio, ORF$start_bio), strand = ORF$Strand)
ov <- findOverlaps(tis_gr, intr, ignore.strand = FALSE)
pr <- data.table(i = queryHits(ov), tx = intr$tx_id[subjectHits(ov)], gene = intr$gene_id[subjectHits(ov)])
pr <- pr[gene == ORF$gene_id[i] & tx != ORF$tx_id[i]]
ORF[, `:=`(cls_D = oid %in% pr$i, alt_ex = oid %in% pr[tx %in% expr_tx, i])]
.add("class D (start codon inside an intron of another isoform of the gene): %d ORFs; with the alternative isoform in the expressed set: %d", sum(ORF$cls_D), sum(ORF$alt_ex))

# ---- RNA-seq reads at the start codons (the read counting of A2) ----
count_sites <- function(chr, pos, strand, gene) {
  pos_gr <- GRanges(chr, IRanges(pos, pos), strand = strand)
  bam_sl <- names(scanBamHeader(PATH_BAM)[[1]]$targets); stopifnot(all(chr %in% bam_sl)); seqlevels(pos_gr) <- bam_sl
  which_gr <- restrict(reduce(resize(pos_gr, width = 401, fix = "center"), ignore.strand = TRUE), start = 1L)
  param <- ScanBamParam(which = which_gr, what = c("flag", "qname"), mapqFilter = 255,
                        flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE))
  gal <- readGAlignments(PATH_BAM, param = param)
  gal <- gal[!duplicated(paste(mcols(gal)$qname, mcols(gal)$flag, start(gal), cigar(gal)))]
  fl <- mcols(gal)$flag; first_mate <- bitwAnd(fl, 64L) > 0L; aln_strand <- as.character(strand(gal))
  sense <- ifelse(first_mate, ifelse(aln_strand == "+", "-", "+"), aln_strand)
  blk <- grglist(gal); jnc <- junctions(gal); ps <- as.character(strand(pos_gr))
  h <- findOverlaps(pos_gr, blk, ignore.strand = TRUE); okh <- sense[subjectHits(h)] == ps[queryHits(h)]
  R <- tabulate(queryHits(h)[okh], nbins = length(pos_gr))
  gaps <- unlist(jnc, use.names = FALSE); gap_read <- togroup(PartitioningByEnd(jnc))
  hg <- findOverlaps(pos_gr, gaps, ignore.strand = TRUE); okg <- sense[gap_read[subjectHits(hg)]] == ps[queryHits(hg)]
  gk <- paste(as.character(seqnames(gaps))[subjectHits(hg)], start(gaps)[subjectHits(hg)], end(gaps)[subjectHits(hg)])
  ik <- paste(as.character(seqnames(intr)), start(intr), end(intr)); ikg <- paste(ik, intr$gene_id)
  other_gene_gap <- (gk %in% ik) & !(paste(gk, gene[queryHits(hg)]) %in% ikg)
  E <- tabulate(queryHits(hg)[okg & !other_gene_gap], nbins = length(pos_gr))
  data.table(R = R, E = E) }
.add("A34: reading RNA-seq alignments at %s start codons", format(nrow(ORF), big.mark = ","))
rc <- count_sites(ORF$chr, ORF$start_bio, ORF$Strand, ORF$gene_id)
ORF[, `:=`(R_TIS = rc$R, E_TIS = rc$E, cov_TIS = rc$R + rc$E)]; ORF[, S_TIS := fifelse(cov_TIS >= MIN_COV, R_TIS / cov_TIS, NA_real_)]
# check: the same counting at the detected TuORF start codons gives the S_TIS of A2
A2 <- fread(file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.csv")); A2[, chr := sub("_.*$", "", uORF_genome_range)]
A2s <- copy(A2); A2s[, Strand := fifelse(uORF_start_genomic < uORF_stop_genomic, "+", "-")]   # the start codon lies 5' of the stop codon
rcA <- count_sites(A2s$chr, A2s$uORF_start_genomic, A2s$Strand, A2s$gene_id)
okA <- (rcA$R == A2s$R_TIS) & (rcA$E == A2s$E_TIS)
.add("check against A2 (7,423 detected TuORFs): retained and spliced-over read counts at the start codon identical for %d of %d TuORFs", sum(okA), nrow(A2s))

# ---- c_k for class D with an expressed alternative isoform, and s_k ----
LL <- leader_lengths(); leader_tx <- LL[leader_len > 0, transcript_id]
TAB <- lapply(c(CAGE = "CAGE", TSSseq = "TSSseq"), function(src) { a <- fread(sprintf("%s/Table_Transcriptome_Impact_PerTranscript_%s.csv", OUT_DIR, src)); a[grepl("^Nuclear", Category) & transcript_id %in% leader_tx] })
tpm_all <- fread(file.path(OUT_DIR, "Table_Transcriptome_Impact_PerTranscript_CAGE.csv")); tpm <- setNames(tpm_all$TPM, tpm_all$transcript_id)
ex_by <- split(ex2, as.character(ex2$transcript_id)); intr_by <- psetdiff(unlist(range(ex_by)), ex_by)
tg2 <- unique(data.table(tx = as.character(ex2$transcript_id), gene = as.character(ex2$gene_id)))
ORF[, c_share := NA_real_]
for (i in which(ORF$cls_D & ORF$alt_ex & is.finite(ORF$S_TIS))) {
  txs <- intersect(tg2[gene == ORF$gene_id[i], tx], names(ex_by)); if (!length(txs)) next
  pos <- GRanges(ORF$chr[i], IRanges(ORF$start_bio[i], ORF$start_bio[i]))
  car <- txs[vapply(txs, function(t) overlapsAny(pos, ex_by[[t]]), logical(1))]
  alt <- txs[vapply(txs, function(t) length(intr_by[[t]]) > 0 && overlapsAny(pos, intr_by[[t]]), logical(1))]
  a <- sum(tpm[car], na.rm = TRUE); b <- sum(tpm[alt], na.rm = TRUE)
  if (length(alt) && a + b > 0) ORF$c_share[i] <- a / (a + b) }
S1 <- fifelse(is.na(ORF$S_TIS), 1, ORF$S_TIS)
ORF[, s := fifelse(!is.na(c_share) & c_share > 0, pmin(1, S1 / c_share), S1)]
fwrite(ORF[, .(tx_id, gene_id, chr, Strand, start_bio, type, pUP_CAGE, pUP_TSSseq, cls_D, alt_ex, R_TIS, E_TIS, cov_TIS, S_TIS, c_share, s)], file.path(OUT_DIR, "A34_Table_predicted_ORF_s.csv"))
.add("start codons with >= %d reads: %d of %d (%.1f%%); S < 0.90: %d; s < 1: %d", MIN_COV, sum(is.finite(ORF$S_TIS)), nrow(ORF), 100 * mean(is.finite(ORF$S_TIS)), sum(ORF$S_TIS < 0.9, na.rm = TRUE), sum(ORF$s < 1))

# ---- the ladder for the predicted sets ----
SETS <- list(`Predicted uORFs ending in the leader` = "ends in the leader", `Predicted overlapping uORFs` = "overlapping the mORF", `All predicted AUG uORFs` = c("ends in the leader", "overlapping the mORF"))
LAD <- list()
for (src in names(TAB)) { tx <- TAB[[src]]; pool <- sum(tx$TPM); G <- uniqueN(tx$gene_id); pc <- paste0("pUP_", src)
  for (nm in names(SETS)) {
    d <- ORF[type %in% SETS[[nm]], .(p = get(pc) / 100, s), by = .(transcript_id = tx_id)]
    t <- d[, { pm <- max(p); .(p_t = pm, s_t = if (pm > 0) min(1, max(p * s) / pm) else 1) }, by = transcript_id]
    m <- merge(tx[, .(transcript_id, gene_id, TPM)], t, by = "transcript_id", all.x = TRUE); car <- !is.na(m$p_t)
    m[is.na(p_t), `:=`(p_t = 0, s_t = 1)]
    r <- data.table(source = src, set = nm, genes = uniqueN(m[car, gene_id]), pct_genes = 100 * uniqueN(m[car, gene_id]) / G, pct_mRNA = 100 * sum(m$TPM[car]) / pool,
                    pct_TSS = 100 * sum(m$TPM * m$p_t) / pool, pct_TSS_splicing = 100 * sum(m$TPM * m$p_t * m$s_t) / pool,
                    TPM_carriers = sum(m$TPM[car]), TPM_TSS = sum(m$TPM * m$p_t), TPM_TSS_splicing = sum(m$TPM * m$p_t * m$s_t), TPM_pool = pool)
    r[, `:=`(loss_points = pct_TSS - pct_TSS_splicing, loss_pct = 100 * (pct_TSS - pct_TSS_splicing) / pct_TSS)]
    LAD[[length(LAD) + 1]] <- r
    .add("%-7s %-38s genes %5d (%.2f%%) | unweighted %.2f%% | TSS-weighted %.3f%% | TSS + splicing %.3f%% | loss %.3f points (%.1f%% of the TSS-weighted share)",
         src, nm, r$genes, r$pct_genes, r$pct_mRNA, r$pct_TSS, r$pct_TSS_splicing, r$loss_points, r$loss_pct) } }
LAD <- rbindlist(LAD); fwrite(LAD, file.path(OUT_DIR, "A34_Table_ladder.csv"))
a33 <- fread(file.path(OUT_DIR, "A33_Table_ladder.csv"))
chk <- function(src, nm33, nm34) abs(a33[source == src & set == nm33, pct_mRNA_TSS_weighted] - LAD[source == src & set == nm34, pct_TSS]) < 1e-6
.add("check against A33: TSS-weighted shares identical for the leader-ending set (%s) and for all predicted uORFs (%s), CAGE and TSS-seq",
     all(chk("CAGE", "Any predicted AUG uORF (ends in the leader; Figure 3C)", "Predicted uORFs ending in the leader"), chk("TSSseq", "Any predicted AUG uORF (ends in the leader; Figure 3C)", "Predicted uORFs ending in the leader")),
     all(chk("CAGE", "Any predicted AUG uORF, including overlapping uORFs", "All predicted AUG uORFs"), chk("TSSseq", "Any predicted AUG uORF, including overlapping uORFs", "All predicted AUG uORFs")))

# ---- figure in the style of Figure 3C ----
mm <- ORF[, .(p = pUP_CAGE / 100, meas = is.finite(S_TIS)), by = .(transcript_id = tx_id)][, .(p_t = max(p), meas_t = any(meas[p == max(p)])), by = transcript_id]
mm <- merge(mm, TAB$CAGE[, .(transcript_id, TPM)], by = "transcript_id"); MEAS <- 100 * mm[meas_t == TRUE, sum(TPM * p_t)] / mm[, sum(TPM * p_t)]
.add("splicing factor measured (>= %d reads) at the most included start codon for %.1f%% of the TSS-weighted mRNA of all predicted uORFs (CAGE)", MIN_COV, MEAS)
source("revision_modules/fig_compact.R")
L <- LAD[source == "CAGE"]; L2 <- LAD[source == "TSSseq"]
SN <- names(SETS); cols <- c(`Predicted uORFs ending in the leader` = "#BDBDBD", `Predicted overlapping uORFs` = "#FDAE6B", `All predicted AUG uORFs` = "#636363")
steps <- c("% of genes\nwith the ORF", "% of nuclear mRNA\n(unweighted)", "% of nuclear mRNA\nafter TSS weighting", "% after TSS weighting\nand leader splicing")
D <- rbindlist(lapply(SN, function(s) data.table(set = s, step = steps, value = unlist(L[set == s, .(pct_genes, pct_mRNA, pct_TSS, pct_TSS_splicing)]))))
D[, `:=`(set = factor(set, levels = SN), step = factor(step, levels = steps))]
p <- ggplot(D, aes(step, value, fill = set)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", value)), position = position_dodge(0.8), vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = cols, name = NULL) + coord_cartesian(ylim = c(0, 70)) + labs(x = NULL, y = "Percentage",
  title = "Sequence-predicted AUG uORFs: effects of TSS usage and leader splicing",
  subtitle = sprintf("CAGE-seq; with TSS-seq, after TSS weighting and leader splicing: %.1f%%, %.1f%% and %.1f%%", L2$pct_TSS_splicing[1], L2$pct_TSS_splicing[2], L2$pct_TSS_splicing[3])) +
  theme_classic(base_size = 10)
pc <- legend_in(compact_wide(p), 0.995, 0.99)
tab <- data.frame(`Predicted set` = c("Ending in the leader", "Overlapping the mORF", "All"), Genes = format(L$genes, big.mark = ","), `TPM` = format(round(L$TPM_carriers), big.mark = ","), `TPM, TSS` = format(round(L$TPM_TSS), big.mark = ","),
                  `TPM, TSS + splicing` = format(round(L$TPM_TSS_splicing), big.mark = ","), `Splicing loss` = sprintf("%.2f points (%.1f%%)", L$loss_points, L$loss_pct), check.names = FALSE)
tg <- tableGrob(tab, rows = NULL, theme = ttheme_minimal(base_size = 8, base_family = FIG_FONT, core = list(fg_params = list(hjust = 0, x = 0.02)), colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))
cap_txt <- sprintf("Pool: %s nuclear transcripts with a 5' leader (%s genes, %s TPM). Predicted uORFs ending in the leader: AUG ORFs whose stop codon is in the annotated leader. Predicted overlapping uORFs: an AUG in the leader without an in-frame stop codon in the leader, not in frame with the mORF. Leader splicing: retention factor from RNA-seq reads across the start codon, as for the detected uORFs (measured at the most included start codon for %.0f%% of the TSS-weighted mRNA of all predicted uORFs).",
                    format(nrow(TAB$CAGE), big.mark = ","), format(uniqueN(TAB$CAGE$gene_id), big.mark = ","), format(round(sum(TAB$CAGE$TPM)), big.mark = ","), MEAS)
cap <- textGrob(paste(strwrap(cap_txt, width = 135), collapse = "\n"), x = 0.01, hjust = 0, gp = gpar(fontsize = 8, fontfamily = FIG_FONT))
g <- arrangeGrob(pc, tg, cap, ncol = 1, heights = unit(c(3.6, 1.0, 0.6), "in"))
fig_save(file.path(OUT_DIR, "A34_Figure_predicted_uORF_ladder.pdf"), g, h = 5.2)
writeLines(.log, file.path(OUT_DIR, "A34_summary.txt"))
