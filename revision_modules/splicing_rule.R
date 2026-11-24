################################################################################
# splicing_rule.R — the per-ORF splicing retention factor s_k used by A26 (pool shares) and A18 (mRNA-weighted shares).
# Sourced by both modules (V23, 2026-11-08; the V22 rule set s_k = 1 for every class-D ORF whose alternative isoform is
# quantified, which overcorrected: for MS1 (AT5G17920) 12,794 junction reads splice over the start codon and 15 keep it,
# yet the isoform quantification gives the carrier isoform 1,307 TPM and the alternative isoform 398 TPM).
#
# Rule (junction-consistent):
#   * S_k = fraction of RNA-seq reads across the start codon that keep it (A2 for TuORFs, A25 for overlapping uORFs);
#     s_k = 1 where fewer than MIN_COV reads cover the start codon.
#   * For a class-D ORF (start codon inside an intron of another isoform) whose alternative isoform is in the expressed isoform set,
#     the isoform quantification has already split the gene's mRNA between the carrier isoforms (exons contain the start codon) and
#     the alternative isoforms (intron contains it). The junction reads are taken as the measure of that split: the share of the
#     gene's mRNA that keeps the start codon is S_k, the carrier isoforms hold the share c_k = TPM(carriers) / TPM(carriers + alternative),
#     so the factor applied to the carrier TPM is s_k = min(1, S_k / c_k). (If the quantification agrees with the junction reads,
#     S_k = c_k and s_k = 1, which is the V22 case; if it disagrees, the junction reads decide.)
#   * Every other ORF: s_k = S_k (the V21 rule).
# orf_s_table() returns one row per ORF with S_k, c_k and s_k under the three rules (V21, V22, V23) for the record.
#
# V25 (2026-11-09): the factor is assigned per transcript, over the ORFs of the set that the transcript carries (start and stop codon in its
# exons; the same ORFs that define the transcript's p): s_t = max_k(p_k s_k) / max_k(p_k) over the carried ORFs k (orf_carriers(), tx_s()).
# In V23 and V24 the factor was a gene value over every ORF of the gene in the audit tables (gene_s(), kept for the record), so TuORFs that
# no expressed isoform carries (673 of the 7,423 TuORFs; among them 323 class-D TuORFs with an expressed alternative isoform, measured S and
# c_k = 0, for which the code fell back to s_k = S_k although their p is close to 1 and few reads keep the start codon) lowered the gene factor
# of 49 genes, 44 of which have every carried TuORF fully retained (for example AT1G07890: two uncarried class-D TuORFs with p = 0.9999 and
# S about 0.005, one carried TuORF with S = 1, gene factor 0.0071). The division by c_k is applied to the class-D ORFs with c_k > 0 (311 ORFs:
# 289 TuORFs and 22 overlapping uORFs); the 375 with c_k = 0 (323 and 52) have their start codon in no exon of an expressed isoform, so no
# transcript carries them and they never enter a share. min_cov is documented only: A2 and A25 already set S to NA below 10 reads.
# TuORF pool share after TSS weighting and splicing: 8.538% (V23/V24) -> 8.546% (per-transcript factor). Reviewed 2026-11-09 (V25 review).
################################################################################
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer) })
orf_s_table <- function(include_ouORF = TRUE, min_cov = 10L,
                        path_A2 = "output/A2_Table_TuORF_leader_splicing.csv", path_A25 = "output/A25_Table_ouORF_leader_splicing.csv",
                        path_tx = "output/Table_Transcriptome_Impact_PerTranscript_CAGE.csv",
                        path_gtf = "Files/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf") {
  A2 <- fread(path_A2)
  orfs <- A2[, .(gene_id, kind = "TuORF", uid = uORF_genome_range, chr = sub("_.*$", "", uORF_genome_range), tis = uORF_start_genomic, stop = uORF_stop_genomic,   # V25: chr and stop codon position for orf_carriers()
                 p_CAGE = pct_uORF_plus / 100, p_TSSseq = pUP_TSSseq / 100,
                 S = S_TIS, cov = cov_TIS, cls = substr(annot_class, 1, 1), alt_ex = other_TIS_in_intron_ex %in% TRUE)]
  if (include_ouORF) {
    A25 <- fread(path_A25)
    orfs <- rbind(orfs, A25[, .(gene_id, kind = "ouORF", uid = ORF_id2, chr = as.character(chr), tis = start_codon_genomic, stop = fifelse(Strand == "+", ouORF_genome_range_max, ouORF_genome_range_min),
                                p_CAGE = pct_uORF_plus_CAGE / 100, p_TSSseq = pct_uORF_plus_TSSseq / 100,
                                S = S_TIS, cov = cov_TIS, cls = substr(annot_class, 1, 1), alt_ex = other_TIS_in_intron_ex %in% TRUE)]) }
  # carrier share c_k from the expressed isoform set (RSEM TPM of the per-transcript table; exons from the same GTF)
  tx <- fread(path_tx); tpm <- setNames(tx$TPM, tx$transcript_id)
  gtf <- import(path_gtf); ex <- gtf[gtf$type == "exon"]; ex_by <- split(ex, as.character(ex$transcript_id))
  tx_gene <- unique(data.table(tx = as.character(ex$transcript_id), gene = as.character(ex$gene_id)))
  intr_by <- psetdiff(unlist(range(ex_by)), ex_by)
  orfs[, c_share := NA_real_]
  idx <- which(orfs$cls == "D" & orfs$alt_ex & is.finite(orfs$S))
  for (i in idx) {
    txs <- intersect(tx_gene[gene == orfs$gene_id[i], tx], names(ex_by)); if (!length(txs)) next
    chr <- as.character(seqnames(ex_by[[txs[1]]]))[1]; pos <- GRanges(chr, IRanges(orfs$tis[i], orfs$tis[i]))
    car <- txs[vapply(txs, function(t) overlapsAny(pos, ex_by[[t]]), logical(1))]
    alt <- txs[vapply(txs, function(t) length(intr_by[[t]]) > 0 && overlapsAny(pos, intr_by[[t]]), logical(1))]
    a <- sum(tpm[car], na.rm = TRUE); b <- sum(tpm[alt], na.rm = TRUE)
    if (length(alt) && a + b > 0) orfs$c_share[i] <- a / (a + b) }
  S1 <- fifelse(is.na(orfs$S), 1, orfs$S)
  orfs[, `:=`(s_V21 = S1, s_V22 = fifelse(cls == "D" & alt_ex, 1, S1), s = fifelse(!is.na(c_share) & c_share > 0, pmin(1, S1 / c_share), S1), s_measured = is.finite(S))]
  orfs[]
}
gene_s <- function(orfs, pcol = "p_CAGE", scol = "s") {   # V23/V24 rule, kept for the record: s for a gene = max_k(p_k s_k) / max_k p_k over every ORF of the gene in the audit tables
  orfs[, { p <- get(pcol); s <- get(scol); pm <- max(p, na.rm = TRUE)
    .(s_gene = if (is.finite(pm) && pm > 0) min(1, max(p * s, na.rm = TRUE) / pm) else 1, s_measured = any(s_measured), n_orfs = .N) }, by = gene_id] }
# V25 (2026-11-09): the ORFs each expressed isoform carries, and the per-transcript factor over them (used by A26 and A18).
orf_carriers <- function(orfs, path_gtf = "Files/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf") {
  # one row per (ORF, transcript of the same gene in the expressed isoform set whose exons contain both the start codon and the stop codon)
  gtf <- import(path_gtf); ex <- gtf[gtf$type == "exon"]
  exr <- GRanges(as.character(seqnames(ex)), ranges(ex)); ex_tx <- as.character(ex$transcript_id); ex_gene <- as.character(ex$gene_id)
  hits <- function(pos) { h <- findOverlaps(GRanges(orfs$chr, IRanges(pos, pos)), exr, ignore.strand = TRUE)
    unique(data.table(i = queryHits(h), transcript_id = ex_tx[subjectHits(h)], gene = ex_gene[subjectHits(h)])) }
  m <- merge(hits(orfs$tis), hits(orfs$stop), by = c("i", "transcript_id", "gene"))
  m <- m[gene == orfs$gene_id[i]]
  unique(m[, .(uid = orfs$uid[i], gene_id = gene, transcript_id)]) }
tx_s <- function(carriers, orfs, pcol = "p_CAGE") {
  # s for a transcript = max_k(p_k s_k) / max_k p_k over the ORFs k of `orfs` that the transcript carries (its best effective inclusion);
  # the V21 and V22 variants are returned for the record; transcripts that carry no ORF of `orfs` are absent (the caller sets s = 1 for them)
  d <- merge(carriers, orfs[, c("uid", "gene_id", pcol, "s", "s_V21", "s_V22", "s_measured", "c_share"), with = FALSE], by = c("uid", "gene_id"))
  d[, { p <- get(pcol); pm <- max(p, na.rm = TRUE); f <- function(sv) if (is.finite(pm) && pm > 0) min(1, max(p * sv, na.rm = TRUE) / pm) else 1
    .(s_tx = f(s), s_tx_V21 = f(s_V21), s_tx_V22 = f(s_V22), p_max = pm, s_measured = any(s_measured), s_quantified = any(!is.na(c_share)), n_orfs = .N) }, by = transcript_id] }

# V23 (2026-11-08): one leader definition per transcript for numerator and denominator of the pool shares (A26, A18): a transcript of the
# expressed isoform set has a 5' leader when exonic bases lie upstream of its CDS start (the expressed-isoform GTF has no UTR features).
leader_lengths <- function(path_gtf = "Files/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf") {
  gtf <- import(path_gtf); ex <- as.data.table(gtf[gtf$type == "exon"]); cds <- as.data.table(gtf[gtf$type == "CDS"])
  ex[, strand := as.character(strand)]; cds[, strand := as.character(strand)]
  cs <- cds[, .(cds_start = if (strand[1] == "+") min(start) else max(end), str = strand[1]), by = transcript_id]
  ed <- merge(ex[, .(transcript_id, start, end)], cs, by = "transcript_id")
  ed[, up := fifelse(str == "+", pmax(0, pmin(end, cds_start - 1) - start + 1), pmax(0, end - pmax(start, cds_start + 1) + 1))]
  ed[, .(leader_len = sum(up)), by = transcript_id] }
