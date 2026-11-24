################################################################################
# A33 (2026-11-19, Larry: "Can you also add the numbers if including overlapping uORFs for sequenced-predicted uORFs?")
# The "any predicted AUG uORF" row of Figure 3C counts only AUG ORFs whose stop codon lies in the annotated leader (canonical script:
# ORFik::findMapORFs on the 5' UTR sequence). Here the same ladder is computed again, and then with the predicted overlapping uORFs added:
# AUGs in the leader whose ORF has no in-frame stop codon in the leader and is not in frame with the mORF (it runs into the mORF in another
# frame). AUGs in frame with the mORF without a stop in the leader are N-terminal extensions of the mORF, not uORFs, and are not counted.
#   Pool, carriers, % of genes, % of mRNA and the TSS weighting are those of A26: nuclear transcripts with a leader from the per-transcript
#   table of the canonical script (output/Table_Transcriptome_Impact_PerTranscript_{CAGE,TSSseq}.csv); p = highest %uORF+/100 among the
#   transcript's ORFs of the set. %uORF+ as in the canonical script: signal at or upstream of the start codon / signal in the window
#   (annotated leader + 100 nt upstream), on genomic coordinates (.make_tss_windows, .compute_pUP_for_targets).
#   Check: the leader-ending set must give the canonical n_pred_uORFs and Max_pUP_pred of every transcript and the A26 values.
# Outputs: output/A33_Table_predicted_uORFs_with_overlapping.csv (one row per predicted AUG ORF), output/A33_Table_ladder.csv, output/A33_summary.txt
# Run from Weak_uORF_submission/:  Rscript revision_modules/A33_predicted_uORFs_with_overlapping.R
################################################################################
suppressPackageStartupMessages({ library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(Biostrings); library(Rsamtools); library(data.table); library(ORFik) })
options(warn = 1)
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
source("revision_modules/splicing_rule.R")   # leader_lengths()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF2 <- file.path(DATA_DIR, "Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf")
PATH_FASTA <- file.path(DATA_DIR, "TAIR10_chr_all_2.fas")
PATHS <- c(CAGE = file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt"), TSSseq = file.path(DATA_DIR, "wt_TSS_Nielsen_2019.txt"))
STOPS <- c("TAA", "TAG", "TGA")

# ---- pool and leaders (as A26) ----
LL <- leader_lengths(); leader_tx <- LL[leader_len > 0, transcript_id]
TAB <- lapply(c(CAGE = "CAGE", TSSseq = "TSSseq"), function(src) { a <- fread(sprintf("%s/Table_Transcriptome_Impact_PerTranscript_%s.csv", OUT_DIR, src)); a[grepl("^Nuclear", Category) & transcript_id %in% leader_tx] })
pool_tx <- TAB$CAGE$transcript_id; stopifnot(setequal(pool_tx, TAB$TSSseq$transcript_id))
.add("pool: %s nuclear transcripts with a leader in %s genes", format(length(pool_tx), big.mark = ","), format(uniqueN(TAB$CAGE$gene_id), big.mark = ","))

# ---- leader sequences (canonical: fiveUTRsByTranscript of the expressed-isoform GTF) ----
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF2, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF2, format = "gtf")
utr <- fiveUTRsByTranscript(txdb, use.names = TRUE); utr <- utr[elementNROWS(utr) > 0]
fa <- FaFile(PATH_FASTA); fa_lv <- seqlevels(scanFaIndex(fa))
lv <- seqlevels(utr); map <- setNames(lv, lv); bad <- !lv %in% fa_lv
if (any(bad)) { alt <- paste0("Chr", sub("^chr", "", lv[bad], ignore.case = TRUE)); map[bad] <- ifelse(alt %in% fa_lv, alt, lv[bad]) }
utr <- renameSeqlevels(utr, map)
utr <- utr[names(utr) %in% pool_tx]
.add("leaders with an annotated 5' UTR in the pool: %s of %s", format(length(utr), big.mark = ","), format(length(pool_tx), big.mark = ","))
seqs <- extractTranscriptSeqs(fa, utr)

# ---- (1) ORFs that end in the leader: the canonical call ----
orf_leader <- ORFik::findMapORFs(utr, seqs, startCodon = "ATG", longestORF = FALSE, groupByTx = FALSE, minimumLength = 0)
rng <- unlist(range(orf_leader)); st <- as.character(strand(rng))
A <- data.table(tx_id = sub("_[0-9]+$", "", names(orf_leader)), chr = as.character(seqnames(rng)), Strand = st,
                start_bio = ifelse(st == "+", start(rng), end(rng)), type = "ends in the leader")
# ---- (2) AUGs without an in-frame stop in the leader ----
atg <- vmatchPattern("ATG", seqs); stp <- lapply(STOPS, function(s) vmatchPattern(s, seqs))
B <- rbindlist(lapply(seq_along(seqs), function(i) {
  a <- start(atg[[i]]); if (!length(a)) return(NULL)
  L <- width(seqs)[i]; s_all <- sort(unlist(lapply(stp, function(x) start(x[[i]]))))
  has_stop <- vapply(a, function(k) any(s_all > k & (s_all - k) %% 3L == 0L), logical(1))
  k <- a[!has_stop]; if (!length(k)) return(NULL)
  data.table(tx_id = names(seqs)[i], tpos = k, inframe = ((L + 1L) - k) %% 3L == 0L) }))
gBl <- pmapFromTranscripts(IRanges(B$tpos, width = 1L), utr[B$tx_id]); gBl <- gBl[width(gBl) > 0L]   # one range per exon is returned; keep the exon that holds the base
stopifnot(all(elementNROWS(gBl) == 1L)); gB <- unlist(gBl, use.names = FALSE)
B[, `:=`(chr = as.character(seqnames(gB)), Strand = as.character(strand(gB)), start_bio = start(gB))]
B[, type := fifelse(inframe, "in frame with the mORF (N-terminal extension, not counted)", "overlapping the mORF")]
ORF <- rbind(A, B[, .(tx_id, chr, Strand, start_bio, type)])
ORF[, gene_id := substr(tx_id, 1, 9)]
.add("predicted AUG ORFs in the pool leaders: %s ending in the leader; %s overlapping the mORF; %s in frame with the mORF (not counted)",
     format(sum(ORF$type == "ends in the leader"), big.mark = ","), format(sum(ORF$type == "overlapping the mORF"), big.mark = ","), format(sum(grepl("^in frame", ORF$type)), big.mark = ","))

# ---- %uORF+ with the canonical windows and formula ----
ul <- unlist(utr, use.names = FALSE); txv <- rep(names(utr), elementNROWS(utr))
plus <- as.character(strand(ul)) == "+"; key <- ifelse(plus, start(ul), -end(ul)); o <- order(txv, key); first <- o[!duplicated(txv[o])]
ul[first] <- resize(ul[first], width = width(ul[first]) + 100L, fix = "end"); ul <- restrict(ul, start = 1L)
WIN <- split(ul, txv)
pUP_for <- function(path) {
  raw <- fread(path, header = FALSE, col.names = c("count", "chr", "pos", "strand"))
  raw[, chr := gsub("^chr|^Chr", "", chr, ignore.case = TRUE)]
  gr <- GRanges(raw$chr, IRanges(raw$pos, width = 1L), strand = raw$strand, count = raw$count)
  base <- function(x) sub("^chr", "", x, ignore.case = TRUE); wl <- seqlevels(WIN); m_ <- setNames(wl, base(wl))   # file chr "1" -> window level (e.g. "Chr1")
  gr <- gr[base(as.character(seqnames(gr))) %in% names(m_)]; gr <- keepSeqlevels(gr, unique(as.character(seqnames(gr))), pruning.mode = "coarse"); gr <- renameSeqlevels(gr, setNames(m_[base(seqlevels(gr))], seqlevels(gr)))
  keep <- intersect(seqlevels(gr), seqlevels(WIN)); gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse"); W <- keepSeqlevels(WIN, keep, pruning.mode = "coarse")
  ov <- findOverlaps(gr, W)
  S <- data.table(tx_id = names(W)[subjectHits(ov)], coord = start(gr)[queryHits(ov)], count = gr$count[queryHits(ov)])
  tot <- S[, .(tc = sum(count)), by = tx_id]
  Q <- ORF[, .(tx_id, Strand, start_bio)][, row := .I]
  J <- S[Q, on = .(tx_id), allow.cartesian = TRUE]
  R <- J[, .(cb = sum(count[ifelse(Strand == "+", coord <= start_bio, coord >= start_bio)], na.rm = TRUE)), by = .(row, tx_id)]
  R <- merge(R, tot, by = "tx_id", all.x = TRUE)
  setorder(R, row); round(ifelse(is.na(R$tc) | R$tc == 0, 0, 100 * R$cb / R$tc), 3) }
for (src in names(PATHS)) set(ORF, j = paste0("pUP_", src), value = pUP_for(PATHS[[src]]))
fwrite(ORF, file.path(OUT_DIR, "A33_Table_predicted_uORFs_with_overlapping.csv"))

# ---- check against the canonical per-transcript table, then the ladder ----
LAD <- list()
for (src in names(PATHS)) {
  tx <- copy(TAB[[src]]); pool <- sum(tx$TPM); G <- uniqueN(tx$gene_id); pc <- paste0("pUP_", src)
  agg <- function(sel) ORF[sel, .(n = .N, p = max(get(pc))), by = .(transcript_id = tx_id)]
  a1 <- agg(ORF$type == "ends in the leader"); a2 <- agg(ORF$type %in% c("ends in the leader", "overlapping the mORF")); a3 <- agg(ORF$type == "overlapping the mORF")
  chk <- merge(tx[, .(transcript_id, n_pred_uORFs, Max_pUP_pred)], a1, by = "transcript_id", all.x = TRUE)
  chk[is.na(n), `:=`(n = 0L, p = 0)]
  .add("%-7s check against the canonical table: n_pred_uORFs equal for %d of %d transcripts; |Max_pUP_pred - p| <= 0.001 for %d of %d", src,
       sum(chk$n_pred_uORFs == chk$n), nrow(chk), sum(abs(chk$Max_pUP_pred - chk$p) <= 0.001), nrow(chk))
  for (nm in c("Any predicted AUG uORF (ends in the leader; Figure 3C)", "Any predicted AUG uORF, including overlapping uORFs", "Predicted overlapping uORFs only (any transcript)")) {
    a <- switch(substr(nm, 1, 30), "Any predicted AUG uORF (ends i" = a1, "Any predicted AUG uORF, includ" = a2, a3)
    d <- merge(tx[, .(transcript_id, gene_id, TPM)], a, by = "transcript_id", all.x = TRUE); d[is.na(n), `:=`(n = 0L, p = 0)]
    r <- data.table(source = src, set = nm, genes = uniqueN(d[n > 0, gene_id]), genes_total = G, pct_genes = 100 * uniqueN(d[n > 0, gene_id]) / G,
                    carriers = sum(d$n > 0), pct_mRNA_unweighted = 100 * d[n > 0, sum(TPM)] / pool, pct_mRNA_TSS_weighted = 100 * d[, sum(TPM * p / 100)] / pool)
    LAD[[length(LAD) + 1]] <- r
    .add("%-7s %-55s genes %5d of %d (%.2f%%) | carriers %5d tx | unweighted %.2f%% | TSS-weighted %.3f%%", src, nm, r$genes, G, r$pct_genes, r$carriers, r$pct_mRNA_unweighted, r$pct_mRNA_TSS_weighted) } }
LAD <- rbindlist(LAD); fwrite(LAD, file.path(OUT_DIR, "A33_Table_ladder.csv"))
# detected overlapping uORFs among the predicted overlapping ones (start codons)
ou <- as.data.table(openxlsx::read.xlsx(file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx"))); ou[, chr := gsub("^chr|^Chr", "", as.character(chr), ignore.case = TRUE)]
ou[, sc := ifelse(Strand == "+", pmin(as.integer(ORF_gstart), as.integer(ORF_gstop)), pmax(as.integer(ORF_gstart), as.integer(ORF_gstop)))]
k_ou <- unique(paste(ou$gene_id, ou$Strand, ou$sc)); ORF[, k := paste(gene_id, Strand, start_bio)]
.add("detected overlapping uORF start codons (%d): among the predicted overlapping uORFs %d, among the leader-ending ORFs %d, in neither %d (not on a pool leader or another isoform)",
     length(k_ou), sum(k_ou %in% ORF[type == "overlapping the mORF", k]), sum(k_ou %in% ORF[type == "ends in the leader", k]), sum(!k_ou %in% ORF$k))
writeLines(.log, file.path(OUT_DIR, "A33_summary.txt"))
