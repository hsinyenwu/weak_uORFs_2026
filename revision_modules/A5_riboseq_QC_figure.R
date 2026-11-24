################################################################################
# A5 — Ribo-seq quality-control figure (Reviewer 2, §2.3)                       #
#   "read-length distribution, reading-frame distribution per read length,     #
#    metagene profiles at start/stop codons of CDSs and of TuORFs"              #
#                                                                               #
# Panels                                                                        #
#  A  read-length distribution of nuclear primary alignments                    #
#     (random 10% subsample of the six merged Col-0 libraries; the full-library  #
#     distribution from Ribo-seQC is used as a cross-check)                      #
#  B  codon position of the read 5' end, per read length, for reads whose 5' end #
#     falls inside an annotated CDS (first/last 15 nt excluded); the offsets     #
#     used in this study (24-28 nt -> 8-12 nt) and the resulting in-frame share  #
#  C  P-site metagene at mORF start and stop codons (the paper's P-site table,   #
#     i.e. the P-sites behind uIF, Fig S12 and Fig S13), coloured by frame       #
#  D  the same at the start and stop codons of the TuORFs of Table S1A of at     #
#     least 2 aa: in a 1-aa uORF (AUG directly followed by the stop codon) the   #
#     start codon is also the last sense codon, so initiation and termination    #
#     signals fall on one codon (Wu et al. 2024, Plant Cell; Larry, 2026-10-29)  #
#  E  per-ORF in-frame fraction of P-sites, TuORFs of >= 2 aa vs mORFs (>= 10    #
#     P-sites; 1-aa uORFs left out as in D)                                     #
#  Termination rule (Larry, 2026-10-29; Wu et al. 2024 Plant Cell, Fig. 9B-E):   #
#     P-sites at the second nucleotide of the last sense codon (frame +1; the    #
#     1-nt 3' shift of terminating-ribosome footprints) count as in frame in the #
#     percentages of C, D and E; the bars keep the frame colours (-2 stays blue); #
#     the plain-frame values are kept in the tables (*_raw).                     #
# Transcripts: one isoform per gene = the isoforms used in the paper (TuORF      #
# source isoforms, Table S1A; No-TuORF controls = highest-expressed isoform).    #
# Outputs (output/A5_*): per-length table, per-ORF table, metagene tables,       #
# summary txt, figure pdf.                                                       #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A5_riboseq_QC_figure.R
# The BAM path can be overridden with the environment variable A5_BAM.          #
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomicAlignments); library(Rsamtools)
  library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(gridExtra); library(grid)
})
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF   <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_PS    <- file.path(DATA_DIR, "CTRL_expressed_P_sites_sort_count")   # count, chr, pos, strand (the paper's P-site table)
PATH_S1    <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_RSEM_RNA  <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_ISO   <- file.path(DATA_DIR, "RSEM_isoforms_TPM_6samples.tsv")
PATH_5UTR  <- file.path(DATA_DIR, "genes_w_5UTR.tsv")
PATH_BAM   <- Sys.getenv("A5_BAM", "/Users/wu/Desktop/Desktop_2026_July/RiboseQC_CTRL/ribo_CTRL_merged_10p.bam")
PATH_RQC   <- Sys.getenv("A5_RIBOSEQC", "/Users/wu/Desktop/Desktop_2026_July/RiboseQC_CTRL/ribo_CTRL_merged.bam_results_RiboseQC")
OFFSETS <- c(`24` = 8L, `25` = 9L, `26` = 10L, `27` = 11L, `28` = 12L)   # the P-site rule of this study (RiboTaper settings)
LEN_SHOW <- 20:35; TRIM <- 15L; MIN_PS <- 10L; FLANK <- 30L; INSIDE <- 60L
for (f in c(PATH_GTF, PATH_PS, PATH_S1, PATH_OUORF, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_ISO, PATH_5UTR, PATH_BAM)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
big <- function(x) trimws(format(x, big.mark = ","))

################################################################################
# 1) Annotation and the transcripts used in the paper
################################################################################
.add("A5: TxDb")
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
ex  <- exonsBy(txdb, by = "tx", use.names = TRUE); cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
utr <- fiveUTRsByTranscript(txdb, use.names = TRUE)

S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)") %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"), chr = as.character(chr),
  g_min = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\1", uORF_genome_range)), g_max = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\2", uORF_genome_range)))
ou  <- read.xlsx(PATH_OUORF)
# No-TuORF controls exactly as in Fig 1C: RNA TPM > 1, annotated 5' UTR, no TuORF and no ouORF; highest-expressed isoform
rna <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); ribo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
TE_tbl <- data.frame(gene_id = rna$gene_id, RNA_TPM = rna$TPM, Ribo_TPM = ribo$TPM[match(rna$gene_id, ribo$gene_id)]) %>%
  mutate(TE = Ribo_TPM / RNA_TPM) %>% filter(gene_id %in% read.delim(PATH_5UTR)$gene_id, RNA_TPM > 1, is.finite(TE))
upstream_translated <- unique(c(S1A$gene_id, as.character(ou$gene_id)))
iso <- fread(PATH_ISO); iso[, mTPM := rowMeans(.SD), .SDcols = patterns("^TPM_")]
iso <- iso[gene_id %in% setdiff(TE_tbl$gene_id, upstream_translated) & transcript_id %in% names(utr)]
setorder(iso, gene_id, -mTPM); ctrl_tx <- iso[!duplicated(gene_id), .(gene_id, tx_id = transcript_id)]
tu_tx  <- S1A %>% distinct(gene_id, .keep_all = TRUE) %>% select(gene_id, tx_id)
rep_tx <- unique(c(S1A$tx_id, ctrl_tx$tx_id))
rep_tx <- rep_tx[rep_tx %in% names(ex) & rep_tx %in% names(cds)]
.add(sprintf("  transcripts: %d TuORF source isoforms (%d genes) + %d No-TuORF control isoforms = %d one-per-gene transcripts with a CDS",
             length(unique(S1A$tx_id)), length(unique(S1A$gene_id)), nrow(ctrl_tx), length(rep_tx)))
utr5_len <- setNames(rep(0L, length(rep_tx)), rep_tx); u <- sum(width(utr[intersect(rep_tx, names(utr))])); utr5_len[names(u)] <- as.integer(u)
cds_len  <- as.integer(sum(width(cds[rep_tx]))); names(cds_len) <- rep_tx
info <- data.table(tx = rep_tx, cds_s = utr5_len[rep_tx] + 1L, cds_e = utr5_len[rep_tx] + cds_len[rep_tx], cds_len = cds_len[rep_tx])

################################################################################
# 2) Read lengths and 5'-end codon position from the BAM (10% subsample)
################################################################################
.add("A5: BAM pass (", basename(PATH_BAM), ")")
bf <- BamFile(PATH_BAM); sl <- seqlengths(seqinfo(bf))
flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
cds_rep <- cds[rep_tx]; cl_rep <- sum(width(cds_rep))
len_tab <- list(); fr_tab <- list()
for (ch in names(sl)) {
  ga <- readGAlignments(bf, param = ScanBamParam(flag = flag, which = GRanges(ch, IRanges(1, sl[[ch]]))))
  L <- qwidth(ga); len_tab[[ch]] <- data.table(chr = ch, len = L)[, .N, by = .(chr, len)]
  if (ch %in% c("1", "2", "3", "4", "5")) {
    fe <- resize(granges(ga), width = 1, fix = "start")                     # 5' end, strand-aware
    m <- mapToTranscripts(fe, cds_rep, ignore.strand = FALSE)              # position within the CDS of a one-per-gene transcript
    m <- m[!duplicated(mcols(m)$xHits)]
    p <- start(m); cl <- cl_rep[as.character(seqnames(m))]; keep <- p > TRIM & p <= cl - TRIM
    fr_tab[[ch]] <- data.table(len = L[mcols(m)$xHits][keep], f5 = (p[keep] - 1L) %% 3L)[, .N, by = .(len, f5)]
  }
  .add(sprintf("    chr %s: %s primary alignments", ch, big(length(ga)))); rm(ga); invisible(gc())
}
len_tab <- rbindlist(len_tab)[, .(N = sum(N)), by = .(chr, len)]; fr_tab <- rbindlist(fr_tab)[, .(N = sum(N)), by = .(len, f5)]
nuc <- len_tab[chr %in% c("1", "2", "3", "4", "5"), .(N = sum(N)), by = len][order(len)]; nuc[, share := N / sum(N)]
n_nuc <- sum(nuc$N)
.add(sprintf("  nuclear primary alignments in the subsample: %s; 24-28 nt = %s; 28 nt = %s; 27 nt = %s; 29 nt = %s",
             big(n_nuc), pct(nuc[len %in% 24:28, sum(share)]), pct(nuc[len == 28, share]), pct(nuc[len == 27, share]), pct(nuc[len == 29, share])))
# frame per length: codon position of the 5' end; P-site frame after the study's offsets for 24-28 nt
fw <- dcast(fr_tab, len ~ f5, value.var = "N", fill = 0); setnames(fw, c("len", "pos1", "pos2", "pos3"))
fw[, tot := pos1 + pos2 + pos3]; fw[, `:=`(s1 = pos1 / tot, s2 = pos2 / tot, s3 = pos3 / tot)]
fw[, max_share := pmax(s1, s2, s3)]; fw[, offset_used := NA_integer_]; fw[as.character(len) %in% names(OFFSETS), offset_used := OFFSETS[as.character(len)]]
# a 5' end at codon position k (0-based) with offset o puts the P-site at codon position (k + o) mod 3; in frame = that equals 0
fw[, inframe_after_offset := NA_real_]
for (i in which(!is.na(fw$offset_used))) { o <- fw$offset_used[i]; k0 <- (3L - (o %% 3L)) %% 3L   # the 5'-end position that becomes frame 0
  fw$inframe_after_offset[i] <- c(fw$s1[i], fw$s2[i], fw$s3[i])[k0 + 1L] }
fw[, share_of_nuclear_reads := nuc$share[match(len, nuc$len)]]
.add("  5'-end codon position per read length (CDS-internal 5' ends):")
for (i in seq_len(nrow(fw))) if (fw$len[i] %in% LEN_SHOW) .add(sprintf("    %2d nt: n=%9s  pos1 %5.1f%%  pos2 %5.1f%%  pos3 %5.1f%%  max %5.1f%%  %s", fw$len[i], big(fw$tot[i]), 100 * fw$s1[i], 100 * fw$s2[i], 100 * fw$s3[i], 100 * fw$max_share[i],
  if (is.na(fw$offset_used[i])) "" else sprintf("offset %d -> in frame %.1f%%", fw$offset_used[i], 100 * fw$inframe_after_offset[i])))
used <- fw[!is.na(offset_used)]
.add(sprintf("  24-28 nt pooled: %.1f%% of CDS-internal P-sites in frame after the offsets (per length %s)", 100 * sum(used$inframe_after_offset * used$tot) / sum(used$tot),
             paste(sprintf("%d nt %.0f%%", used$len, 100 * used$inframe_after_offset), collapse = ", ")))
# cross-check against the Ribo-seQC run on the full merged BAM
if (file.exists(PATH_RQC)) { e <- new.env(); load(PATH_RQC, envir = e); rld <- e$res_all$read_stats$rld["nucl", ]
  full <- data.table(len = as.integer(sub("reads_", "", names(rld))), N_full = as.numeric(rld)); full[, share_full := N_full / sum(N_full)]
  cmp <- merge(nuc, full, by = "len"); .add(sprintf("  cross-check vs Ribo-seQC on the full merged BAM (%s nuclear reads): r = %.4f between the length distributions; 28 nt = %s (full) vs %s (subsample)",
    big(sum(full$N_full)), cor(cmp$share, cmp$share_full), pct(full[len == 28, share_full]), pct(nuc[len == 28, share])))
  fw[, share_of_nuclear_reads_full := full$share_full[match(len, full$len)]] }

################################################################################
# 3) P-site metagenes from the paper's P-site table (transcript coordinates)
################################################################################
.add("A5: transcript tiling + P-site join")
tile_tx <- function(grl) { ul <- unlist(grl, use.names = FALSE)
  d <- data.table(tx = rep(names(grl), elementNROWS(grl)), chr = as.character(seqnames(ul)), strand = as.character(strand(ul)), s = start(ul), e = end(ul))
  d <- d[, .(pos = s:e), by = .(tx, chr, strand, s, e)][, .(tx, chr, strand, pos)]
  d[, key := ifelse(strand == "+", pos, -pos)]; setorder(d, tx, key); d[, txpos := seq_len(.N), by = tx]; d[, key := NULL]; d }
P <- tile_tx(ex[rep_tx])
ps <- fread(PATH_PS, header = FALSE, col.names = c("count", "chr", "pos", "strand")); ps[, chr := as.character(chr)]
setkey(ps, chr, strand, pos); setkey(P, chr, strand, pos); P[, n := 0]; P[ps, n := i.count]
P[info, on = "tx", `:=`(cds_s = i.cds_s, cds_e = i.cds_e)]
# sanity: the genomic position at the CDS start in transcript coordinates must be the annotated CDS start
chk <- P[txpos == cds_s, .(tx, chr, strand, pos)]; g1 <- ifelse(as.character(unlist(runValue(strand(cds_rep)))) == "+", min(start(cds_rep)), max(end(cds_rep))); names(g1) <- names(cds_rep)
.add(sprintf("  check: transcript-coordinate CDS start = annotated CDS start in %.2f%% of %d transcripts", 100 * mean(chk$pos == g1[chk$tx]), nrow(chk)))
.add(sprintf("  transcript positions: %s; P-sites on these transcripts: %s of %s in the table (%s)", big(nrow(P)), big(sum(P$n)), big(sum(ps$count)), pct(sum(P$n) / sum(ps$count))))
P[, rel_s := txpos - cds_s]; P[, rel_e := txpos - (cds_e - 2L)]                      # 0 = first nt of the start / stop codon
# Termination rule (Larry, 2026-10-29; Wu et al. 2024 Plant Cell, Fig. 9B-E and Discussion): the footprint of a terminating ribosome is 1 nt shorter at
# the 5' end and 1 nt longer at the 3' end than that of an elongating ribosome, so the fixed 5'-end offset puts its P-site on the second nucleotide of
# the last sense codon (frame +1, rel_e == -2). Those P-sites are counted as in frame (fr = 0) in the in-frame percentages of C and D and in the per-ORF
# fractions (E); the metagene bars themselves keep the plain frame colours (the -2 bar is blue). The plain frame is kept alongside as *_raw for comparison.
term_frame <- function(rel_s, rel_e) ifelse(rel_e == -2L, 0L, as.integer(rel_s %% 3L))
P[, fr := term_frame(rel_s, rel_e)]; P[, `:=`(o_s = cds_s, o_e = cds_e, oid = tx)]
mg_of <- function(D, set_name) {   # metagene by position and plain frame (rel %% 3): the bars keep the frame colours, so the termination peak at -2 stays blue (Larry, 2026-10-29); the termination rule applies to the in-frame percentages only
  a <- D[rel_s >= -FLANK & rel_s <= INSIDE & txpos <= o_e]; b <- D[rel_e >= -INSIDE & rel_e <= FLANK & txpos >= o_s]
  mg <- rbind(merge(a[, .(n = sum(n)), by = .(rel = rel_s, frame = rel_s %% 3L)], a[, .(n_orf = uniqueN(oid)), by = .(rel = rel_s)], by = "rel")[, end := "start codon"],
              merge(b[, .(n = sum(n)), by = .(rel = rel_e, frame = rel_e %% 3L)], b[, .(n_orf = uniqueN(oid)), by = .(rel = rel_e)], by = "rel")[, end := "stop codon"])
  mg[, frame := factor(frame, levels = 0:2)]; mg[, set := set_name]; setorder(mg, end, rel, frame); mg[] }
mg_cds <- mg_of(P, "mORF")
in_cds <- P[txpos >= cds_s & txpos <= cds_e - 3L]                                    # mORF without the stop codon (as the uIF denominator)
orf_m <- in_cds[, .(tot = sum(n), f0 = sum(n[fr == 0]), f1 = sum(n[fr == 1]), f2 = sum(n[fr == 2]), f0_raw = sum(n[rel_s %% 3 == 0]), term = sum(n[rel_e == -2L])), by = tx][, `:=`(set = "mORF", id = tx)]
pool_m <- in_cds[, .(all = sum(n), f0 = sum(n[fr == 0]), f0_raw = sum(n[rel_s %% 3 == 0]), term = sum(n[rel_e == -2L]))]
pool_m_int <- in_cds[rel_s >= TRIM & txpos <= cds_e - TRIM, .(all = sum(n), f0 = sum(n[fr == 0]))]
.add(sprintf("  mORF P-sites in frame (stop codon excluded; termination rule): %s of %s (%s; without the rule %s; P-sites at the second nucleotide of the last sense codon = %s of all); CDS-internal (first/last %d nt excluded): %s", big(pool_m$f0), big(pool_m$all), pct(pool_m$f0 / pool_m$all), pct(pool_m$f0_raw / pool_m$all), pct(pool_m$term / pool_m$all), TRIM, pct(pool_m_int$f0 / pool_m_int$all)))

# TuORFs: start / last nucleotide in transcript coordinates
ex_u <- ex[unique(S1A$tx_id)]
gs <- GRanges(S1A$chr, IRanges(ifelse(S1A$Strand == "+", S1A$g_min, S1A$g_max), width = 1), strand = S1A$Strand)
ge <- GRanges(S1A$chr, IRanges(ifelse(S1A$Strand == "+", S1A$g_max, S1A$g_min), width = 1), strand = S1A$Strand)
map1 <- function(g) { m <- mapToTranscripts(g, ex_u, ignore.strand = FALSE); m <- m[names(ex_u)[mcols(m)$transcriptsHits] == S1A$tx_id[mcols(m)$xHits]]
  out <- rep(NA_integer_, length(g)); out[mcols(m)$xHits] <- start(m); out }
S1A$u_s <- map1(gs); S1A$u_e <- map1(ge); S1A$L_tx <- S1A$u_e - S1A$u_s + 1L
.add(sprintf("  TuORFs mapped to transcript coordinates: %d of %d; transcript length = uORF_length_nt + 3 in %.2f%%", sum(!is.na(S1A$L_tx)), nrow(S1A), 100 * mean(S1A$L_tx == S1A$uORF_length_nt + 3L, na.rm = TRUE)))
U <- as.data.table(S1A[!is.na(S1A$L_tx), c("uid", "gene_id", "tx_id", "u_s", "u_e", "L_tx", "uORF_length_nt", "uIF", "pct_uORF_plus", "detection_method")])
setnames(U, "tx_id", "tx")
U_D <- U[uORF_length_nt > 3L]   # Figure S1 D: 1-aa uORFs left out (start codon = last sense codon)
.add(sprintf("  Figure S1 D: %d 1-aa uORFs left out (start codon = last sense codon); %d TuORFs of >= 2 aa in the metagene", nrow(U) - nrow(U_D), nrow(U_D)))
W <- U[, .(txpos = (u_s - FLANK):(u_e + FLANK)), by = .(uid, tx, u_s, u_e, L_tx)]
setkey(P, tx, txpos); W <- merge(W, P[, .(tx, txpos, n)], by = c("tx", "txpos"), all.x = FALSE)
W[, rel_s := txpos - u_s]; W[, rel_e := txpos - (u_e - 2L)]
W[, fr := term_frame(rel_s, rel_e)]; W[, `:=`(o_s = u_s, o_e = u_e, oid = uid)]   # termination rule as for the mORFs
WD <- W[uid %in% U_D$uid]
mg_tu <- mg_of(WD, "TuORF")
# which positions does the uIF count cover? test "whole ORF" and "ORF without the stop codon" against Table S1A uORF_P_sites
tot_all  <- W[txpos >= u_s & txpos <= u_e,      .(tot = sum(n)), by = uid]
tot_nost <- W[txpos >= u_s & txpos <= u_e - 3L, .(tot = sum(n)), by = uid]
s1 <- as.data.table(S1A)[, .(uid, uORF_P_sites)]
m_all <- merge(tot_all, s1, by = "uid"); m_nost <- merge(tot_nost, s1, by = "uid")
.add(sprintf("  check vs Table S1A uORF_P_sites: whole ORF identical in %.1f%%, ORF without the stop codon identical in %.1f%% of %d TuORFs", 100 * mean(m_all$tot == m_all$uORF_P_sites), 100 * mean(m_nost$tot == m_nost$uORF_P_sites), nrow(m_all)))
EXCL_STOP <- mean(m_nost$tot == m_nost$uORF_P_sites) > mean(m_all$tot == m_all$uORF_P_sites)
.add(sprintf("  -> per-ORF frame counts use the %s (the uIF definition)", if (EXCL_STOP) "ORF without the stop codon" else "whole ORF"))
in_u <- if (EXCL_STOP) W[txpos >= u_s & txpos <= u_e - 3L] else W[txpos >= u_s & txpos <= u_e]
orf_u <- in_u[, .(tot = sum(n), f0 = sum(n[fr == 0]), f1 = sum(n[fr == 1]), f2 = sum(n[fr == 2]), f0_raw = sum(n[rel_s %% 3 == 0]), term = sum(n[rel_e == -2L])), by = .(uid, tx)][, `:=`(set = "TuORF", id = uid)]
pool_u_all <- in_u[, .(all = sum(n), f0 = sum(n[fr == 0]), f0_raw = sum(n[rel_s %% 3 == 0]))]
pool_u <- in_u[uid %in% U_D$uid, .(all = sum(n), f0 = sum(n[fr == 0]), f0_raw = sum(n[rel_s %% 3 == 0]), term = sum(n[rel_e == -2L]))]   # the panel D set
.add(sprintf("  TuORF P-sites in frame (termination rule): %s of %s (%s; without the rule %s) over %d TuORFs (%d TuORFs of Table S1A without transcript coordinates omitted)", big(pool_u_all$f0), big(pool_u_all$all), pct(pool_u_all$f0 / pool_u_all$all), pct(pool_u_all$f0_raw / pool_u_all$all), nrow(U), nrow(S1A) - nrow(U)))
.add(sprintf("  TuORF P-sites in frame, TuORFs of >= 2 aa (panel D; termination rule): %s of %s (%s; without the rule %s; P-sites at the second nucleotide of the last sense codon = %s = %s of all) over %d TuORFs", big(pool_u$f0), big(pool_u$all), pct(pool_u$f0 / pool_u$all), pct(pool_u$f0_raw / pool_u$all), big(pool_u$term), pct(pool_u$term / pool_u$all), nrow(U_D)))
# per-ORF in-frame fraction
orf <- rbind(orf_m[, .(set, id, tx, tot, f0, f1, f2, f0_raw, term)], orf_u[, .(set, id, tx, tot, f0, f1, f2, f0_raw, term)]); orf[, inframe := f0 / tot]; orf[, inframe_raw := f0_raw / tot]   # inframe = termination rule; inframe_raw = plain frame
orf[tot == 0, c("inframe", "inframe_raw") := NA_real_]   # V41 (2026-11-12): 0 / 0 is not defined; NA gives an empty cell in the table (NaN was shown as #NUM! in Excel); the csv file is unchanged
ok <- orf[tot >= MIN_PS]
one_aa <- U[uORF_length_nt == 3L, uid]
ok_E <- ok[!(set == "TuORF" & id %in% one_aa)]   # Figure S1 E: 1-aa uORFs left out as in D (Larry, 2026-10-29)
orf[, in_figure_S1E := tot >= MIN_PS & !(set == "TuORF" & id %in% one_aa)]
for (s in c("mORF", "TuORF")) { z <- ok_E[set == s]; .add(sprintf("  per-ORF in-frame fraction, %s with >= %d P-sites (n = %s): median %.2f, IQR %.2f-%.2f; >= 0.5 in %s; >= 0.6 in %s (without the termination rule: median %.2f, >= 0.5 in %s)", s, MIN_PS, big(nrow(z)), median(z$inframe), quantile(z$inframe, .25), quantile(z$inframe, .75), pct(mean(z$inframe >= 0.5)), pct(mean(z$inframe >= 0.6)), median(z$inframe_raw), pct(mean(z$inframe_raw >= 0.5)))) }
z1 <- ok[set == "TuORF" & id %in% U[uORF_length_nt == 3L, uid]]; z2 <- ok[set == "TuORF" & !(id %in% U[uORF_length_nt == 3L, uid])]   # record the effect of the 1-aa uORFs on panel E
.add(sprintf("  panel E with and without 1-aa uORFs: TuORFs of >= 2 aa with >= %d P-sites n = %s, median in-frame %.2f, >= 0.5 in %s; 1-aa uORFs n = %s, median %.2f, >= 0.5 in %s", MIN_PS, big(nrow(z2)), median(z2$inframe), pct(mean(z2$inframe >= 0.5)), big(nrow(z1)), median(z1$inframe), pct(mean(z1$inframe >= 0.5))))
ku <- ok_E[set == "TuORF"]; ku <- merge(ku, U[, .(uid, uIF, pct_uORF_plus, uORF_length_nt, detection_method)], by.x = "id", by.y = "uid")
.add(sprintf("  TuORF in-frame fraction vs uIF: Spearman rho = %.2f; vs %%uORF+: %.2f; vs length: %.2f (n = %d; without the termination rule %.2f / %.2f / %.2f)", cor(ku$inframe, ku$uIF, method = "spearman", use = "complete.obs"), cor(ku$inframe, ku$pct_uORF_plus, method = "spearman", use = "complete.obs"), cor(ku$inframe, ku$uORF_length_nt, method = "spearman"), nrow(ku), cor(ku$inframe_raw, ku$uIF, method = "spearman", use = "complete.obs"), cor(ku$inframe_raw, ku$pct_uORF_plus, method = "spearman", use = "complete.obs"), cor(ku$inframe_raw, ku$uORF_length_nt, method = "spearman")))
for (dm in sort(unique(ku$detection_method))) { z <- ku[detection_method == dm]; .add(sprintf("    %s: n = %d, median in-frame %.2f", dm, nrow(z), median(z$inframe))) }
# start/stop peak sizes
peak <- function(mg) { s <- mg[end == "start codon"]; e <- mg[end == "stop codon"]; term <- e[rel >= -3 & rel <= 2, .(n = sum(n)), by = rel]; tp <- term[which.max(n)]
  c(start0 = s[rel == 0, sum(n)], body_f0 = mean(s[rel >= 3 & rel <= INSIDE & rel %% 3 == 0 & frame == "0", n]), body_f12 = mean(s[rel >= 3 & rel <= INSIDE & rel %% 3 != 0 & frame != "0", n]),   # in-frame positions vs the other positions; termination-rule P-sites (frame 0 at a +1 position) in neither
    term_pos = tp$rel, term_n = tp$n, body_end_f0 = mean(e[rel >= -INSIDE & rel <= -6 & rel %% 3 == 0 & frame == "0", n]), stop_codon = e[rel >= 0 & rel <= 2, sum(n)], before = e[rel >= -60 & rel <= -1, sum(n)]) }
for (nm in c("mORF", "TuORF")) { pk <- peak(if (nm == "mORF") mg_cds else mg_tu)
  .add(sprintf("  %s metagene: start-codon P-sites = %.1fx the mean in-frame position of the next %d nt (frame 0 : other frames inside the ORF = %.1f : 1); termination peak at position %+d relative to the first nucleotide of the stop codon (%.1fx the mean in-frame position of the preceding %d nt); the three stop-codon positions carry %.1f%% of the P-sites of the preceding 60 nt",
               nm, pk["start0"] / pk["body_f0"], INSIDE, pk["body_f0"] / pk["body_f12"], as.integer(pk["term_pos"]), pk["term_n"] / pk["body_end_f0"], INSIDE, 100 * pk["stop_codon"] / pk["before"])) }

################################################################################
# 4) Tables
################################################################################
fw_out <- fw[order(len), .(read_length = len, n_CDS_internal_5p_ends = tot, share_5p_end_codon_pos1 = round(s1, 4), share_codon_pos2 = round(s2, 4), share_codon_pos3 = round(s3, 4), max_share = round(max_share, 4),
                           offset_used_in_study = offset_used, inframe_after_offset = round(inframe_after_offset, 4), share_of_nuclear_reads_subsample = round(share_of_nuclear_reads, 5),
                           share_of_nuclear_reads_full = if ("share_of_nuclear_reads_full" %in% names(fw)) round(share_of_nuclear_reads_full, 5) else NA_real_)]
fwrite(fw_out, file.path(OUT_DIR, "A5_Table_read_length_frame.csv"))
fwrite(rbind(mg_cds, mg_tu)[order(set, end, rel), .(set, end, rel_position = rel, frame, P_sites = n, n_ORFs = n_orf)], file.path(OUT_DIR, "A5_Table_metagenes.csv"))
fwrite(orf[order(set, -tot), .(set, id, tx, P_sites = tot, frame0 = f0, frame1 = f1, frame2 = f2, frame1_at_last_sense_codon_counted_in_frame0 = term, inframe = round(inframe, 4), inframe_without_termination_rule = round(inframe_raw, 4), in_figure_S1E)], file.path(OUT_DIR, "A5_Table_inframe_per_ORF.csv"))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S2 — Ribo-seq quality metrics behind Figure S1.",
  "Libraries: the six Col-0 seedling Ribo-seq libraries of Wu et al. 2024 (GSE183264), merged; primary alignments.",
  sprintf("read_length_frame: random 10%% subsample of the alignments (%s nuclear primary alignments). Codon position of the read 5' end for 5' ends inside the coding sequence (CDS) of the TuORF source isoforms and one isoform per No-TuORF control gene (%d transcripts; first/last %d nt of each CDS excluded). offset_used_in_study / inframe_after_offset: the P-site rule of this study (24-28 nt, offsets 8-12 nt) and the share of P-sites that fall in frame 0 with it. share_of_nuclear_reads_full: Ribo-seQC on all alignments.", big(n_nuc), length(rep_tx), TRIM),
  sprintf("metagenes: P-sites of the study's P-site table (24-28 nt reads, offsets 8-12 nt, all six libraries) summed per position relative to the first nucleotide of the start or stop codon; frame relative to the ORF (plain frame, as plotted); the P-sites at the second nucleotide of the last sense codon (frame +1, the blue bar at -2: the footprints of terminating ribosomes are shifted 1 nt towards the 3' end, Wu et al. 2024) are counted as frame 0 in the in-frame percentages of Figure S1C-D and in inframe_per_ORF; mORF = %d one-per-gene transcripts, TuORF = %d TuORFs of Table S1A of at least 2 aa (1-aa uORFs left out because their start codon is also the last sense codon). Positions are inside the ORF plus a %d-nt flank; n_ORFs = ORFs contributing to a position. inframe_per_ORF lists all TuORFs (in_figure_S1E marks those in panel E).", length(rep_tx), nrow(U_D), FLANK),
  sprintf("inframe_per_ORF: P-sites by frame inside each mORF and each TuORF (stop codon excluded); frame0 includes the P-sites at the second nucleotide of the last sense codon (also listed in frame1_at_last_sense_codon_counted_in_frame0); inframe = frame0 / all; inframe_without_termination_rule = plain frame 0 / all (Figure S1E uses inframe for ORFs with >= %d P-sites, TuORFs of at least 2 aa: in_figure_S1E). Both fractions are NA (empty cell) for the %d ORFs without P-sites.", MIN_PS, sum(orf$tot == 0))
)), colNames = FALSE)
addWorksheet(wb, "read_length_frame"); writeData(wb, "read_length_frame", fw_out)
addWorksheet(wb, "metagenes"); writeData(wb, "metagenes", rbind(mg_cds, mg_tu)[order(set, end, rel), .(set, end, rel_position = rel, frame, P_sites = n, n_ORFs = n_orf)])
addWorksheet(wb, "inframe_per_ORF"); writeData(wb, "inframe_per_ORF", orf[order(set, -tot), .(set, id, tx, P_sites = tot, frame0 = f0, frame1 = f1, frame2 = f2, frame1_at_last_sense_codon_counted_in_frame0 = term, inframe = round(inframe, 4), inframe_without_termination_rule = round(inframe_raw, 4), in_figure_S1E)])
saveWorkbook(wb, file.path(OUT_DIR, "A5_Table_riboseq_QC.xlsx"), overwrite = TRUE)

################################################################################
# 5) Figure
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "top", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.5))
frame_cols <- c(`0` = "#D7191C", `1` = "#2C7BB6", `2` = "#1A9641")
# A
da <- nuc[len %in% LEN_SHOW]; da[, used := ifelse(len %in% 24:28, "used (24–28 nt)", "not used")]
pA <- ggplot(da, aes(factor(len), 100 * share, fill = used)) + geom_col(width = 0.8) + scale_fill_manual(values = c(`used (24–28 nt)` = "grey25", `not used` = "grey70"), name = NULL) +
  labs(x = "Read length (nt)", y = "% of nuclear reads", title = "A  Read-length distribution", subtitle = sprintf("10%% subsample of the six merged libraries (n = %.1f million nuclear primary alignments; 24–28 nt = %s)", n_nuc / 1e6, pct(nuc[len %in% 24:28, sum(share)]))) + th
# B
db <- melt(fw[len %in% LEN_SHOW, .(len, `1` = s1, `2` = s2, `3` = s3)], id.vars = "len", variable.name = "cpos", value.name = "share")
lab <- fw[len %in% LEN_SHOW][, lbl := ifelse(is.na(offset_used), sprintf("%.0f%%", 100 * max_share), sprintf("+%d\n%.0f%%", offset_used, 100 * inframe_after_offset))]
pB <- ggplot(db, aes(factor(len), 100 * share, fill = cpos)) + geom_col(width = 0.8) + scale_fill_manual(values = c(`1` = "#D7191C", `2` = "#2C7BB6", `3` = "#1A9641"), name = "codon position of the read 5′ end") +
  geom_text(data = lab, aes(factor(len), 103, label = gsub("\n", " ", lbl)), inherit.aes = FALSE, size = 2.1, angle = 90, hjust = 0, vjust = 0.5) + scale_y_continuous(limits = c(0, 150), breaks = seq(0, 100, 25)) +
  annotate("rect", xmin = which(LEN_SHOW == 24) - 0.5, xmax = which(LEN_SHOW == 28) + 0.5, ymin = -1, ymax = 101, fill = NA, colour = "black", linewidth = 0.4) +
  labs(x = "Read length (nt)", y = "% of CDS-internal 5′ ends", title = "B  Reading frame by read length", subtitle = "labels: offset used here (nt) and the resulting in-frame percentage; other lengths: percentage of the dominant position") + th
# C / D
mg_plot <- function(mg, ttl, sub) { mg <- copy(mg); mg[, end := factor(end, levels = c("start codon", "stop codon"))]
  ggplot(mg, aes(rel, n / 1e3, fill = frame)) + geom_col(width = 0.9) + facet_wrap(~end, scales = "free_x") + scale_fill_manual(values = frame_cols, name = "frame relative to the ORF", labels = c("0 (in frame)", "1", "2")) +
    geom_vline(xintercept = -0.5, linetype = 3, linewidth = 0.3) + labs(x = "Position of the P-site relative to the first nucleotide of the codon (nt)", y = "P-sites (thousands)", title = ttl, subtitle = sub) + th + theme(panel.grid.major.x = element_line()) }
pC <- mg_plot(mg_cds, "C  P-site metagene at mORF start and stop codons", sprintf("%s transcripts (TuORF source isoforms and one isoform per No-TuORF gene); P-sites of this study (24–28-nt reads, offsets 8–12 nt); the blue bar at position −2 of the stop-codon panel (frame +1 of the last sense codon) comes from terminating ribosomes, whose footprints are shifted by 1 nt (Wu et al. 2024), and is counted as in frame in the percentage: %s of mORF P-sites in frame", big(length(rep_tx)), pct(pool_m$f0 / pool_m$all)))
pD <- mg_plot(mg_tu, "D  P-site metagene at TuORF start and stop codons", sprintf("%s TuORFs of at least 2 aa (Table S1A; 1-aa uORFs left out: their start codon is also the last sense codon); positions inside the TuORF plus %d-nt flanks (the flank after the stop codon includes mORF starts); the blue bar at −2 is counted as in frame as in C; %s of their P-sites in frame", big(nrow(U_D)), FLANK, pct(pool_u$f0 / pool_u$all)))
# E
de <- ok_E[, .(set = factor(set, levels = c("TuORF", "mORF")), inframe)]
medE <- de[, .(m = median(inframe), n = .N), by = set]
pE <- ggplot(de, aes(inframe, colour = set)) + stat_ecdf(linewidth = 0.7) + geom_vline(xintercept = 1/3, linetype = 3, linewidth = 0.3) + scale_colour_manual(values = c(TuORF = "#D7191C", mORF = "grey30"), name = NULL, labels = setNames(sprintf("%s (n = %s; median %.2f)", medE$set, big(medE$n), medE$m), as.character(medE$set))) +   # named labels: each label follows its key whatever the row order of medE
  labs(x = "In-frame fraction of P-sites per ORF", y = "Cumulative fraction\nof ORFs", title = "E  Per-ORF in-frame fraction", subtitle = sprintf("ORFs with ≥ %d P-sites (stop codon excluded; the P-sites of the blue bar at −2 in C–D counted as in frame), TuORFs of at least 2 aa\ndotted line = 1/3 (no periodicity)", MIN_PS)) + th + theme(legend.direction = "vertical")
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
lay <- rbind(c(1, 2), c(3, 3), c(4, 4), c(5, 5))
# V30 (2026-11-09, Larry): the legends of A and E inside the plots (B, C and D keep theirs above the panel: bars and labels fill these panels); aligned plotting areas (fig_grid)
cS1 <- c(lapply(list(pA, pB), compact), lapply(list(pC, pD, pE), compact_wide)); cS1[[1]] <- legend_in(cS1[[1]], 0.03, 0.98); cS1[[5]] <- legend_in(cS1[[5]], 0.03, 0.95)
g <- fig_grid(cS1, design = "AB\nCC\nDD\nEE", heights = c(1.1, 1, 1, 0.85), top = fig_title("Figure S1. Ribo-seq quality metrics of the libraries used in this study"))
fig_save(file.path(OUT_DIR, "A5_Figure_riboseq_QC.pdf"), g, h = 9)
.add(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A5_summary.txt"))
