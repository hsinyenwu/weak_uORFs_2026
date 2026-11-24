################################################################################
# A21 — Are %uORF+ and the mORF percentage of Supplemental Fig. S11 comparable?  #
# (Larry, 2026-11-03: "the calculation for the uORF and the mORF is not          #
#  identical now; they were both done from 5'UTR upstream 100 to the end of CDS, #
#  but now they are not the same")                                               #
# Since the A1 window fix, %uORF+ counts CAGE signal in the TuORF isoform's      #
#   5' leader + 100 nt upstream (the "leader window"); the mORF percentage of    #
#   Supplemental Fig. S11 counts signal over the whole transcript + 100 nt (the  #
#   "full window"), because within the leader window every mORF start codon is  #
#   at 100% by construction. This script measures what the mismatch does:        #
#   1) reproduces S11A (one transcript per protein-coding gene, full window)     #
#      and counts the genes whose window has no CAGE signal at all (scored 0%);  #
#   2) for the single-TuORF cohort computes, on the same isoform, %uORF+ with    #
#      the leader window (manuscript), %uORF+ with the full window (submitted    #
#      definition) and the mORF percentage with the full window; reports how     #
#      often the leader-window %uORF+ exceeds the full-window mORF percentage    #
#      (impossible under one window), the share of signal inside the CDS/3' UTR, #
#      and the gap between the two start codons under the common full window.   #
# Inputs: Files/Araport11+CTRL_20181206.gtf, Files/wt_R123_genome_RiboPlotR.txt,  #
#   output/Table_S1_all_TuORFs.xlsx (S1A), output/Table_S2_single_TuORF_genes.xlsx#
# Outputs (output/): A21_mORF_window_check.csv (cohort, per gene),               #
#   A21_summary.txt. Results only; NOT part of the manuscript.                   #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A21_mORF_window_check.R
################################################################################
suppressPackageStartupMessages({ library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(data.table); library(openxlsx) })
options(warn = 1); t0 <- Sys.time(); .log <- character(0); .add <- function(...) { m <- paste0(...); cat(m, "\n"); .log <<- c(.log, m) }
DATA_DIR <- "Files"; OUT_DIR <- "output"; UP <- 100L
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf"); PATH_CAGE <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt")
stopifnot(file.exists(PATH_GTF), file.exists(PATH_CAGE))

txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else makeTxDbFromGFF(PATH_GTF, format = "gtf")
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE); cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
tx_all <- transcripts(txdb, use.names = TRUE, columns = c("tx_name", "gene_id"))
txinfo <- data.table(tx_id = tx_all$tx_name, gene_id = as.character(unlist(tx_all$gene_id)))
d <- fread(PATH_CAGE, header = FALSE, select = 1:4, col.names = c("count", "chr", "pos", "strand")); d[, chr := sub("^[Cc]hr", "", as.character(chr))]
cage <- GRanges(d$chr, IRanges(d$pos, width = 1), strand = d$strand, count = as.numeric(d$count)); cage <- cage[seqnames(cage) %in% seqlevels(txdb)]
.add(sprintf("A21: TxDb %d transcripts, %d with a 5' leader, %d with a CDS; CAGE positions %d; %.1f min", nrow(txinfo), length(utr5), length(cds), length(cage), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# windows (same definitions as .make_tss_windows in the canonical script)
win_full <- function(tx_ids) { tx <- tx_all[tx_all$tx_name %in% tx_ids]; names(tx) <- tx$tx_name; trim(suppressWarnings(promoters(tx, upstream = UP, downstream = width(tx)))) }
win_leader <- function(tx_ids) { u <- utr5[names(utr5) %in% tx_ids]; ul <- unlist(u, use.names = TRUE); k <- names(ul)
  dd <- data.table(i = seq_along(ul), k = k, pos5 = ifelse(as.character(strand(ul)) == "+", start(ul), -end(ul))); first <- dd[order(k, pos5), .SD[1], by = k]$i
  ul[first] <- resize(ul[first], width = width(ul[first]) + UP, fix = "end"); restrict(ul, start = 1L) }
# share of the window's CAGE signal at or upstream of a start codon (first nucleotide), per transcript; 0 when no signal (as in the pipeline)
share_up <- function(win, starts) {   # starts: data.table(tx_id, s0, Strand)
  ov <- findOverlaps(cage, win); dt <- data.table(tx_id = names(win)[subjectHits(ov)], pos = start(cage)[queryHits(ov)], count = cage$count[queryHits(ov)])
  dt <- merge(dt, starts, by = "tx_id")
  r <- dt[, .(total = sum(count), up = sum(count[if (Strand[1] == "+") pos <= s0[1] else pos >= s0[1]])), by = tx_id]
  r <- merge(starts[, .(tx_id)], r, by = "tx_id", all.x = TRUE); r[is.na(total), `:=`(total = 0, up = 0)]; r[, pct := ifelse(total > 0, 100 * up / total, 0)]; r }
cds_start <- function(tx_ids) { c1 <- cds[names(cds) %in% tx_ids]; dd <- as.data.table(c1)[, .(Strand = as.character(strand[1]), s0 = if (as.character(strand[1]) == "+") min(start) else max(end)), by = .(tx_id = group_name)]; dd }
grp <- function(p) cut(p, c(-Inf, 30, 80, Inf), labels = c("0~30", "30~80", "80~100"), right = FALSE)

# 1) Supplemental Fig. S11A as in the canonical script: one transcript per protein-coding gene (first tx_id with both a CDS and a 5' leader), full window
valid <- intersect(names(cds), names(utr5)); tmap <- txinfo[tx_id %in% valid][order(tx_id)][, .SD[1], by = gene_id]
st_m <- cds_start(tmap$tx_id); m_all <- share_up(win_full(tmap$tx_id), st_m)
.add(sprintf("\n== S11A (all protein-coding genes with a CDS and a 5' leader, one transcript each, full window): n = %d; groups 0~30 / 30~80 / 80~100 = %d / %d / %d (%.1f / %.1f / %.1f%%)",
             nrow(m_all), sum(grp(m_all$pct) == "0~30"), sum(grp(m_all$pct) == "30~80"), sum(grp(m_all$pct) == "80~100"), 100 * mean(grp(m_all$pct) == "0~30"), 100 * mean(grp(m_all$pct) == "30~80"), 100 * mean(grp(m_all$pct) == "80~100")))
.add(sprintf("   genes with NO CAGE signal in the window (scored 0%%): %d = %.1f%% of all, %d of the %d genes in the 0~30 group; genes with signal but 0%% upstream: %d; genes with signal < 10 tags: %d",
             sum(m_all$total == 0), 100 * mean(m_all$total == 0), sum(m_all$total == 0 & grp(m_all$pct) == "0~30"), sum(grp(m_all$pct) == "0~30"), sum(m_all$total > 0 & m_all$pct == 0), sum(m_all$total > 0 & m_all$total < 10)))
m_sig <- m_all[total > 0]
.add(sprintf("   genes WITH signal (n = %d): groups 0~30 / 30~80 / 80~100 = %.1f / %.1f / %.1f%%; median share of signal downstream of the mORF start (inside CDS / 3' UTR) = %.1f%%; > 20%% in %.1f%% of genes; at 100%% upstream: %.1f%%",
             nrow(m_sig), 100 * mean(grp(m_sig$pct) == "0~30"), 100 * mean(grp(m_sig$pct) == "30~80"), 100 * mean(grp(m_sig$pct) == "80~100"), median(100 - m_sig$pct), 100 * mean(100 - m_sig$pct > 20), 100 * mean(m_sig$pct >= 99.95)))

# 2) single-TuORF cohort: the same isoform under both windows
S1A <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1A (CAGE)"))
S2 <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx"), sheet = 1))
U <- unique(S1A[gene_id %in% S2$gene_id, .(gene_id, tx_id, Strand, s0 = as.integer(uORF_start_genomic), pUP_table = as.numeric(pct_uORF_plus))])
U <- U[tx_id %in% names(utr5) & tx_id %in% names(cds)]; U <- U[!duplicated(gene_id)]
.add(sprintf("\n== single-TuORF cohort: %d genes in Table S2, %d with their TuORF isoform in the TxDb (one TuORF each)", uniqueN(S2$gene_id), nrow(U)))
st_u <- U[, .(tx_id, s0, Strand)]; st_c <- cds_start(U$tx_id)
u_leader <- share_up(win_leader(U$tx_id), st_u); u_full <- share_up(win_full(U$tx_id), st_u); m_full <- share_up(win_full(U$tx_id), st_c)
C <- merge(U, setnames(u_leader[, .(tx_id, total, pct)], c("tx_id", "signal_leader", "uORF_leader")), by = "tx_id")
C <- merge(C, setnames(u_full[, .(tx_id, total, pct)], c("tx_id", "signal_full", "uORF_full")), by = "tx_id")
C <- merge(C, setnames(m_full[, .(tx_id, pct)], c("tx_id", "mORF_full")), by = "tx_id")
C[, `:=`(cds_internal_share = 100 - mORF_full, gap_full = mORF_full - uORF_full, leader_exceeds_mORF = uORF_leader > mORF_full + 1e-6)]
.add(sprintf("   check against Table S1A: max |leader-window %%uORF+ - table| = %.3f points (%d differ by > 0.5)", max(abs(C$uORF_leader - C$pUP_table)), sum(abs(C$uORF_leader - C$pUP_table) > 0.5)))
.add(sprintf("   signal inside the CDS / 3' UTR (full window minus leader window): median %.1f%% of a gene's full-window signal; > 20%% in %d genes (%.1f%%); > 50%% in %d (%.1f%%)",
             median(C$cds_internal_share), sum(C$cds_internal_share > 20), 100 * mean(C$cds_internal_share > 20), sum(C$cds_internal_share > 50), 100 * mean(C$cds_internal_share > 50)))
.add(sprintf("   %%uORF+ leader window (manuscript) vs full window (submitted): median %.1f vs %.1f; groups 0~30 / 30~80 / 80~100: %.1f / %.1f / %.1f%% vs %.1f / %.1f / %.1f%%; %d genes (%.1f%%) change group",
             median(C$uORF_leader), median(C$uORF_full), 100 * mean(grp(C$uORF_leader) == "0~30"), 100 * mean(grp(C$uORF_leader) == "30~80"), 100 * mean(grp(C$uORF_leader) == "80~100"),
             100 * mean(grp(C$uORF_full) == "0~30"), 100 * mean(grp(C$uORF_full) == "30~80"), 100 * mean(grp(C$uORF_full) == "80~100"), sum(grp(C$uORF_leader) != grp(C$uORF_full)), 100 * mean(grp(C$uORF_leader) != grp(C$uORF_full))))
.add(sprintf("   mORF percentage (full window) in the cohort: median %.1f%%; groups 0~30 / 30~80 / 80~100: %.1f / %.1f / %.1f%%", median(C$mORF_full), 100 * mean(grp(C$mORF_full) == "0~30"), 100 * mean(grp(C$mORF_full) == "30~80"), 100 * mean(grp(C$mORF_full) == "80~100")))
.add(sprintf("   MISMATCH: leader-window %%uORF+ EXCEEDS the full-window mORF percentage in %d of %d genes (%.1f%%) (impossible under one window); by > 10 points in %d, > 20 points in %d; %d genes at 100%% uORF+ have an mORF percentage below 80%%",
             sum(C$leader_exceeds_mORF), nrow(C), 100 * mean(C$leader_exceeds_mORF), sum(C$uORF_leader - C$mORF_full > 10), sum(C$uORF_leader - C$mORF_full > 20), sum(C$uORF_leader >= 99.95 & C$mORF_full < 80)))
.add(sprintf("   under ONE window (full, both start codons): mORF %% >= uORF %% in all %d genes (check: %s); gap = signal between the two start codons: median %.1f points, > 20 points in %.1f%% of genes, 0 in %.1f%%",
             nrow(C), all(C$gap_full >= -1e-6), median(C$gap_full), 100 * mean(C$gap_full > 20), 100 * mean(abs(C$gap_full) < 1e-6)))
.add(sprintf("   within-gene Spearman: leader %%uORF+ vs full-window mORF %%: %.3f; full-window %%uORF+ vs mORF %%: %.3f", cor(C$uORF_leader, C$mORF_full, method = "spearman"), cor(C$uORF_full, C$mORF_full, method = "spearman")))
fwrite(C[order(gene_id)], file.path(OUT_DIR, "A21_mORF_window_check.csv"))
.add(sprintf("\nA21 done in %.1f min; per-gene table output/A21_mORF_window_check.csv", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A21_summary.txt"))
