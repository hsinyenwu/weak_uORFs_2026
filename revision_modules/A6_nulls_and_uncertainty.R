################################################################################
# A6 — Null models and uncertainty for the TSS effects (Reviewer 2, §2.4)       #
#                                                                               #
# (a) Fig S9: the distal-minus-proximal %uORF+ difference in multi-TuORF genes  #
#     against two position-matched baselines:                                   #
#       null 1 = the same two positions shifted together to every possible      #
#                place in the same leader (spacing kept; exact enumeration);     #
#       null 2 = positions at the same distances from the mORF start in         #
#                No-TuORF control leaders.                                       #
# (b) Uncertainty of %uORF+: Jeffreys 95% intervals (and a binomial bootstrap)  #
#     on effective CAGE tag counts (the CAGE file is a sum of three TPM-scaled   #
#     replicates; one tag = one of three unit values), the same on the integer  #
#     TSS-seq counts, agreement between the two sources, depth dependence of    #
#     the "<100%" and "<90%" calls, and depth-standardised estimates.           #
# (c) uORF length vs exclusion (Fig 5C-E): a random-placement positional null   #
#     (an ORF of the same length placed at every possible position of the same  #
#     leader; exclusion = %uORF+ <= 30), stratification by the distance from    #
#     the uORF start codon to the mORF, and logistic models with and without    #
#     that distance.                                                            #
# Outputs (output/A6_*): per-gene tables, model table, summary txt, figure pdf.  #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A6_nulls_and_uncertainty.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(gridExtra); library(grid)
})
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
set.seed(20261028)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF   <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_CAGE  <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt")
PATH_TSS   <- file.path(DATA_DIR, "wt_TSS_Nielsen_2019.txt")
PATH_S1    <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"); PATH_S2 <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_RSEM_RNA  <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_ISO   <- file.path(DATA_DIR, "RSEM_isoforms_TPM_6samples.tsv")
PATH_5UTR  <- file.path(DATA_DIR, "genes_w_5UTR.tsv")
UPSTREAM <- 100L; EXCL <- 30; N_CTRL <- 50L; B <- 2000L; STOP_OVERLAP <- 4L
for (f in c(PATH_GTF, PATH_CAGE, PATH_TSS, PATH_S1, PATH_S2, PATH_OUORF, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_ISO, PATH_5UTR)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x); big <- function(x) trimws(format(x, big.mark = ","))
med <- function(x) median(x[is.finite(x)])

################################################################################
# 1) Annotation, tables, gene sets, leader windows with inclusion fractions
################################################################################
.add("A6: TxDb")
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
ex <- exonsBy(txdb, by = "tx", use.names = TRUE); utr <- fiveUTRsByTranscript(txdb, use.names = TRUE)
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)") %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"), chr = as.character(chr))
S1B <- read.xlsx(PATH_S1, sheet = "Table S1B (TSS-seq)") %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"))
S2  <- read.xlsx(PATH_S2, sheet = 1) %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"), chr = as.character(chr))
ou  <- read.xlsx(PATH_OUORF)
rna <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); ribo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
TE_tbl <- data.frame(gene_id = rna$gene_id, RNA_TPM = rna$TPM, Ribo_TPM = ribo$TPM[match(rna$gene_id, ribo$gene_id)]) %>%
  mutate(TE = Ribo_TPM / RNA_TPM) %>% filter(gene_id %in% read.delim(PATH_5UTR)$gene_id, RNA_TPM > 1, is.finite(TE))
upstream_translated <- unique(c(S1A$gene_id, as.character(ou$gene_id)))
iso <- fread(PATH_ISO); iso[, mTPM := rowMeans(.SD), .SDcols = patterns("^TPM_")]
iso <- iso[gene_id %in% setdiff(TE_tbl$gene_id, upstream_translated) & transcript_id %in% names(utr)]
setorder(iso, gene_id, -mTPM); ctrl_tx <- iso[!duplicated(gene_id), .(gene_id, tx_id = transcript_id)]
tx_use <- unique(c(S1A$tx_id, S2$tx_id, ctrl_tx$tx_id)); tx_use <- tx_use[tx_use %in% names(utr)]
.add(sprintf("  transcripts with leader windows: %d (TuORF isoforms %d, No-TuORF controls %d)", length(tx_use), length(unique(S1A$tx_id)), nrow(ctrl_tx)))

# windows = annotated 5' UTR exons + 100 nt upstream (the %uORF+ window), tiled per nucleotide, ordered 5'->3' (exon space)
win <- utr[tx_use]; ul <- unlist(win, use.names = FALSE); txv <- rep(names(win), elementNROWS(win))
plus <- as.character(strand(ul)) == "+"; key <- ifelse(plus, start(ul), -end(ul)); o <- order(txv, key); idx <- o[!duplicated(txv[o])]
ul$tx <- txv; ul$annot <- TRUE
up_piece <- ul[idx]; up_piece <- flank(up_piece, width = UPSTREAM, start = TRUE); up_piece <- restrict(up_piece, start = 1L); up_piece$annot <- FALSE
tile <- function(gr) { d <- data.table(tx = gr$tx, chr = as.character(seqnames(gr)), strand = as.character(strand(gr)), annot = gr$annot, s = start(gr), e = end(gr)); d[, .(pos = s:e), by = .(tx, chr, strand, annot, s, e)][, .(tx, chr, strand, annot, pos)] }
P <- unique(rbind(tile(ul), tile(up_piece)))
P[, key := ifelse(strand == "+", pos, -pos)]; setorder(P, tx, key); P[, idx := seq_len(.N), by = tx]; P[, key := NULL]
# CAGE (TPM-scaled sum of 3 replicates) and TSS-seq (integer counts)
cg <- fread(PATH_CAGE, header = FALSE, col.names = c("count", "chr", "pos", "strand")); cg[, chr := gsub("^chr", "", chr, ignore.case = TRUE)]
ts <- fread(PATH_TSS,  header = FALSE, col.names = c("count", "chr", "pos", "strand")); ts[, chr := gsub("^chr", "", chr, ignore.case = TRUE)]
setkey(cg, chr, strand, pos); setkey(ts, chr, strand, pos); setkey(P, chr, strand, pos)
P[, cage := 0]; P[cg, cage := i.count]; P[, tss := 0]; P[ts, tss := i.count]
setorder(P, tx, idx)
P[, `:=`(cum_c = cumsum(cage), tot_c = sum(cage), cum_t = cumsum(tss), tot_t = sum(tss)), by = tx]
P[, f_c := ifelse(tot_c > 0, cum_c / tot_c, NA_real_)]; P[, f_t := ifelse(tot_t > 0, cum_t / tot_t, NA_real_)]
# effective CAGE tags: the file is a sum of three TPM-scaled replicates; the three smallest values are the single-tag units
units <- sort(unique(cg$count))[1:3]; u_mean <- mean(units)
.add(sprintf("  CAGE single-tag units: %s (TPM per tag in the three replicates); effective tags = signal / %.5f", paste(sprintf("%.5f", units), collapse = ", "), u_mean))
TXI <- P[, .(first_annot = min(idx[annot]), last_annot = max(idx[annot]), leader_len = sum(annot), idx_cds = max(idx[annot]) + 1L, tot_c = tot_c[1], tot_t = tot_t[1]), by = tx]
.add(sprintf("  windows tiled: %s positions; %d transcripts without CAGE in the window, %d without TSS-seq", big(nrow(P)), sum(TXI$tot_c == 0), sum(TXI$tot_t == 0)))
# TuORF starts -> exon-space index and inclusion fractions
start_idx <- function(d) { q <- data.table(tx = d$tx_id, chr = d$chr, strand = d$Strand, pos = d$uORF_start_genomic); setkey(P, tx, chr, strand, pos); P[q, .(idx, f_c, f_t, cum_c, tot_c, cum_t, tot_t)] }
s1 <- start_idx(S1A); S1A$idx <- s1$idx; S1A$f_c <- s1$f_c; S1A$f_t <- s1$f_t
.add(sprintf("  check: inclusion fraction at the TuORF start vs Table S1A %%uORF+ (CAGE): r = %.4f (n = %d); vs Table S1B (TSS-seq): r = %.4f", cor(100 * S1A$f_c, S1A$pct_uORF_plus, use = "complete.obs"), sum(!is.na(S1A$f_c)), cor(100 * S1A$f_t, S1B$pct_uORF_plus[match(S1A$uid, S1B$uid)], use = "complete.obs")))

################################################################################
# 2) (a) Figure S9 baseline: distal vs proximal %uORF+ in multi-TuORF genes
################################################################################
.add("A6a: Fig S3 position-matched baselines")
R <- S1A %>% mutate(uIF = ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, NA_real_)) %>% filter(is.finite(uIF), uIF > 0) %>%   # unrounded uIF (2026-11-01)
  group_by(gene_id, uORF_start_genomic, uORF_stop_genomic) %>% slice_max(uIF, n = 1, with_ties = FALSE) %>% ungroup() %>%   # V23: de-duplication restored (a comment had disabled it; no duplicates today)
  group_by(gene_id) %>% arrange(ifelse(Strand == "+", uORF_start_genomic, -uORF_start_genomic), .by_group = TRUE) %>% mutate(rank = row_number(), N = n()) %>% ungroup()
pairs <- bind_rows(lapply(2:4, function(k) { r <- R %>% filter(N == k, rank %in% c(1L, k))
  d <- r %>% filter(rank == 1L) %>% select(gene_id, dist_uid = uid, dist_tx = tx_id, dist_pos = uORF_start_genomic, dist_pUP = pct_uORF_plus, chr, Strand)
  p <- r %>% filter(rank == k) %>% select(gene_id, prox_uid = uid, prox_tx = tx_id, prox_pos = uORF_start_genomic, prox_pUP = pct_uORF_plus)
  inner_join(d, p, by = "gene_id") %>% mutate(N = k) }))
pairs$delta_obs <- pairs$prox_pUP - pairs$dist_pUP
for (k in 2:4) { z <- pairs[pairs$N == k, ]; .add(sprintf("  %d-TuORF genes (n = %d): median %%uORF+ distal %.1f vs proximal %.1f; median difference %.1f points; distal < proximal in %s", k, nrow(z), med(z$dist_pUP), med(z$prox_pUP), med(z$delta_obs), pct(mean(z$delta_obs > 0)))) }
# exon-space indices of both starts on one transcript of the gene (try the proximal TuORF's isoform, then the distal one's)
setkey(P, tx, chr, strand, pos)
get_idx <- function(tx_, chr_, strand_, pos_) { q <- data.table(tx = tx_, chr = chr_, strand = strand_, pos = as.integer(pos_)); P[q, idx] }
pairs <- as.data.frame(pairs)
pairs$tx <- pairs$prox_tx; pairs$i_prox <- get_idx(pairs$prox_tx, pairs$chr, pairs$Strand, pairs$prox_pos); pairs$i_dist <- get_idx(pairs$prox_tx, pairs$chr, pairs$Strand, pairs$dist_pos)
alt <- is.na(pairs$i_dist) | is.na(pairs$i_prox)
pairs$tx[alt] <- pairs$dist_tx[alt]; pairs$i_prox[alt] <- get_idx(pairs$dist_tx[alt], pairs$chr[alt], pairs$Strand[alt], pairs$prox_pos[alt]); pairs$i_dist[alt] <- get_idx(pairs$dist_tx[alt], pairs$chr[alt], pairs$Strand[alt], pairs$dist_pos[alt])
pairs <- pairs %>% left_join(TXI %>% select(tx, first_annot, idx_cds, leader_len, tot_c), by = "tx") %>% mutate(spacing = i_prox - i_dist, d_prox = idx_cds - i_prox, d_dist = idx_cds - i_dist)
ok <- !is.na(pairs$i_prox) & !is.na(pairs$i_dist) & pairs$spacing > 0 & pairs$tot_c > 0
.add(sprintf("  pairs with both start codons on one isoform's leader and CAGE present: %d of %d", sum(ok), nrow(pairs)))
# null 1: shift the pair (spacing kept) to every valid position of the same leader; exact expectation and percentile of the observed difference
Pc <- P[, .(tx, idx, annot, f_c)]; setkey(Pc, tx, idx)
lookup_f <- function(tx_, idx_) { q <- data.table(tx = rep_len(tx_, length(idx_)), idx = as.integer(idx_)); Pc[q, f_c] }
null1 <- function(tx_, i_dist, i_prox, first_annot, idx_cds) { s <- i_prox - i_dist; if (idx_cds - 1L - s < first_annot) return(c(NA, NA, NA)); cand <- first_annot:(idx_cds - 1L - s)
  fd <- lookup_f(tx_, cand); fp <- lookup_f(tx_, cand + s); dl <- 100 * (fp - fd); obs <- 100 * (lookup_f(tx_, i_prox) - lookup_f(tx_, i_dist))
  c(mean(dl, na.rm = TRUE), mean(dl >= obs, na.rm = TRUE), length(cand)) }
n1 <- t(mapply(null1, pairs$tx[ok], pairs$i_dist[ok], pairs$i_prox[ok], pairs$first_annot[ok], pairs$idx_cds[ok]))
pairs$null1_mean <- NA_real_; pairs$null1_pct_ge_obs <- NA_real_; pairs$null1_n_positions <- NA_integer_
pairs$null1_mean[ok] <- n1[, 1]; pairs$null1_pct_ge_obs[ok] <- n1[, 2]; pairs$null1_n_positions[ok] <- as.integer(n1[, 3])
pairs$delta_obs_tx <- NA_real_; pairs$delta_obs_tx[ok] <- 100 * (lookup_f(pairs$tx[ok], pairs$i_prox[ok]) - lookup_f(pairs$tx[ok], pairs$i_dist[ok]))
# null 2: the same distances from the mORF start in No-TuORF control leaders (N_CTRL random controls per pair, leader long enough, CAGE present)
ctrl_i <- TXI[tx %in% ctrl_tx$tx_id & tot_c > 0]
null2 <- function(dd_, dp_) { cc <- ctrl_i[leader_len >= dd_ + 1L]; if (nrow(cc) < 5) return(NA_real_); cc <- cc[sample(.N, min(N_CTRL, .N))]
  fd <- lookup_f(cc$tx, cc$idx_cds - dd_); fp <- lookup_f(cc$tx, cc$idx_cds - dp_); mean(100 * (fp - fd), na.rm = TRUE) }
pairs$null2_mean <- NA_real_; pairs$null2_mean[ok] <- mapply(null2, pairs$d_dist[ok], pairs$d_prox[ok])
s3_tab <- bind_rows(lapply(2:4, function(k) { z <- pairs[pairs$N == k & ok, ]
  w1 <- wilcox.test(z$delta_obs_tx, z$null1_mean, paired = TRUE); w2 <- wilcox.test(z$delta_obs_tx, z$null2_mean, paired = TRUE)
  data.frame(N_TuORFs = k, n_pairs = nrow(z), median_dist_pUP = med(z$dist_pUP), median_prox_pUP = med(z$prox_pUP), median_delta_obs = med(z$delta_obs_tx), median_null1 = med(z$null1_mean), median_null2 = med(z$null2_mean),
             pct_obs_gt_null1 = 100 * mean(z$delta_obs_tx > z$null1_mean, na.rm = TRUE), pct_obs_gt_null2 = 100 * mean(z$delta_obs_tx > z$null2_mean, na.rm = TRUE), pct_obs_in_top5_of_null1 = 100 * mean(z$null1_pct_ge_obs <= 0.05, na.rm = TRUE), q25_delta_obs = unname(quantile(z$delta_obs_tx, .25)), q75_delta_obs = unname(quantile(z$delta_obs_tx, .75)), q25_null1 = unname(quantile(z$null1_mean, .25, na.rm = TRUE)), q75_null1 = unname(quantile(z$null1_mean, .75, na.rm = TRUE)), pct_delta_obs_zero = 100 * mean(z$delta_obs_tx <= 0), pct_delta_obs_gt20 = 100 * mean(z$delta_obs_tx > 20),
             p_wilcoxon_vs_null1 = w1$p.value, p_wilcoxon_vs_null2 = w2$p.value, median_spacing_nt = med(z$spacing), median_d_prox = med(z$d_prox), median_leader_len = med(z$leader_len)) }))
for (i in seq_len(nrow(s3_tab))) with(s3_tab[i, ], .add(sprintf("  %d-TuORF genes (n = %d): observed difference median %.1f points | null 1 (same leader, spacing kept) %.1f | null 2 (No-TuORF leaders, same distances to the mORF) %.1f | observed > null 1 in %.1f%% of genes (p = %.1e), > null 2 in %.1f%% (p = %.1e) | observed in the top 5%% of its own null in %.1f%% | median spacing %d nt, proximal start %d nt from the mORF, leader %d nt",
  N_TuORFs, n_pairs, median_delta_obs, median_null1, median_null2, pct_obs_gt_null1, p_wilcoxon_vs_null1, pct_obs_gt_null2, p_wilcoxon_vs_null2, pct_obs_in_top5_of_null1, as.integer(median_spacing_nt), as.integer(median_d_prox), as.integer(median_leader_len))))
for (i in seq_len(nrow(s3_tab))) with(s3_tab[i, ], .add(sprintf("    %d-TuORF genes: observed difference IQR %.1f-%.1f (null 1 IQR %.1f-%.1f); difference <= 0 in %.1f%% of genes, > 20 points in %.1f%%", N_TuORFs, q25_delta_obs, q75_delta_obs, q25_null1, q75_null1, pct_delta_obs_zero, pct_delta_obs_gt20)))

################################################################################
# 3) (b) Uncertainty of %uORF+ in the single-TuORF cohort
################################################################################
.add("A6b: uncertainty of %uORF+")
s2i <- start_idx(S2)
U <- data.table(gene_id = S2$gene_id, uid = S2$uid, tx = S2$tx_id, pUP = S2$pct_uORF_plus, idx = s2i$idx, f_c = s2i$f_c, cum_c = s2i$cum_c, tot_c = s2i$tot_c, cum_t = s2i$cum_t, tot_t = s2i$tot_t)
U[, pUP_TSS := S1B$pct_uORF_plus[match(uid, S1B$uid)]]
U[, `:=`(n_eff = round(tot_c / u_mean), k_eff = pmin(round(cum_c / u_mean), round(tot_c / u_mean)))]
jeff <- function(k, n) { lo <- qbeta(0.025, k + 0.5, n - k + 0.5); hi <- qbeta(0.975, k + 0.5, n - k + 0.5); lo[!is.na(k) & k == 0] <- 0; hi[!is.na(k) & k == n] <- 1; lo[is.na(n) | n == 0] <- NA; hi[is.na(n) | n == 0] <- NA; list(lo = 100 * lo, hi = 100 * hi) }   # Jeffreys interval, boundary-modified (k = 0 -> lower 0; k = n -> upper 100)
j <- jeff(U$k_eff, U$n_eff); U[, `:=`(ci_lo = j$lo, ci_hi = j$hi)]
jt <- jeff(U$cum_t, U$tot_t); U[, `:=`(ci_lo_TSS = jt$lo, ci_hi_TSS = jt$hi)]
# binomial bootstrap (percentile) as a check of the Jeffreys interval
bt <- t(mapply(function(k, n) if (!is.na(n) && n > 0) 100 * quantile(rbinom(B, n, k / n) / n, c(0.025, 0.975)) else c(NA_real_, NA_real_), U$k_eff, U$n_eff)); U[, `:=`(boot_lo = bt[, 1], boot_hi = bt[, 2])]
.add(sprintf("  single-TuORF cohort n = %d (%d with the start codon on the isoform leader window); effective CAGE tags per window: median %d (IQR %d-%d); TSS-seq tags: median %d (IQR %d-%d)", nrow(U), sum(!is.na(U$idx)), as.integer(med(U$n_eff)), as.integer(quantile(U$n_eff, .25, na.rm = TRUE)), as.integer(quantile(U$n_eff, .75, na.rm = TRUE)), as.integer(med(U$tot_t)), as.integer(quantile(U$tot_t, .25, na.rm = TRUE)), as.integer(quantile(U$tot_t, .75, na.rm = TRUE))))
.add(sprintf("  Jeffreys vs bootstrap 95%% limits: median |difference| %.2f (lower) / %.2f (upper) points", med(abs(U$ci_lo - U$boot_lo)), med(abs(U$ci_hi - U$boot_hi))))
.add(sprintf("  CI width (CAGE): median %.1f points (IQR %.1f-%.1f); genes with width > 20 points: %s", med(U$ci_hi - U$ci_lo), quantile(U$ci_hi - U$ci_lo, .25, na.rm = TRUE), quantile(U$ci_hi - U$ci_lo, .75, na.rm = TRUE), pct(mean(U$ci_hi - U$ci_lo > 20, na.rm = TRUE))))
.add(sprintf("  point estimates: below 100%% uORF+ %s (CAGE), %s (TSS-seq), %s (both); below 90%%: %s / %s / %s", pct(mean(U$pUP < 100)), pct(mean(U$pUP_TSS < 100, na.rm = TRUE)), pct(mean(U$pUP < 100 & U$pUP_TSS < 100, na.rm = TRUE)), pct(mean(U$pUP < 90)), pct(mean(U$pUP_TSS < 90, na.rm = TRUE)), pct(mean(U$pUP < 90 & U$pUP_TSS < 90, na.rm = TRUE))))
.add(sprintf("  upper 95%% limit below 100%% (confidently affected): %s (CAGE), %s (TSS-seq), %s (both); upper limit below 90%%: %s / %s / %s", pct(mean(U$ci_hi < 100, na.rm = TRUE)), pct(mean(U$ci_hi_TSS < 100, na.rm = TRUE)), pct(mean(U$ci_hi < 100 & U$ci_hi_TSS < 100, na.rm = TRUE)), pct(mean(U$ci_hi < 90, na.rm = TRUE)), pct(mean(U$ci_hi_TSS < 90, na.rm = TRUE)), pct(mean(U$ci_hi < 90 & U$ci_hi_TSS < 90, na.rm = TRUE))))
.add(sprintf("  interval includes 100%% (not distinguishable from full inclusion): %s (CAGE), %s (TSS-seq)", pct(mean(U$ci_hi >= 100 - 1e-9, na.rm = TRUE)), pct(mean(U$ci_hi_TSS >= 100 - 1e-9, na.rm = TRUE))))
# genes at exactly 100%: what downstream fraction could have been missed?
g100 <- U[pUP == 100]
.add(sprintf("  genes at 100%% uORF+ (n = %d): median %d effective tags; the 95%% upper bound on an undetected downstream fraction (100 - lower limit) has median %.1f points (IQR %.1f-%.1f); %s of them have TSS-seq %%uORF+ below 100", nrow(g100), as.integer(med(g100$n_eff)), med(100 - g100$ci_lo), quantile(100 - g100$ci_lo, .25, na.rm = TRUE), quantile(100 - g100$ci_lo, .75, na.rm = TRUE), pct(mean(g100$pUP_TSS < 100, na.rm = TRUE))))
# depth dependence
U[, depth_q := cut(n_eff, quantile(n_eff, 0:4 / 4, na.rm = TRUE), include.lowest = TRUE, labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))]
dq <- U[!is.na(depth_q), .(n = .N, median_tags = as.integer(median(n_eff)), pct_below_100 = 100 * mean(pUP < 100), pct_below_90 = 100 * mean(pUP < 90), pct_le_20 = 100 * mean(pUP <= 20), pct_ci_hi_below_100 = 100 * mean(ci_hi < 100, na.rm = TRUE), pct_ci_hi_below_90 = 100 * mean(ci_hi < 90, na.rm = TRUE), pct_le_30 = 100 * mean(pUP <= 30)), by = depth_q][order(depth_q)]
for (i in seq_len(nrow(dq))) with(dq[i], .add(sprintf("  depth %s: n = %d, median %d tags; below 100%%: %.1f%%; below 90%%: %.1f%%; <= 30%%: %.1f%%; upper limit below 100%%: %.1f%%, below 90%%: %.1f%%", depth_q, n, median_tags, pct_below_100, pct_below_90, pct_le_30, pct_ci_hi_below_100, pct_ci_hi_below_90)))
# depth-standardised estimates: downsample every window to a fixed number of effective tags
ds <- function(n_fix, reps = 200) { z <- U[!is.na(n_eff) & n_eff >= n_fix]; p <- z$k_eff / z$n_eff
  m <- replicate(reps, { k <- rbinom(nrow(z), n_fix, p); c(below100 = mean(k < n_fix), below90 = mean(k / n_fix < 0.9), le30 = mean(k / n_fix <= 0.3), le20 = mean(k / n_fix <= 0.2)) })
  c(n_genes = nrow(z), rowMeans(m)) }
ds_tab <- do.call(rbind, lapply(c(20, 50, 100, 200), function(nf) c(tags = nf, ds(nf)))); ds_tab <- as.data.frame(ds_tab)
for (i in seq_len(nrow(ds_tab))) with(ds_tab[i, ], .add(sprintf("  downsampled to %d tags (n = %d genes with at least that depth): below 100%% %.1f%%, below 90%% %.1f%%, <= 30%% %.1f%%, <= 20%% %.1f%%", tags, n_genes, 100 * below100, 100 * below90, 100 * le30, 100 * le20)))

# (b2) the same thresholds for ALL TuORFs (Table S1A/S1B), per TuORF and per gene
.add("A6b2: all TuORFs (Table S1A), per TuORF and per gene")
A <- data.table(gene_id = S1A$gene_id, uid = S1A$uid, tx = S1A$tx_id, pUP = S1A$pct_uORF_plus, pUP_TSS = S1B$pct_uORF_plus[match(S1A$uid, S1B$uid)], n_in_gene = S1A$n_TuORFs_in_gene,
                idx = s1$idx, cum_c = s1$cum_c, tot_c = s1$tot_c, cum_t = s1$cum_t, tot_t = s1$tot_t)
A[, `:=`(n_eff = round(tot_c / u_mean), k_eff = pmin(round(cum_c / u_mean), round(tot_c / u_mean)))]
jA <- jeff(A$k_eff, A$n_eff); A[, `:=`(ci_lo = jA$lo, ci_hi = jA$hi)]; jAt <- jeff(A$cum_t, A$tot_t); A[, `:=`(ci_lo_TSS = jAt$lo, ci_hi_TSS = jAt$hi)]
.add(sprintf("  all TuORFs (n = %d in %d genes; %d in multi-TuORF genes): below 100%% %s (CAGE) / %s (TSS-seq) / %s (both); below 90%% %s / %s / %s; <= 30%% %s / %s; upper 95%% limit below 90%% %s (CAGE); effective tags median %d",
             nrow(A), uniqueN(A$gene_id), sum(A$n_in_gene > 1), pct(mean(A$pUP < 100)), pct(mean(A$pUP_TSS < 100, na.rm = TRUE)), pct(mean(A$pUP < 100 & A$pUP_TSS < 100, na.rm = TRUE)),
             pct(mean(A$pUP < 90)), pct(mean(A$pUP_TSS < 90, na.rm = TRUE)), pct(mean(A$pUP < 90 & A$pUP_TSS < 90, na.rm = TRUE)), pct(mean(A$pUP <= 30)), pct(mean(A$pUP_TSS <= 30, na.rm = TRUE)), pct(mean(A$ci_hi < 90, na.rm = TRUE)), as.integer(med(A$n_eff))))
for (k in c("1", "2", "3", ">=4")) { z <- if (k == ">=4") A[n_in_gene >= 4] else A[n_in_gene == as.integer(k)]
  .add(sprintf("    TuORFs in genes with %s TuORF(s) (n = %d): below 100%% %s, below 90%% %s, <= 30%% %s (CAGE)", k, nrow(z), pct(mean(z$pUP < 100)), pct(mean(z$pUP < 90)), pct(mean(z$pUP <= 30)))) }
G <- A[, .(n = .N, min_pUP = min(pUP), max_pUP = max(pUP), min_TSS = suppressWarnings(min(pUP_TSS, na.rm = TRUE)), max_TSS = suppressWarnings(max(pUP_TSS, na.rm = TRUE))), by = gene_id]
G[!is.finite(min_TSS), min_TSS := NA]; G[!is.finite(max_TSS), max_TSS := NA]
.add(sprintf("  all TuORF genes (n = %d; %d with one TuORF, %d with several): at least one TuORF below 90%% (= the 5'-most): %s (CAGE) / %s (TSS-seq); every TuORF below 90%% (= the mORF-proximal one): %s / %s; at least one below 100%%: %s / %s; every TuORF below 100%%: %s / %s; at least one <= 30%%: %s; every TuORF <= 30%%: %s",
             nrow(G), sum(G$n == 1), sum(G$n > 1), pct(mean(G$min_pUP < 90)), pct(mean(G$min_TSS < 90, na.rm = TRUE)), pct(mean(G$max_pUP < 90)), pct(mean(G$max_TSS < 90, na.rm = TRUE)), pct(mean(G$min_pUP < 100)), pct(mean(G$min_TSS < 100, na.rm = TRUE)), pct(mean(G$max_pUP < 100)), pct(mean(G$max_TSS < 100, na.rm = TRUE)), pct(mean(G$min_pUP <= 30)), pct(mean(G$max_pUP <= 30))))
dsA <- function(n_fix, reps = 200) { z <- A[!is.na(n_eff) & n_eff >= n_fix]; p <- z$k_eff / z$n_eff
  m <- replicate(reps, { k <- rbinom(nrow(z), n_fix, p); c(below100 = mean(k < n_fix), below90 = mean(k / n_fix < 0.9), le30 = mean(k / n_fix <= 0.3), le20 = mean(k / n_fix <= 0.2)) }); c(n_TuORFs = nrow(z), rowMeans(m)) }
dsA_tab <- as.data.frame(do.call(rbind, lapply(c(20, 50, 100, 200), function(nf) c(tags = nf, dsA(nf)))))
for (i in seq_len(nrow(dsA_tab))) with(dsA_tab[i, ], .add(sprintf("  all TuORFs downsampled to %d tags (n = %d): below 100%% %.1f%%, below 90%% %.1f%%, <= 30%% %.1f%%, <= 20%% %.1f%%", tags, n_TuORFs, 100 * below100, 100 * below90, 100 * le30, 100 * le20)))
all_tab <- A[, .(uid, gene_id, tx, n_TuORFs_in_gene = n_in_gene, pUP_CAGE = pUP, eff_tags_CAGE = n_eff, ci95_lo = round(ci_lo, 2), ci95_hi = round(ci_hi, 2), pUP_TSSseq = pUP_TSS, tags_TSSseq = tot_t, ci95_lo_TSSseq = round(ci_lo_TSS, 2), ci95_hi_TSSseq = round(ci_hi_TSS, 2))]
fwrite(all_tab, file.path(OUT_DIR, "A6_Table_pUP_uncertainty_all_TuORFs.csv")); fwrite(dsA_tab, file.path(OUT_DIR, "A6_Table_pUP_downsampled_all_TuORFs.csv"))

################################################################################
# 4) (c) uORF length vs exclusion: positional null and adjusted models
################################################################################
.add("A6c: uORF length vs exclusion — positional null")
L <- data.table(gene_id = S2$gene_id, uid = S2$uid, tx = S2$tx_id, len_aa = S2$uORF_length_aa, len_nt = S2$uORF_length_nt + 3L, pUP = S2$pct_uORF_plus, uIF = S2$uIF, idx_start = s2i$idx, f_c = s2i$f_c)
L <- merge(L, TXI[, .(tx, first_annot, idx_cds, leader_len, tot_c)], by = "tx")
L[, `:=`(d_start = idx_cds - idx_start, pos_from_TSS = idx_start - first_annot)]; L[, d_stop := d_start - len_nt]
L[, excluded := pUP <= EXCL]
L[, bin := cut(len_aa, c(0, 1, 10, 20, 30, Inf), labels = c("1", "2-10", "11-20", "21-30", ">30"))]
.add(sprintf("  cohort n = %d; start codon on the isoform's leader window: %d; CAGE present: %d; uORF start to mORF start: median %d nt (IQR %d-%d); leader %d nt (IQR %d-%d)", nrow(L), sum(!is.na(L$idx_start)), sum(L$tot_c > 0), as.integer(med(L$d_start)), as.integer(quantile(L$d_start, .25, na.rm = TRUE)), as.integer(quantile(L$d_start, .75, na.rm = TRUE)), as.integer(med(L$leader_len)), as.integer(quantile(L$leader_len, .25)), as.integer(quantile(L$leader_len, .75))))
.add(sprintf("  Spearman rho: length vs start-to-mORF distance %.2f; length vs leader length %.2f; start-to-mORF distance vs %%uORF+ %.2f; length vs %%uORF+ %.2f", cor(L$len_aa, L$d_start, method = "spearman", use = "complete.obs"), cor(L$len_aa, L$leader_len, method = "spearman"), cor(L$d_start, L$pUP, method = "spearman", use = "complete.obs"), cor(L$len_aa, L$pUP, method = "spearman")))
# random-placement null: an ORF of the same length at every valid position of the same leader (stop codon may overlap the mORF start by <= STOP_OVERLAP nt)
pl_null <- function(tx_, len_, first_annot, idx_cds) { hi <- min(idx_cds - 1L, idx_cds - 1L + STOP_OVERLAP - len_ + 1L); if (hi < first_annot) return(c(NA, NA, NA)); cand <- first_annot:hi
  f <- lookup_f(tx_, cand); c(mean(100 * f <= EXCL, na.rm = TRUE), mean(100 * f, na.rm = TRUE), length(cand)) }
okL <- !is.na(L$idx_start) & L$tot_c > 0
pn <- t(mapply(pl_null, L$tx[okL], L$len_nt[okL], L$first_annot[okL], L$idx_cds[okL]))
L[, `:=`(exp_excl = NA_real_, exp_pUP = NA_real_, n_positions = NA_integer_)]; L$exp_excl[okL] <- pn[, 1]; L$exp_pUP[okL] <- pn[, 2]; L$n_positions[okL] <- as.integer(pn[, 3])
# also the null restricted to placements with the same distance class to the mORF? (not needed: the regression below adjusts for d_start)
boot_diff <- function(obs, expct, B = 2000) { n <- length(obs); d <- replicate(B, { i <- sample.int(n, n, TRUE); mean(obs[i]) - mean(expct[i]) }); quantile(d, c(0.025, 0.975)) }
len_tab <- L[okL & !is.na(exp_excl), { ci <- boot_diff(as.numeric(excluded), exp_excl)
  .(n = .N, observed_pct_excluded = 100 * mean(excluded), expected_pct_excluded = 100 * mean(exp_excl), excess_points = 100 * (mean(excluded) - mean(exp_excl)), excess_ci_lo = 100 * ci[1], excess_ci_hi = 100 * ci[2], ratio_obs_exp = mean(excluded) / mean(exp_excl),
    median_len_aa = median(len_aa), median_d_start = as.integer(median(d_start)), median_positions = as.integer(median(n_positions))) }, by = bin][order(bin)]
for (i in seq_len(nrow(len_tab))) with(len_tab[i], .add(sprintf("  length %s aa (n = %d): observed excluded %.1f%% | expected from random placement %.1f%% | excess %+.1f points (95%% CI %+.1f to %+.1f) | ratio %.2f | start-to-mORF median %d nt | placements %d", bin, n, observed_pct_excluded, expected_pct_excluded, excess_points, excess_ci_lo, excess_ci_hi, ratio_obs_exp, median_d_start, median_positions)))
# stratification by start-to-mORF distance
L[, d_q := cut(d_start, quantile(d_start, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE, dig.lab = 5)]
L[, len_cls := cut(len_aa, c(0, 10, 20, Inf), labels = c("1-10 aa", "11-20 aa", ">20 aa"))]
strat <- L[!is.na(d_q), .(n = .N, pct_excluded = 100 * mean(excluded), median_d_start = as.integer(median(d_start))), by = .(d_q, len_cls)][order(d_q, len_cls)]
.add("  exclusion by start-to-mORF distance quintile x length class:")
for (q in levels(L$d_q)) { z <- strat[d_q == q]; .add(sprintf("    %-12s (median %4d nt): %s", q, z$median_d_start[1], paste(sprintf("%s %.0f%% (n=%d)", z$len_cls, z$pct_excluded, z$n), collapse = " | "))) }
# logistic models
# V49 (2026-11-17, Larry): all TuORFs, the 1-aa uORFs included. They have a term of their own (one_aa), because they are excluded more often than the trend of
# the longer uORFs predicts; 1-aa uORFs can stall ribosomes and repress strongly. Length is modelled per doubling as before.
Lm <- L[!is.na(d_start) & tot_c > 0]; Lm[, one_aa := as.integer(len_aa == 1)]
m0 <- glm(excluded ~ log2(len_aa) + one_aa, data = Lm, family = binomial)
m1 <- glm(excluded ~ log2(d_start) + log2(leader_len), data = Lm, family = binomial)
m2 <- glm(excluded ~ log2(d_start) + log2(leader_len) + log2(len_aa) + one_aa, data = Lm, family = binomial)
m3 <- glm(excluded ~ log2(pos_from_TSS + 1) + log2(leader_len) + log2(len_aa) + one_aa, data = Lm, family = binomial)
orx <- function(m, term) { cf <- coef(summary(m)); ci <- suppressMessages(confint.default(m)); c(OR = exp(cf[term, 1]), lo = exp(ci[term, 1]), hi = exp(ci[term, 2]), p = cf[term, 4]) }
T1 <- "1-aa uORF (vs the length trend)"
mod_tab <- rbind(data.frame(model = "excluded ~ length + 1-aa", term = "length (per doubling)", t(orx(m0, "log2(len_aa)"))),
                 data.frame(model = "excluded ~ length + 1-aa", term = T1, t(orx(m0, "one_aa"))),
                 data.frame(model = "excluded ~ start-to-mORF distance + leader length", term = "distance (per doubling)", t(orx(m1, "log2(d_start)"))),
                 data.frame(model = "excluded ~ start-to-mORF distance + leader length", term = "leader length (per doubling)", t(orx(m1, "log2(leader_len)"))),
                 data.frame(model = "excluded ~ distance + leader length + length + 1-aa", term = "length (per doubling)", t(orx(m2, "log2(len_aa)"))),
                 data.frame(model = "excluded ~ distance + leader length + length + 1-aa", term = T1, t(orx(m2, "one_aa"))),
                 data.frame(model = "excluded ~ distance + leader length + length + 1-aa", term = "distance (per doubling)", t(orx(m2, "log2(d_start)"))),
                 data.frame(model = "excluded ~ distance + leader length + length + 1-aa", term = "leader length (per doubling)", t(orx(m2, "log2(leader_len)"))),
                 data.frame(model = "excluded ~ distance from the annotated 5′ end + leader length + length + 1-aa", term = "length (per doubling)", t(orx(m3, "log2(len_aa)"))),
                 data.frame(model = "excluded ~ distance from the annotated 5′ end + leader length + length + 1-aa", term = T1, t(orx(m3, "one_aa"))),
                 data.frame(model = "excluded ~ distance from the annotated 5′ end + leader length + length + 1-aa", term = "distance from the annotated 5′ end (per doubling)", t(orx(m3, "log2(pos_from_TSS + 1)"))))
lrt <- anova(m1, m2, test = "Chisq")
.add(sprintf("  logistic models (n = %d, all TuORFs; the 1-aa uORFs with a term of their own): OR per doubling of length %.2f (95%% CI %.2f-%.2f), 1-aa %.2f (%.2f-%.2f); after adjusting for start-to-mORF distance and leader length: length %.2f (%.2f-%.2f), p = %.1e, 1-aa %.2f (%.2f-%.2f), p = %.1e; distance OR per doubling %.2f (%.2f-%.2f); likelihood-ratio test for adding length and the 1-aa term: p = %.1e; model with distance from the annotated 5′ end instead: length OR %.2f (%.2f-%.2f), 1-aa %.2f (%.2f-%.2f)",
             nrow(Lm), mod_tab$OR[1], mod_tab$lo[1], mod_tab$hi[1], mod_tab$OR[2], mod_tab$lo[2], mod_tab$hi[2], mod_tab$OR[5], mod_tab$lo[5], mod_tab$hi[5], mod_tab$p[5], mod_tab$OR[6], mod_tab$lo[6], mod_tab$hi[6], mod_tab$p[6],
             mod_tab$OR[7], mod_tab$lo[7], mod_tab$hi[7], lrt$`Pr(>Chi)`[2], mod_tab$OR[9], mod_tab$lo[9], mod_tab$hi[9], mod_tab$OR[10], mod_tab$lo[10], mod_tab$hi[10]))
# how much of the length gradient does position explain? predicted exclusion by bin from the distance-only model
Lm[, pred_m1 := predict(m1, type = "response")]
pm <- Lm[, .(n = .N, observed = 100 * mean(excluded), predicted_from_position_only = 100 * mean(pred_m1)), by = bin][order(bin)]
for (i in seq_len(nrow(pm))) with(pm[i], .add(sprintf("  length %s aa: observed %.1f%% excluded vs %.1f%% predicted by the start-to-mORF distance + leader length model (n = %d)", bin, observed, predicted_from_position_only, n)))

################################################################################
# 5) Tables
################################################################################
fwrite(pairs %>% select(gene_id, N, dist_uid, prox_uid, tx, dist_pUP, prox_pUP, delta_obs, delta_obs_tx, spacing, d_dist, d_prox, leader_len, null1_mean, null1_pct_ge_obs, null1_n_positions, null2_mean), file.path(OUT_DIR, "A6_Table_S3_pairs_with_nulls.csv"))
fwrite(s3_tab, file.path(OUT_DIR, "A6_Table_S3_summary.csv"))
fwrite(U[, .(gene_id, uid, tx, pUP_CAGE = pUP, eff_tags_CAGE = n_eff, eff_tags_upstream = k_eff, ci95_lo = round(ci_lo, 2), ci95_hi = round(ci_hi, 2), boot_lo = round(boot_lo, 2), boot_hi = round(boot_hi, 2), pUP_TSSseq = pUP_TSS, tags_TSSseq = tot_t, ci95_lo_TSSseq = round(ci_lo_TSS, 2), ci95_hi_TSSseq = round(ci_hi_TSS, 2), depth_quartile = depth_q)], file.path(OUT_DIR, "A6_Table_pUP_uncertainty_per_gene.csv"))
fwrite(dq, file.path(OUT_DIR, "A6_Table_pUP_by_depth.csv")); fwrite(ds_tab, file.path(OUT_DIR, "A6_Table_pUP_downsampled.csv"))
fwrite(L[, .(gene_id, uid, tx, len_aa, len_nt, bin, pUP, excluded, d_start, d_stop, pos_from_TSS, leader_len, n_positions, expected_excluded = round(exp_excl, 4), expected_pUP = round(exp_pUP, 2))], file.path(OUT_DIR, "A6_Table_length_positional_null_per_gene.csv"))
fwrite(len_tab, file.path(OUT_DIR, "A6_Table_length_positional_null.csv")); fwrite(strat, file.path(OUT_DIR, "A6_Table_length_by_distance_strata.csv")); fwrite(mod_tab, file.path(OUT_DIR, "A6_Table_logistic_models.csv")); fwrite(pm, file.path(OUT_DIR, "A6_Table_length_predicted_from_position.csv"))
# V49 (2026-11-17): the two %uORF+ calls by sequencing depth. Each read of a library is kept with probability f (1% to 100% of the reads);
# all single-TuORF genes, CAGE-seq (effective tags) and TSS-seq; mean of 100 replicates. Own seed, so that the random numbers of the analyses above are unchanged.
.seed_keep <- .Random.seed; set.seed(20261117)
FRAC <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
LIB <- c(`CAGE-seq` = sum(cg$count) / u_mean, `TSS-seq` = sum(ts$count))   # reads in the two files (CAGE-seq as effective tags)
sat_in <- list(`CAGE-seq` = U[!is.na(n_eff) & n_eff > 0, .(n = n_eff, k = k_eff)], `TSS-seq` = U[!is.na(tot_t) & tot_t > 0, .(n = tot_t, k = cum_t)])
sat_tab <- rbindlist(lapply(names(sat_in), function(m) { z <- sat_in[[m]]
  rbindlist(lapply(FRAC, function(f) { r <- rowMeans(replicate(100, { ku <- rbinom(nrow(z), z$k, f); kd <- rbinom(nrow(z), z$n - z$k, f); ok <- ku + kd > 0
                                                                     c(sum(ok), mean(kd[ok] > 0), mean(ku[ok] / (ku[ok] + kd[ok]) < 0.9)) }))
    data.table(method = m, fraction_of_reads = f, million_reads = f * LIB[[m]] / 1e6, genes_with_reads = r[1], pct_below_100 = 100 * r[2], pct_below_90 = 100 * r[3]) })) }))
.Random.seed <- .seed_keep
for (m in names(sat_in)) with(sat_tab[method == m], .add(sprintf("  library sampling, %s (%.1f million reads): below 100%% uORF+ %.1f / %.1f / %.1f%%, below 90%% %.1f / %.1f / %.1f%% with 10%% / 50%% / 100%% of the reads",
  m, LIB[[m]] / 1e6, pct_below_100[fraction_of_reads == 0.1], pct_below_100[fraction_of_reads == 0.5], pct_below_100[fraction_of_reads == 1], pct_below_90[fraction_of_reads == 0.1], pct_below_90[fraction_of_reads == 0.5], pct_below_90[fraction_of_reads == 1])))
fwrite(sat_tab, file.path(OUT_DIR, "A6_Table_pUP_library_fraction.csv"))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S4 — Null models and uncertainty for the TSS effects (Figure S3).",
  "multiTuORF_pairs: genes with exactly 2, 3 or 4 TuORFs (uIF > 0; the Supplemental Fig. S9 set); distal = 5'-most, proximal = mORF-proximal TuORF; %uORF+ from CAGE (pUP in the column names). delta_obs_tx = proximal minus distal %uORF+ on the common isoform; delta_obs = the same difference with the Table S1A values of the two TuORFs, each on its listed isoform; null1_mean = mean difference when the pair is shifted (spacing kept) to every position of the same leader; null1_pct_ge_obs = fraction (0 to 1) of those placements with a difference at least as large as observed; null2_mean = mean difference at the same distances from the mORF start in No-TuORF leaders (50 random controls). N = number of TuORFs of the gene; spacing = nt between the two start codons; d_dist and d_prox = nt from the start codon of the distal and of the proximal TuORF to the mORF start; null1_n_positions = number of placements of the pair. delta_obs_tx and the null columns are NA (empty cell) when the two TuORFs have the same start codon (spacing = 0), when a start codon is not in the leader window of the common isoform, or when that window has no CAGE signal. The null1 columns are also NA when the spacing is at least the length of the annotated leader (leader_len), so that the pair has no placement in it.",
  sprintf("pUP_uncertainty: single-TuORF cohort; effective CAGE tags = CAGE signal / %.5f (mean single-tag unit of the three TPM-scaled replicates); ci95_lo and ci95_hi = limits of the 95%% Jeffreys interval on the effective tag counts (a binomial percentile bootstrap gave limits that differ by a median of 0.13 points for the lower and 0.06 points for the upper limit; not listed); TSS-seq columns use the integer TSS-seq counts and Table S1B. Without downsampling, %s of %s genes (%.0f%%) have less than 90%% uORF+ with CAGE-seq and %s of %s (%.0f%%) with TSS-seq (the %d genes without TSS-seq reads count as 0%% uORF+).",
          u_mean, format(sum(U$pUP < 90), big.mark = ","), format(nrow(U), big.mark = ","), 100 * mean(U$pUP < 90), format(sum(U$pUP_TSS < 90, na.rm = TRUE), big.mark = ","), format(sum(!is.na(U$pUP_TSS)), big.mark = ","), 100 * mean(U$pUP_TSS < 90, na.rm = TRUE), sum(U$tot_t == 0, na.rm = TRUE)),   # V41 (2026-11-12): the 57% and 54% of the Results are printed here
  "pUP_by_depth / pUP_downsampled: genes below 100% / 90% / at or below 30% / at or below 20% uORF+ by depth quartile (percentages of genes) and after downsampling every window to a fixed tag number (fractions of genes, 0 to 1); pct_ci_hi_below_100 and pct_ci_hi_below_90 = percentage of genes whose upper limit (ci95_hi) is below 100% or 90%, computed before the limits were rounded to two decimals.",
  "pUP_all_TuORFs / pUP_downsampled_allTuORFs: the same intervals and downsampling for all 7,423 TuORFs of Table S1A (per TuORF).",
  "pUP_library_fraction: single-TuORF genes; each read of the CAGE-seq or TSS-seq library kept at random with the probability fraction_of_reads (mean of 100 replicates); genes_with_reads = genes with at least one read; pct_below_100 and pct_below_90 = percentages of these genes below 100% and below 90% uORF+ (Supplemental Fig. S3B).",
  "multiTuORF_summary: medians of the observed and expected proximal-minus-distal differences by number of TuORFs (its columns with delta_obs in the name are computed from delta_obs_tx). length_null_per_gene: one row per TuORF; excluded = TRUE when %uORF+ <= 30; expected_excluded = fraction (0 to 1) of the random placements with %uORF+ <= 30; d_start = nt from the uORF start codon to the mORF start; pos_from_TSS = nt from the annotated 5' end of the leader to the start codon; n_positions = number of placements; bin = length class (aa). length_pred_position: exclusion predicted by the distance and leader-length model, by length class.",
  sprintf("length_null: random-placement null — an ORF of the same length placed at every valid position of the same leader (stop codon may overlap the mORF start by up to %d nt); the table gives, by length class (bin), the observed and the expected percentage of uORFs with %%uORF+ <= %d; excess_points = observed minus expected percentage; excess_ci_lo and excess_ci_hi = limits of its 95%% interval (bootstrap over genes, 2,000 resamples); median_d_start = median nt from the uORF start codon to the mORF start; median_positions = median number of placements. length_strata: exclusion by quintile of the start-to-mORF distance (d_q) and length class (len_cls). logistic_models: odds ratios (OR) per doubling; lo and hi = limits of the 95%% confidence interval.", STOP_OVERLAP, EXCL))), colNames = FALSE)
addWorksheet(wb, "multiTuORF_summary"); writeData(wb, "multiTuORF_summary", s3_tab)
addWorksheet(wb, "multiTuORF_pairs"); writeData(wb, "multiTuORF_pairs", pairs %>% select(gene_id, N, dist_uid, prox_uid, tx, dist_pUP, prox_pUP, delta_obs, delta_obs_tx, spacing, d_dist, d_prox, leader_len, null1_mean, null1_pct_ge_obs, null1_n_positions, null2_mean))
addWorksheet(wb, "pUP_uncertainty"); writeData(wb, "pUP_uncertainty", U[, .(gene_id, uid, tx, pUP_CAGE = pUP, eff_tags_CAGE = n_eff, eff_tags_upstream = k_eff, ci95_lo = round(ci_lo, 2), ci95_hi = round(ci_hi, 2), pUP_TSSseq = pUP_TSS, tags_TSSseq = tot_t, ci95_lo_TSSseq = round(ci_lo_TSS, 2), ci95_hi_TSSseq = round(ci_hi_TSS, 2), depth_quartile = depth_q)])
addWorksheet(wb, "pUP_by_depth"); writeData(wb, "pUP_by_depth", dq); addWorksheet(wb, "pUP_downsampled"); writeData(wb, "pUP_downsampled", ds_tab)
addWorksheet(wb, "pUP_all_TuORFs"); writeData(wb, "pUP_all_TuORFs", all_tab); addWorksheet(wb, "pUP_downsampled_allTuORFs"); writeData(wb, "pUP_downsampled_allTuORFs", dsA_tab)
addWorksheet(wb, "pUP_library_fraction"); writeData(wb, "pUP_library_fraction", sat_tab)
addWorksheet(wb, "length_null"); writeData(wb, "length_null", len_tab); addWorksheet(wb, "length_null_per_gene"); writeData(wb, "length_null_per_gene", L[, .(gene_id, uid, tx, len_aa, bin, pUP, excluded, d_start, pos_from_TSS, leader_len, n_positions, expected_excluded = round(exp_excl, 4))])
addWorksheet(wb, "length_strata"); writeData(wb, "length_strata", strat); addWorksheet(wb, "logistic_models"); writeData(wb, "logistic_models", mod_tab); addWorksheet(wb, "length_pred_position"); writeData(wb, "length_pred_position", pm)
saveWorkbook(wb, file.path(OUT_DIR, "A6_Table_nulls_and_uncertainty.xlsx"), overwrite = TRUE)

################################################################################
# 6) Figure
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3))
# panel C: baselines for the genes with two TuORFs (drawn as pA; the size of the earlier versions); the two nulls are explained in the schematic at its right (pNull)
lab_n1 <- "null 1"; lab_n2 <- "null 2"
z2 <- pairs[pairs$N == 2 & ok, ]
da <- rbind(data.frame(what = "observed", d = z2$delta_obs_tx), data.frame(what = lab_n1, d = z2$null1_mean), data.frame(what = lab_n2, d = z2$null2_mean))
da$what <- factor(da$what, levels = c("observed", lab_n1, lab_n2))
pA <- ggplot(da[!is.na(da$d), ], aes(d, colour = what)) + stat_ecdf(linewidth = 0.7) + scale_colour_manual(values = c("black", "#D7191C", "#2C7BB6"), name = NULL) + guides(colour = guide_legend(nrow = 3)) +
  labs(x = "proximal minus distal %uORF+ (percentage points)", y = "Fraction\nof genes", title = "C  Distal vs proximal TuORF", subtitle = sprintf("genes with two TuORFs (n = %s)\nmedians: observed %.0f, null 1: %.0f, null 2: %.0f points", big(nrow(z2)), med(z2$delta_obs_tx), med(z2$null1_mean), med(z2$null2_mean))) + th + theme(legend.direction = "vertical")
# panel A: pUP with intervals (drawn as pB)
db <- U[order(pUP)]; db[, rank := seq_len(.N) / .N]
pB <- ggplot(db, aes(x = rank)) + geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "grey70") + geom_line(aes(y = pUP), colour = "black", linewidth = 0.5) + geom_hline(yintercept = c(90, 100), linetype = 3, linewidth = 0.3) +
  labs(x = "Single-TuORF genes ranked by %uORF+", y = "%uORF+ (CAGE),\n95% interval", title = "A  %uORF+ and its 95% interval", subtitle = sprintf("Jeffreys interval on effective CAGE tags (median %d per gene); upper limit below 100%%: %s of genes, below 90%%: %s", as.integer(med(U$n_eff)), pct(mean(U$ci_hi < 100, na.rm = TRUE)), pct(mean(U$ci_hi < 90, na.rm = TRUE)))) + th
# panel B: the two calls by sequencing depth (pSat; the depth quartiles are in Table S4, sheet pUP_by_depth)
ds_plot <- melt(sat_tab[, .(method, million_reads, fraction_of_reads, `below 100%` = pct_below_100, `below 90%` = pct_below_90)], id.vars = c("method", "million_reads", "fraction_of_reads"), variable.name = "call", value.name = "pct")
pSat <- ggplot(ds_plot, aes(million_reads, pct, colour = method, linetype = call, group = interaction(method, call))) + geom_line(linewidth = 0.6) + geom_point(size = 1) +
  geom_point(data = ds_plot[fraction_of_reads == 1], shape = 21, fill = "white", size = 2, stroke = 0.7) +
  scale_colour_manual(values = c(`CAGE-seq` = "black", `TSS-seq` = "#2C7BB6"), name = NULL) + scale_linetype_manual(values = c(`below 100%` = "solid", `below 90%` = "22"), name = NULL) +
  scale_y_continuous(limits = c(45, 108), breaks = c(60, 80, 100)) +   # V57: room for the keys above the curves (the panel is a little lower)
  labs(x = "Reads sampled from the library (million)", y = "% of single-TuORF\ngenes", title = "B  %uORF+ calls by read depth", subtitle = "reads sampled at random; open circles: all reads") + th
# D: length vs exclusion observed vs random placement
dd <- melt(merge(len_tab[, .(bin, observed = observed_pct_excluded, expected = expected_pct_excluded)], pm[, .(bin, predicted = predicted_from_position_only)], by = "bin", all.x = TRUE), id.vars = "bin", variable.name = "what", value.name = "pct")
dd[, what := factor(what, levels = c("observed", "expected", "predicted"))]   # V49: short keys; the schematic at the right of D (pExp) explains them
# V57 (2026-11-24, Larry: "The legend and the right side of the panel D is too simple and cannot be understood. Add more clear explaination in both."): keys and subtitle say what each bar is
pD <- ggplot(dd[!is.na(pct)], aes(bin, pct, fill = what)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(observed = "#D7191C", expected = "grey60", predicted = "#2C7BB6"), labels = c(observed = "observed", expected = "expected for all possible positions", predicted = "predicted from the real position"), name = NULL) + guides(fill = guide_legend(nrow = 3)) +
  geom_text(data = len_tab, aes(bin, pmax(observed_pct_excluded, expected_pct_excluded) + 3, label = sprintf("%+.0f", excess_points)), inherit.aes = FALSE, size = 2.6) +
  labs(x = "uORF length (aa)", y = sprintf("%% excluded\n(%%uORF+ ≤ %d)", EXCL), title = "D  uORF length vs exclusion", subtitle = sprintf("%s of the %s single-TuORF genes; gray and blue are higher for longer uORFs (from 2–10 aa to more than 30 aa) because a longer uORF starts further from the mORF start codon, so fewer transcripts contain it; numbers above the bars: red minus gray (percentage points)", big(nrow(L[okL & !is.na(exp_excl)])), big(nrow(L)))) + th + theme(legend.direction = "vertical", legend.text = element_text(size = 6.5))
# room for the labels of D (the strata by distance are in Table S4, sheet length_strata)
pD <- pD + scale_x_discrete(labels = function(x) sub("-", "–", x, fixed = TRUE)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.45)))   # room for the difference labels above the bars and for the keys (V57: longer keys)
# panel E: odds ratios
# V49 (2026-11-17, Larry: "Yes, simplify E"): two models in four rows, the model with uORF length and the 1-aa term, before and after adjustment for position;
# the model with position only and the model with the distance from the annotated 5' end stay in Table S4 (sheet logistic_models)
M_BEFORE <- "excluded ~ length + 1-aa"; M_ADJ <- "excluded ~ distance + leader length + length + 1-aa"
mf <- as.data.table(mod_tab)[model %in% c(M_BEFORE, M_ADJ)]
mf[, row := factor(c(`length (per doubling)` = "uORF length (per doubling)", `1-aa uORF (vs the length trend)` = "1-aa uORF (compared with the length trend)",
                     `distance (per doubling)` = "distance from the uORF start to the mORF (per doubling)", `leader length (per doubling)` = "leader length (per doubling)")[term],
                   levels = rev(c("uORF length (per doubling)", "1-aa uORF (compared with the length trend)", "distance from the uORF start to the mORF (per doubling)", "leader length (per doubling)")))]
mf[, fit := factor(ifelse(model == M_BEFORE, "before adjustment", "adjusted for position"), levels = c("adjusted for position", "before adjustment"))]   # this order puts the model before adjustment above the adjusted one in each row
pF <- ggplot(mf, aes(OR, row, colour = fit)) + geom_vline(xintercept = 1, linetype = 3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", position = position_dodge(0.8)) + geom_point(size = 2, position = position_dodge(0.8)) +   # V57: 0.8 (was 0.55), so that the two labels of a row do not overlap
  geom_text(aes(x = hi, label = sprintf("%.2f", OR), vjust = ifelse(fit == "before adjustment", 0.15, 0.85)), hjust = -0.25, size = 2.8, position = position_dodge(0.8), show.legend = FALSE) +   # V57: the gray label above its line, the red label below its line
  scale_y_discrete(expand = expansion(add = c(0.75, 0.9))) +   # V57: room above the top row for its gray label
  scale_x_log10(expand = expansion(mult = c(0.05, 0.12))) + scale_colour_manual(values = c(`before adjustment` = "grey45", `adjusted for position` = "#D7191C"), breaks = c("before adjustment", "adjusted for position"), name = NULL) +
  labs(x = "Odds ratio for exclusion (95% confidence interval, log scale)", y = NULL, title = "E  Odds of exclusion before and after adjustment for position",
       subtitle = sprintf("%s of the %s single-TuORF genes; both models have uORF length and a term for the 1-aa uORFs; position means the distance from the uORF start codon to the mORF start codon and the leader length", big(nrow(Lm)), big(nrow(L)))) + th
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
# panel C, right: how the observed difference and the two nulls are calculated (schematic with text, 8 pt; V49, 2026-11-17)
PD <- 28; PP <- 52; XM <- 78; TXT <- FIG_PT_MIN / .pt; RED <- "#D7191C"; BLU <- "#2C7BB6"
sch_leader <- function(y, x0, pk_x = NULL, pk_h = NULL) list(
  annotate("segment", x = x0, xend = XM, y = y, yend = y, linewidth = 0.45, colour = "grey30"),
  if (length(pk_x)) annotate("segment", x = pk_x, xend = pk_x, y = y, yend = y + pk_h, linewidth = 0.8, colour = "grey60"),
  annotate("rect", xmin = XM, xmax = XM + 14, ymin = y - 2, ymax = y + 2, fill = "grey85", colour = "grey40", linewidth = 0.3))
note <- function(y, lab) annotate("text", x = 0, y = y, label = lab, hjust = 0, vjust = 1, size = TXT, family = FIG_FONT, lineheight = 0.95)
pk1 <- c(6, 10, 36, 44, 64); ph1 <- c(5, 3, 4, 6, 3)
y2 <- c(17, 10, 3); x2 <- c(10, 20, 4)
pNull <- ggplot() + scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + scale_y_continuous(limits = c(-26, 100), expand = c(0, 0)) + theme_void() + coord_cartesian(clip = "off") +
  annotate("rect", xmin = PD, xmax = PP, ymin = 86, ymax = 92, fill = "grey88") + sch_leader(88, 2, pk1, ph1) +
  annotate("segment", x = c(PD, PP), xend = c(PD, PP), y = 85, yend = 92, linewidth = 0.8, colour = "black") +
  annotate("text", x = c(PD, PP, XM + 7), y = 94, label = c("distal", "proximal", "mORF"), vjust = 0, size = TXT, family = FIG_FONT) +
  note(81, "Observed (black): % of the CAGE signal (gray bars)\nbetween the distal and the proximal start codons") +
  sch_leader(50, 2, pk1, ph1) +
  annotate("segment", x = c(4, PD, PP), xend = c(4, PD, PP) + (PP - PD), y = c(58, 61, 64), yend = c(58, 61, 64), linewidth = 0.7, colour = RED) +
  annotate("segment", x = c(4, PD, PD, PP, PP, PP + PP - PD), xend = c(4, PD, PD, PP, PP, PP + PP - PD), y = rep(c(58, 61, 64), each = 2) - 1.3, yend = rep(c(58, 61, 64), each = 2) + 1.3, linewidth = 0.7, colour = RED) +
  note(45, "null 1 (red): the pair, with its spacing kept,\nmoved to every position of the same leader;\nmean of the differences") +
  annotate("segment", x = c(PD, PP), xend = c(PD, PP), y = 0, yend = 21, linetype = 2, linewidth = 0.4, colour = BLU) +
  annotate("rect", xmin = PD, xmax = PP, ymin = y2 - 1.5, ymax = y2 + 1.5, fill = BLU, alpha = 0.2) +
  sch_leader(y2[1], x2[1], c(14, 40, 60), c(4, 5, 3)) + sch_leader(y2[2], x2[2], c(24, 48, 70), c(3, 4, 5)) + sch_leader(y2[3], x2[3], c(8, 33, 56), c(4, 3, 4)) +
  annotate("segment", x = rep(c(PD, PP), 3), xend = rep(c(PD, PP), 3), y = rep(y2, each = 2) - 2, yend = rep(y2, each = 2) + 2, linewidth = 0.8, colour = BLU) +
  note(-4, sprintf("null 2 (blue): positions at the same distances\nfrom the mORF start in %d random No-TuORF\nleaders; mean of the differences", N_CTRL)) +
  labs(title = "How the curves of C are calculated") +
  theme(plot.title = element_text(family = FIG_FONT, size = FIG_PT_MIN + 2, face = "bold"), plot.title.position = "plot", plot.margin = margin(4, 10, 6, 8))
# panel D, right: how the three bars of D are calculated (schematic with text, 8 pt; V49, 2026-11-17)
# V57 (2026-11-24, Larry: "The legend and the right side of the panel D is too simple and cannot be understood. Add more clear explaination in both."):
# one drawing and one paragraph for each bar, and a last paragraph on why red is lower than gray (why gray and blue rise with length: subtitle of D)
# y: about 0.85 pt per unit (the cell is about 2.8 in high; text lines are about 10 pt apart); at most about 72 characters per line in the 60% width
xs_g <- c(4, 17, 30, 43, 56); ys_g <- c(165.5, 169.5, 165.5, 169.5, 165.5); ph_pt <- c(8, 5, 7, 10, 5); ph_g <- c(5, 3, 4, 6, 3)
lab_on <- function(x, y, lab, hjust = 0.5) annotate("label", x = x, y = y, label = lab, hjust = hjust, size = TXT, family = FIG_FONT, fill = "white", border.colour = "white", label.padding = unit(0.6, "pt"))
pExp <- ggplot() + scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + scale_y_continuous(limits = c(-6, 200), expand = c(0, 0)) + theme_void() + coord_cartesian(clip = "off") +
  sch_leader(192, 2, pk1, ph_pt) + annotate("rect", xmin = 40, xmax = 58, ymin = 190.2, ymax = 193.8, fill = RED) +
  note(186, "Red, observed: percentage of TuORFs with %uORF+ ≤ 30 (excluded)") +
  sch_leader(155, 2, pk1, ph_g) +
  annotate("rect", xmin = xs_g, xmax = xs_g + 18, ymin = ys_g - 1.5, ymax = ys_g + 1.5, fill = "grey60") +
  note(150, "Gray, expected for all possible positions: an ORF of the same length is\nplaced at every nucleotide position of the leader where it fits. The bar\nis the percentage of positions with %uORF+ ≤ 30 (calculated from the CAGE\nsignal, gray ticks), averaged over the TuORFs of the group.") +
  sch_leader(97, 2) +   # no CAGE ticks: the model uses only the two distances
  annotate("rect", xmin = 22, xmax = 40, ymin = 95.2, ymax = 98.8, fill = BLU) +
  annotate("segment", x = c(22, 2), xend = c(XM, XM), y = c(87.5, 76), yend = c(87.5, 76), linewidth = 0.45, colour = "grey20", arrow = arrow(length = unit(1.1, "mm"), ends = "both", type = "closed")) +
  lab_on((22 + XM) / 2, 87.5, "distance to the mORF start codon") + lab_on((2 + XM) / 2, 76, "leader length") +
  note(68, "Blue, predicted from the real position: a logistic model with only these\ntwo values, fitted to all TuORFs of D, gives each TuORF a probability of\nexclusion. The bar is the mean probability. Over all groups, blue and red\nhave the same mean; red above blue means more exclusion than predicted.") +
  note(17, "Red is lower than gray in every group because a uORF that is present in\nfew mRNAs is less likely to be detected by Ribo-seq.") +
  labs(title = "How the three bars of D are calculated") +
  theme(plot.title = element_text(family = FIG_FONT, size = FIG_PT_MIN + 2, face = "bold"), plot.title.position = "plot", plot.margin = margin(4, 10, 6, 8))
# V30 (2026-11-09, Larry): legends inside the plots or at the side; the panels of a row share their plotting area (fig_grid)
# V49 (2026-11-17): panels in the order of their first citation in the text: A = interval (pB), B = the two calls by sequencing depth (pSat),
# C = genes with two TuORFs (pA, half width as before) with the schematic of the two nulls at its right (pNull, not aligned with C), D = length (half width, short keys
# inside) with the explanation of "expected" and "predicted" at its right (pExp), E = two logistic models (pF)
cS6 <- c(lapply(list(pB, pSat, pA), compact), list(patchwork::free(pNull), compact(pD), patchwork::free(pExp), compact_wide(pF)))
cS6[[2]] <- legend_in(cS6[[2]], 0.02, 0.98) + theme(legend.box = "horizontal")
cS6[[3]] <- legend_in(cS6[[3]], 0.99, 0.03); cS6[[5]] <- legend_in(cS6[[5]], 0.02, 0.98); cS6[[7]] <- legend_in(cS6[[7]], 0.01, 0.98)
g <- fig_grid(cS6, design = "AAAAABBBBB\nCCCCCDDDDD\nEEEEFFFFFF\nGGGGGGGGGG", heights = c(0.95, 1.4, 1.75, 1.0),   # V57: D 40% and its explanation 60% of the width, more height for this row
              top = fig_title("Figure S3. Null models and uncertainty for the TSS effects"))   # C, D and E need the full width
fig_save(file.path(OUT_DIR, "A6_Figure_nulls_and_uncertainty.pdf"), g, h = 9)
.add(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A6_summary.txt"))
