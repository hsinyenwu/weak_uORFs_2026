################################################################################
# A7 — ATG frequency across the TSS with composition controls (Reviewer 2, §2.5) #
#                                                                               #
# Fig 3A-B show ATG frequency around the 50%-CAGE position (the TSS reference    #
# point). The reviewer notes that the TSS also marks the intergenic/promoter ->  #
# exon composition shift (Fig S11), so a drop in ATG cannot by itself be         #
# attributed to selection; and that the TuORF vs No-TuORF comparison carries    #
# ascertainment bias (TuORF genes must contain a translated ATG uORF).           #
# Here, for the same transcript windows as Fig 3 (one transcript per gene):      #
#  * ATG versus the five trinucleotides of identical base composition          #
#    (AGT, GAT, GTA, TAG, TGA), position by position;                            #
#  * ATG observed / expected, where the expectation at each position comes from  #
#    the local mononucleotide (zero-order) or dinucleotide (first-order Markov)   #
#    composition of the same sequences;                                          #
#  * each gene set (No-TuORF genes; single-TuORF genes) against its own          #
#    expectation, upstream (-60..-1) versus downstream (0..+58) of the TSS point. #
# Outputs (output/A7_*): per-position tables, summary txt, figure pdf.           #
# V47 (2026-11-16): the figure has panels A-F.                                  #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A7_ATG_composition_controls.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(Rsamtools); library(Biostrings); library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(gridExtra); library(grid)
})
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
set.seed(20261028)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF   <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_CAGE  <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt")
PATH_FASTA <- file.path(DATA_DIR, "TAIR10_chr_all_2.fas")
PATH_S1    <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"); PATH_S2 <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
UPSTREAM <- 100L; HALF <- 60L; B <- 500L
MATCHED <- c("AGT", "GAT", "GTA", "TAG", "TGA")           # same base composition as ATG
for (f in c(PATH_GTF, PATH_CAGE, PATH_FASTA, PATH_S1, PATH_S2, PATH_OUORF)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x); big <- function(x) trimws(format(x, big.mark = ","))

################################################################################
# 1) Windows (1-exon 5' UTR + 100 nt), CAGE 50% point, sequences  (as Fig 3A-B)
################################################################################
.add("A7: TxDb, windows, sequences")
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
utr <- fiveUTRsByTranscript(txdb, use.names = TRUE); utr1 <- utr[elementNROWS(utr) == 1]
w <- unlist(utr1, use.names = TRUE); w <- resize(w, width = width(w) + UPSTREAM, fix = "end")   # strand-aware 5' extension by 100 nt
fa <- FaFile(PATH_FASTA); idx <- scanFaIndex(fa); fa_names <- as.character(seqnames(idx))
lv <- seqlevels(w); map <- setNames(lv, lv)
if (!all(lv %in% fa_names)) { alt <- paste0("Chr", lv); ok <- alt %in% fa_names; map[ok] <- alt[ok] }
keep_lv <- lv[map[lv] %in% fa_names]; w <- GenomeInfoDb::keepSeqlevels(w, keep_lv, pruning.mode = "coarse"); seqlevels(w) <- unname(map[seqlevels(w)])
seqlengths(w) <- seqlengths(idx)[seqlevels(w)]; w0 <- w; w <- trim(w); w <- w[width(w) == width(w0)[match(names(w), names(w0))]]
seqs <- getSeq(fa, w); names(seqs) <- names(w)
.add(sprintf("  transcripts with a single-exon 5' UTR: %s; windows (UTR + %d nt) retrieved from the genome: %s", big(length(utr1)), UPSTREAM, big(length(seqs))))
# CAGE -> window coordinates (single exon: a direct offset), 50% point
cg <- fread(PATH_CAGE, header = FALSE, col.names = c("count", "chr", "pos", "strand")); cg[, chr := gsub("^chr", "", chr, ignore.case = TRUE)]
cg_chr <- unname(map[cg$chr]); cg <- cg[!is.na(cg_chr)]; cg[, chr := cg_chr[!is.na(cg_chr)]]
cgr <- GRanges(cg$chr, IRanges(cg$pos, cg$pos), strand = cg$strand, count = cg$count)
ov <- findOverlaps(cgr, w, ignore.strand = FALSE)
D <- data.table(tx = names(w)[subjectHits(ov)], gpos = start(cgr)[queryHits(ov)], count = cgr$count[queryHits(ov)], s = start(w)[subjectHits(ov)], e = end(w)[subjectHits(ov)], str = as.character(strand(w))[subjectHits(ov)])
D[, idx := ifelse(str == "+", gpos - s + 1L, e - gpos + 1L)]
D <- D[, .(count = sum(count)), by = .(tx, idx)]; setorder(D, tx, idx); D[, `:=`(cum = cumsum(count), tot = sum(count)), by = tx]; D <- D[tot > 0.1]
P50 <- D[cum / tot >= 0.5, .SD[1], by = tx][, .(tx, idx50 = idx, tot)]
P50[, wlen := width(w)[match(tx, names(w))]]
P50 <- P50[idx50 >= HALF + 1L & (wlen - idx50) >= HALF + 2L]                    # the +-60 nt window must fit (as in Fig 3)
# gene sets: single-TuORF genes (Fig 3 'TuORF' group) and No-TuORF genes (no detected upstream translation at all)
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)"); S2 <- read.xlsx(PATH_S2, sheet = 1); ou <- read.xlsx(PATH_OUORF)
gid <- function(tx) sub("\\.[0-9]+$", "", tx)
P50[, gene := gid(tx)]; upstream_translated <- unique(c(S1A$gene_id, as.character(ou$gene_id)))
P50[, set := fifelse(gene %in% S2$gene_id, "single-TuORF genes", fifelse(gene %in% upstream_translated, "drop", "No-TuORF genes"))]
P50all <- copy(P50)                                                                            # 2026-11-19: all genes, for panels G-H (section 2b)
P50 <- P50[set != "drop"]
n_tx_all <- nrow(P50); setorder(P50, gene, -tot); P50 <- P50[!duplicated(gene)]                  # one transcript per gene (most CAGE)
.add(sprintf("  transcripts with a 50%% point and a fitting window: %s -> one per gene: %s (No-TuORF %s, single-TuORF %s)", big(n_tx_all), big(nrow(P50)), big(sum(P50$set == "No-TuORF genes")), big(sum(P50$set == "single-TuORF genes"))))
win <- subseq(seqs[P50$tx], start = P50$idx50 - HALF, end = P50$idx50 + HALF + 2L)       # 123 nt: positions -60..+62 (0 = the 50% point)
M <- do.call(rbind, strsplit(as.character(win), "")); storage.mode(M) <- "character"
REL <- (-HALF):(HALF)                                                                         # trinucleotide start positions -60..+60 -> columns 1..121
tri_at <- function(Mx, j) paste0(Mx[, j], Mx[, j + 1L], Mx[, j + 2L])

################################################################################
# 2) Per-position frequencies, composition-matched trinucleotides, expectations
################################################################################
.add("A7: per-position frequencies and expectations")
freq_tables <- function(Mx) {
  n <- nrow(Mx); out <- vector("list", length(REL))
  for (k in seq_along(REL)) { j <- k; tri <- tri_at(Mx, j)
    pA <- mean(Mx[, j] == "A"); pT <- mean(Mx[, j + 1L] == "T"); pG <- mean(Mx[, j + 2L] == "G")
    cAT <- sum(Mx[, j] == "A" & Mx[, j + 1L] == "T"); cTG <- sum(Mx[, j + 1L] == "T" & Mx[, j + 2L] == "G"); nA <- sum(Mx[, j] == "A"); nT <- sum(Mx[, j + 1L] == "T")
    e1 <- pA * (if (nA > 0) cAT / nA else 0) * (if (nT > 0) cTG / nT else 0)                 # first-order Markov expectation for ATG at this position
    e0 <- pA * pT * pG                                                                           # zero-order expectation
    tb <- table(factor(tri, levels = c("ATG", MATCHED)))
    out[[k]] <- data.table(rel = REL[k], n = n, ATG = unname(tb["ATG"]) / n, matched_mean = mean(tb[MATCHED]) / n, AGT = unname(tb["AGT"]) / n, GAT = unname(tb["GAT"]) / n, GTA = unname(tb["GTA"]) / n, TAG = unname(tb["TAG"]) / n, TGA = unname(tb["TGA"]) / n,
                           exp0 = e0, exp1 = e1, A = pA, C = mean(Mx[, j] == "C"), G = mean(Mx[, j] == "G"), T = mean(Mx[, j] == "T")) }
  rbindlist(out) }
summ <- function(ft) { up <- ft[rel <= -1]; dn <- ft[rel >= 0 & rel <= HALF - 2L]
  c(ATG_up = mean(up$ATG), ATG_dn = mean(dn$ATG), matched_up = mean(up$matched_mean), matched_dn = mean(dn$matched_mean),
    oe1_up = mean(up$ATG) / mean(up$exp1), oe1_dn = mean(dn$ATG) / mean(dn$exp1), oe0_up = mean(up$ATG) / mean(up$exp0), oe0_dn = mean(dn$ATG) / mean(dn$exp0),
    AGT_r = mean(dn$AGT) / mean(up$AGT), GAT_r = mean(dn$GAT) / mean(up$GAT), GTA_r = mean(dn$GTA) / mean(up$GTA), TAG_r = mean(dn$TAG) / mean(up$TAG), TGA_r = mean(dn$TGA) / mean(up$TGA),
    A_up = mean(up$A), A_dn = mean(dn$A), T_up = mean(up$T), T_dn = mean(dn$T), G_up = mean(up$G), G_dn = mean(dn$G), C_up = mean(up$C), C_dn = mean(dn$C)) }
derived <- function(v) c(ATG_ratio = unname(v["ATG_dn"] / v["ATG_up"]), matched_ratio = unname(v["matched_dn"] / v["matched_up"]), excess = unname((v["ATG_dn"] / v["ATG_up"]) / (v["matched_dn"] / v["matched_up"])),
                         oe1_up = unname(v["oe1_up"]), oe1_dn = unname(v["oe1_dn"]), oe1_change = unname(v["oe1_dn"] / v["oe1_up"]), oe0_up = unname(v["oe0_up"]), oe0_dn = unname(v["oe0_dn"]), oe0_change = unname(v["oe0_dn"] / v["oe0_up"]))
FT <- list(); SM <- list(); CI <- list()
for (s in c("No-TuORF genes", "single-TuORF genes")) {
  Mx <- M[P50$set == s, , drop = FALSE]; ft <- freq_tables(Mx); ft[, set := s]; FT[[s]] <- ft
  v <- summ(ft); d <- derived(v); SM[[s]] <- c(n = nrow(Mx), v, d)
  bs <- replicate(B, { i <- sample.int(nrow(Mx), nrow(Mx), TRUE); derived(summ(freq_tables(Mx[i, , drop = FALSE]))) })
  CI[[s]] <- apply(bs, 1, quantile, c(0.025, 0.975))
  .add(sprintf("  %s (n = %s genes):", s, big(nrow(Mx))))
  .add(sprintf("    ATG per position: upstream %.2f%% of genes, downstream %.2f%% -> downstream/upstream ratio %.2f (95%% CI %.2f-%.2f)", 100 * v["ATG_up"], 100 * v["ATG_dn"], d["ATG_ratio"], CI[[s]][1, "ATG_ratio"], CI[[s]][2, "ATG_ratio"]))
  .add(sprintf("    composition-matched trinucleotides (mean of %s): upstream %.2f%%, downstream %.2f%% -> ratio %.2f (%.2f-%.2f); individual ratios AGT %.2f, GAT %.2f, GTA %.2f, TAG %.2f, TGA %.2f", paste(MATCHED, collapse = "/"), 100 * v["matched_up"], 100 * v["matched_dn"], d["matched_ratio"], CI[[s]][1, "matched_ratio"], CI[[s]][2, "matched_ratio"], v["AGT_r"], v["GAT_r"], v["GTA_r"], v["TAG_r"], v["TGA_r"]))
  .add(sprintf("    ATG drop relative to the matched trinucleotides: %.2f (%.2f-%.2f)  [1 = ATG falls as much as its composition-matched trinucleotides]", d["excess"], CI[[s]][1, "excess"], CI[[s]][2, "excess"]))
  .add(sprintf("    observed/expected ATG (first-order Markov, local dinucleotides): upstream %.2f, downstream %.2f (%.2f-%.2f); change %.2f (%.2f-%.2f)", d["oe1_up"], d["oe1_dn"], CI[[s]][1, "oe1_dn"], CI[[s]][2, "oe1_dn"], d["oe1_change"], CI[[s]][1, "oe1_change"], CI[[s]][2, "oe1_change"]))
  .add(sprintf("    observed/expected ATG (zero-order, local mononucleotides): upstream %.2f, downstream %.2f; change %.2f (%.2f-%.2f)", d["oe0_up"], d["oe0_dn"], d["oe0_change"], CI[[s]][1, "oe0_change"], CI[[s]][2, "oe0_change"]))
  .add(sprintf("    base composition upstream -> downstream: A %.1f -> %.1f%%, T %.1f -> %.1f%%, G %.1f -> %.1f%%, C %.1f -> %.1f%%", 100 * v["A_up"], 100 * v["A_dn"], 100 * v["T_up"], 100 * v["T_dn"], 100 * v["G_up"], 100 * v["G_dn"], 100 * v["C_up"], 100 * v["C_dn"]))
}
FT <- rbindlist(FT)
# No-TuORF vs single-TuORF: the downstream ATG excess of TuORF genes relative to their own expectation (ascertainment: they must carry a translated ATG uORF)
.add(sprintf("  downstream ATG observed/expected (first-order): No-TuORF %.2f vs single-TuORF %.2f; the TuORF set is ascertained on carrying a translated AUG uORF, so its downstream excess is expected by construction", SM[["No-TuORF genes"]]["oe1_dn"], SM[["single-TuORF genes"]]["oe1_dn"]))
# position-resolved: where does the ATG o/e drop? (10-nt bins)
FT[, bin := cut(rel, seq(-60, 60, 10), include.lowest = TRUE, right = FALSE)]
binned <- FT[, .(ATG = mean(ATG), matched = mean(matched_mean), exp1 = mean(exp1), exp0 = mean(exp0), oe1 = mean(ATG) / mean(exp1), oe0 = mean(ATG) / mean(exp0)), by = .(set, bin)][order(set, bin)]
.add("  ATG observed/expected (first-order) by 10-nt bin, No-TuORF genes: ", paste(sprintf("%s %.2f", binned[set == "No-TuORF genes", as.character(bin)], binned[set == "No-TuORF genes", oe1]), collapse = " | "))
.add("  ATG observed/expected (first-order) by 10-nt bin, single-TuORF genes: ", paste(sprintf("%s %.2f", binned[set == "single-TuORF genes", as.character(bin)], binned[set == "single-TuORF genes", oe1]), collapse = " | "))

################################################################################
# 2b) Two further tests (Larry, 2026-11-19: "Could we add (1) the ATG vs TAC frequency analysis and (2) the non overlapping vs overlapping ATG
#     frequency to support our point further?"; the strand test uses CAT, ATG read on the other strand; exploratory version with all tests:
#     exploratory_20261118_ATG_selection/). Same windows, upstream (-60..-1) and downstream (0..+58) ranges and the same five permutations as above.
#  G  strand: CAT is ATG on the other DNA strand (ribosomes read only the mRNA strand); its permutations are ACT, ATC, CTA, TAC and TCA.
#  H  ORF type: the ORF read from the trinucleotide in its own frame up to the mORF start codon either has an in-frame stop codon in the leader or
#     has none; without a stop it overlaps the mORF (other frame) or is in frame with it. Index = (downstream/upstream of ATG in the class) /
#     (the same for the five permutations in the same class). In-frame ATGs without a stop are not shown: gene models usually start the CDS at the
#     first in-frame ATG, so these ATGs are rare in annotated leaders by definition (counts in the table).
#  Gene sets: No-TuORF genes (as in A-F) and all genes with a 50% point and a fitting window, whatever their uORF status.
################################################################################
.add("A7: strand and ORF-type tests (panels G-H)")
PERM_CAT <- c("ACT", "ATC", "CTA", "TAC", "TCA"); STOPS <- c("TAA", "TAG", "TGA")
setorder(P50all, gene, -tot); P50all <- P50all[!duplicated(gene)]
gsA <- as.character(seqs[P50all$tx]); nA <- length(gsA); npos <- length(REL)
TRI <- matrix(NA_character_, nA, npos); ORFT <- matrix(NA_character_, nA, npos)
for (i in seq_len(nA)) {
  s_ <- gsA[i]; L <- nchar(s_); ch <- strsplit(s_, "", fixed = TRUE)[[1]]
  tri_all <- paste0(ch[1:(L - 2)], ch[2:(L - 1)], ch[3:L]); k <- P50all$idx50[i] + REL
  TRI[i, ] <- tri_all[k]; st <- which(tri_all %in% STOPS)
  for (fr in 0:2) { kk <- which(k %% 3L == fr); if (!length(kk)) next
    sf <- st[st %% 3L == fr]
    nxt <- if (length(sf)) sf[findInterval(k[kk], sf) + 1L] else rep(NA_integer_, length(kk))   # first in-frame stop after the start codon
    nxt[!is.na(nxt) & nxt <= k[kk]] <- NA
    ORFT[i, kk] <- ifelse(!is.na(nxt), "stop codon in the leader", ifelse(((L + 1L) - k[kk]) %% 3L == 0L, "in frame with the mORF", "overlapping the mORF")) } }
iu <- which(REL <= -1); idn <- which(REL >= 0 & REL <= HALF - 2L)
ORFc <- c("stop codon in the leader", "overlapping the mORF", "in frame with the mORF")
IND <- list(ATG = matrix(TRI == "ATG", nA, npos), pATG = matrix(TRI %in% MATCHED, nA, npos), CAT = matrix(TRI == "CAT", nA, npos), pCAT = matrix(TRI %in% PERM_CAT, nA, npos))
for (cc in ORFc) { IND[[paste0("ATG|", cc)]] <- IND$ATG & ORFT == cc; IND[[paste0("pATG|", cc)]] <- IND$pATG & ORFT == cc }
gh_stats <- function(rows) {
  r <- function(nm, k = 1) { sh <- colMeans(IND[[nm]][rows, , drop = FALSE]) / k; mean(sh[idn]) / mean(sh[iu]) }
  out <- c(ATG = r("ATG"), ATG_perm = r("pATG", 5), CAT = r("CAT"), CAT_perm = r("pCAT", 5))
  for (cc in ORFc) out[paste0("index|", cc)] <- r(paste0("ATG|", cc)) / r(paste0("pATG|", cc), 5)
  out["overlap_vs_stop"] <- out["index|overlapping the mORF"] / out["index|stop codon in the leader"]; out }
GH_SETS <- list(`No-TuORF genes` = which(P50all$set == "No-TuORF genes"), `all genes` = seq_len(nA))
GH <- list(); GHN <- list()
for (s in names(GH_SETS)) { rows <- GH_SETS[[s]]; est <- gh_stats(rows)
  bs <- replicate(B, gh_stats(sample(rows, length(rows), TRUE))); ci <- apply(bs, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
  GH[[s]] <- data.table(set = s, n_genes = length(rows), stat = names(est), estimate = unname(est), ci_lo = ci[1, ], ci_hi = ci[2, ])
  GHN[[s]] <- rbindlist(lapply(c("ATG", "CAT", paste0("ATG|", ORFc)), function(nm) data.table(set = s, what = nm, upstream = sum(IND[[nm]][rows, iu]), downstream = sum(IND[[nm]][rows, idn])))) }
GH <- rbindlist(GH); GHN <- rbindlist(GHN)
gv <- function(s, st) GH[set == s & stat == st]
stopifnot(abs(gv("No-TuORF genes", "ATG")$estimate - SM[["No-TuORF genes"]]["ATG_ratio"]) < 1e-9, abs(gv("No-TuORF genes", "ATG_perm")$estimate - SM[["No-TuORF genes"]]["matched_ratio"]) < 1e-9)   # same as panels A-F
for (s in names(GH_SETS)) {
  f <- function(st) { x <- gv(s, st); sprintf("%.2f (%.2f-%.2f)", x$estimate, x$ci_lo, x$ci_hi) }
  .add(sprintf("  %s (n = %s): downstream/upstream ATG %s, its permutations %s; CAT %s, its permutations %s", s, big(length(GH_SETS[[s]])), f("ATG"), f("ATG_perm"), f("CAT"), f("CAT_perm")))
  .add(sprintf("    ATG index (ATG drop / drop of the permutations in the same class): stop codon in the leader %s, overlapping the mORF %s, in frame with the mORF %s; overlapping / stop %s",
               f("index|stop codon in the leader"), f("index|overlapping the mORF"), f("index|in frame with the mORF"), f("overlap_vs_stop")))
  .add("    ATG counts upstream/downstream: ", paste(sprintf("%s %d/%d", GHN[set == s]$what, GHN[set == s]$upstream, GHN[set == s]$downstream), collapse = "; ")) }

################################################################################
# 3) Tables
################################################################################
sum_tab <- rbindlist(lapply(names(SM), function(s) { v <- SM[[s]]; ci <- CI[[s]]
  data.table(set = s, n_genes = v["n"], ATG_pct_upstream = 100 * v["ATG_up"], ATG_pct_downstream = 100 * v["ATG_dn"], ATG_down_up_ratio = v["ATG_ratio"], ATG_ratio_ci_lo = ci[1, "ATG_ratio"], ATG_ratio_ci_hi = ci[2, "ATG_ratio"],
             matched_pct_upstream = 100 * v["matched_up"], matched_pct_downstream = 100 * v["matched_dn"], matched_down_up_ratio = v["matched_ratio"], matched_ratio_ci_lo = ci[1, "matched_ratio"], matched_ratio_ci_hi = ci[2, "matched_ratio"],
             ATG_vs_matched = v["excess"], ATG_vs_matched_ci_lo = ci[1, "excess"], ATG_vs_matched_ci_hi = ci[2, "excess"],
             obs_exp1_upstream = v["oe1_up"], obs_exp1_downstream = v["oe1_dn"], obs_exp1_dn_ci_lo = ci[1, "oe1_dn"], obs_exp1_dn_ci_hi = ci[2, "oe1_dn"], obs_exp1_change = v["oe1_change"], obs_exp1_change_ci_lo = ci[1, "oe1_change"], obs_exp1_change_ci_hi = ci[2, "oe1_change"],
             obs_exp0_upstream = v["oe0_up"], obs_exp0_downstream = v["oe0_dn"], obs_exp0_change = v["oe0_change"], AGT_ratio = v["AGT_r"], GAT_ratio = v["GAT_r"], GTA_ratio = v["GTA_r"], TAG_ratio = v["TAG_r"], TGA_ratio = v["TGA_r"],
             A_up = v["A_up"], A_dn = v["A_dn"], T_up = v["T_up"], T_dn = v["T_dn"], G_up = v["G_up"], G_dn = v["G_dn"], C_up = v["C_up"], C_dn = v["C_dn"]) }))
fwrite(FT[, .(set, rel, n, ATG, matched_mean, AGT, GAT, GTA, TAG, TGA, exp0, exp1, A, C, G, T)], file.path(OUT_DIR, "A7_Table_ATG_per_position.csv"))
fwrite(sum_tab, file.path(OUT_DIR, "A7_Table_ATG_summary.csv")); fwrite(binned, file.path(OUT_DIR, "A7_Table_ATG_binned.csv"))
fwrite(P50[, .(gene, tx, set, idx50, window_len = wlen, cage_total = tot)], file.path(OUT_DIR, "A7_Table_genes_used.csv"))
fwrite(GH, file.path(OUT_DIR, "A7_Table_ATG_strand_ORF_tests.csv")); fwrite(GHN, file.path(OUT_DIR, "A7_Table_ATG_strand_ORF_counts.csv"))   # 2026-11-19, panels G-H
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S7 — ATG frequency across the TSS reference point with composition controls (Figure S12).",
  sprintf("Windows as in Figure 3A-B: transcripts with a single-exon 5' leader, leader + %d nt upstream, 50%%-cumulative-CAGE position = 0, positions -60..+60 (trinucleotide start); one transcript per gene (most CAGE signal). No-TuORF genes = no detected TuORF or overlapping uORF; single-TuORF genes = the Table S3 cohort.", UPSTREAM),
  "per_position: share of genes with the given trinucleotide starting at each position; matched_mean = mean of AGT, GAT, GTA, TAG, TGA (same base composition as ATG); exp0 = expected ATG share from the local mononucleotide frequencies at the three positions; exp1 = expected share from the local dinucleotide frequencies (first-order Markov: P(A) P(T|A) P(G|T)); A/C/G/T = base composition at the position.",
  "summary: upstream = positions -60..-1, downstream = 0..+58; ratios are downstream/upstream; ATG_vs_matched = ATG ratio divided by the matched-trinucleotide ratio; obs_exp = mean observed ATG share / mean expected share (obs_exp1 with exp1, obs_exp0 with exp0); obs_exp1_change and obs_exp0_change = the downstream value divided by the upstream value; _up and _dn in column names = upstream and downstream; 95% confidence intervals by bootstrap over genes (500 resamples; columns ending in ci_lo and ci_hi).",
  "binned: the same by 10-nt bins (oe1 and oe0 = ATG share divided by exp1 and by exp0). genes_used: the transcripts and 50% points (idx50 = position of the 50% point in the window, counted from its 5' end; window_len = length of the window in nt; cage_total = CAGE signal in the window).",
  "strand_ORF_tests (Figure S12G-H; 2026-11-19): downstream/upstream frequency (same windows and ranges) for ATG, the mean of its five permutations (ATG_perm), CAT (ATG on the other strand) and the mean of its permutations ACT, ATC, CTA, TAC and TCA (CAT_perm); index|class = (downstream/upstream of ATG starting an ORF of that class) / (the same for the five permutations at positions of the same class), where the ORF read in the frame of the trinucleotide up to the mORF start codon has a stop codon in the leader, overlaps the mORF (no in-frame stop, other frame) or is in frame with the mORF (no stop; rare by definition, because gene models usually start the CDS at the first in-frame ATG); 95% bootstrap intervals over genes (500 resamples). Gene sets: No-TuORF genes and all genes with a 50% point and a fitting window, whatever their uORF status. strand_ORF_counts: the numbers of ATG and CAT occurrences behind these values.")), colNames = FALSE)
addWorksheet(wb, "summary"); writeData(wb, "summary", sum_tab); addWorksheet(wb, "per_position"); writeData(wb, "per_position", FT[, .(set, rel, n, ATG, matched_mean, AGT, GAT, GTA, TAG, TGA, exp0, exp1, A, C, G, T)])
addWorksheet(wb, "binned"); writeData(wb, "binned", binned); addWorksheet(wb, "genes_used"); writeData(wb, "genes_used", P50[, .(gene, tx, set, idx50, window_len = wlen, cage_total = tot)])
addWorksheet(wb, "strand_ORF_tests"); writeData(wb, "strand_ORF_tests", GH); addWorksheet(wb, "strand_ORF_counts"); writeData(wb, "strand_ORF_counts", GHN)   # 2026-11-19, panels G-H
saveWorkbook(wb, file.path(OUT_DIR, "A7_Table_ATG_composition_controls.xlsx"), overwrite = TRUE)

################################################################################
# 4) Figure
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3))
sm <- function(x, k = 5) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
FTs <- copy(FT); FTs[, `:=`(ATG_s = sm(ATG), matched_s = sm(matched_mean), exp1_s = sm(exp1), exp0_s = sm(exp0)), by = set]
# A/B: ATG vs matched trinucleotides
pAB <- function(s, lab) { d <- FTs[set == s]; dl <- melt(d[, .(rel, ATG = 100 * ATG_s, `composition-matched trinucleotides (mean)` = 100 * matched_s, `expected from local dinucleotides` = 100 * exp1_s)], id.vars = "rel", variable.name = "what", value.name = "pct")
  ggplot(dl[!is.na(pct)], aes(rel, pct, colour = what)) + geom_vline(xintercept = 0, linetype = 3) + geom_line(linewidth = 0.7) + scale_colour_manual(values = c("#D7191C", "grey40", "#2C7BB6"), name = NULL) + guides(colour = guide_legend(nrow = 3)) +
    labs(x = "Position relative to the 50%-CAGE TSS (nt)", y = "% of genes with the trinucleotide at the position", title = lab, subtitle = sprintf("%s (n = %s); 5-nt running mean", s, big(SM[[s]]["n"]))) + th + theme(legend.direction = "vertical") }
pA <- pAB("No-TuORF genes", "A  ATG and matched trinucleotides"); pB <- pAB("single-TuORF genes", "B  The same in single-TuORF genes")
# C: obs/exp across positions (binned)
pC <- ggplot(binned, aes(bin, oe1, colour = set, group = set)) + geom_hline(yintercept = 1, linetype = 3) + geom_vline(xintercept = 6.5, linetype = 3) + geom_line(linewidth = 0.6) + geom_point(size = 1.8) + scale_colour_manual(values = c(`No-TuORF genes` = "grey20", `single-TuORF genes` = "#D7191C"), name = NULL) +
  labs(x = "Position relative to the 50%-CAGE TSS (10-nt bins)", y = "Observed / expected\nATG", title = "C  ATG observed / expected") + th + theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1))
# D: down/up ratios for ATG and each matched trinucleotide
dd <- rbindlist(lapply(names(SM), function(s) { v <- SM[[s]]; data.table(set = s, tri = c("ATG", MATCHED), ratio = c(v["ATG_ratio"], v["AGT_r"], v["GAT_r"], v["GTA_r"], v["TAG_r"], v["TGA_r"])) }))
dd[, tri := factor(tri, levels = c("ATG", MATCHED))]
pD <- ggplot(dd, aes(tri, ratio, fill = set)) + geom_hline(yintercept = 1, linetype = 3) + geom_col(position = position_dodge(0.8), width = 0.75) + scale_fill_manual(values = c(`No-TuORF genes` = "grey40", `single-TuORF genes` = "#D7191C"), name = NULL) +
  labs(x = "Trinucleotide (ATG and its base-composition matches)", y = "Downstream /\nupstream", title = "D  ATG versus matched trinucleotides") + th
# E: obs/exp up vs down with CIs
de <- rbindlist(lapply(names(SM), function(s) { v <- SM[[s]]; ci <- CI[[s]]; data.table(set = s, side = c("upstream", "downstream"), oe = c(v["oe1_up"], v["oe1_dn"]), lo = c(NA, ci[1, "oe1_dn"]), hi = c(NA, ci[2, "oe1_dn"])) }))
de[, side := factor(side, levels = c("upstream", "downstream"))]
pE <- ggplot(de, aes(side, oe, fill = set)) + geom_hline(yintercept = 1, linetype = 3) + geom_col(position = position_dodge(0.8), width = 0.75) + geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.2, na.rm = TRUE) +
  scale_fill_manual(values = c(`No-TuORF genes` = "grey40", `single-TuORF genes` = "#D7191C"), name = NULL) + labs(x = NULL, y = "Observed / expected\nATG", title = "E  Observed / expected, up and down") + th
# F: base composition
dfc <- melt(FTs[, .(set, rel, A, C, G, T)], id.vars = c("set", "rel"), variable.name = "base", value.name = "frac")
pF <- ggplot(dfc, aes(rel, 100 * frac, colour = base)) + geom_vline(xintercept = 0, linetype = 3) + geom_line(linewidth = 0.5) + facet_wrap(~set) + theme(panel.spacing.x = unit(10, "pt")) + scale_colour_manual(values = c(A = "#1A9641", C = "#2C7BB6", G = "#FDAE61", T = "#D7191C"), name = NULL) +
  labs(x = "Position relative to the 50%-CAGE TSS (nt)", y = "% of genes", title = "F  Base composition") + th
# G-H (2026-11-19): strand and ORF type
colGH <- c(`No-TuORF genes` = "grey40", `all genes` = "black")
dG <- GH[stat %in% c("ATG", "ATG_perm", "CAT", "CAT_perm")][, what := factor(c(ATG = "ATG", ATG_perm = "ATG\npermutations", CAT = "CAT", CAT_perm = "CAT\npermutations")[stat],
       levels = c("ATG", "ATG\npermutations", "CAT", "CAT\npermutations"))]
pG <- ggplot(dG, aes(what, estimate, colour = set)) + geom_hline(yintercept = 1, linetype = 3) + geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(0.5), size = 0.25) +
  scale_colour_manual(values = colGH, name = NULL) + labs(x = NULL, y = "Downstream /\nupstream", title = "G  ATG versus CAT") + th
dH <- GH[stat %in% paste0("index|", ORFc[1:2])][, what := factor(sub("index\\|", "", stat), levels = ORFc[1:2], labels = c("stop codon\nin the leader", "overlapping\nthe mORF"))]
pH <- ggplot(dH, aes(what, estimate, colour = set)) + geom_hline(yintercept = 1, linetype = 3) + geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(0.5), size = 0.25) +
  scale_colour_manual(values = colGH, name = NULL) + labs(x = NULL, y = "ATG drop / drop\nof permutations", title = "H  ATG drop by ORF type") + th
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
# V30 (2026-11-09, Larry): the legends of C-F inside the plots (D, E and F: the y axis runs a little higher to make room); A and B keep the legend above the panel (no free room inside); aligned plotting areas (fig_grid)
cS14 <- lapply(list(pA, pB, pC, pD + expand_limits(y = 1.62), pE + expand_limits(y = 1.32), pF + expand_limits(y = 88) + guides(colour = guide_legend(ncol = 2)), pG + expand_limits(y = c(0, 1.25)), pH + expand_limits(y = c(0, 1.25))), compact)   # G-H added 2026-11-19
cS14[[3]] <- legend_in(cS14[[3]], 0.99, 0.37, just = c(1, 0)); cS14[[4]] <- legend_in(cS14[[4]], 0.99, 0.99); cS14[[5]] <- legend_in(cS14[[5]], 0.03, 0.99); cS14[[6]] <- legend_in(cS14[[6]], 0.995, 0.99)
cS14[[7]] <- legend_in(cS14[[7]], 0.99, 0.03, just = c(1, 0)); cS14[[8]] <- legend_in(cS14[[8]], 0.02, 0.03, just = c(0, 0))
g <- fig_grid(cS14, ncol = 2, heights = c(1.25, 1, 1, 0.9), top = fig_title("Figure S12. ATG frequency across the TSS with composition controls"))   # 2026-11-19: row G-H added; the subtitles of C-F (stated in the legend) removed to make room
fig_save(file.path(OUT_DIR, "A7_Figure_ATG_composition_controls.pdf"), g, h = 9)
.add(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A7_summary.txt"))
