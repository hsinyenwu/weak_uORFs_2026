################################################################################
# A25 (2026-11-06, Larry) — Leader-splicing audit for translated overlapping uORFs (ouORFs)
#
# Same question and method as A2, applied to the 965 CiPS ouORFs (Files/CiPS_ouORF_all_Final.xlsx). An ouORF starts in the
# 5' leader and ends inside the main ORF, so only its start codon can be removed by leader splicing; the stop codon lies in the
# coding sequence and is not audited.
#   (1) ANNOTATION (Araport11+CTRL isoforms): does the start codon lie inside an intron of another isoform of the same gene
#       (class D)? Does the gene carry leader introns elsewhere (B)? None (A)? Introns inside the ouORF body (C1, C2) are
#       recorded for completeness; for an ouORF the body includes coding-sequence introns, so these classes do not mean that
#       the start codon is affected.
#   (2) RNA-seq JUNCTIONS at the start codon (same merged control RNA-seq as A2, unique sense reads): R = reads whose aligned
#       block covers the first nucleotide of the start codon, E = reads whose splice gap spans it; S_TIS = R / (R + E) when
#       R + E >= MIN_COV; spliced over when S_TIS < S_CLEAN. Effective inclusion = %uORF+ x S_TIS, with %uORF+ per ouORF from
#       output/Table_ouORFs_pUP.csv (CAGE-seq and TSS-seq, canonical script).
# Outputs (OUT_DIR, prefix A25_): per-ouORF table (csv + xlsx with a column legend), summary txt.
#
# Run from Weak_uORF_submission/:  Rscript revision_modules/A25_ouORF_leader_splicing_audit.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomicAlignments); library(Rsamtools)
  library(data.table); library(dplyr); library(openxlsx); library(rtracklayer)
})
options(stringsAsFactors = FALSE, warn = 1)
t0 <- Sys.time()

DATA_DIR <- "Files"
OUT_DIR  <- "output"
PATH_GTF  <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_GTF2 <- file.path(DATA_DIR, "Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf")
PATH_BAM  <- "/Users/wu/Desktop/CTRL_v1/RNA_CTRL_merged.bam"          # STAR 2.6.0c, D1/D2/D3/D5/D6/D7 merged, paired, dUTP (read 2 = sense)
PATH_OU   <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_PUP  <- file.path(OUT_DIR, "Table_ouORFs_pUP.csv")
PATH_S1   <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
PATH_A2   <- file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.csv")
MIN_COV   <- 10      # reads needed at the start codon before S is interpreted
S_CLEAN   <- 0.90    # S >= 0.90 = no appreciable splicing over the start codon
stopifnot(file.exists(PATH_GTF), file.exists(PATH_GTF2), file.exists(PATH_BAM), file.exists(PATH_OU), file.exists(PATH_PUP), file.exists(PATH_S1), file.exists(PATH_A2))
dir.create(OUT_DIR, showWarnings = FALSE)
.log <- character(0)
.add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
pct <- function(x, n = length(x)) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * sum(x, na.rm = TRUE) / n)

################################################################################
# 1) Annotation: introns of coding isoforms, CDS starts, leader introns
################################################################################
.add("A25: building TxDb from ", PATH_GTF)
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else
  GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
tx_all  <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx2gene <- setNames(as.character(unlist(mcols(tx_all)$gene_id)), mcols(tx_all)$tx_name)
cds_by  <- cdsBy(txdb, by = "tx", use.names = TRUE)
intr_by <- intronsByTranscript(txdb, use.names = TRUE)
intr_by <- intr_by[names(intr_by) %in% names(cds_by)]          # coding isoforms only
cds_start <- vapply(seq_along(cds_by), function(i) {
  g <- cds_by[[i]]; if (as.character(strand(g)[1]) == "+") min(start(g)) else max(end(g))
}, numeric(1)); names(cds_start) <- names(cds_by)
intr <- unlist(intr_by, use.names = TRUE)
intr$tx_id   <- names(intr); names(intr) <- NULL
intr$gene_id <- unname(tx2gene[intr$tx_id])
intr$cds_start <- unname(cds_start[intr$tx_id])
intr$leader  <- ifelse(as.character(strand(intr)) == "+", end(intr) < intr$cds_start, start(intr) > intr$cds_start)
expr_tx <- tryCatch({ g2 <- rtracklayer::import(PATH_GTF2); unique(as.character(mcols(g2)$transcript_id)) }, error = function(e) character(0))
.add(sprintf("  introns in coding isoforms: %d | in 5' leaders: %d | expressed-isoform GTF: %d transcripts", length(intr), sum(intr$leader), length(expr_tx)))

################################################################################
# 2) ouORFs and their %uORF+
################################################################################
o <- read.xlsx(PATH_OU, sheet = 1)
pup <- fread(PATH_PUP); pup[, chr := as.character(chr)]
O <- o %>% transmute(gene_id, tx_id, ORF_id2, chr = as.character(chr), Strand,
                     TIS   = as.integer(ORF_gstart),                      # first nt of the start codon (transcript orientation on both strands)
                     g_min = as.integer(pmin(ORF_gstart, ORF_gstop)), g_max = as.integer(pmax(ORF_gstart, ORF_gstop)),
                     AA = as.integer(AA), startcodon_identical_to_nuORF = as.logical(startcodon_identical_to_nuORF))
stopifnot(all(ifelse(O$Strand == "+", O$TIS == O$g_min, O$TIS == O$g_max)))
O <- O %>% left_join(as.data.frame(pup) %>% transmute(chr, g_min = g_start, g_max = g_end, Strand, uorf_gene, pct_uORF_plus = pUP_CAGE, pUP_TSSseq = pUP_TSS),
                     by = c("chr", "g_min", "g_max", "Strand"))
stopifnot(!any(is.na(O$pct_uORF_plus)), all(O$gene_id == O$uorf_gene)); O$uorf_gene <- NULL
O$uid <- paste(O$gene_id, paste(O$chr, O$g_min, O$g_max, sep = "_"), sep = "|")
O$n_isoform_calls <- as.integer(table(O$uid)[O$uid])            # the same ORF (same start and stop codon) called on two isoforms of the gene
n_calls <- nrow(O); O <- O[!duplicated(O$uid), ]                 # one row per unique ouORF, as the %uORF+ table and the TuORF audit (A2) do
.add(sprintf("  ouORF calls: %d -> unique ouORFs (start and stop codon): %d (%d called on two isoforms)", n_calls, nrow(O), sum(O$n_isoform_calls > 1)))
stopifnot(!any(duplicated(O$uid)))
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)")
O$start_shared_with_TuORF <- paste(O$chr, O$TIS, O$Strand) %in% paste(as.character(S1A$chr), S1A$uORF_start_genomic, S1A$Strand)   # same start codon listed as a TuORF on another isoform
O$gene_has_TuORF <- O$gene_id %in% S1A$gene_id
.add(sprintf("  ouORFs: %d in %d genes | genes that also have a TuORF: %d | ouORFs whose start codon is also a listed TuORF start codon (isoform-dependent stop): %d",
             nrow(O), n_distinct(O$gene_id), n_distinct(O$gene_id[O$gene_has_TuORF]), sum(O$start_shared_with_TuORF)))

# --- 2a) annotation flags per ouORF -------------------------------------------
O_gr <- GRanges(O$chr, IRanges(O$g_min, O$g_max), strand = O$Strand)
ov <- findOverlaps(O_gr, intr, ignore.strand = FALSE)
pair <- data.table(i = queryHits(ov), j = subjectHits(ov))
pair[, `:=`(same_gene = O$gene_id[i] == intr$gene_id[j],
            same_tx   = O$tx_id[i] == intr$tx_id[j],
            tis_inside = O$TIS[i] >= start(intr)[j] & O$TIS[i] <= end(intr)[j],
            within_body = start(intr)[j] > O$g_min[i] & end(intr)[j] < O$g_max[i],
            other_expressed = intr$tx_id[j] %in% expr_tx)]
pair <- pair[same_gene == TRUE]
flag_by_o <- pair[, .(
  src_intron_in_body     = any(same_tx & within_body),               # ouORF body spans an intron of its own isoform (start codon kept)
  other_TIS_in_intron    = any(!same_tx & tis_inside),               # another isoform splices the start codon out
  other_TIS_in_intron_ex = any(!same_tx & tis_inside & other_expressed),
  other_body_overlap     = any(!same_tx & !tis_inside),              # another isoform's intron overlaps the body (start codon kept)
  other_body_overlap_ex  = any(!same_tx & !tis_inside & other_expressed)
), by = i]
for (k in c("src_intron_in_body", "other_TIS_in_intron", "other_TIS_in_intron_ex", "other_body_overlap", "other_body_overlap_ex")) O[[k]] <- FALSE
O[flag_by_o$i, c("src_intron_in_body", "other_TIS_in_intron", "other_TIS_in_intron_ex", "other_body_overlap", "other_body_overlap_ex")] <-
  flag_by_o[, .(src_intron_in_body, other_TIS_in_intron, other_TIS_in_intron_ex, other_body_overlap, other_body_overlap_ex)]
leader_genes <- unique(intr$gene_id[intr$leader]); leader_tx <- unique(intr$tx_id[intr$leader])
O$gene_has_leader_intron <- O$gene_id %in% leader_genes
O$src_has_leader_intron  <- O$tx_id %in% leader_tx
O$n_isoforms <- as.integer(table(tx2gene[names(cds_by)])[O$gene_id]); O$n_isoforms[is.na(O$n_isoforms)] <- 1L
O$annot_class <- with(O, case_when(
  other_TIS_in_intron ~ "D. start codon inside an intron of another isoform",
  src_intron_in_body  ~ "C1. ouORF body spans an intron of its own isoform (start codon kept)",
  other_body_overlap  ~ "C2. intron of another isoform overlaps the ouORF body (start codon kept)",
  gene_has_leader_intron ~ "B. leader intron(s) in the gene, start codon untouched",
  TRUE ~ "A. no annotated leader intron in any isoform"))
O$annot_class_ex <- with(O, case_when(other_TIS_in_intron_ex ~ "D", src_intron_in_body ~ "C1", other_body_overlap_ex ~ "C2", gene_has_leader_intron ~ "B", TRUE ~ "A"))
O$annot_TIS_affected <- O$other_TIS_in_intron                       # only class D removes the start codon

################################################################################
# 3) RNA-seq junction evidence at the start codon
################################################################################
.add("A25: reading RNA-seq alignments at ", nrow(O), " start codons from ", basename(PATH_BAM))
pos_gr <- GRanges(O$chr, IRanges(O$TIS, O$TIS), strand = O$Strand, uid = O$uid)
bam_sl <- names(scanBamHeader(PATH_BAM)[[1]]$targets)
stopifnot(all(as.character(seqnames(pos_gr)) %in% bam_sl))
seqlevels(pos_gr) <- bam_sl
which_gr <- reduce(resize(pos_gr, width = 401, fix = "center"), ignore.strand = TRUE)   # merge nearby sites so no read is read twice
which_gr <- restrict(which_gr, start = 1L)
param <- ScanBamParam(which = which_gr, what = c("flag", "qname"), mapqFilter = 255,
                      flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE))
gal <- readGAlignments(PATH_BAM, param = param)
gal <- gal[!duplicated(paste(mcols(gal)$qname, mcols(gal)$flag, start(gal), cigar(gal)))]   # an alignment fetched by two adjacent windows is kept once
.add(sprintf("  alignments read: %d (unique mappers; %.1f%% spliced)", length(gal), 100 * mean(grepl("N", cigar(gal)))))
fl <- mcols(gal)$flag
first_mate <- bitwAnd(fl, 64L) > 0L
aln_strand <- as.character(strand(gal))
sense <- ifelse(first_mate, ifelse(aln_strand == "+", "-", "+"), aln_strand)   # dUTP: read 2 is sense, read 1 antisense
blk <- grglist(gal); jnc <- junctions(gal)
pos_strand <- as.character(strand(pos_gr))
cnt <- function(h) { ok <- sense[subjectHits(h)] == pos_strand[queryHits(h)]; tabulate(queryHits(h)[ok], nbins = length(pos_gr)) }
O$R_TIS <- cnt(findOverlaps(pos_gr, blk, ignore.strand = TRUE))
# V23 (2026-11-08): splice gaps that are annotated introns of another gene only are not splicing over this ORF (as in A2)
gaps <- unlist(jnc, use.names = FALSE); gap_read <- togroup(PartitioningByEnd(jnc))
hg <- findOverlaps(pos_gr, gaps, ignore.strand = TRUE); ok_g <- sense[gap_read[subjectHits(hg)]] == pos_strand[queryHits(hg)]
gk <- paste(as.character(seqnames(gaps))[subjectHits(hg)], start(gaps)[subjectHits(hg)], end(gaps)[subjectHits(hg)])
ik <- paste(as.character(seqnames(intr)), start(intr), end(intr)); ikg <- paste(ik, intr$gene_id)
other_gene_gap <- (gk %in% ik) & !(paste(gk, O$gene_id[queryHits(hg)]) %in% ikg)
O$E_TIS <- tabulate(queryHits(hg)[ok_g & !other_gene_gap], nbins = length(pos_gr))
.add(sprintf("  splice gaps spanning a start codon: %d sense gaps; %d are annotated introns of another gene only and were not counted (%d ORFs affected)", sum(ok_g), sum(ok_g & other_gene_gap), length(unique(queryHits(hg)[ok_g & other_gene_gap]))))
O$AS_TIS <- { h <- findOverlaps(pos_gr, blk, ignore.strand = TRUE); ok <- sense[subjectHits(h)] != pos_strand[queryHits(h)]; tabulate(queryHits(h)[ok], nbins = length(pos_gr)) }
rm(gal, blk, jnc); invisible(gc())
O <- O %>% mutate(
  cov_TIS = R_TIS + E_TIS,
  S_TIS = ifelse(cov_TIS >= MIN_COV, R_TIS / cov_TIS, NA_real_),
  junction_status = case_when(
    cov_TIS < MIN_COV ~ "unassessable (<10 reads at the start codon)",
    S_TIS < S_CLEAN ~ "spliced over (S < 0.90 at the start codon)",
    TRUE ~ "retained (S >= 0.90 at the start codon)"),
  TIS_spliced = !is.na(S_TIS) & S_TIS < S_CLEAN,
  pUP_eff_CAGE   = ifelse(is.na(S_TIS), NA_real_, pct_uORF_plus * S_TIS),
  pUP_eff_TSSseq = ifelse(is.na(S_TIS), NA_real_, pUP_TSSseq * S_TIS),
  proxy_status = case_when(
    annot_TIS_affected | TIS_spliced ~ "start codon affected by splicing",
    junction_status == "unassessable (<10 reads at the start codon)" ~ "annotation-clean, junctions unassessable",
    TRUE ~ "clean (annotation + junctions)"))

################################################################################
# 4) Summaries, with the TuORF audit (A2) for comparison
################################################################################
A2 <- fread(PATH_A2)
.add(""); .add("=== 4a. Annotation view (all ouORFs; Araport11+CTRL isoforms) ===")
tabA <- O %>% count(annot_class) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabA))) .add(sprintf("  %-78s %5d (%4.1f%%)", tabA$annot_class[k], tabA$n[k], tabA$pct[k]))
.add(sprintf("  ouORFs in genes with >=1 leader intron in any isoform: %s | source isoform itself has a leader intron: %s", pct(O$gene_has_leader_intron), pct(O$src_has_leader_intron)))
.add(sprintf("  start codon inside an intron of another isoform (D): %s; counting only alternative isoforms in the expressed 2-isoform set: %s  [TuORFs, A2: %s]",
             pct(O$annot_TIS_affected), pct(O$annot_class_ex == "D"), pct(grepl("^D", A2$annot_class), nrow(A2))))
.add(""); .add("=== 4b. RNA-seq junction view at the start codon (control RNA-seq, unique sense reads) ===")
.add(sprintf("  coverage at the start codon: median %d reads; sites with >= %d reads: %s  [TuORFs: median %d; %s]",
             median(O$cov_TIS), MIN_COV, pct(O$cov_TIS >= MIN_COV), median(A2$cov_TIS), pct(A2$cov_TIS >= MIN_COV, nrow(A2))))
.add(sprintf("  antisense reads discarded: median %d (%.1f%% of all reads at the sites)", median(O$AS_TIS), 100 * sum(O$AS_TIS) / sum(O$AS_TIS + O$cov_TIS)))
hasT <- !is.na(O$S_TIS); hasA <- !is.na(A2$S_TIS)
.add(sprintf("  start codons with >= %d reads (n = %d): no spliced read at all %s | S >= 0.90 %s | S < 0.50 %s | median S = %.3f", MIN_COV, sum(hasT),
             pct(O$E_TIS[hasT] == 0, sum(hasT)), pct(O$S_TIS[hasT] >= S_CLEAN, sum(hasT)), pct(O$S_TIS[hasT] < 0.5, sum(hasT)), median(O$S_TIS[hasT])))
.add(sprintf("  [TuORFs, A2 (n = %d): no spliced read %s | S >= 0.90 %s | S < 0.50 %s]", sum(hasA),
             pct(A2$E_TIS[hasA] == 0, sum(hasA)), pct(A2$S_TIS[hasA] >= S_CLEAN, sum(hasA)), pct(A2$S_TIS[hasA] < 0.5, sum(hasA))))
tabJ <- O %>% count(junction_status) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabJ))) .add(sprintf("  %-52s %5d (%4.1f%%)", tabJ$junction_status[k], tabJ$n[k], tabJ$pct[k]))
xt <- table(substr(O$annot_class, 1, 2), O$junction_status)
.add("  junction status by annotation class (rows = class; columns = retained / spliced over / unassessable):")
for (r in rownames(xt)) .add(sprintf("    %-4s retained =%d  spliced o=%d  unassessa=%d", r, xt[r, grep("^retained", colnames(xt))], xt[r, grep("^spliced", colnames(xt))], xt[r, grep("^unassess", colnames(xt))]))
.add(sprintf("  ouORFs with annotation-clean leaders (A/B) but junction evidence of splicing over the start codon: %d", sum(O$TIS_spliced & substr(O$annot_class, 1, 1) %in% c("A", "B"))))
.add(""); .add("=== 4c. Scope of the TSS proxy ===")
tabP <- O %>% count(proxy_status) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabP))) .add(sprintf("  ouORFs: %-48s %5d (%4.1f%%)", tabP$proxy_status[k], tabP$n[k], tabP$pct[k]))
.add(sprintf("  genes with >=1 ouORF whose start codon is affected by splicing: %s of %d ouORF genes",
             pct(tapply(O$proxy_status == "start codon affected by splicing", O$gene_id, any), n_distinct(O$gene_id)), n_distinct(O$gene_id)))
.add(""); .add("=== 4d. Effective inclusion (CAGE %uORF+ x S at the start codon; ouORFs with >= 10 reads) ===")
.add(sprintf("  mean %%uORF+ %.1f -> mean effective inclusion %.1f (n = %d); ouORFs with >= 50%% uORF+ %s -> with effective inclusion >= 50%% %s",
             mean(O$pct_uORF_plus[hasT]), mean(O$pUP_eff_CAGE[hasT]), sum(hasT), pct(O$pct_uORF_plus[hasT] >= 50, sum(hasT)), pct(O$pUP_eff_CAGE[hasT] >= 50, sum(hasT))))
.add(sprintf("  ouORFs where splicing removes the start codon from >= 50%% of the transcripts (S < 0.50): %s of %d covered | TSS selection removes the ouORF from >= 50%% (%%uORF+ <= 50): %s of %d (all ouORFs)",
             pct(O$S_TIS[hasT] < 0.5, sum(hasT)), sum(hasT), pct(O$pct_uORF_plus <= 50), nrow(O)))

################################################################################
# 5) Tables
################################################################################
out <- O %>% select(gene_id, tx_id, ORF_id2, n_isoform_calls, chr, Strand, start_codon_genomic = TIS, ouORF_genome_range_min = g_min, ouORF_genome_range_max = g_max, AA,
                    pct_uORF_plus_CAGE = pct_uORF_plus, pct_uORF_plus_TSSseq = pUP_TSSseq, gene_has_TuORF, start_shared_with_TuORF,
                    n_isoforms, gene_has_leader_intron, src_has_leader_intron, other_TIS_in_intron, other_TIS_in_intron_ex, src_intron_in_body, other_body_overlap,
                    annot_class, annot_class_expressed_isoforms = annot_class_ex, annot_TIS_affected,
                    R_TIS, E_TIS, AS_TIS, cov_TIS, S_TIS, TIS_spliced, junction_status, pUP_eff_CAGE, pUP_eff_TSSseq, proxy_status)
fwrite(out, file.path(OUT_DIR, "A25_Table_ouORF_leader_splicing.csv"))
source("revision_modules/A25_tables.R")   # V41 (2026-11-12): the shipped tables are written by their own file, which can be run alone from the saved csv (Rscript revision_modules/A25_tables.R)
################################################################################
# 7) Figure (V21, 2026-11-06, Larry: ouORF results in their own Results section): overlapping uORFs versus TuORFs
#    A: %uORF+ (CAGE) cumulative distributions; B: annotation classes; C: S at the start codon (covered ORFs); D: shares of ORFs
#    that TSS selection or leader splicing removes from more than 0%, at least 10%, 50% and 80% of the transcripts (A19 rule).
################################################################################
source("revision_modules/A25_figure.R")   # V22 (2026-11-07): the figure code lives in its own file so it can be redrawn alone (A25_FIGURE_ONLY=1 Rscript revision_modules/A25_figure.R) without the 83-min BAM scan
.add(sprintf("A25 done in %.1f min: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(OUT_DIR, "A25_Table_ouORF_leader_splicing.csv")))
writeLines(c(sprintf("A25 ouORF leader-splicing audit — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), file.path(OUT_DIR, "A25_summary.txt"))
