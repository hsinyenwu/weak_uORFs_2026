################################################################################
# A2 — Leader-splicing audit for TuORFs (Reviewer 2, §2.1)                     #
#                                                                              #
# Question: how often does splicing of the 5' leader, rather than TSS choice,  #
# decide whether a TuORF is present in the mature mRNA?  Two complementary     #
# views per TuORF:                                                             #
#   (1) ANNOTATION: do any Araport11+CTRL isoforms of the gene carry an intron #
#       in the 5' leader, and does such an intron remove the TuORF start codon #
#       (TIS) or alter the TuORF body/stop?                                    #
#   (2) RNA-seq JUNCTIONS (same seedlings, control RNA-seq, 6 libraries        #
#       merged): at the TIS and at the stop codon, reads whose aligned block   #
#       covers the position (R, retained) vs reads whose splice gap spans it   #
#       (E, spliced over).  S = R / (R + E).  Effective inclusion =            #
#       %uORF+ x S_TIS.                                                        #
# Outputs (OUT_DIR, prefix A2_): per-TuORF table (xlsx + csv), summary txt,    #
# sensitivity table for the single-TuORF cohort, one multi-panel PDF.          #
#                                                                              #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A2_leader_splicing_audit.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomicAlignments); library(Rsamtools)
  library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(rtracklayer)
})
options(stringsAsFactors = FALSE, warn = 1)
t0 <- Sys.time()

DATA_DIR <- "Files"
OUT_DIR  <- "output"
PATH_GTF  <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_GTF2 <- file.path(DATA_DIR, "Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf")
PATH_BAM  <- "/Users/wu/Desktop/CTRL_v1/RNA_CTRL_merged.bam"          # STAR 2.6.0c, D1/D2/D3/D5/D6/D7 merged, paired, dUTP (read 2 = sense)
PATH_S1   <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
PATH_S2   <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
MIN_COV   <- 10      # reads needed at a position before S is interpreted
S_CLEAN   <- 0.90    # S >= 0.90 = no appreciable splicing over the position
stopifnot(file.exists(PATH_GTF), file.exists(PATH_BAM), file.exists(PATH_S1), file.exists(PATH_S2))
dir.create(OUT_DIR, showWarnings = FALSE)
.log <- character(0)
.add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }

################################################################################
# 1) Annotation: introns, CDS starts, leader introns, per gene
################################################################################
.add("A2: building TxDb from ", PATH_GTF)
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
.add(sprintf("  introns in coding isoforms: %d | in 5' leaders: %d (%.1f%%) | coding isoforms with >=1 leader intron: %d of %d | genes: %d of %d",
             length(intr), sum(intr$leader), 100 * mean(intr$leader),
             length(unique(intr$tx_id[intr$leader])), length(cds_by),
             length(unique(intr$gene_id[intr$leader])), length(unique(tx2gene[names(cds_by)]))))

# expressed-isoform set (the 2-isoform GTF used for the transcriptome-impact module)
expr_tx <- tryCatch({ g2 <- rtracklayer::import(PATH_GTF2); unique(as.character(mcols(g2)$transcript_id)) }, error = function(e) character(0))
.add(sprintf("  expressed-isoform GTF: %d transcripts", length(expr_tx)))

################################################################################
# 2) TuORFs (Table S1A = all 7,423 unique TuORFs; S1B for the TSS-seq %uORF+)
################################################################################
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)")
S1B <- read.xlsx(PATH_S1, sheet = "Table S1B (TSS-seq)")
S2  <- read.xlsx(PATH_S2, sheet = 1)
stopifnot(nrow(S1A) == nrow(S1B))
S1A$pUP_TSSseq <- S1B$pct_uORF_plus[match(paste(S1A$gene_id, S1A$uORF_genome_range), paste(S1B$gene_id, S1B$uORF_genome_range))]
U <- S1A %>%
  mutate(chr = as.character(chr),
         g_min = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\1", uORF_genome_range)),
         g_max = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\2", uORF_genome_range)),
         TIS  = as.integer(uORF_start_genomic),           # first nt of the start codon
         STOP = as.integer(uORF_stop_genomic),            # last nt of the stop codon
         uIF_raw = uORF_TL / mORF_TL,
         uid = paste(gene_id, uORF_genome_range, sep = "|"))
stopifnot(!any(duplicated(U$uid)))
.add(sprintf("  TuORFs: %d in %d genes | single-TuORF cohort (Table S2): %d genes", nrow(U), n_distinct(U$gene_id), nrow(S2)))

# --- 2a) annotation flags per TuORF ------------------------------------------
U_gr <- GRanges(U$chr, IRanges(U$g_min, U$g_max), strand = U$Strand, uid = U$uid, gene_id = U$gene_id, tx_id = U$tx_id, TIS = U$TIS)
# intron-TuORF pairs within the same gene
ov <- findOverlaps(U_gr, intr, ignore.strand = FALSE)
pair <- data.table(i = queryHits(ov), j = subjectHits(ov))
pair[, `:=`(same_gene = U$gene_id[i] == intr$gene_id[j],
            same_tx   = U$tx_id[i] == intr$tx_id[j],
            tis_inside = U$TIS[i] >= start(intr)[j] & U$TIS[i] <= end(intr)[j],
            within_body = start(intr)[j] > U$g_min[i] & end(intr)[j] < U$g_max[i],
            other_expressed = intr$tx_id[j] %in% expr_tx)]
pair <- pair[same_gene == TRUE]
flag_by_u <- pair[, .(
  src_intron_in_body     = any(same_tx & within_body),               # TuORF spans an intron of its own isoform
  other_TIS_in_intron    = any(!same_tx & tis_inside),               # another isoform splices the start codon out
  other_TIS_in_intron_ex = any(!same_tx & tis_inside & other_expressed),
  other_body_overlap     = any(!same_tx & !tis_inside),              # another isoform's intron overlaps the body/stop (TIS kept)
  other_body_overlap_ex  = any(!same_tx & !tis_inside & other_expressed)
), by = i]
U$src_intron_in_body <- FALSE; U$other_TIS_in_intron <- FALSE; U$other_TIS_in_intron_ex <- FALSE; U$other_body_overlap <- FALSE; U$other_body_overlap_ex <- FALSE
U[flag_by_u$i, c("src_intron_in_body", "other_TIS_in_intron", "other_TIS_in_intron_ex", "other_body_overlap", "other_body_overlap_ex")] <-
  flag_by_u[, .(src_intron_in_body, other_TIS_in_intron, other_TIS_in_intron_ex, other_body_overlap, other_body_overlap_ex)]

# leader introns anywhere in the gene (any isoform) and in the source isoform upstream of the TuORF
leader_genes <- unique(intr$gene_id[intr$leader])
leader_tx    <- unique(intr$tx_id[intr$leader])
U$gene_has_leader_intron <- U$gene_id %in% leader_genes
U$src_has_leader_intron  <- U$tx_id %in% leader_tx
U$n_isoforms <- as.integer(table(tx2gene[names(cds_by)])[U$gene_id]); U$n_isoforms[is.na(U$n_isoforms)] <- 1L

U$annot_class <- with(U, case_when(
  other_TIS_in_intron ~ "D. start codon inside an intron of another isoform",
  src_intron_in_body  ~ "C1. TuORF spans an intron of its own isoform",
  other_body_overlap  ~ "C2. intron of another isoform alters the TuORF body/stop (start codon kept)",
  gene_has_leader_intron ~ "B. leader intron(s) in the gene, TuORF untouched",
  TRUE ~ "A. no annotated leader intron in any isoform"))
U$annot_class_ex <- with(U, case_when(          # same, counting only alternative isoforms in the expressed set
  other_TIS_in_intron_ex ~ "D",
  src_intron_in_body ~ "C1",
  other_body_overlap_ex ~ "C2",
  gene_has_leader_intron ~ "B",
  TRUE ~ "A"))
U$annot_affected <- U$annot_class %in% c("D. start codon inside an intron of another isoform",
                                         "C1. TuORF spans an intron of its own isoform",
                                         "C2. intron of another isoform alters the TuORF body/stop (start codon kept)")

################################################################################
# 3) RNA-seq junction evidence at the TIS and the stop codon
################################################################################
.add("A2: reading RNA-seq alignments at ", 2 * nrow(U), " positions from ", basename(PATH_BAM))
pos_gr <- c(GRanges(U$chr, IRanges(U$TIS, U$TIS), strand = U$Strand, uid = U$uid, site = "TIS"),
            GRanges(U$chr, IRanges(U$STOP, U$STOP), strand = U$Strand, uid = U$uid, site = "STOP"))
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
blk <- grglist(gal)          # aligned blocks
jnc <- junctions(gal)        # splice gaps (N)
pos_strand <- as.character(strand(pos_gr))
cnt <- function(h) { ok <- sense[subjectHits(h)] == pos_strand[queryHits(h)]; tabulate(queryHits(h)[ok], nbins = length(pos_gr)) }
R_cnt <- cnt(findOverlaps(pos_gr, blk, ignore.strand = TRUE))
# V23 (2026-11-08): a splice gap that is an annotated intron of ANOTHER gene (and not an intron of the TuORF's gene) is read-through of that gene,
# not splicing of this one (six TuORFs sit in expressed introns of neighbouring genes); such gaps are left out of E
gaps <- unlist(jnc, use.names = FALSE); gap_read <- togroup(PartitioningByEnd(jnc))
hg <- findOverlaps(pos_gr, gaps, ignore.strand = TRUE); ok_g <- sense[gap_read[subjectHits(hg)]] == pos_strand[queryHits(hg)]
gk <- paste(as.character(seqnames(gaps))[subjectHits(hg)], start(gaps)[subjectHits(hg)], end(gaps)[subjectHits(hg)])
ik <- paste(as.character(seqnames(intr)), start(intr), end(intr)); ikg <- paste(ik, intr$gene_id)
own_gene <- U$gene_id[(queryHits(hg) - 1L) %% nrow(U) + 1L]
other_gene_gap <- (gk %in% ik) & !(paste(gk, own_gene) %in% ikg)
E_cnt <- tabulate(queryHits(hg)[ok_g & !other_gene_gap], nbins = length(pos_gr))
.add(sprintf("  splice gaps spanning a site: %d sense gaps; %d (%.1f%%) are annotated introns of another gene only and were not counted as splicing over the TuORF (%d sites affected)", sum(ok_g), sum(ok_g & other_gene_gap), 100 * mean(other_gene_gap[ok_g]), length(unique(queryHits(hg)[ok_g & other_gene_gap]))))
# antisense coverage for the record (how much signal was discarded by the strand rule)
A_cnt <- { h <- findOverlaps(pos_gr, blk, ignore.strand = TRUE); ok <- sense[subjectHits(h)] != pos_strand[queryHits(h)]; tabulate(queryHits(h)[ok], nbins = length(pos_gr)) }
rm(gal, blk, jnc); invisible(gc())
site <- pos_gr$site
U$R_TIS <- R_cnt[site == "TIS"]; U$E_TIS <- E_cnt[site == "TIS"]; U$AS_TIS <- A_cnt[site == "TIS"]
U$R_STOP <- R_cnt[site == "STOP"]; U$E_STOP <- E_cnt[site == "STOP"]; U$AS_STOP <- A_cnt[site == "STOP"]
U <- U %>% mutate(
  cov_TIS = R_TIS + E_TIS, cov_STOP = R_STOP + E_STOP,
  S_TIS  = ifelse(cov_TIS  >= MIN_COV, R_TIS  / cov_TIS,  NA_real_),
  S_STOP = ifelse(cov_STOP >= MIN_COV, R_STOP / cov_STOP, NA_real_),
  S_min  = pmin(S_TIS, S_STOP, na.rm = FALSE),
  junction_status = case_when(
    cov_TIS < MIN_COV & cov_STOP < MIN_COV ~ "unassessable (<10 reads at both sites)",
    (!is.na(S_TIS) & S_TIS < S_CLEAN) | (!is.na(S_STOP) & S_STOP < S_CLEAN) ~ "spliced over (S < 0.90 at TIS or stop)",
    TRUE ~ "retained (S >= 0.90 at every covered site)"),
  TIS_spliced = !is.na(S_TIS) & S_TIS < S_CLEAN,
  STOP_spliced = !is.na(S_STOP) & S_STOP < S_CLEAN,
  pUP_eff_CAGE = ifelse(is.na(S_TIS), NA_real_, pct_uORF_plus * S_TIS),
  pUP_eff_TSSseq = ifelse(is.na(S_TIS), NA_real_, pUP_TSSseq * S_TIS),
  proxy_status = case_when(
    annot_affected | junction_status == "spliced over (S < 0.90 at TIS or stop)" ~ "affected by splicing",
    junction_status == "unassessable (<10 reads at both sites)" ~ "annotation-clean, junctions unassessable",
    TRUE ~ "clean (annotation + junctions)"))

################################################################################
# 4) Summaries
################################################################################
pct <- function(x, n = length(x)) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * sum(x, na.rm = TRUE) / n)
.add(""); .add("=== 4a. Annotation view (all unique TuORFs; Araport11+CTRL isoforms) ===")
tabA <- U %>% count(annot_class) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabA))) .add(sprintf("  %-78s %5d (%4.1f%%)", tabA$annot_class[k], tabA$n[k], tabA$pct[k]))
.add(sprintf("  TuORFs in genes with >=1 leader intron in any isoform: %s | source isoform itself has a leader intron: %s",
             pct(U$gene_has_leader_intron), pct(U$src_has_leader_intron)))
.add(sprintf("  annotation-affected TuORFs (C1 + C2 + D): %s; counting only alternative isoforms in the expressed 2-isoform set: %s",
             pct(U$annot_affected), pct(U$annot_class_ex %in% c("C1", "C2", "D"))))
.add(sprintf("  TIS retained but body/stop altered by annotated splicing (C1 + C2): %s", pct(U$annot_class %in% c("C1. TuORF spans an intron of its own isoform", "C2. intron of another isoform alters the TuORF body/stop (start codon kept)"))))

.add(""); .add("=== 4b. RNA-seq junction view (control RNA-seq, unique sense reads) ===")
.add(sprintf("  coverage at the TIS: median %d reads; sites with >= %d reads: %s | at the stop: median %d; >= %d: %s",
             median(U$cov_TIS), MIN_COV, pct(U$cov_TIS >= MIN_COV), median(U$cov_STOP), MIN_COV, pct(U$cov_STOP >= MIN_COV)))
.add(sprintf("  antisense reads discarded at the TIS: median %d (%.1f%% of all reads at covered sites)", median(U$AS_TIS), 100 * sum(U$AS_TIS) / sum(U$AS_TIS + U$cov_TIS)))
hasT <- !is.na(U$S_TIS)
.add(sprintf("  TIS with >= %d reads (n = %d): no spliced read at all %s | S >= 0.90 %s | S < 0.50 %s | median S = %.3f",
             MIN_COV, sum(hasT), pct(U$E_TIS[hasT] == 0, sum(hasT)), pct(U$S_TIS[hasT] >= S_CLEAN, sum(hasT)), pct(U$S_TIS[hasT] < 0.5, sum(hasT)), median(U$S_TIS[hasT])))
hasS <- !is.na(U$S_STOP)
.add(sprintf("  stop with >= %d reads (n = %d): no spliced read %s | S >= 0.90 %s | S < 0.50 %s",
             MIN_COV, sum(hasS), pct(U$E_STOP[hasS] == 0, sum(hasS)), pct(U$S_STOP[hasS] >= S_CLEAN, sum(hasS)), pct(U$S_STOP[hasS] < 0.5, sum(hasS))))
tabJ <- U %>% count(junction_status) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabJ))) .add(sprintf("  %-50s %5d (%4.1f%%)", tabJ$junction_status[k], tabJ$n[k], tabJ$pct[k]))
.add("  junction status by annotation class (rows = class; columns = retained / spliced over / unassessable):")
xt <- table(substr(U$annot_class, 1, 2), U$junction_status)
for (r in rownames(xt)) .add(sprintf("    %-4s %s", r, paste(sprintf("%s=%d", substr(colnames(xt), 1, 9), xt[r, ]), collapse = "  ")))
.add(sprintf("  TuORFs with annotation-clean leaders (A/B) but junction evidence of splicing over the TIS or stop: %d (unannotated or CTRL-only events)",
             sum(!U$annot_affected & U$junction_status == "spliced over (S < 0.90 at TIS or stop)")))

.add(""); .add("=== 4c. Scope of the TSS proxy ===")
tabP <- U %>% count(proxy_status) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabP))) .add(sprintf("  TuORFs: %-45s %5d (%4.1f%%)", tabP$proxy_status[k], tabP$n[k], tabP$pct[k]))
gene_status <- U %>% group_by(gene_id) %>% summarise(st = if (any(proxy_status == "affected by splicing")) "affected by splicing" else if (all(proxy_status == "clean (annotation + junctions)")) "clean (annotation + junctions)" else "annotation-clean, junctions unassessable", .groups = "drop")
tabG <- gene_status %>% count(st) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabG))) .add(sprintf("  genes:  %-45s %5d (%4.1f%%)", tabG$st[k], tabG$n[k], tabG$pct[k]))
.add(sprintf("  => TSS proxy not contradicted by splicing (clean + unassessable): %.1f%% of TuORFs, %.1f%% of TuORF genes; strictly clean: %.1f%% of TuORFs, %.1f%% of genes",
             100 * mean(U$proxy_status != "affected by splicing"), 100 * mean(gene_status$st != "affected by splicing"),
             100 * mean(U$proxy_status == "clean (annotation + junctions)"), 100 * mean(gene_status$st == "clean (annotation + junctions)")))

# effective inclusion
hasE <- !is.na(U$pUP_eff_CAGE)
.add(""); .add("=== 4d. Effective inclusion = %uORF+ (CAGE) x S_TIS, TuORFs with >= 10 reads at the TIS ===")
.add(sprintf("  n = %d | median %%uORF+ %.1f -> median effective %.1f | TuORFs moving from >= 80%% to < 80%%: %d | from > 0 to effective < 30%%: %d | Pearson r(%%uORF+, effective) = %.3f",
             sum(hasE), median(U$pct_uORF_plus[hasE]), median(U$pUP_eff_CAGE[hasE]),
             sum(U$pct_uORF_plus[hasE] >= 80 & U$pUP_eff_CAGE[hasE] < 80), sum(U$pct_uORF_plus[hasE] >= 30 & U$pUP_eff_CAGE[hasE] < 30),
             cor(U$pct_uORF_plus[hasE], U$pUP_eff_CAGE[hasE])))

################################################################################
# 5) Sensitivity: single-TuORF cohort (Table S3) on the clean set
################################################################################
C <- S2 %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|")) %>%
  left_join(U %>% select(uid, annot_class, annot_affected, junction_status, proxy_status, S_TIS, S_STOP, cov_TIS, pUP_eff_CAGE, pUP_eff_TSSseq), by = "uid") %>%
  mutate(uIF_raw = uORF_TL / mORF_TL, TE = mORF_TE)
stopifnot(!any(is.na(C$proxy_status)))
stats_set <- function(d, pcol = "pct_uORF_plus", label) {
  d <- d %>% filter(!is.na(.data[[pcol]]), is.finite(uIF_raw), uIF_raw > 0)
  x <- d[[pcol]]
  dT <- d %>% filter(is.finite(TE), TE > 0)
  grp <- cut(x, c(-Inf, 30, 80, Inf), right = FALSE, labels = c("<30", "30-80", "80-100"))
  medTE <- tapply(dT$TE, cut(dT[[pcol]], c(-Inf, 30, 80, Inf), right = FALSE), median)
  data.frame(set = label, n = nrow(d), n_TE = nrow(dT),
             pct_100 = round(100 * mean(x >= 100), 1), pct_lt90 = round(100 * mean(x < 90), 1),
             grp_lt30 = round(100 * mean(grp == "<30"), 1), grp_30_80 = round(100 * mean(grp == "30-80"), 1), grp_80_100 = round(100 * mean(grp == "80-100"), 1),
             r_pearson_pUP_log2uIF = round(cor(x, log2(d$uIF_raw)), 3), rho_pUP_uIF = round(cor(x, d$uIF_raw, method = "spearman"), 3),
             r_pearson_pUP_log2TE = round(cor(dT[[pcol]], log2(dT$TE)), 3), rho_pUP_TE = round(cor(dT[[pcol]], dT$TE, method = "spearman"), 3),
             medTE_lt30 = round(medTE[1], 3), medTE_30_80 = round(medTE[2], 3), medTE_80_100 = round(medTE[3], 3),
             uIF_strong = sum(d$uIF_raw > 1), uIF_medium = sum(d$uIF_raw > 0.2 & d$uIF_raw <= 1), uIF_weak = sum(d$uIF_raw <= 0.2))
}
sens <- bind_rows(
  stats_set(C, "pct_uORF_plus", "all single-TuORF genes (Table S3)"),
  stats_set(C %>% filter(proxy_status == "clean (annotation + junctions)"), "pct_uORF_plus", "clean: no annotated leader-splicing conflict, S >= 0.90 at covered sites"),
  stats_set(C %>% filter(proxy_status != "affected by splicing"), "pct_uORF_plus", "not affected (clean + junction-unassessable)"),
  stats_set(C %>% filter(!annot_affected), "pct_uORF_plus", "annotation-clean only (classes A/B)"),
  stats_set(C %>% filter(proxy_status == "affected by splicing"), "pct_uORF_plus", "affected by splicing (annotation or junctions)"),
  stats_set(C %>% filter(!is.na(pUP_eff_CAGE)), "pct_uORF_plus", "subset with >= 10 RNA-seq reads at the TIS, plain %uORF+ (reference for the next row)"),
  stats_set(C, "pUP_eff_CAGE", "same subset, %uORF+ replaced by effective inclusion (%uORF+ x S_TIS)"),
  stats_set(C %>% filter(proxy_status == "clean (annotation + junctions)"), "pUP_eff_CAGE", "clean, effective inclusion")
)
.add(""); .add("=== 5. Single-TuORF cohort: sensitivity to leader splicing (CAGE %uORF+) ===")
for (k in seq_len(nrow(sens))) .add(sprintf("  %-85s n=%4d | 100%%: %4.1f%% | <90%%: %4.1f%% | groups %4.1f/%4.1f/%4.1f | r(pUP,log2 uIF)=%.3f | r(pUP,log2 TE)=%.3f | median TE %.2f/%.2f/%.2f",
                                            sens$set[k], sens$n[k], sens$pct_100[k], sens$pct_lt90[k], sens$grp_lt30[k], sens$grp_30_80[k], sens$grp_80_100[k],
                                            sens$r_pearson_pUP_log2uIF[k], sens$r_pearson_pUP_log2TE[k], sens$medTE_lt30[k], sens$medTE_30_80[k], sens$medTE_80_100[k]))
tabC <- C %>% count(proxy_status) %>% mutate(pct = round(100 * n / sum(n), 1))
for (k in seq_len(nrow(tabC))) .add(sprintf("  cohort: %-45s %5d (%4.1f%%)", tabC$proxy_status[k], tabC$n[k], tabC$pct[k]))

# Fig 5 ladder (uORF length vs exclusion) on all vs clean, single-TuORF cohort, non-overlapping uORFs
lad <- function(d, label) d %>% filter(ORF_type == "uORF", uORF_length_aa >= 2) %>%
  mutate(bin = cut(uORF_length_aa, c(1, 10, 20, 30, Inf), labels = c("2-10", "11-20", "21-30", ">30"))) %>%
  group_by(bin) %>% summarise(n = n(), pct_excluded = round(100 * mean(pct_uORF_plus <= 30), 1), .groups = "drop") %>% mutate(set = label)
ladder <- bind_rows(lad(C, "all"), lad(C %>% filter(proxy_status == "clean (annotation + junctions)"), "clean"))
.add(""); .add("=== 5b. uORF length vs exclusion (%uORF+ <= 30), single-TuORF cohort ===")
for (s in unique(ladder$set)) { l <- ladder %>% filter(set == s); .add(sprintf("  %-6s %s", s, paste(sprintf("%s aa: %.1f%% (n=%d)", l$bin, l$pct_excluded, l$n), collapse = " | "))) }

# Fig S9 (S3 as submitted; 2-TuORF genes): both TuORFs clean
two <- U %>% filter(is.finite(uIF_raw), uIF_raw > 0) %>% group_by(gene_id) %>% filter(n() == 2) %>%   # same uIF > 0 filter as the multi-uORF figure module (Figure S9 of the revision since V39, S3 as submitted)
  arrange(ifelse(Strand == "+", g_min, -g_min), .by_group = TRUE) %>% mutate(rank = row_number(), both_clean = all(proxy_status == "clean (annotation + junctions)")) %>% ungroup()
.add(""); .add("=== 5c. Genes with two TuORFs: median %uORF+ of the first vs second TuORF (CAGE) ===")
.add(sprintf("  all 2-TuORF genes (n=%d): %.0f%% vs %.0f%% | both TuORFs clean (n=%d): %.0f%% vs %.0f%%",
             n_distinct(two$gene_id), median(two$pct_uORF_plus[two$rank == 1]), median(two$pct_uORF_plus[two$rank == 2]),
             n_distinct(two$gene_id[two$both_clean]), median(two$pct_uORF_plus[two$rank == 1 & two$both_clean]), median(two$pct_uORF_plus[two$rank == 2 & two$both_clean])))

################################################################################
# 6) Outputs
################################################################################
legend <- data.frame(Column = c("gene_id", "tx_id", "uORF_genome_range", "uORF_start_genomic", "uORF_stop_genomic", "detection_method", "n_TuORFs_in_gene",
                                "pct_uORF_plus", "pUP_TSSseq", "uIF", "mORF_TE", "n_isoforms", "gene_has_leader_intron", "src_has_leader_intron",
                                "src_intron_in_body", "other_TIS_in_intron", "other_TIS_in_intron_ex", "other_body_overlap", "other_body_overlap_ex",
                                "annot_class", "annot_class_ex", "R_TIS", "E_TIS", "AS_TIS", "cov_TIS", "S_TIS", "R_STOP", "E_STOP", "AS_STOP", "cov_STOP", "S_STOP",
                                "junction_status", "pUP_eff_CAGE", "pUP_eff_TSSseq", "proxy_status"),
                     Description = c("Gene", "Source isoform of the TuORF call", "chr_start_end of the TuORF", "First nt of the start codon (genomic)", "Last nt of the stop codon (genomic)",
                                     "RiboTaper / CiPS / Both", "TuORFs detected in the gene", "%uORF+ from CAGE-seq (Table S1A)", "%uORF+ from TSS-seq (Table S1B)", "uORF intensity factor (relative ribosome occupancy; Table S1A)", "mORF TE (Table S1A)",
                                     "Coding isoforms of the gene in Araport11+CTRL", "Any coding isoform of the gene has an intron entirely upstream of the start of its coding sequence",
                                     "The source isoform has a leader intron", "An intron of the source isoform lies inside the TuORF (TuORF spans a splice junction)",
                                     "The start codon (translation initiation site, TIS) lies inside an intron of another isoform of the gene (that isoform lacks the start codon)", "Same, counting only alternative isoforms in the expressed 2-isoform set",
                                     "An intron of another isoform overlaps the TuORF without containing the start codon (TIS kept, body/stop altered)", "Same, expressed isoforms only",
                                     "Annotation class (D > C1 > C2 > B > A)", "Annotation class using expressed alternative isoforms only",
                                     "RNA-seq sense reads whose aligned block covers the first nt of the start codon (retained)", "RNA-seq sense reads whose splice gap spans that position (spliced over); gaps that are annotated introns of another gene only are not counted",
                                     "Antisense reads covering the position (not used)", "R_TIS + E_TIS", "R_TIS / (R_TIS + E_TIS); NA if fewer than 10 reads",
                                     "RNA-seq sense reads whose aligned block covers the last nt of the stop codon (retained)", "RNA-seq sense reads whose splice gap spans that position (spliced over; same rule as E_TIS)",
                                     "Antisense reads covering the position (not used)", "R_STOP + E_STOP", "R_STOP / (R_STOP + E_STOP); NA if fewer than 10 reads",
                                     "retained / spliced over (S < 0.90 at TIS or stop) / unassessable (< 10 reads at both)", "%uORF+ (CAGE) x S_TIS; NA where S_TIS is NA (fewer than 10 reads)", "%uORF+ (TSS-seq) x S_TIS; NA where S_TIS is NA",
                                     "clean (annotation + junctions) / annotation-clean, junctions unassessable / affected by splicing"))
out_cols <- c(legend$Column)
Uout <- U %>% select(all_of(out_cols)) %>% arrange(gene_id, uORF_start_genomic)
write.xlsx(list("A2 TuORF leader splicing" = Uout, "Column legend" = legend), file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.xlsx"), overwrite = TRUE)
fwrite(Uout, file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.csv"))
fwrite(sens, file.path(OUT_DIR, "A2_sensitivity_single_TuORF_cohort.csv"))
fwrite(ladder, file.path(OUT_DIR, "A2_length_ladder_all_vs_clean.csv"))
fwrite(C %>% select(gene_id, uORF_genome_range, pct_uORF_plus, pUP_eff_CAGE, pUP_eff_TSSseq, uIF, mORF_TE, annot_class, junction_status, proxy_status, S_TIS, S_STOP, cov_TIS),
       file.path(OUT_DIR, "A2_single_TuORF_cohort_splicing_flags.csv"))

# --- figure ---
cls_short <- function(x) sub("\\..*$", "", x)
U$cls <- factor(cls_short(U$annot_class), levels = c("A", "B", "C2", "C1", "D"))
C$cls <- factor(cls_short(C$annot_class), levels = c("A", "B", "C2", "C1", "D"))
pal_cls <- c(A = "#4D4D4D", B = "#8C8C8C", C2 = "#F4A582", C1 = "#D6604D", D = "#B2182B")
dfA <- bind_rows(U %>% transmute(set = "all TuORFs", cls), C %>% transmute(set = "single-TuORF genes", cls))
pA <- ggplot(dfA %>% count(set, cls) %>% group_by(set) %>% mutate(p = 100 * n / sum(n)), aes(x = set, y = p, fill = cls)) +
  geom_col(width = 0.6, colour = "white") + scale_fill_manual(values = pal_cls, name = "Annotation class",
    labels = c(A = "A: no leader intron", B = "B: leader intron, TuORF untouched", C2 = "C2: other isoform alters body/stop", C1 = "C1: TuORF spans own intron", D = "D: start codon in other isoform's intron")) +
  theme_classic() + labs(title = "A", x = NULL, y = "% of TuORFs") + scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) + theme(plot.title = element_text(face = "bold", size = 16), legend.position = "bottom", legend.direction = "vertical")
dB <- U %>% filter(!is.na(S_TIS))
pB <- ggplot(dB, aes(x = S_TIS, colour = cls)) + stat_ecdf(geom = "step", linewidth = 1) + scale_colour_manual(values = pal_cls, name = "Class") +
  geom_vline(xintercept = S_CLEAN, linetype = "dashed") + theme_classic() +
  labs(title = "B", x = "Fraction of RNA-seq reads retaining the start codon (S_TIS)", y = "Fraction of TuORFs") + theme(plot.title = element_text(face = "bold", size = 16), legend.position = c(0.25, 0.65))
dC <- U %>% filter(!is.na(pUP_eff_CAGE))
pC <- ggplot(dC, aes(x = pct_uORF_plus, y = pUP_eff_CAGE, colour = cls)) + geom_point(alpha = 0.5, size = 0.9) + scale_colour_manual(values = pal_cls, name = "Class") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + theme_classic() +
  labs(title = "C", x = "%uORF+ (CAGE)", y = "Effective inclusion = %uORF+ x S_TIS") + theme(plot.title = element_text(face = "bold", size = 16), legend.position = "bottom") + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1)))
dD <- C %>% filter(is.finite(uIF_raw), uIF_raw > 0) %>% mutate(set = ifelse(proxy_status == "clean (annotation + junctions)", "clean", ifelse(proxy_status == "affected by splicing", "affected", "unassessable")))
pD <- ggplot(dD %>% filter(set != "unassessable"), aes(x = pct_uORF_plus, y = log2(uIF_raw), colour = set)) + geom_point(alpha = 0.45, size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) + scale_colour_manual(values = c(clean = "#2166AC", affected = "#B2182B"), name = NULL) + theme_classic() +
  labs(title = "D", x = "%uORF+ (CAGE)", y = bquote(log[2]~"(uIF)")) + theme(plot.title = element_text(face = "bold", size = 16), legend.position = c(0.15, 0.9))
dE <- dD %>% filter(is.finite(TE), TE > 0, set == "clean") %>% mutate(grp = cut(pct_uORF_plus, c(-Inf, 30, 80, Inf), right = FALSE, labels = c("0~30", "30~80", "80~100")))
pE <- ggplot(dE %>% filter(TE <= 3.5), aes(x = TE, colour = grp)) + stat_ecdf(geom = "step", linewidth = 1.1) + coord_cartesian(xlim = c(0, 3.5)) + theme_classic() +
  scale_colour_manual(values = c("80~100" = "#F8766D", "30~80" = "#00BA38", "0~30" = "#619CFF"), name = "%uORF+ (clean set)") +   # V23: the colours of Figure 2F
  labs(title = "E", x = "mORF TE", y = "Fraction") + theme(plot.title = element_text(face = "bold", size = 16), legend.position = c(0.7, 0.4))
pF <- ggplot(ladder, aes(x = bin, y = pct_excluded, fill = set)) + geom_col(position = position_dodge(0.7), width = 0.65) +
  scale_fill_manual(values = c(all = "#8C8C8C", clean = "#2166AC"), name = NULL) + theme_classic() +
  labs(title = "F", x = "TuORF length (aa)", y = "% TuORFs with %uORF+ <= 30") + theme(plot.title = element_text(face = "bold", size = 16), legend.position = c(0.2, 0.9))
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
# V30 (2026-11-09, Larry): legends inside the plots (A: in free room made at the right of the bars); the panels of a row share their plotting area (fig_grid)
cS9 <- suppressMessages(lapply(list(pA + scale_x_discrete(labels = function(x) gsub(" ", "\n", x), expand = expansion(add = c(0.6, 2.9))), pB, pC + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))), pD, pE, pF), compact))
cS9 <- list(legend_in(cS9[[1]], 1, 0.5, just = c(1, 0.5), title = TRUE, wrap = wrap_w(24)), legend_in(cS9[[2]], 0.5, 0.08, just = c(0.5, 0), title = TRUE), legend_in(cS9[[3]], 0.03, 0.99, title = TRUE),
            legend_in(cS9[[4]], 0.03, 0.99), legend_in(cS9[[5]], 0.99, 0.03, title = TRUE), legend_in(cS9[[6]], 0.03, 0.99))
fig <- fig_grid(cS9, ncol = 2, top = fig_title("Figure S7. Leader splicing and TSS-based estimates of TuORF inclusion"), h = 8)
fig_save(file.path(OUT_DIR, "A2_Figure_leader_splicing_audit.pdf"), fig, h = 8)

.add(""); .add(sprintf("A2 done in %.1f min. Outputs: A2_Table_TuORF_leader_splicing.xlsx/.csv, A2_sensitivity_single_TuORF_cohort.csv, A2_length_ladder_all_vs_clean.csv, A2_single_TuORF_cohort_splicing_flags.csv, A2_Figure_leader_splicing_audit.pdf, A2_summary.txt",
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(c(sprintf("A2 leader-splicing audit — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), file.path(OUT_DIR, "A2_summary.txt"))
