################################################################################
# A12 — Uncapped mRNA (PARE of Nagarajan et al. 2019) versus uIF and mORF TE     #
# Repeat of A11 (GMUCT of Yu et al. 2021) with an independent degradome data    #
# set: poly(A)+ PARE of two-week-old Col-0 seedlings from Nagarajan, Kukulich,   #
# von Hagel & Green 2019 Nucleic Acids Res 47:9216-9230 (GEO GSE119706).         #
#                                                                                #
# PARE: ~/Desktop/Desktop_2026_July/PARE_analysis/BAMtoBED_2ends/WT_polyA+a.tsv  #
#   (GEO WT poly(A)+ 5'P PARE replicate 2, ATH237, SRR8429733) and WT_polyA+b.tsv#
#   (replicate 1, ATH362). Four columns, no header: count, chromosome, 1-based   #
#   genomic position of the read 5' end (both strands; never add 1), strand.    #
# RNA-seq denominator: WT poly(A)+ RNA-seq replicate 1 of the same study        #
#   (SRR7812609; STAR 2.7.9a / TAIR10 BAM in ~/Desktop/Desktop_2026_July/xrn4/    #
#   STAR). The local replicate-2 BAM is truncated (samtools quickcheck fails),   #
#   so the denominator is a single replicate. Sensitivity denominator: our own  #
#   RNA-seq TPM (RSEM), as in A11.                                               #
# Counting: PARE 5' ends and uniquely mapped RNA-seq reads (MAPQ 255) within the #
#   exons of a gene (union over the Araport11 isoforms of Files/Araport11+CTRL_  #
#   20181206.gtf, AT identifiers only; PARE strand-aware, RNA-seq strandedness   #
#   inferred from the data); positions / reads overlapping two genes discarded.  #
# RPM over the reads assigned to nuclear genes (chromosomes 1-5); uncapped ratio #
#   U = mean PARE RPM of the two libraries / RNA-seq RPM (as A11; = 2 to the     #
#   power of the "proportion uncapped" of Yu et al. 2021).                       #
# Cohorts, statistics and panels as A11 (Figure S19): single-TuORF genes (Table #
#   S3, uIF classes of Figure 1C) and the No-TuORF genes of Figure 1C. Panel E   #
#   compares the PARE ratio with the GMUCT ratio of A11 gene by gene.            #
# Outputs (output/): A12_Figure_PARE_uncapped_vs_uIF.pdf, A12_Table_PARE_uncapped#
#   .xlsx, A12_Table_genes.csv, A12_gene_counts.csv (all genes), A12_summary.txt, #
#   A12_run.log                                                                  #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A12_PARE_uncapped_vs_uIF.R
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(openxlsx); library(gridExtra); library(grid); library(multcompView)
  library(GenomicRanges); library(GenomicAlignments); library(Rsamtools); library(GenomeInfoDb); library(SummarizedExperiment)
})
options(stringsAsFactors = FALSE, warn = 1)
set.seed(12)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PARE_DIR <- path.expand("~/Desktop/Desktop_2026_July/PARE_analysis/BAMtoBED_2ends")
PARE     <- c(PARE_a = "WT_polyA+a.tsv", PARE_b = "WT_polyA+b.tsv")      # a = GEO replicate 2 (ATH237), b = GEO replicate 1 (ATH362)
PATH_BAM <- path.expand("~/Desktop/Desktop_2026_July/xrn4/STAR/WT_PolyA+_rep1Aligned.sortedByCoord.out.bam")
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_S1  <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
S2_CANDIDATES <- c(Sys.getenv("A12_S2", ""), "/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS3_single_TuORF_genes_revised.xlsx", "/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS2_single_TuORF_genes_revised.xlsx")
PATH_S2  <- S2_CANDIDATES[nzchar(S2_CANDIDATES) & file.exists(S2_CANDIDATES)][1]
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_RSEM_RNA  <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_5UTR <- file.path(DATA_DIR, "genes_w_5UTR.tsv")
PATH_HL   <- file.path(DATA_DIR, "RNA-half-lives-5EU.tsv")
PATH_A11  <- file.path(OUT_DIR, "A11_Table_genes.csv")                   # GMUCT-based uncapped ratios of A11, for the cross-data-set comparison
for (f in c(file.path(PARE_DIR, PARE), PATH_BAM, paste0(PATH_BAM, ".bai"), PATH_GTF, PATH_S1, PATH_S2, PATH_OUORF, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_5UTR, PATH_HL, PATH_A11)) stopifnot(file.exists(f))
MIN_RNA <- 10L            # RNA-seq reads required for a gene's uncapped ratio (one replicate of 24 M mapped reads)
CUTOFFS <- c(5L, 10L, 20L, 50L); NBOOT <- 1000L
CLASSES <- c("High uIF (> 1)", "Medium uIF (0.2 to 1)", "Low uIF (0.2 or less)", "No TuORF")
PUP_CLASSES <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
big <- function(x) trimws(format(x, big.mark = ",")); f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

################################################################################
# 1) Exon models per gene (union over isoforms)
################################################################################
.add("A12: exon models per gene from ", PATH_GTF)
gtf <- fread(PATH_GTF, header = FALSE, sep = "\t", quote = "", select = c(1, 3, 4, 5, 7, 9), col.names = c("chr", "type", "start", "end", "strand", "attr"))
ex <- gtf[type == "exon", .(chr = as.character(chr), start, end, strand, gene = sub('.*gene_id "([^"]+)".*', "\\1", attr))]
ex <- ex[grepl("^AT[1-5CM]G", gene)]                                     # Araport11 genes only; the CTRL (MSTRG) assemblies are left out so that they do not make positions ambiguous
exg <- reduce(split(GRanges(ex$chr, IRanges(ex$start, ex$end), strand = ex$strand), ex$gene))
gene_len <- sum(width(exg)); nuclear <- grepl("^AT[1-5]G", names(exg))
.add(sprintf("  %s Araport11 genes (%s nuclear), %s exonic bases in total", big(length(exg)), big(sum(nuclear)), big(sum(gene_len))))

################################################################################
# 2) PARE 5' ends per gene
################################################################################
count_pare <- function(path) {
  d <- fread(path, header = FALSE, col.names = c("n", "chr", "pos", "strand"), colClasses = c("integer", "character", "integer", "character"))
  gr <- GRanges(d$chr, IRanges(d$pos, width = 1), strand = d$strand)      # 1-based 5' end on both strands, used as is
  h <- findOverlaps(gr, exg)                                              # strand-aware: PARE reads are sense to the mRNA
  nh <- tabulate(queryHits(h), nbins = length(gr)); keep <- nh[queryHits(h)] == 1L
  cnt <- tapply(as.numeric(d$n[queryHits(h)[keep]]), names(exg)[subjectHits(h)[keep]], sum)
  out <- setNames(numeric(length(exg)), names(exg)); out[names(cnt)] <- cnt
  list(counts = out, positions = nrow(d), total = sum(as.numeric(d$n)), nuclear_chr = sum(as.numeric(d$n[d$chr %in% as.character(1:5)])),
       assigned = sum(as.numeric(d$n[nh == 1L])), ambiguous = sum(as.numeric(d$n[nh > 1L])))
}
pare <- lapply(file.path(PARE_DIR, PARE), count_pare); names(pare) <- names(PARE)
for (j in names(pare)) .add(sprintf("  %s (%s): %s distinct 5' ends, %s reads (%s on chromosomes 1-5); %s (%.1f%%) assigned to the exons of one gene on the same strand, %s (%.1f%%) overlapping two genes discarded",
  j, PARE[[j]], big(pare[[j]]$positions), big(pare[[j]]$total), big(pare[[j]]$nuclear_chr), big(pare[[j]]$assigned), 100 * pare[[j]]$assigned / pare[[j]]$total,
  big(pare[[j]]$ambiguous), 100 * pare[[j]]$ambiguous / pare[[j]]$total))

################################################################################
# 3) RNA-seq reads per gene (strandedness inferred from the data)
################################################################################
infer_strandedness <- function(bam, exg, n_genes = 3000L) {                 # share of unique reads on the gene strand over isolated nuclear genes
  gr <- unlist(range(exg)); gr <- gr[countOverlaps(gr, gr, ignore.strand = TRUE) == 1L & grepl("^AT[1-5]G", names(gr)) & width(gr) >= 1500]
  gr <- gr[sample(length(gr), min(n_genes, length(gr)))]
  rd <- readGAlignments(bam, param = ScanBamParam(which = reduce(gr, ignore.strand = TRUE), mapqFilter = 255L, flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)))
  h <- findOverlaps(rd, gr, ignore.strand = TRUE)
  c(same = mean(as.character(strand(rd))[queryHits(h)] == as.character(strand(gr))[subjectHits(h)]), reads = length(h))
}
st <- infer_strandedness(PATH_BAM, exg)
lib_type <- if (st[["same"]] >= 0.85) "stranded, read = mRNA strand" else if (st[["same"]] <= 0.15) "stranded, read = reverse of mRNA" else "unstranded"
.add(sprintf("  RNA-seq %s: %.1f%% of %s unique reads over isolated genes lie on the gene strand -> counted as %s", basename(PATH_BAM), 100 * st[["same"]], big(st[["reads"]]), lib_type))
idx <- idxstatsBam(PATH_BAM)
.add(sprintf("  alignments in the BAM: %s (%s on chromosomes 1-5)", big(sum(idx$mapped)), big(sum(idx$mapped[idx$seqnames %in% as.character(1:5)]))))
bf <- BamFile(PATH_BAM, yieldSize = 2e6)
sbp <- ScanBamParam(mapqFilter = 255L, flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE))
se <- if (grepl("reverse", lib_type)) {
  summarizeOverlaps(exg, bf, mode = "Union", singleEnd = TRUE, ignore.strand = FALSE, inter.feature = TRUE, param = sbp, preprocess.reads = invertStrand)
} else {
  summarizeOverlaps(exg, bf, mode = "Union", singleEnd = TRUE, ignore.strand = (lib_type == "unstranded"), inter.feature = TRUE, param = sbp)
}
rna_counts <- as.numeric(assay(se)[, 1]); names(rna_counts) <- rownames(se)
.add(sprintf("  RNA-seq reads assigned to the exons of one gene (unique, MAPQ 255): %s, of which %s to nuclear genes (%.1f min so far)", big(sum(rna_counts)), big(sum(rna_counts[nuclear])), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

################################################################################
# 4) RPM and uncapped ratio U
################################################################################
cnt <- data.table(gene_id = names(exg), nuclear = nuclear, exon_length = gene_len, PARE_a = pare$PARE_a$counts, PARE_b = pare$PARE_b$counts, RNA_rep1 = rna_counts[names(exg)])
libsize <- c(PARE_a = sum(cnt[nuclear == TRUE, PARE_a]), PARE_b = sum(cnt[nuclear == TRUE, PARE_b]), RNA_rep1 = sum(cnt[nuclear == TRUE, RNA_rep1]))   # reads assigned to nuclear genes
cnt[, `:=`(rpm_PARE_a = PARE_a / libsize[["PARE_a"]] * 1e6, rpm_PARE_b = PARE_b / libsize[["PARE_b"]] * 1e6, RNA_rpm = RNA_rep1 / libsize[["RNA_rep1"]] * 1e6)]
cnt[, `:=`(PARE_reads = PARE_a + PARE_b, RNA_reads = RNA_rep1, PARE_rpm = (rpm_PARE_a + rpm_PARE_b) / 2)]
cnt[, U := PARE_rpm / RNA_rpm]
cnt[, U_a := ifelse(RNA_rep1 > 0 & PARE_a > 0, rpm_PARE_a / RNA_rpm, NA_real_)]
cnt[, U_b := ifelse(RNA_rep1 > 0 & PARE_b > 0, rpm_PARE_b / RNA_rpm, NA_real_)]
nuc <- cnt[nuclear == TRUE]
.add(sprintf("  nuclear-gene reads used for RPM: PARE a %s, PARE b %s, RNA-seq %s", big(libsize[["PARE_a"]]), big(libsize[["PARE_b"]]), big(libsize[["RNA_rep1"]])))
.add(sprintf("  nuclear genes with >= %d RNA-seq reads: %s of %s; agreement of the two PARE libraries among them (Spearman): PARE RPM %.3f, uncapped ratio %.3f",
             MIN_RNA, big(nuc[RNA_reads >= MIN_RNA, .N]), big(nrow(nuc)), nuc[RNA_reads >= MIN_RNA, cor(rpm_PARE_a, rpm_PARE_b, method = "spearman")],
             nuc[RNA_reads >= MIN_RNA & is.finite(U_a) & is.finite(U_b), cor(U_a, U_b, method = "spearman")]))
fwrite(cnt, file.path(OUT_DIR, "A12_gene_counts.csv"))

################################################################################
# 5) Cohorts: single-TuORF genes (Table S3) and the No-TuORF controls of Figure 1C
################################################################################
S2 <- as.data.table(read.xlsx(PATH_S2, 1))
S2 <- S2[, .(gene_id, uIF = ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, uIF), uIF_Group, pUP = pct_uORF_plus, pUP_Group, RNA_TPM, Ribo_TPM, TE = mORF_TE)]   # unrounded uIF (2026-11-01)
S2[, class := fcase(grepl("^(High|Strong)", uIF_Group), CLASSES[1], grepl("^Medium", uIF_Group), CLASSES[2], grepl("^(Low|Weak)", uIF_Group), CLASSES[3])]
S2[, pclass := fcase(pUP_Group == "0~30", PUP_CLASSES[1], pUP_Group == "30~80", PUP_CLASSES[2], pUP_Group == "80~100", PUP_CLASSES[3])]
stopifnot(!anyNA(S2$class), !anyNA(S2$pclass))
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)"); ou <- read.xlsx(PATH_OUORF)
rna <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); ribo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
ctrl <- data.table(gene_id = rna$gene_id, RNA_TPM = rna$TPM, Ribo_TPM = ribo$TPM[match(rna$gene_id, ribo$gene_id)])
ctrl[, TE := Ribo_TPM / RNA_TPM]
ctrl <- ctrl[gene_id %in% read.delim(PATH_5UTR)$gene_id & RNA_TPM > 1 & is.finite(TE) & !gene_id %in% unique(c(S1A$gene_id, as.character(ou$gene_id)))]
ctrl[, `:=`(uIF = NA_real_, uIF_Group = NA_character_, pUP = NA_real_, pUP_Group = NA_character_, class = CLASSES[4], pclass = PUP_CLASSES[4])]
G <- rbind(S2[, cohort := "single-TuORF"], ctrl[, cohort := "No TuORF"], fill = TRUE)
hl <- fread(PATH_HL); setnames(hl, 1, "gene_id"); hl <- hl[, .(gene_id, half_life = Mean)]
G <- merge(G, hl, by = "gene_id", all.x = TRUE)
a11 <- fread(PATH_A11, select = c("gene_id", "uncapped_ratio", "in_uncapped_ratio_set", "GMUCT_RPM", "RNA_RPM")); setnames(a11, c("gene_id", "U_GMUCT", "in_GMUCT_set", "GMUCT_rpm", "RNA_rpm_yu"))
a11[!(in_GMUCT_set %in% TRUE) | !is.finite(U_GMUCT), `:=`(U_GMUCT = NA_real_, GMUCT_rpm = NA_real_, RNA_rpm_yu = NA_real_)]   # GMUCT values only for genes in the A11 set
G <- merge(G, a11[, .(gene_id, U_GMUCT, GMUCT_rpm, RNA_rpm_yu)], by = "gene_id", all.x = TRUE)
G <- merge(G, cnt[, .(gene_id, PARE_a, PARE_b, RNA_rep1, PARE_reads, RNA_reads, rpm_PARE_a, rpm_PARE_b, PARE_rpm, RNA_rpm, U, U_a, U_b)], by = "gene_id", all.x = TRUE)
G[, U_labRNA := PARE_rpm / RNA_TPM]                                      # sensitivity: our own RNA-seq as denominator (CDS-only RSEM quantification; see note)
G[, U_GMUCT_nag := GMUCT_rpm / RNA_rpm]                                  # GMUCT (Yu et al.) over the Nagarajan RNA-seq: the same denominator as the PARE ratio, isolating the agreement of the two degradome signals
G[, class := factor(class, CLASSES)]; G[, pclass := factor(pclass, PUP_CLASSES)]
.add(sprintf("  cohorts: %s single-TuORF genes (Table S2; %d / %d / %d high / medium / low uIF) and %s No-TuORF control genes; %s have counts, %s a half-life, %s a GMUCT ratio (A11)",
             big(sum(G$cohort == "single-TuORF")), sum(G$class == CLASSES[1]), sum(G$class == CLASSES[2]), sum(G$class == CLASSES[3]),
             big(sum(G$cohort == "No TuORF")), big(sum(!is.na(G$RNA_reads))), big(sum(!is.na(G$half_life))), big(sum(is.finite(G$U_GMUCT)))))
D <- G[RNA_reads >= MIN_RNA & PARE_reads >= 1]
.add(sprintf("  genes with >= %d RNA-seq reads and >= 1 PARE read (uncapped-ratio set): %s single-TuORF (%d / %d / %d) + %s No TuORF; median RNA-seq reads %s (single-TuORF) / %s (No TuORF), median PARE reads %s / %s",
             MIN_RNA, big(sum(D$cohort == "single-TuORF")), sum(D$class == CLASSES[1]), sum(D$class == CLASSES[2]), sum(D$class == CLASSES[3]), big(sum(D$cohort == "No TuORF")),
             big(median(D[cohort == "single-TuORF", RNA_reads])), big(median(D[cohort == "No TuORF", RNA_reads])), big(median(D[cohort == "single-TuORF", PARE_reads])), big(median(D[cohort == "No TuORF", PARE_reads]))))
for (cl in CLASSES) .add(sprintf("    %s: %d of %d genes kept (%.0f%%); median RNA TPM (our data) %.1f", cl, sum(D$class == cl), sum(G$class == cl), 100 * sum(D$class == cl) / sum(G$class == cl), median(D[class == cl, RNA_TPM], na.rm = TRUE)))

################################################################################
# 6) Statistics helpers (as A11)
################################################################################
sp_ci <- function(x, y, B = NBOOT) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(c(rho = NA, lo = NA, hi = NA, n = n))
  r <- cor(x, y, method = "spearman"); rx <- rank(x); ry <- rank(y)
  b <- replicate(B, { i <- sample.int(n, n, TRUE); cor(rx[i], ry[i]) })
  c(rho = r, lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)), n = n)
}
psp <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z); rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  rxy <- cor(rx, ry); rxz <- cor(rx, rz); ryz <- cor(ry, rz); (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}
fold_ci <- function(a, b, B = NBOOT) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; f <- median(a) / median(b)
  bb <- replicate(B, median(sample(a, length(a), TRUE)) / median(sample(b, length(b), TRUE)))
  c(fold = f, lo = unname(quantile(bb, 0.025)), hi = unname(quantile(bb, 0.975)))
}
ks_letters <- function(d, grp, val) {
  d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]), ]; lv <- levels(droplevels(d[[grp]])); d[[grp]] <- factor(d[[grp]], lv)
  cmb <- combn(lv, 2, simplify = FALSE)
  pw <- rbindlist(lapply(cmb, function(g) { a <- d[[val]][d[[grp]] == g[1]]; b <- d[[val]][d[[grp]] == g[2]]
    data.table(group1 = g[1], group2 = g[2], n1 = length(a), n2 = length(b), KS_p = suppressWarnings(ks.test(a, b)$p.value),
               Wilcoxon_p = wilcox.test(a, b)$p.value, median_fold_1_over_2 = median(a) / median(b)) }))
  ids <- setNames(paste0("G", seq_along(lv)), lv); p <- pw$KS_p; names(p) <- paste(ids[pw$group1], ids[pw$group2], sep = "-")
  L <- multcompLetters(p, compare = "<", threshold = 0.05, Letters = letters)$Letters
  L <- L[ids]; names(L) <- lv
  ch <- unique(unlist(strsplit(unname(L), "")))
  mp <- setNames(letters[seq_along(ch)], ch); L[] <- vapply(unname(L), function(s) paste(mp[strsplit(s, "")[[1]]], collapse = ""), "")
  n <- table(d[[grp]])[lv]
  list(pairs = pw, letters = L, labels = setNames(sprintf("%s (n = %s) %s", lv, trimws(big(as.integer(n))), L), lv), medians = tapply(d[[val]], d[[grp]], median)[lv])
}

################################################################################
# 7) Correlations (single-TuORF genes with an uncapped ratio)
################################################################################
S <- D[cohort == "single-TuORF"]
cors <- list()
addc <- function(label, x, y, extra = "") { v <- sp_ci(x, y); cors[[length(cors) + 1]] <<- data.table(comparison = label, note = extra, n = v[["n"]], rho = v[["rho"]], ci_lo = v[["lo"]], ci_hi = v[["hi"]]); v }
v1  <- addc("uncapped ratio vs uIF", S$U, S$uIF, "PARE (both libraries) / RNA-seq replicate 1 of Nagarajan et al. 2019; panel A")
v2  <- addc("uncapped ratio vs mORF TE", S$U, S$TE, "panel C")
v3  <- addc("uncapped ratio vs %uORF+", S$U, S$pUP)
v4  <- addc("uncapped ratio vs RNA-seq RPM (Nagarajan et al.)", S$U, S$RNA_rpm, "expression dependence of the ratio")
v5  <- addc("uncapped ratio vs RNA TPM (this study)", S$U, S$RNA_TPM, "expression dependence of the ratio")
v6  <- addc("uncapped ratio (our RNA-seq as denominator) vs uIF", S$U_labRNA, S$uIF, "sensitivity; TE shares this denominator, so TE is not tested with it")
v7  <- addc("uncapped ratio vs mORF TE, No-TuORF genes", D[cohort == "No TuORF", U], D[cohort == "No TuORF", TE], "panel F; does decapping track TE without a uORF")
v8  <- addc("uncapped ratio vs half-life, all genes with both", D$U, D$half_life, "validation of the ratio as a decay readout")
v9  <- addc("uncapped ratio (PARE) vs uncapped ratio (GMUCT, A11), all genes in both sets", D$U, D$U_GMUCT, "panel E")
v10 <- addc("uncapped ratio (PARE) vs uncapped ratio (GMUCT, A11), single-TuORF genes", S$U, S$U_GMUCT)
v9b <- addc("uncapped ratio (PARE) vs GMUCT RPM / Nagarajan RNA-seq RPM (same denominator), all genes in both sets", D$U, D$U_GMUCT_nag, "panel E; isolates the agreement of the two degradome signals")
v9c <- addc("PARE RPM vs GMUCT RPM (the two numerators), all genes in both sets", D$PARE_rpm, D$GMUCT_rpm)
v9d <- addc("RNA-seq RPM, Nagarajan et al. vs Yu et al. (the two denominators), all genes in both sets", D$RNA_rpm, D$RNA_rpm_yu)
v11 <- addc("uncapped ratio, library a only, vs uIF", S$U_a, S$uIF, "WT_polyA+a = GEO replicate 2")
v12 <- addc("uncapped ratio, library b only, vs uIF", S$U_b, S$uIF, "WT_polyA+b = GEO replicate 1")
v13 <- addc("uncapped ratio, library a only, vs mORF TE", S$U_a, S$TE)
v14 <- addc("uncapped ratio, library b only, vs mORF TE", S$U_b, S$TE)
v15 <- addc("uIF vs mORF TE, same genes", S$uIF, S$TE, "the Supplemental Fig. S2A relationship within the uncapped-ratio set")
v16 <- addc("half-life vs uncapped ratio, single-TuORF genes", S$half_life, S$U)
cor_tab <- rbindlist(cors)
pr_log <- function(x, y) { ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0; cor(log10(x[ok]), log10(y[ok])) }   # V41 (2026-11-12): Pearson r on log10 values without the genes that have a value of 0 (log10(0) is not finite; one cell of the table was NaN and showed #NUM! in Excel)
cor_tab[, pearson_r_log := NA_real_]
cor_tab[comparison == "uncapped ratio vs uIF", pearson_r_log := pr_log(S$U, S$uIF)]
cor_tab[comparison == "uncapped ratio vs mORF TE", pearson_r_log := pr_log(S$U, S$TE)]
pr_n0 <- function(x, y) { ok <- is.finite(x) & is.finite(y); c(all = sum(ok & !(x > 0 & y > 0)), y0 = sum(ok & y <= 0)) }   # genes that the Pearson r leaves out (all) and those of them whose second value is 0 (y0); the note of the row says how many (column n is the n of rho)
for (cmp in list(list("uncapped ratio vs uIF", S$U, S$uIF, "uIF"), list("uncapped ratio vs mORF TE", S$U, S$TE, "mORF TE"))) { k <- pr_n0(cmp[[2]], cmp[[3]])
  if (k[["all"]] > 0) cor_tab[comparison == cmp[[1]], note := paste0(note, ifelse(nzchar(note), "; ", ""), sprintf("Pearson r without %d gene%s with a %s of 0", k[["all"]], ifelse(k[["all"]] == 1, "", "s"), ifelse(k[["y0"]] == k[["all"]], cmp[[4]], "value")))] }
part <- data.table(comparison = c("uncapped ratio vs uIF | RNA-seq RPM (Nagarajan et al.)", "uncapped ratio vs uIF | RNA TPM (this study)",
                                  "uncapped ratio vs mORF TE | RNA-seq RPM (Nagarajan et al.)", "uncapped ratio vs mORF TE | RNA TPM (this study)",
                                  "uncapped ratio vs uIF | mORF TE", "uncapped ratio vs mORF TE | uIF"),
                   partial_rho = c(psp(S$U, S$uIF, S$RNA_rpm), psp(S$U, S$uIF, S$RNA_TPM), psp(S$U, S$TE, S$RNA_rpm), psp(S$U, S$TE, S$RNA_TPM), psp(S$U, S$uIF, S$TE), psp(S$U, S$TE, S$uIF)))
.add("  Spearman correlations (single-TuORF genes unless stated):")
for (i in seq_len(nrow(cor_tab))) .add(sprintf("    %-78s rho = %6.3f (95%% CI %.3f to %.3f), n = %s", cor_tab$comparison[i], cor_tab$rho[i], cor_tab$ci_lo[i], cor_tab$ci_hi[i], big(cor_tab$n[i])))
for (i in seq_len(nrow(part))) .add(sprintf("    partial %-70s rho = %6.3f", part$comparison[i], part$partial_rho[i]))
sens <- rbindlist(lapply(CUTOFFS, function(k) { s <- G[cohort == "single-TuORF" & RNA_reads >= k & PARE_reads >= 1]
  rbind(data.table(RNA_reads_min = k, PARE = "both libraries", n = nrow(s), rho_uIF = cor(s$U, s$uIF, method = "spearman"), rho_TE = cor(s$U, s$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s$U, s$pUP, method = "spearman")),
        { s1 <- s[is.finite(U_a)]; data.table(RNA_reads_min = k, PARE = "library a only", n = nrow(s1), rho_uIF = cor(s1$U_a, s1$uIF, method = "spearman"), rho_TE = cor(s1$U_a, s1$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s1$U_a, s1$pUP, method = "spearman")) },
        { s2 <- s[is.finite(U_b)]; data.table(RNA_reads_min = k, PARE = "library b only", n = nrow(s2), rho_uIF = cor(s2$U_b, s2$uIF, method = "spearman"), rho_TE = cor(s2$U_b, s2$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s2$U_b, s2$pUP, method = "spearman")) }) }))
.add("  sensitivity to the RNA-seq read cutoff and to the PARE library:")
for (i in seq_len(nrow(sens))) .add(sprintf("    >= %2d reads, %-15s n = %5d: rho(U, uIF) = %6.3f, rho(U, TE) = %6.3f, rho(U, %%uORF+) = %6.3f", sens$RNA_reads_min[i], sens$PARE[i], sens$n[i], sens$rho_uIF[i], sens$rho_TE[i], sens$rho_pUP[i]))

################################################################################
# 8) Group comparisons (Figure 1C style)
################################################################################
grp_summary <- function(d, grp, val, lv_ref) {
  k <- ks_letters(as.data.frame(d), grp, val); lv <- names(k$letters)
  tab <- rbindlist(lapply(lv, function(g) { x <- d[[val]][d[[grp]] == g & is.finite(d[[val]])]
    fc <- if (g != lv_ref) fold_ci(x, d[[val]][d[[grp]] == lv_ref & is.finite(d[[val]])]) else c(fold = 1, lo = NA, hi = NA)
    data.table(group = g, n = length(x), median = median(x), q25 = unname(quantile(x, .25)), q75 = unname(quantile(x, .75)),
               fold_vs_ref = fc[["fold"]], fold_lo = fc[["lo"]], fold_hi = fc[["hi"]], KS_letter = k$letters[[g]]) }))
  list(table = tab, pairs = k$pairs, labels = k$labels)
}
gU  <- grp_summary(D, "class", "U", CLASSES[4])                              # uncapped ratio by uIF class
gUp <- grp_summary(D, "pclass", "U", PUP_CLASSES[4])                         # uncapped ratio by %uORF+ group
gUa <- grp_summary(D[is.finite(U_a)], "class", "U_a", CLASSES[4])            # library a only
gUb <- grp_summary(D[is.finite(U_b)], "class", "U_b", CLASSES[4])            # library b only
gUl <- grp_summary(D, "class", "U_labRNA", CLASSES[4])                       # our RNA-seq as denominator
gG  <- grp_summary(D[is.finite(U_GMUCT)], "class", "U_GMUCT", CLASSES[4])    # GMUCT ratio (A11) on the genes present in both sets
gT  <- grp_summary(D, "class", "TE", CLASSES[4])                             # TE in the same set
gR  <- grp_summary(D, "class", "RNA_rpm", CLASSES[4])                        # expression by class in the set
show <- function(nm, g, unit = "") { .add("  ", nm, ":"); for (i in seq_len(nrow(g$table))) .add(sprintf("    %-24s n = %5d, median %s%s, %.2fx vs No TuORF (%.2f to %.2f), letter %s", g$table$group[i], g$table$n[i], f2(g$table$median[i], 3), unit, g$table$fold_vs_ref[i], g$table$fold_lo[i], g$table$fold_hi[i], g$table$KS_letter[i]))
  for (i in seq_len(nrow(g$pairs))) .add(sprintf("      %s vs %s: KS p = %.2g, Wilcoxon p = %.2g, median fold %.2f", g$pairs$group1[i], g$pairs$group2[i], g$pairs$KS_p[i], g$pairs$Wilcoxon_p[i], g$pairs$median_fold_1_over_2[i])) }
show("uncapped ratio (PARE) by uIF class", gU); show("uncapped ratio (PARE) by %uORF+ group", gUp)
show("uncapped ratio, library a only, by uIF class", gUa); show("uncapped ratio, library b only, by uIF class", gUb)
show("uncapped ratio with our RNA-seq as denominator, by uIF class", gUl); show("GMUCT uncapped ratio (A11) by uIF class, genes in both sets", gG)
show("mORF TE by uIF class, uncapped-ratio set", gT); show("RNA-seq RPM (Nagarajan et al.) by uIF class, uncapped-ratio set", gR)
D[, rna_tertile := cut(RNA_rpm, quantile(RNA_rpm, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low RNA", "mid RNA", "high RNA"))]
strat <- D[, .(n = .N, median_U = median(U)), by = .(rna_tertile, class)][order(rna_tertile, class)]
strat_w <- dcast(strat, class ~ rna_tertile, value.var = "median_U")
.add("  median uncapped ratio by class within RNA-seq abundance tertiles:"); for (i in seq_len(nrow(strat_w))) .add(sprintf("    %-24s low %.3f, mid %.3f, high %.3f", strat_w$class[i], strat_w$`low RNA`[i], strat_w$`mid RNA`[i], strat_w$`high RNA`[i]))

################################################################################
# 9) Figure
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(0.3, "cm"),
                                      legend.text = element_text(size = 7), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3), legend.margin = margin(0, 0, 0, 0))
pal4 <- c("#F8766D", "#00BA38", "#619CFF", "grey40")      # V38 (2026-11-11): the colours of Supplemental Fig. S2B (High, Medium, Low uIF; No TuORF)
pal_pup <- c("#619CFF", "#00BA38", "#F8766D", "grey40")   # V38: the colours of Figure 2F (0-30%, 30-80%, 80-100% uORF+; No TuORF). Before, the 0-30% group had the red of the 80-100% group of Figure 2F.
lab_plain <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
scat <- function(d, x, y, xl, yl, title, sub, v) {
  d <- d[is.finite(d[[x]]) & is.finite(d[[y]]) & d[[x]] > 0 & d[[y]] > 0, ]
  ggplot(d, aes(x = .data[[x]], y = .data[[y]])) + geom_point(alpha = 0.25, size = 0.7, colour = "grey30") +
    geom_smooth(method = "lm", se = FALSE, colour = "#D7301F", linewidth = 0.7, formula = y ~ x) + scale_x_log10(labels = lab_plain) + scale_y_log10(labels = lab_plain) +
    labs(title = title, subtitle = sprintf("%s; Spearman rho = %.2f (%.2f to %.2f), n = %s", sub, v[["rho"]], v[["lo"]], v[["hi"]], trimws(big(v[["n"]]))), x = xl, y = yl) + th   # statistics in the subtitle, never over the points
}
ecdfp <- function(d, grp, val, g, xl, title, sub, pal = pal4, xlim = NULL) {
  d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]) & d[[val]] > 0, ]; d$lab <- factor(g$labels[as.character(d[[grp]])], g$labels)
  p <- ggplot(d, aes(x = .data[[val]], colour = lab)) + stat_ecdf(geom = "step", linewidth = 0.8) + scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(title = title, subtitle = sub, x = xl, y = "Fraction of genes") + th + guides(colour = guide_legend(ncol = 1)) + scale_x_log10(labels = lab_plain)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}
U_lab <- "Uncapped ratio (PARE / RNA-seq, reads per million)"
xl <- unname(quantile(D$U, c(0.005, 0.995)))                                   # ECDF view: central 99% of the set
clip_txt <- sprintf("x axis clipped to %s to %s", signif(xl[1], 2), signif(xl[2], 2))
pA <- scat(S, "uIF", "U", "uIF", "Uncapped ratio", "A", sprintf("Single-TuORF genes, >= %d RNA-seq reads (Nagarajan et al. 2019); %d with uIF = 0 not drawn", MIN_RNA, sum(S$uIF == 0)), v1)
pB <- ecdfp(D, "class", "U", gU, U_lab, "B", paste0("Letters: pairwise Kolmogorov–Smirnov tests, p < 0.05 (as Fig. S2B); ", clip_txt), xlim = xl)
pC <- scat(S, "TE", "U", "mORF TE", "Uncapped ratio", "C", "Same genes as A", v2)
pD <- ecdfp(D, "pclass", "U", gUp, U_lab, "D", paste0("By %uORF+ group (single-TuORF genes) and No-TuORF genes; ", clip_txt), pal = pal_pup, xlim = xl)
pE <- scat(D[is.finite(U_GMUCT_nag)], "U_GMUCT_nag", "U", "GMUCT (Yu et al. 2021) / RNA-seq (Nagarajan et al. 2019), reads per million", U_lab, "E",
           sprintf("Genes in both sets, same denominator for both ratios (own RNA-seq each: rho = %.2f)", v9[["rho"]]), v9b)
pF <- scat(D[cohort == "No TuORF"], "TE", "U", "mORF TE", U_lab, "F", "No-TuORF genes", v7)
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
g <- arrangeGrob(patchwork::patchworkGrob(patchwork::wrap_plots(lapply(list(pA, pB, pC, pD, pE, pF), function(p) compact(p) + ggplot2::theme(panel.heights = grid::unit(1.25, "in"), panel.widths = grid::unit(2.0, "in"))), ncol = 2)), top = fig_title("Figure S20. Uncapped mRNA measured by PARE versus uIF and mORF TE"))   # V21 (2026-11-07, Larry): six plotting areas of equal size, panels aligned
fig_save(file.path(OUT_DIR, "A12_Figure_PARE_uncapped_vs_uIF.pdf"), g, h = 9)

################################################################################
# 10) Tables and summary
################################################################################
genes_out <- G[order(cohort != "single-TuORF", class, gene_id), .(gene_id, cohort, uIF_class = class, pct_uORF_plus_group = pclass, uIF, pct_uORF_plus = pUP, mORF_TE = TE, RNA_TPM_this_study = RNA_TPM,
  PARE_reads_a = PARE_a, PARE_reads_b = PARE_b, RNA_reads_rep1 = RNA_rep1, PARE_RPM = PARE_rpm, RNA_RPM = RNA_rpm,
  uncapped_ratio = U, uncapped_ratio_a = U_a, uncapped_ratio_b = U_b, uncapped_ratio_ourRNA = U_labRNA, uncapped_ratio_GMUCT_A11 = U_GMUCT, uncapped_ratio_GMUCT_over_Nagarajan_RNA = U_GMUCT_nag,
  in_uncapped_ratio_set = RNA_reads >= MIN_RNA & PARE_reads >= 1, half_life_h = half_life)]
fwrite(genes_out, file.path(OUT_DIR, "A12_Table_genes.csv"))
libstats <- data.table(library = c("PARE a (WT_polyA+a; GEO replicate 2, ATH237)", "PARE b (WT_polyA+b; GEO replicate 1, ATH362)", "RNA-seq replicate 1 (SRR7812609; unique reads)"),
  reads_total = c(pare$PARE_a$total, pare$PARE_b$total, sum(idx$mapped)), reads_chr1to5 = c(pare$PARE_a$nuclear_chr, pare$PARE_b$nuclear_chr, sum(idx$mapped[idx$seqnames %in% as.character(1:5)])),
  assigned_to_one_gene = c(pare$PARE_a$assigned, pare$PARE_b$assigned, sum(rna_counts)), assigned_nuclear_genes = c(libsize[["PARE_a"]], libsize[["PARE_b"]], libsize[["RNA_rep1"]]),
  ambiguous_two_genes = c(pare$PARE_a$ambiguous, pare$PARE_b$ambiguous, NA), note = c("", "", sprintf("%.1f%% of reads over isolated genes on the gene strand; counted as %s", 100 * st[["same"]], lib_type)))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S11. Uncapped mRNA (PARE) versus uIF and mORF TE (Figure S20): repeat of the GMUCT analysis (Figure S19 / Table S10) with the poly(A)+ PARE libraries of Nagarajan et al. 2019 Nucleic Acids Res 47:9216-9230 (GEO GSE119706; two-week-old Col-0 seedlings).",
  "PARE: WT_polyA+a (GEO replicate 2, ATH237, SRR8429733) and WT_polyA+b (replicate 1, ATH362); read 5' ends counted within the exons of each Araport11 gene (union over isoforms, same strand; 5' ends overlapping two genes discarded).",
  "RNA-seq denominator: WT poly(A)+ RNA-seq replicate 1 of the same study (SRR7812609), uniquely mapped reads (STAR, MAPQ 255) counted the same way (one replicate; the alignment of the second replicate was incomplete).",
  "Reads per million (RPM) use the reads assigned to nuclear genes (chromosomes 1-5) of each library; the uncapped ratio of a gene is the mean PARE RPM of the two libraries divided by the RNA-seq RPM (= 2 to the power of the 'proportion uncapped' of Yu et al. 2021).",
  sprintf("Genes enter the uncapped-ratio analyses with at least %d RNA-seq reads and at least one PARE read; cutoffs of 5, 20 and 50 reads and single PARE libraries are in 'cutoff_sensitivity' (rho_uIF, rho_TE and rho_pUP = Spearman rho of the uncapped ratio with uIF, mORF TE and %%uORF+). In the sheet 'genes', uncapped_ratio is NA (empty cell) for genes without RNA-seq reads; uncapped_ratio_a and uncapped_ratio_b are NA when the gene has no RNA-seq read or no read in that PARE library; uncapped_ratio_ourRNA = PARE RPM divided by the RNA TPM of this study; uncapped_ratio_GMUCT_A11 = the uncapped ratio of Table S10 (GMUCT), given for the genes of its uncapped-ratio set; uncapped_ratio_GMUCT_over_Nagarajan_RNA = GMUCT RPM divided by the RNA-seq RPM of this table, for the same genes. 'A11' in a row label means the GMUCT analysis of Table S10.", MIN_RNA),
  "Cohorts: single-TuORF genes of Table S3 with their uIF classes (Supplemental Fig. S2B) and %uORF+ groups; No-TuORF genes as in Supplemental Fig. S2B (annotated 5' leader, RNA TPM > 1, TE defined, no TuORF and no overlapping uORF).",
  "Sheets: correlations (Spearman rho with bootstrap 95% intervals over genes, 1,000 resamples; the column pearson_r_log gives the Pearson r of the log10 values for some rows, and genes with a value of 0 are not used for it), partial_correlations (rank-based partial correlations), cutoff_sensitivity, group_summaries (median, quartiles, median fold versus No TuORF with bootstrap 95% interval (1,000 resamples; fold_lo and fold_hi), letters from pairwise Kolmogorov–Smirnov (KS) tests at p < 0.05 as in Supplemental Fig. S2B), pairwise_tests, expression_strata (median uncapped ratio by class within RNA-seq abundance tertiles), library_stats, genes (per-gene values; the GMUCT ratio of Table S10 is included for comparison).",
  "Code: revision_modules/A12_PARE_uncapped_vs_uIF.R in the Supplemental Code.")), colNames = FALSE)
addWorksheet(wb, "correlations"); writeData(wb, "correlations", cor_tab)
addWorksheet(wb, "partial_correlations"); writeData(wb, "partial_correlations", part)
addWorksheet(wb, "cutoff_sensitivity"); writeData(wb, "cutoff_sensitivity", sens)
gs <- rbind(cbind(measure = "uncapped ratio (PARE)", grouping = "uIF class", gU$table), cbind(measure = "uncapped ratio (PARE)", grouping = "%uORF+ group", gUp$table),
            cbind(measure = "uncapped ratio, library a only", grouping = "uIF class", gUa$table), cbind(measure = "uncapped ratio, library b only", grouping = "uIF class", gUb$table),
            cbind(measure = "uncapped ratio, our RNA-seq as denominator", grouping = "uIF class", gUl$table), cbind(measure = "uncapped ratio (GMUCT, A11), genes in both sets", grouping = "uIF class", gG$table),
            cbind(measure = "mORF TE (uncapped-ratio set)", grouping = "uIF class", gT$table), cbind(measure = "RNA-seq RPM, Nagarajan et al. (uncapped-ratio set)", grouping = "uIF class", gR$table))
addWorksheet(wb, "group_summaries"); writeData(wb, "group_summaries", gs)
pt <- rbind(cbind(measure = "uncapped ratio (PARE)", grouping = "uIF class", gU$pairs), cbind(measure = "uncapped ratio (PARE)", grouping = "%uORF+ group", gUp$pairs),
            cbind(measure = "uncapped ratio, library a only", grouping = "uIF class", gUa$pairs), cbind(measure = "uncapped ratio, library b only", grouping = "uIF class", gUb$pairs),
            cbind(measure = "uncapped ratio, our RNA-seq as denominator", grouping = "uIF class", gUl$pairs), cbind(measure = "uncapped ratio (GMUCT, A11), genes in both sets", grouping = "uIF class", gG$pairs),
            cbind(measure = "mORF TE (uncapped-ratio set)", grouping = "uIF class", gT$pairs), cbind(measure = "RNA-seq RPM, Nagarajan et al. (uncapped-ratio set)", grouping = "uIF class", gR$pairs))
addWorksheet(wb, "pairwise_tests"); writeData(wb, "pairwise_tests", pt)
addWorksheet(wb, "expression_strata"); writeData(wb, "expression_strata", strat)
addWorksheet(wb, "library_stats"); writeData(wb, "library_stats", libstats)
genes_x <- copy(genes_out); genes_x[!is.finite(uncapped_ratio), uncapped_ratio := NA_real_]   # V41 (2026-11-12): the ratio is not defined without RNA-seq reads (Inf or NaN, shown as #NUM! in Excel); the sheet gets an empty cell there, the csv file is unchanged
addWorksheet(wb, "genes"); writeData(wb, "genes", genes_x)
saveWorkbook(wb, file.path(OUT_DIR, "A12_Table_PARE_uncapped.xlsx"), overwrite = TRUE)
.add(sprintf("A12 done in %.1f min; outputs: output/A12_Figure_PARE_uncapped_vs_uIF.pdf, A12_Table_PARE_uncapped.xlsx, A12_Table_genes.csv, A12_gene_counts.csv, A12_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A12_summary.txt"))
